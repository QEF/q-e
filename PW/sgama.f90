!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine sgama (nrot, nat, s, sname, t_rev, at, bg, tau, ityp, nsym,&
     nr1, nr2, nr3, irt, ftau, npk, nks, xk, wk, invsym, minus_q, xq, &
     modenum, time_reversal, magnetic_sym, m_loc)
  !-----------------------------------------------------------------------
  !
  !     This routine performs the following tasks:
  !     1)  It finds the point group of the crystal, by eliminating the
  !         symmetries of the Bravais lattice which are not allowed
  !         by the atomic positions.
  !     1a) If xq.ne.0 it restricts the symmetries to those of the small
  !         group of q. In this case the small group of q is determined
  !         seeking all sym.op. such that Sq=q+G (and Sq=-q+G is also
  !         considered) when iswitch=-2, while only sym.op. such that
  !         Sq=q (exactly, without G) when iswitch=-3.
  !     1b) if iswitch.eq.-4 it keep only the symmetries which send
  !         a mode in itself. The mode is given by modenum
  !     2)  It finds the special points in the irreducible wedge of the
  !         true point group (or small group of q) of the crystal starting
  !         from the points in the irreducible wedge of the point group
  !         of the Bravais lattice.
  !     3)  It checks if the point group has the inversion symmetry.
  !
  !     This routine is mainly the driver of separate routines which
  !     perform each single task.
  !
  !     Modified by SdG to include the "small group of q" stuff for the
  !     linear-response preparation run.
  !
#include "f_defs.h"
  USE kinds, only : DP
  implicit none
  !
  !    First the I/O variables
  !

  integer :: nrot, nat, s (3, 3, 48), ityp (nat), nsym, nr1, nr2, &
       nr3, irt (48, nat), ftau (3, 48), npk, nks, modenum
  ! input: number of symmetries of the original
  ! input: number of atoms in the cell
  ! input: matrices of the symmetry operations
  ! input: type of each atom
  ! output: real number of symmetry operations
  !  input: dimensions of the fft mesh
  ! output: for each atom gives the rotated ato
  ! output: fractionary translation of each sym
  ! input: maximum number of k points
  ! input-output: starting and ending number of
  ! input: the mode to be computed
  ! input: main switch of the program; used whe
  !        xq<>0 to restrict the small group of

  real(DP) :: at (3, 3), bg (3, 3), tau (3, nat), xk (3, npk), &
       wk (npk), xq (3), m_loc(3,nat)
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  ! input: coordinates of atomic positions
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  ! input: coordinates of a q-point
  logical, intent(in) :: time_reversal, magnetic_sym
  ! time_reversal=true : use time-reversal symmetry (q=>-q)
  ! magnetic_sym =true : find symmetries that leave magnetization unchanged
  !
  logical, intent(out) :: invsym, minus_q
  ! output: if true the crystal has inversion
  ! output: if true a symmetry sends q->-q+G
  character :: sname (48) * 45
  ! input: name of the rotation part of each symmetry operation
  !
  !    And then the local variables
  !

  real(DP), allocatable :: rtau (:,:,:)
  ! direct translations of each point
  integer :: table (48, 48), irot, jrot, ipol, jpol, invs (3, 3, 48) &
       , irg (48), temp, na
  ! multiplication table of the group
  ! counter over the rotations
  ! counter over the rotations
  ! counter over the polarizations
  ! counter over the polarizations
  ! contains the inverse of each rotation
  ! gives the correspondence of symmetry
  ! operations forming a n-th coset
  ! auxilary variable
  ! counter on atoms

  integer :: t_rev(48)
  ! for magnetic symmetries: if 1 there is time reversal operation
  logical :: sym (48)
  ! if true the corresponding operation is a symmetry operation

  allocate(rtau (3, 48, nat))
  !
  !    Here we find the true symmetries of the crystal
  !
  IF ( magnetic_sym ) THEN
     CALL sgam_at_mag (nrot, s, nat, tau, ityp, at, bg, &
                  nr1, nr2, nr3, sym, irt, ftau, m_loc, sname, t_rev)
  ELSE
     CALL sgam_at (nrot, s, nat, tau, ityp, at, bg, nr1, nr2, nr3, sym, &
       irt, ftau)
  END IF
  !
  !    If xq.ne.(0,0,0) this is a preparatory run for a linear response
  !    calculation at xq. The relevant point group is therefore only the
  !    small group of q. Here we exclude from the list the symmetries
  !    that do not belong to it
  !
  call smallg_q (xq, modenum, at, bg, nrot, s, ftau, sym, minus_q)
  !
  IF ( .not. time_reversal ) THEN
     !
     minus_q=.false.
!     IF ( ABS(DOT_PRODUCT(xq,xq)) > 1.0D-07 ) CALL errore ('sgama', &
!          'phonon not implemented with non collinear magnetism', 1)
     ! If somebody wants to implement phonon calculations in non 
     ! collinear magnetic case he/she has to pay attention to the
     ! fact that in non collinear case the symmetry k -> -k is not
     ! always allowed as in collinear case. Adriano
  ENDIF
  !
  if (modenum .ne. 0) then
     call sgam_ph (at, bg, nrot, s, irt, tau, rtau, nat, sym)
     call mode_group (modenum, xq, at, bg, nat, nrot, s, irt, rtau, &
          sym, minus_q)
  endif
  !
  !    We compute the multiplication table of the group
  !
  call multable (nrot, s, table)
  !
  !   And we set the matrices of the inverse
  !
  call inverse_s (nrot, s, table, invs)
  !
  !    Find the coset in the point group of the Bravais lattice
  !
  call coset (nrot, table, sym, nsym, irg)
  !
  !    here we set the k-points in the irreducible wedge of the point grou
  !    of the crystal
  !
  call irrek (npk, nks, xk, wk, at, bg, nrot, invs, nsym, irg, minus_q)
  !
  ! copy symm. operations in sequential order so that
  ! s(i,j,irot) , irot <= nsym          are the sym.ops. of the crystal
  !               nsym+1 < irot <= nrot are the sym.ops. of the lattice
  !
  jrot = 0
  do irot = 1, nrot
     if (sym (irot) ) then
        jrot = jrot + 1
        do ipol = 1, 3
           do jpol = 1, 3
              temp = s (ipol, jpol, jrot)
              s (ipol, jpol, jrot) = s (ipol, jpol, irot)
              s (ipol, jpol, irot) = temp
           enddo
           ftau (ipol, jrot) = ftau (ipol, irot)
        enddo
        do na = 1, nat
           irt (jrot, na) = irt (irot, na)
        enddo
        sname (jrot) = sname (irot)
        t_rev(jrot) = t_rev(irot)
     endif
  enddo
  if (jrot.ne.nsym) call errore ('sgama', 'unexpected', 1)
  !
  ! Sets to zero the first matrix that is not a symmetry of the crystal.
  ! This will be used by d3toten program.
  !
  if (nrot.lt.48) then
     do ipol = 1, 3
        do jpol = 1, 3
           s (ipol, jpol, nrot + 1) = 0
        enddo
     enddo
  endif
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = .true.
  irot = nsym / 2 + 1
  do ipol = 1, 3
     do jpol = 1, 3
        invsym = invsym.and.s (ipol, jpol, irot) .eq. - s (ipol, jpol, 1)
     enddo

  enddo

  deallocate (rtau)
  return

end subroutine sgama
!-----------------------------------------------------------------------

