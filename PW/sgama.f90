  !!
! Copyright (C) 2001-2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine sgama (nrot, nat, s, sname, t_rev, at, bg, tau, ityp, nsym,&
     nr1, nr2, nr3, irt, ftau, invsym, minus_q, xq, &
     modenum, time_reversal, magnetic_sym, m_loc)
  !-----------------------------------------------------------------------
  !
  !     This routine performs the following tasks:
  !     1)  It finds the point group of the crystal, by eliminating the
  !         symmetries of the Bravais lattice which are not allowed
  !         by the atomic positions.
  !     2)  If xq.ne.0 it restricts the symmetries to those of the small
  !         group of q. In this case the small group of q is determined
  !         seeking all sym.op. such that Sq=q+G (and Sq=-q+G is also
  !         considered) 
  !     3)  if modenum.ne.0 keep only symmetries which send mode
  !         "modenum" into itself
  !     4)  It checks if the point group has the inversion symmetry.
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
  integer, intent(in) :: nrot, nat, ityp (nat), nr1, nr2, nr3,  modenum
  real(DP), intent(in) :: at (3,3), bg (3,3), tau (3,nat), xq (3), m_loc(3,nat)
  logical, intent(in) :: time_reversal, magnetic_sym
  !
  character(len=45), intent(inout) :: sname (48)
  ! name of the rotation part of each symmetry operation
  integer, intent(inout) :: s(3,3,48)
  !
  integer, intent(out) :: nsym, irt (48, nat), ftau (3, 48)
  logical, intent(out) :: invsym, minus_q
  ! minus_q : if true a symmetry sends q->-q+G
  !
  real(DP), allocatable :: rtau (:,:,:)
  ! direct translations of each point
  integer :: stemp(3,3), irot, jrot, ipol, jpol, na
  ! counters
  integer :: t_rev(48)
  ! for magnetic symmetries: if 1 there is time reversal operation
  logical :: sym (48)
  ! if true the corresponding operation is a symmetry operation
  !
  !    Here we find the true symmetries of the crystal
  !
  CALL sgam_at (nrot, s, nat, tau, ityp, at, bg, nr1, nr2, nr3, &
                sym, irt, ftau)
  IF ( magnetic_sym ) &
     CALL sgam_at_mag (nrot, s, nat, bg, irt, m_loc, sname, sym, t_rev)
  !
  !    If xq.ne.(0,0,0) this is a preparatory run for a linear response
  !    calculation at xq. The relevant point group is therefore only the
  !    small group of q. Here we exclude from the list the symmetries
  !    that do not belong to it
  !
  call smallg_q (xq, modenum, at, bg, nrot, s, ftau, sym, minus_q)
  !
  IF ( .not. time_reversal ) minus_q = .false.
  !
  ! If somebody wants to implement phonon calculations in non 
  ! collinear magnetic case he/she has to pay attention to the
  ! fact that in non collinear case the symmetry k -> -k is not
  ! always allowed as in collinear case. Adriano
  !
  if (modenum /= 0) then
     allocate(rtau (3, 48, nat))
     call sgam_ph (at, bg, nrot, s, irt, tau, rtau, nat, sym)
     call mode_group (modenum, xq, at, bg, nat, nrot, s, irt, rtau, &
          sym, minus_q)
     deallocate (rtau)
  endif
  !
  ! copy symm. operations in sequential order so that
  ! s(i,j,irot) , irot <= nsym          are the sym.ops. of the crystal
  !               nsym+1 < irot <= nrot are the sym.ops. of the lattice
  jrot = 0
  do irot = 1, nrot
     if (sym (irot) ) then
        jrot = jrot + 1
        stemp = s(:,:,jrot)
        s (:,:, jrot) = s (:,:, irot)
        s (:,:, irot) = stemp
        ftau (:, jrot) = ftau (:, irot)
        irt (jrot,1:nat) = irt (irot,1:nat)
        sname (jrot) = sname (irot)
        t_rev (jrot) = t_rev(irot)
     endif
  enddo
  nsym = jrot
  sym (1:nsym) = .true.
  sym (nsym+1:nrot) = .false.
  !
  ! Sets to zero the first matrix that is not a symmetry of the crystal.
  ! This will be used by d3toten program.
  !
  if (nrot < 48) s(:,:, nrot+1) = 0
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  irot = nsym/2+1
  invsym = ALL ( s(:,:,irot) == -s(:,:,1) )
  !
  return
  !
end subroutine sgama
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
subroutine irreducible_BZ (nrot, s, nsym, at, bg, npk, nks, xk, wk, minus_q)
  !-----------------------------------------------------------------------
  !
  !     This routine finds the special points in the irreducible wedge of the
  !     true point group (or small group of q) of the crystal, starting
  !     from the points in the irreducible wedge of the point group
  !     of the Bravais lattice.
  !
  USE kinds, only : DP
  implicit none
  !
  integer,  intent(in) :: nrot, nsym, npk, s(3,3,48)
  real(DP), intent(in) :: at (3,3), bg (3,3)
  logical,  intent(in) :: minus_q
  integer,  intent(inout) :: nks
  real(DP), intent(inout) :: xk (3, npk), wk (npk)
  !
  integer :: table (48, 48), invs (3, 3, 48), irg (48)
  ! table: multiplication table of the group
  ! invs : contains the inverse of each rotation
  ! irg  : gives the correspondence of symmetry operations forming a n-th coset
  logical :: sym(48)
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
  sym(1:nsym) = .true.
  sym(nsym+1:)= .false.
  call coset (nrot, table, sym, nsym, irg)
  !
  !    here we set the k-points in the irreducible wedge of the point grou
  !    of the crystal
  !
  call irrek (npk, nks, xk, wk, at, bg, nrot, invs, nsym, irg, minus_q)
  !
  return
  !
end subroutine irreducible_BZ 
!-----------------------------------------------------------------------

