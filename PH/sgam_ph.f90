!
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
  integer :: t_rev(48)
  ! for magnetic symmetries: if 1 there is time reversal operation
  logical :: sym (48)
  ! if true the corresponding operation is a symmetry operation
  integer, external :: copy_sym
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
  !    Here we re-order all rotations in such a way that true sym.ops
  !    are the first nsym; rotations that are not sym.ops. follow
  !
  nsym = copy_sym ( nrot, sym, s, sname, ftau, nat, irt, t_rev )
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = ALL ( s(:,:,nsym/2+1) == -s(:,:,1) )
  !
  return
  !
end subroutine sgama
!
!-----------------------------------------------------------------------
subroutine sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the vector rtau which contains for each
  !     atom and each rotation the vector S\tau_a - \tau_b, where
  !     b is the rotated a atom, given by the array irt. These rtau are
  !     non zero only if fractional translations are present.
  !
#include "f_defs.h"
  USE kinds
  implicit none
  !
  !     first the dummy variables
  !
  integer, intent(in) :: nsym, s (3, 3, 48), nat, irt (48, nat)
  ! nsym: number of symmetries of the point group
  ! s:    matrices of symmetry operations
  ! nat : number of atoms in the unit cell
  ! irt(n,m) = transformed of atom m for symmetry n
  real(DP), intent(in) :: at (3, 3), bg (3, 3), tau (3, nat)
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  ! tau: coordinates of the atoms
  logical, intent(in) :: sym (nsym)
  ! sym(n)=.true. if operation n is a symmetry
  real(DP), intent(out):: rtau (3, 48, nat)
  ! rtau: the direct translations
  !
  !    here the local variables
  !
  integer :: na, nb, isym, ipol
  ! counters on: atoms, symmetry operations, polarization
  real(DP) , allocatable :: xau (:,:)
  real(DP) :: ft (3)
  !
  allocate (xau(3,nat))    
  !
  !   compute the atomic coordinates in crystal axis, xau
  !
  do na = 1, nat
     do ipol = 1, 3
        xau (ipol, na) = bg (1, ipol) * tau (1, na) + &
                         bg (2, ipol) * tau (2, na) + &
                         bg (3, ipol) * tau (3, na)
     enddo
  enddo
  !
  !    for each symmetry operation, compute the atomic coordinates
  !    of the rotated atom, ft, and calculate rtau = Stau'-tau
  !
  do isym = 1, nsym
     if (sym (isym) ) then
        do na = 1, nat
           nb = irt (isym, na)
           do ipol = 1, 3
              ft (ipol) = s (1, ipol, isym) * xau (1, na) + &
                          s (2, ipol, isym) * xau (2, na) + &
                          s (3, ipol, isym) * xau (3, na) - xau (ipol, nb)
           enddo
           do ipol = 1, 3
              rtau (ipol, isym, na) = at (ipol, 1) * ft (1) + &
                                      at (ipol, 2) * ft (2) + &
                                      at (ipol, 3) * ft (3)
           enddo
        enddo
     endif
  enddo
  !
  !    deallocate workspace
  !
  deallocate(xau)
  return
end subroutine sgam_ph
!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
subroutine smallg_q (xq, modenum, at, bg, nrot, s, ftau, sym, minus_q)
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave q unchanged.
  ! Furthermore it checks if one of the above matrices send q --> -q+G.
  ! In this case minus_q is set true.
  !
  !  input-output variables
  !
  USE kinds
  implicit none

  real(DP) :: bg (3, 3), at (3, 3), xq (3)
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  ! input: the q point of the crystal

  integer :: s (3, 3, 48), nrot, ftau (3, 48), modenum
  ! input: the symmetry matrices
  ! input: number of symmetry operations
  ! input: fft grid dimension (units for ftau)
  ! input: fractionary translation of each symmetr
  ! input: main switch of the program, used for
  !        q<>0 to restrict the small group of q
  !        to operation such that Sq=q (exactly,
  !        without G vectors) when iswitch = -3.
  logical :: sym (48), minus_q
  ! input-output: .true. if symm. op. S q = q + G
  ! output: .true. if there is an op. sym.: S q = - q + G
  !
  !  local variables
  !

  real(DP) :: aq (3), raq (3), zero (3)
  ! q vector in crystal basis
  ! the rotated of the q vector
  ! the zero vector

  integer :: irot, ipol, jpol
  ! counter on symmetry op.
  ! counter on polarizations
  ! counter on polarizations

  logical :: eqvect
  ! logical function, check if two vectors are equa
  !
  ! return immediately (with minus_q=.true.) if xq=(0,0,0)
  !
  minus_q = .true.
  if ( (xq (1) == 0.d0) .and. (xq (2) == 0.d0) .and. (xq (3) == 0.d0) ) &
       return
  !
  !   Set to zero some variables
  !
  minus_q = .false.
  zero(:) = 0.d0
  !
  !   Transform xq to the crystal basis
  !
  aq = xq
  call cryst_to_cart (1, aq, at, - 1)
  !
  !   Test all symmetries to see if this operation send Sq in q+G or in -q+G
  !
  do irot = 1, nrot
     if (.not.sym (irot) ) goto 100
     raq(:) = 0.d0
     do ipol = 1, 3
        do jpol = 1, 3
           raq(ipol) = raq(ipol) + DBLE( s(ipol,jpol,irot) ) * aq( jpol)
        enddo
     enddo
     sym (irot) = eqvect (raq, aq, zero)
     !
     !  if "iswitch.le.-3" (modenum.ne.0) S must be such that Sq=q exactly !
     !
     if (modenum.ne.0 .and. sym(irot) ) then
        do ipol = 1, 3
           sym(irot) = sym(irot) .and. (abs(raq(ipol)-aq(ipol)) < 1.0d-5)
        enddo
     endif
     if (sym (irot) .and..not.minus_q) then
        ! l'istruzione "originale" in kreductor era la seguente...
        !         if (.not. minus_q) then
        raq = - raq
        minus_q = eqvect (raq, aq, zero)
     endif
100  continue
  enddo
  !
  !  if "iswitch.le.-3" (modenum.ne.0) time reversal symmetry is not included !
  !
  if (modenum.ne.0) minus_q = .false.
  !
  return
end subroutine smallg_q

