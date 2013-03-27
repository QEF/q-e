!
! Copyright (C) 2008-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
  SUBROUTINE set_small_group_of_q(nsymq, invsymq, minus_q)
!-----------------------------------------------------------------------
!
!  This routine is a driver that sets the small group of q. It rotates
!  the matrices s so that the first nsymq elements are the small group
!  of q, and tells to the colling code is minus_q is true.
!  It deals also with the case modenum /= 0
!
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE ions_base, ONLY : nat, tau
  USE symm_base, ONLY : s, nsym, ftau, irt, time_reversal
  USE control_flags, ONLY : modenum
  USE qpoint, ONLY : xq
  USE symm_base, ONLY : copy_sym, d1, d2, d3, inverse_s, s_axis_to_cart
  USE paw_variables, ONLY : okpaw
  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: nsymq
  LOGICAL, INTENT(INOUT) :: minus_q, invsymq
  !
  REAL(DP), ALLOCATABLE :: rtau(:,:,:)

  LOGICAL :: sym(48)

  sym(1:nsym)=.true.
  call smallg_q (xq, modenum, at, bg, nsym, s, ftau, sym, minus_q)
  IF ( .not. time_reversal ) minus_q = .false.
  IF (modenum /= 0) THEN
!
!   in this case remove also the symmetries that do not send the mode
!   in itself
!
     ALLOCATE(rtau (3, 48, nat))
     CALL sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
     CALL mode_group (modenum, xq, at, bg, nat, nsym, s, irt, minus_q, &
                                                            rtau, sym)
     DEALLOCATE(rtau)
  ENDIF
  nsymq = copy_sym ( nsym, sym )
  
  call inverse_s ( )
  !
  ! check if inversion (I) is a symmetry. If so, there should be nsymq/2
  ! symmetries without inversion, followed by nsymq/2 with inversion
  ! Since identity is always s(:,:,1), inversion should be s(:,:,1+nsymq/2)
  !
  invsymq = ALL ( s(:,:,nsymq/2+1) == -s(:,:,1) )
  !
  !  Since the order of the s matrices is changed we need to recalculate:
  !
  call s_axis_to_cart ( )
  IF (okpaw) CALL d_matrix(d1,d2,d3)

  RETURN
  END SUBROUTINE set_small_group_of_q
!
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
  USE kinds, ONLY : DP
  implicit none

  real(DP), parameter :: accep = 1.e-5_dp

  real(DP), intent(in) :: bg (3, 3), at (3, 3), xq (3)
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  ! input: the q point of the crystal

  integer, intent(in) :: s (3, 3, 48), nrot, ftau (3, 48), modenum
  ! input: the symmetry matrices
  ! input: number of symmetry operations
  ! input: fft grid dimension (units for ftau)
  ! input: fractionary translation of each symmetr
  ! input: main switch of the program, used for
  !        q<>0 to restrict the small group of q
  !        to operation such that Sq=q (exactly,
  !        without G vectors) when iswitch = -3.
  logical, intent(inout) :: sym (48), minus_q
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
     sym (irot) = eqvect (raq, aq, zero, accep)
     !
     !  if "iswitch.le.-3" (modenum.ne.0) S must be such that Sq=q exactly !
     !
     if (modenum.ne.0 .and. sym(irot) ) then
        do ipol = 1, 3
           sym(irot) = sym(irot) .and. (abs(raq(ipol)-aq(ipol)) < 1.0d-5)
        enddo
     endif
!     if (.not.minus_q) then
     if (sym(irot).and..not.minus_q) then
        raq = - raq
        minus_q = eqvect (raq, aq, zero, accep)
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

