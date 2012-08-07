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

