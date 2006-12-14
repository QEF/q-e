!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------------
SUBROUTINE reset_k_points()
  !-----------------------------------------------------------------------------
  !
  USE klist,             ONLY : nks
  USE ktetra,            ONLY : nk1_   => nk1, &
                                nk2_   => nk2, &
                                nk3_   => nk3, &
                                k1_    => k1,  &
                                k2_    => k2,  &
                                k3_    => k3
  USE klist,             ONLY : lxkcry, &
                                xk_    => xk, &
                                wk_    => wk
  USE input_parameters,  ONLY : k_points, xk, wk, nk1, nk2, nk3, &
                                k1, k2, k3, nkstot
  !
  IMPLICIT NONE
  !
  ! 
  IF ( k_points == 'automatic' ) THEN
    !
    ! ... automatic generation of k-points
    !
    lxkcry = .FALSE.
    nks  = 0
    nk1_ = nk1
    nk2_ = nk2
    nk3_ = nk3
    k1_  = k1
    k2_  = k2
    k3_  = k3
    !
  ELSE IF ( k_points == 'tpiba' ) THEN
    !
    ! ... input k-points are in 2pi/a units
    !
    lxkcry = .FALSE.
    nks = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    !
  ELSE IF ( k_points == 'crystal' ) THEN
    !
    ! ... input k-points are in crystal (reciprocal lattice) axis
    !
    lxkcry = .TRUE.
    nks = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    !
  ELSE IF ( k_points == 'gamma' ) THEN
    !
    ! ... Only Gamma (k=0) is used
    !
    lxkcry = .FALSE.
    nks = 1
    xk_(:,1) = 0.0d0
    wk_(1)   = 1.0d0
    !
  ELSE
    !
    ! ... default: input k-points are in 2pi/a units
    !
    lxkcry = .FALSE.
    nks  = nkstot
    xk_(:,1:nks) = xk(:,1:nks)
    wk_(1:nks)   = wk(1:nks)
    !
  END IF  
  !
END SUBROUTINE reset_k_points
