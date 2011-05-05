!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------
SUBROUTINE reset_k_points()
  !-----------------------------------------------------------------------------
  ! ... Copy input data for k-points into internal variables
  ! ... Used both at startup and in neb/strings calculations
  !
  USE start_k,          ONLY : nk1_   => nk1, &
                                nk2_   => nk2, &
                                nk3_   => nk3, &
                                k1_    => k1,  &
                                k2_    => k2,  &
                                k3_    => k3
  USE klist,             ONLY : lxkcry, &
                                xk_     => xk, &
                                wk_     => wk, &
                                nkstot_ => nkstot
  USE input_parameters,  ONLY : k_points, xk, wk, nk1, nk2, nk3, &
                                k1, k2, k3, nkstot
  !
  IMPLICIT NONE
  !
  nk1_ = 0
  nk2_ = 0
  nk3_ = 0
  k1_  = 0
  k2_  = 0
  k3_  = 0
  lxkcry = .FALSE.
  !
  IF ( k_points == 'automatic' ) THEN
    !
    ! ... automatic generation of k-points
    !
    nkstot_  = 0
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
    nkstot_ = nkstot
    xk_(:,1:nkstot_) = xk(:,1:nkstot_)
    wk_(1:nkstot_)   = wk(1:nkstot_)
    !
  ELSE IF ( k_points == 'crystal' ) THEN
    !
    ! ... input k-points are in crystal (reciprocal lattice) axis
    !
    lxkcry = .TRUE.
    nkstot_ = nkstot
    xk_(:,1:nkstot_) = xk(:,1:nkstot_)
    wk_(1:nkstot_)   = wk(1:nkstot_)
    !
  ELSE IF ( k_points == 'gamma' ) THEN
    !
    ! ... Only Gamma (k=0) is used
    !
    nkstot_ = 1
    xk_(:,1) = 0.0d0
    wk_(1)   = 1.0d0
    !
  ELSE
    !
    ! ... default: input k-points are in 2pi/a units
    !
    nkstot_  = nkstot
    xk_(:,1:nkstot_) = xk(:,1:nkstot_)
    wk_(1:nkstot_)   = wk(1:nkstot_)
    !
  END IF  
  ! 
END SUBROUTINE reset_k_points
