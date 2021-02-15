!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE uspp_data
  !
  !! These parameters are needed with the US pseudopotentials.
  !
  USE upf_kinds,      ONLY : DP
  !
  SAVE
  PRIVATE
  !
  PUBLIC :: nqxq, nqx, dq, spline_ps
  PUBLIC :: qrad,   tab,   tab_at,   tab_d2y
  PUBLIC :: qrad_d, tab_d, tab_at_d, tab_d2y_d
  PUBLIC :: deallocate_uspp_data
  !
  INTEGER :: nqxq
  !! size of interpolation table
  INTEGER :: nqx 
  !! number of interpolation points
  REAL(DP), PARAMETER:: dq = 0.01D0
  !! space between points in the pseudopotential tab.
  REAL(DP), ALLOCATABLE :: qrad(:,:,:,:)
  !! radial FT of Q functions
  REAL(DP), ALLOCATABLE :: tab(:,:,:)
  !! interpolation table for PPs
  REAL(DP), ALLOCATABLE :: tab_at(:,:,:)
  !! interpolation table for atomic wfc
  LOGICAL :: spline_ps = .FALSE.
  REAL(DP), ALLOCATABLE :: tab_d2y(:,:,:)
  !! for cubic splines
  !
  ! GPUs vars
  !
  REAL(DP), ALLOCATABLE :: qrad_d(:,:,:,:)
  REAL(DP), ALLOCATABLE :: tab_d(:,:,:)
  REAL(DP), ALLOCATABLE :: tab_at_d(:,:,:)
  REAL(DP), ALLOCATABLE :: tab_d2y_d(:,:,:)
  !   
#if defined(__CUDA)
  attributes (DEVICE) :: qrad_d, tab_d, tab_at_d, tab_d2y_d
#endif
  !
contains
  !
  subroutine deallocate_uspp_data()
     IMPLICIT NONE
     IF( ALLOCATED( qrad ) )      DEALLOCATE( qrad )
     IF( ALLOCATED( tab ) )       DEALLOCATE( tab )
     IF( ALLOCATED( tab_at ) )    DEALLOCATE( tab_at )
     IF( ALLOCATED( tab_d2y ) )   DEALLOCATE( tab_d2y )
     !
     IF( ALLOCATED( qrad_d ) )    DEALLOCATE( qrad_d )
     IF( ALLOCATED( tab_d ) )     DEALLOCATE( tab_d )
     IF( ALLOCATED( tab_at_d ) )  DEALLOCATE( tab_at_d )
     IF( ALLOCATED( tab_d2y_d ) ) DEALLOCATE( tab_d2y_d )
  end subroutine 
  !
END MODULE uspp_data

