!
! Copyright (C) 2005 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE coarsegrained_vars
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: max_fe_iter    = 10
  INTEGER, PARAMETER :: max_shake_iter = 5
  INTEGER, PARAMETER :: num_acc        = 8
  !
  REAL (KIND=DP), PARAMETER :: fe_step     = 0.2D0
  REAL (KIND=DP), PARAMETER :: fe_grad_thr = 1.D-4
  !
  REAL (KIND=DP), ALLOCATABLE :: dfe_acc(:,:)
  REAL (KIND=DP), ALLOCATABLE :: fe_grad(:)
  REAL (KIND=DP), ALLOCATABLE :: new_target(:)
  REAL (KIND=DP), ALLOCATABLE :: to_target(:)
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_coarsegrained_vars( nconstr )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nconstr
      !
      !
      ALLOCATE( dfe_acc( nconstr, num_acc ) )
      !
      ALLOCATE( fe_grad(    nconstr ) )
      ALLOCATE( new_target( nconstr ) )
      ALLOCATE( to_target(  nconstr ) )
      !
      RETURN
      !
    END SUBROUTINE allocate_coarsegrained_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_coarsegrained_vars()
      !------------------------------------------------------------------------
      !
      IF ( ALLOCATED( dfe_acc ) )     DEALLOCATE( dfe_acc )
      IF ( ALLOCATED( fe_grad ) )     DEALLOCATE( fe_grad )
      IF ( ALLOCATED( new_target ) )  DEALLOCATE( new_target )
      IF ( ALLOCATED( to_target ) )   DEALLOCATE( to_target )
      !
      RETURN
      !
    END SUBROUTINE deallocate_coarsegrained_vars
    !
END MODULE coarsegrained_vars
