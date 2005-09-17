!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
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
  INTEGER, PARAMETER :: max_fe_iter    = 50
  INTEGER, PARAMETER :: max_shake_iter = 5
  !
  REAL(DP), PARAMETER :: fe_step = 0.4D0
  !
  REAL(DP), ALLOCATABLE :: dfe_acc(:)
  REAL(DP), ALLOCATABLE :: fe_grad(:)
  REAL(DP), ALLOCATABLE :: new_target(:)
  REAL(DP), ALLOCATABLE :: to_target(:)
  !
  LOGICAL  :: to_new_target
  !
  ! ... Laio-Parrinello meta-dynamics
  !
  INTEGER, PARAMETER :: max_metadyn_iter = 500
  !
  REAL(DP), ALLOCATABLE :: metadyn_history(:,:)
  !
  REAL(DP), PARAMETER :: A            = 0.01D0
  REAL(DP), PARAMETER :: sigma        = 0.4D0
  REAL(DP), PARAMETER :: sigma_sq     = sigma**2
  REAL(DP), PARAMETER :: two_sigma_sq = 2.D0 * sigma_sq
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_coarsegrained_vars( nconstr, nstep )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,           INTENT(IN) :: nconstr
      INTEGER, OPTIONAL, INTENT(IN) :: nstep
      !
      !
      ALLOCATE( dfe_acc(    nconstr ) )
      ALLOCATE( fe_grad(    nconstr ) )
      ALLOCATE( new_target( nconstr ) )
      ALLOCATE( to_target(  nconstr ) )
      !
      IF ( PRESENT( nstep ) ) ALLOCATE( metadyn_history( nconstr, nstep ) )
      !
      RETURN
      !
    END SUBROUTINE allocate_coarsegrained_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_coarsegrained_vars()
      !------------------------------------------------------------------------
      !
      IF ( ALLOCATED( dfe_acc ) )          DEALLOCATE( dfe_acc )
      IF ( ALLOCATED( fe_grad ) )          DEALLOCATE( fe_grad )
      IF ( ALLOCATED( new_target ) )       DEALLOCATE( new_target )
      IF ( ALLOCATED( to_target ) )        DEALLOCATE( to_target )
      IF ( ALLOCATED( metadyn_history ) )  DEALLOCATE( metadyn_history )
      !
      RETURN
      !
    END SUBROUTINE deallocate_coarsegrained_vars
    !
END MODULE coarsegrained_vars
!
!----------------------------------------------------------------------------
MODULE coarsegrained_base
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: set_target, add_gaussians
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_target()
      !------------------------------------------------------------------------
      !
      USE coarsegrained_vars, ONLY : to_target, to_new_target, max_shake_iter
      USE constraints_module, ONLY : target
      !
      !
      IF ( to_new_target ) THEN
         !
         target(:) = target(:) + to_target(:) / DBLE( max_shake_iter )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE set_target
    !
    !------------------------------------------------------------------------
    SUBROUTINE add_gaussians( iter )
      !------------------------------------------------------------------------
      !
      USE constraints_module, ONLY : nconstr
      USE coarsegrained_vars, ONLY : metadyn_history, fe_grad, dfe_acc
      USE coarsegrained_vars, ONLY : A, sigma_sq, two_sigma_sq
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: iter
      !
      INTEGER               :: i
      REAL(DP), ALLOCATABLE :: delta(:)
      !
      !
      ! ... history dependent term
      !
      ALLOCATE( delta( nconstr ) )
      !
      dfe_acc = 0.D0
      !
      DO i = 1, iter - 1
         !
         delta = metadyn_history(:,i) - metadyn_history(:,iter)
         !
         dfe_acc(:) = dfe_acc(:) + delta(:) * &
                      EXP( - ( delta(:) .dot. delta(:) ) / two_sigma_sq )
         !
      END DO
      !
      fe_grad(:) = fe_grad(:) + A / sigma_sq * dfe_acc(:)
      !
      DEALLOCATE( delta )
      !
      RETURN
      !
    END SUBROUTINE add_gaussians
    !
END MODULE coarsegrained_base
