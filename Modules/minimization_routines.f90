!
! Copyright (C) 2003-2004 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define USE_SMART_STEP
!#define DEBUG_SMART_STEP
!
!--------------------------------------------------------------------------
MODULE minimization_routines
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the optimization of the reaction path (NEB calculations)
  ! ... Written by Carlo Sbraccia ( 04-11-2003 )  
  !
  USE kinds,          ONLY :  DP
  USE constants,      ONLY :  au, eV_to_kelvin, eps32  
  USE neb_variables,  ONLY :  ds, pos, grad, norm_grad
  USE basic_algebra_routines
  !  
  IMPLICIT NONE
  !
  CONTAINS
     !
     ! ... steepest descent 
     !
     !----------------------------------------------------------------------
     SUBROUTINE steepest_descent( index )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( norm_grad(index) > eps32 ) THEN
          !
          pos(:,index) = pos(:,index) - ds(index) * grad(:,index)
          !
       END IF
       !      
       RETURN
       !
     END SUBROUTINE steepest_descent 
     !
     ! ... Molecular Dynamics based algorithms 
     ! ... velocity Verlet and quick min       
     !
     !----------------------------------------------------------------------
     SUBROUTINE velocity_Verlet_first_step( index )
       !---------------------------------------------------------------------- 
       !
       USE neb_variables, ONLY : vel, mass, vel_zeroed
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( vel_zeroed(index) ) THEN
          !
          vel(:,index) = vel(:,index) - 0.5D0 * grad(:,index) / mass(:)
          !
          pos(:,index) = pos(:,index) + vel(:,index)
          !
          vel_zeroed(index) = .FALSE.
          !
       ELSE
          !
          vel(:,index) = vel(:,index) - &
                         0.5D0 * ds(index) * grad(:,index) / mass(:)
          !
          pos(:,index) = pos(:,index) + ds(index) * vel(:,index)
          !
       END IF
       !       
       RETURN
       !
     END SUBROUTINE velocity_Verlet_first_step
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE velocity_Verlet_second_step( index )
       !----------------------------------------------------------------------
       !
       USE neb_variables, ONLY : vel, mass, damp, ldamped_dyn, lmol_dyn
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( ldamped_dyn ) THEN
          !
          vel(:,index) = damp * ( vel(:,index) - &
                                  0.5D0 * ds(index) * grad(:,index) / mass(:) )
          !
       ELSE IF ( lmol_dyn ) THEN
          !       
          vel(:,index) = vel(:,index) - &
                         0.5D0 * ds(index) * grad(:,index) / mass(:)
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE velocity_Verlet_second_step
     !
     !----------------------------------------------------------------------
     SUBROUTINE quick_min_second_step( index )
       !----------------------------------------------------------------------
       !
       USE constants,        ONLY : eps8
       USE neb_variables,    ONLY : pos_old, grad_old, vel, &
                                    mass, dim, vel_zeroed
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN)         :: index
       REAL (KIND=DP)              :: force_versor(dim)
       REAL (KIND=DP)              :: vel_component
       REAL (KIND=DP), ALLOCATABLE :: y(:), s(:)
       !
       !
       vel(:,index) = vel(:,index) - &
                      0.5D0 * ds(index) * grad(:,index) / mass(:)               
       !
       IF ( norm_grad(index) > eps32 ) THEN
          !
          force_versor = - grad(:,index) / norm_grad(index)
          !
          vel_component = ( vel(:,index) .dot. force_versor )
          !
          IF ( vel_component > 0.D0 ) THEN
             ! 
             vel(:,index) = vel_component * force_versor
             !
          ELSE
             !
#if defined (DEBUG_SMART_STEP) && defined (USE_SMART_STEP)
             PRINT '(/5X,"IMAGE = ",I2,"  resetting velocity"/)', index
#endif
             !
             vel(:,index) = 0.D0
             !
#if defined (USE_SMART_STEP)
             !
             vel_zeroed(index) = .TRUE.
             !
             ! ... an approximate newton-raphson step is performed
             !
             IF ( norm( pos_old(:,index) ) > 0.D0 .AND. &
                  norm( grad_old(:,index) ) > 0.D0 ) THEN
                !
                ALLOCATE( y( dim ), s( dim ) )
                !
                y = grad(:,index) - grad_old(:,index)
                s =  pos(:,index) -  pos_old(:,index)
                !
#  if defined (DEBUG_SMART_STEP)
                PRINT '(5X,"projection = ",F7.4)', &
                    ( s .dot. grad(:,index) ) / ( norm( s ) * norm_grad(index) )
#  endif
                !
                IF ( ABS( y .dot. s ) > eps8 ) THEN
                   !
                   grad(:,index) = 2.D0 * s * &
                                ABS( ( s .dot. grad(:,index) ) / ( y .dot. s ) )
                   !
                END IF
                !
#  if defined (DEBUG_SMART_STEP)
                PRINT '(5X,"step length:  ",F12.8)', norm( grad(:,index) )
#  endif
                !
                DEALLOCATE( y, s )
                !
             END IF
             !
#endif
             !
          END IF
          !
       ELSE
          !
          vel(:,index) = 0.D0
          !
       END IF
       !
       ! ... pos_old and grad_old are updated here
       !   
       pos_old(:,index)  = pos(:,index)
       grad_old(:,index) = grad(:,index)       
       !
       RETURN
       !
     END SUBROUTINE quick_min_second_step
     !
     ! ... routine for fixed temperature dynamics
     !
     !---------------------------------------------------------------------- 
     SUBROUTINE thermalization( N_in , N_fin )
       !----------------------------------------------------------------------
       !
       USE neb_variables, ONLY : vel, mass, temp, temp_req, deg_of_freedom
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: N_in, N_fin
       REAL (KIND=DP)      :: local_temp
       INTEGER             :: image
       !
       !
       temp = 0.D0  
       !
       DO image = N_in, N_fin
          !
          local_temp = ( ( mass(:) * vel(:,image) ) .dot. vel(:,image) ) / &
                       REAL( deg_of_freedom ) 
          !
          temp = temp + local_temp
          !
          vel(:,image) = vel(:,image) * SQRT( temp_req / local_temp )
          !
       END DO
       !
       temp = temp / REAL( N_fin - N_in )
       !
#if defined ( _DEBUG ) || ( _DEBUG_THERMALIZATION )
       !
       PRINT *, "temperature before thermalization", temp * au * eV_to_kelvin
       !
       temp = 0.D0  
       !
       DO image = N_in, N_fin
          !
          temp = temp + ( ( mass(:) * vel(:,image) ) .dot. vel(:,image) )
          !
       END DO
       !
       temp = temp / REAL( ( N_fin - N_in ) * deg_of_freedom )      
       !
       PRINT *, "temperature after thermalization", temp * au * eV_to_kelvin
       !
#endif
       !
       RETURN
       !
     END SUBROUTINE thermalization
     !
     ! ... this routine computes the optimal time step using an approximation
     ! ... of the second derivative along the last displacement vector
     !
     !------------------------------------------------------------------------
     FUNCTION optimal_time_step( index )
       !------------------------------------------------------------------------
       !
       USE constants,     ONLY : eps8
       USE neb_variables, ONLY : dim, pos_old, grad_old, istep_neb
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN)        :: index
       REAL(KIND=DP)              :: optimal_time_step
       REAL(KIND=DP), ALLOCATABLE :: y(:), s(:)
       REAL(KIND=DP)              :: projection
       INTEGER, PARAMETER         :: steps_at_fixed_ds = 10
         ! for the first "steps_at_fixed_ds" steps the update of the time
         ! step is disabled
       REAL(KIND=DP), PARAMETER   :: ds_min = 0.1D0, &
                                     ds_max = 4.0D0
         ! minimum and maximum allowed time-steps
       !
       !
       ALLOCATE( y( dim ), s( dim ) )
       !
       y = grad(:,index) - grad_old(:,index)
       s =  pos(:,index) -  pos_old(:,index)
       !
       IF ( istep_neb < steps_at_fixed_ds ) THEN
          !
          optimal_time_step = ds(index)
          !
       ELSE
          !
          IF ( ABS( y .dot. s ) > eps8 ) THEN
             !
             projection = ( s .dot. grad(:,index) ) / &
                          ( norm( s ) * norm_grad(index) )
             !
             PRINT '(5X,"projection = ",F7.4)', projection
             !
             ! ... to avoid numerical instabilities we take the 30 %
             ! ... of the optimal time step
             !
             optimal_time_step = 0.3D0 * ABS( projection ) * &
                                 SQRT( 2.D0 * ( s .dot. s ) / ABS( y .dot. s ) )
             !
          ELSE
             !
             optimal_time_step = ds(index)
             !
          END IF
          !
          optimal_time_step = MAX( ds_min, optimal_time_step )
          optimal_time_step = MIN( ds_max, optimal_time_step )
          !
          PRINT '(5X,"optimal_time_step = ",F12.8)', optimal_time_step
          !
       END IF
       !
       DEALLOCATE( y, s )
       !
     END FUNCTION optimal_time_step
     !
END MODULE minimization_routines
