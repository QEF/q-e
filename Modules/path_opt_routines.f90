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
MODULE path_opt_routines
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the optimization of the reaction path (NEB and SMD calculations)
  ! ... Written by Carlo Sbraccia ( 2003-2004 )  
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : eps32
  USE path_variables, ONLY : ds
  USE path_variables, ONLY : pos, grad, norm_grad, frozen
  USE path_variables, ONLY : ft_pos, ft_grad, norm_ft_grad, ft_frozen
  !
  USE basic_algebra_routines
  !  
  IMPLICIT NONE
  !
  CONTAINS
     !
     ! ... "real space" routines
     !
     !----------------------------------------------------------------------
     SUBROUTINE r_steepest_descent( index )
       !----------------------------------------------------------------------
       !
       ! ... this routine is also used for the langevin dynamics
       !
       USE path_variables, ONLY : llangevin, lang
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( frozen(index) ) RETURN
       !
       pos(:,index) = pos(:,index) - ds * grad(:,index)
       !
       IF ( llangevin ) &
          pos(:,index) = pos(:,index) + SQRT( ds ) * lang(:,index)
       !
       RETURN
       !
     END SUBROUTINE r_steepest_descent 
     !
     ! ... Molecular Dynamics based algorithms 
     ! ... velocity Verlet and quick min       
     !
     !----------------------------------------------------------------------
     SUBROUTINE r_velocity_Verlet_first_step( index )
       !---------------------------------------------------------------------- 
       !
       USE path_variables, ONLY : vel, vel_zeroed, pos_old, grad_old
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( frozen(index) ) RETURN
       !
       ! ... pos_old and grad_old are updated here
       !   
       pos_old(:,index)  = pos(:,index)
       grad_old(:,index) = grad(:,index)
       !
       IF ( vel_zeroed(index) ) THEN
          !
          vel(:,index) = vel(:,index) - 0.5D0 * grad(:,index)
          !
          pos(:,index) = pos(:,index) + vel(:,index)
          !
          vel_zeroed(index) = .FALSE.
          !
       ELSE
          !
          vel(:,index) = vel(:,index) - 0.5D0 * ds * grad(:,index)
          !
          pos(:,index) = pos(:,index) + ds * vel(:,index)
          !
       END IF
       !       
       RETURN
       !
     END SUBROUTINE r_velocity_Verlet_first_step
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE r_velocity_Verlet_second_step( index )
       !----------------------------------------------------------------------
       !
       USE path_variables, ONLY : vel, damp, ldamped_dyn, lmol_dyn
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( frozen(index) ) RETURN
       !
       IF ( ldamped_dyn ) THEN
          !
          vel(:,index) = damp * ( vel(:,index) - &
                                  0.5D0 * ds * grad(:,index) )
          !
       ELSE IF ( lmol_dyn ) THEN
          !       
          vel(:,index) = vel(:,index) - 0.5D0 * ds * grad(:,index)
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE r_velocity_Verlet_second_step
     !
     !----------------------------------------------------------------------
     SUBROUTINE r_quick_min_second_step( index )
       !----------------------------------------------------------------------
       !
       USE constants,       ONLY : eps8
       USE path_variables,  ONLY : pos_old, grad_old, vel, &
                                   dim, vel_zeroed
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN)         :: index
       REAL (KIND=DP)              :: force_versor(dim)
       REAL (KIND=DP)              :: vel_component
       REAL (KIND=DP), ALLOCATABLE :: y(:), s(:)
       !
       !
       IF ( frozen(index) ) RETURN
       !
       vel(:,index) = vel(:,index) - 0.5D0 * ds * grad(:,index)
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
                s = pos(:,index)  -  pos_old(:,index)
                !
#  if defined (DEBUG_SMART_STEP)
                PRINT '(5X,"projection = ",F7.4)', &
                    ( s .dot. grad(:,index) ) / &
                    ( norm( s ) * norm_grad(index) )
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
       RETURN
       !
     END SUBROUTINE r_quick_min_second_step
     !
     ! ... "reciprocal space" routines
     !
     ! ... steepest descent 
     !
     !----------------------------------------------------------------------
     SUBROUTINE ft_steepest_descent( mode )
       !----------------------------------------------------------------------
       !
       ! ... this routine is also used for the langevin dynamics
       !
       USE path_variables, ONLY : llangevin, ft_lang
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: mode
       !
       !
       IF ( ft_frozen(mode) ) RETURN
       !
       ft_pos(:,mode) = ft_pos(:,mode) - ds * ft_grad(:,mode)
       !
       IF ( llangevin ) &
          ft_pos(:,mode) = ft_pos(:,mode) + SQRT( ds ) * ft_lang(:,mode)
       !      
       RETURN
       !
     END SUBROUTINE ft_steepest_descent 
     !
     ! ... Molecular Dynamics based algorithms 
     ! ... velocity Verlet and quick min
     !
     !----------------------------------------------------------------------
     SUBROUTINE ft_velocity_Verlet_first_step( mode )
       !---------------------------------------------------------------------- 
       !
       USE path_variables, ONLY : ft_vel, ft_vel_zeroed, &
                                  ft_pos_old, ft_grad_old
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: mode
       !
       !
       IF ( ft_frozen(mode) ) RETURN
       !
       ! ... ft_pos_old  and  ft_grad_old  are updated here
       !   
       ft_pos_old(:,mode)  = ft_pos(:,mode)
       ft_grad_old(:,mode) = ft_grad(:,mode)
       !
       IF ( ft_vel_zeroed(mode) ) THEN
          !
          ft_vel(:,mode) = ft_vel(:,mode) - 0.5D0 * ft_grad(:,mode)
          !
          ft_pos(:,mode) = ft_pos(:,mode) + ft_vel(:,mode)
          !
          ft_vel_zeroed(mode) = .FALSE.
          !
       ELSE
          !
          ft_vel(:,mode) = ft_vel(:,mode) - &
                           0.5D0 * ds * ft_grad(:,mode)
          !
          ft_pos(:,mode) = ft_pos(:,mode) + ds * ft_vel(:,mode)
          !
       END IF
       !       
       RETURN
       !
     END SUBROUTINE ft_velocity_Verlet_first_step
     !
     !----------------------------------------------------------------------
     SUBROUTINE ft_velocity_Verlet_second_step( mode )
       !----------------------------------------------------------------------
       !
       USE path_variables, ONLY : ft_vel, damp, ldamped_dyn, lmol_dyn
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: mode
       !
       !
       IF ( ft_frozen(mode) ) RETURN
       !
       IF ( ldamped_dyn ) THEN
          !
          ft_vel(:,mode) = damp * ( ft_vel(:,mode) - &
                                    0.5D0 * ds * ft_grad(:,mode) )
          !
       ELSE IF ( lmol_dyn ) THEN
          !       
          ft_vel(:,mode) = ft_vel(:,mode) - 0.5D0 * ds * ft_grad(:,mode)
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE ft_velocity_Verlet_second_step
     !
     !----------------------------------------------------------------------
     SUBROUTINE ft_quick_min_second_step( mode )
       !----------------------------------------------------------------------
       !
       USE constants,      ONLY : eps8
       USE path_variables, ONLY : ft_pos_old, ft_grad_old, &
                                  ft_vel, dim, ft_vel_zeroed
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN)         :: mode
       REAL (KIND=DP)              :: ft_force_versor(dim)
       REAL (KIND=DP)              :: ft_vel_component
       REAL (KIND=DP), ALLOCATABLE :: y(:), s(:)
       !
       !
       IF ( ft_frozen(mode) ) RETURN
       !
       ft_vel(:,mode) = ft_vel(:,mode) - 0.5D0 * ds * ft_grad(:,mode)
       !
       IF ( norm_ft_grad(mode) > eps32 ) THEN
          !
          ft_force_versor = - ft_grad(:,mode) / norm_ft_grad(mode)
          !
          ft_vel_component = ( ft_vel(:,mode) .dot. ft_force_versor )
          !
          IF ( ft_vel_component > 0.D0 ) THEN
             ! 
             ft_vel(:,mode) = ft_vel_component * ft_force_versor
             !
          ELSE
             !
#if defined (DEBUG_SMART_STEP) && defined (USE_SMART_STEP)
             PRINT '(/5X,"MODE = ",I2,"  resetting velocity"/)', mode
#endif
             !
             ft_vel(:,mode) = 0.D0
             !
#if defined (USE_SMART_STEP)
             !
             ft_vel_zeroed(mode) = .TRUE.
             !
             ! ... an approximate newton-raphson step is performed
             !
             IF ( norm( ft_pos_old(:,mode) ) > 0.D0 .AND. &
                  norm( ft_grad_old(:,mode) ) > 0.D0 ) THEN
                !
                ALLOCATE( y( dim ), s( dim ) )
                !
                y = ft_grad(:,mode) - ft_grad_old(:,mode)
                s =  ft_pos(:,mode) -  ft_pos_old(:,mode)
                !
#  if defined (DEBUG_SMART_STEP)
                PRINT '(5X,"projection = ",F7.4)', &
                ( s .dot. ft_grad(:,mode) ) / ( norm( s ) * norm_ft_grad(mode) )
#  endif
                !
                IF ( ABS( y .dot. s ) > eps8 ) THEN
                   !
                   ft_grad(:,mode) = 2.D0 * s * &
                              ABS( ( s .dot. ft_grad(:,mode) ) / ( y .dot. s ) )
                   !
                END IF
                !
#  if defined (DEBUG_SMART_STEP)
                PRINT '(5X,"step length:  ",F12.8)', norm( ft_grad(:,mode) )
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
          ft_vel(:,mode) = 0.D0
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE ft_quick_min_second_step
     !
END MODULE path_opt_routines
