!
! Copyright (C) 2003-2005 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define USE_SMART_STEP
!#define DEBUG_SMART_STEP
!
!--------------------------------------------------------------------------
MODULE path_opt_routines
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the optimization of the reaction path (NEB and SMD calculations)
  !
  ! ... Written by Carlo Sbraccia ( 2003-2005 )
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : eps32
  USE path_variables, ONLY : ds
  USE path_variables, ONLY : pos, grad, norm_grad, frozen
  !
  USE basic_algebra_routines
  !  
  IMPLICIT NONE
  !
  CONTAINS
     !
     !----------------------------------------------------------------------
     SUBROUTINE steepest_descent( index )
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
     END SUBROUTINE steepest_descent 
     !
     ! ... Molecular Dynamics based algorithms 
     ! ... velocity Verlet and quick min       
     !
     !----------------------------------------------------------------------
     SUBROUTINE velocity_Verlet_first_step( index )
       !---------------------------------------------------------------------- 
       !
       USE path_variables, ONLY : vel, vel_zeroed, pos_old, grad_old
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( frozen(index) ) THEN
          !
          vel_zeroed(index) = .FALSE.
          !
          RETURN
          !
       END IF
       !
       ! ... pos_old and grad_old are updated here
       !   
       pos_old(:,index)  = pos(:,index)
       grad_old(:,index) = grad(:,index)
       !
       IF ( vel_zeroed(index) ) THEN
          !
          pos(:,index) = pos(:,index) - 0.5D0 * grad(:,index)
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
     END SUBROUTINE velocity_Verlet_first_step
     !
     !----------------------------------------------------------------------
     SUBROUTINE velocity_Verlet_second_step( index )
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
          vel(:,index) = damp * ( vel(:,index) - 0.5D0 * ds * grad(:,index) )
          !
       ELSE IF ( lmol_dyn ) THEN
          !       
          vel(:,index) = vel(:,index) - 0.5D0 * ds * grad(:,index)
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
       USE constants,       ONLY : eps8
       USE path_variables,  ONLY : pos_old, grad_old, vel, dim, vel_zeroed
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
             vel(:,index) = 0.D0
             !
#if defined (USE_SMART_STEP)
             !
#  if defined (DEBUG_SMART_STEP)
             PRINT '(/5X,"IMAGE = ",I2,"  resetting velocity"/)', index
#  endif             
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
                    ( s .dot. grad(:,index) ) / ( norm( s ) * norm_grad(index) )
#  endif
                !
                IF ( ( y .dot. s ) > eps8 ) THEN
                   !
                   grad(:,index) = 2.D0 * s * &
                                   ( s .dot. grad(:,index) ) / ( y .dot. s )
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
     END SUBROUTINE quick_min_second_step
     !
     ! ... Broyden (rank one) optimisation
     !
     !-----------------------------------------------------------------------
     SUBROUTINE broyden()
       !-----------------------------------------------------------------------
       !
       USE io_files,       ONLY : broy_file, iunbroy
       USE path_variables, ONLY : dim, num_of_images, reset_broyden
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), ALLOCATABLE :: g(:), s(:,:)
       INTEGER                     :: j
       INTEGER                     :: k
       REAL (KIND=DP)              :: s_norm
       LOGICAL                     :: exists
       REAL (KIND=DP), PARAMETER   :: J0           = 4.0D0
       REAL (KIND=DP), PARAMETER   :: step_max     = 1.0D0
       INTEGER,        PARAMETER   :: broyden_ndim = 32
       !
       !
       ALLOCATE( g( dim * num_of_images ) )
       ALLOCATE( s( dim * num_of_images, broyden_ndim ) )
       !
       g(:) = RESHAPE( SOURCE = grad, SHAPE = (/ dim * num_of_images /) )
       !
       IF ( norm( g ) == 0.D0 ) RETURN
       !
       ! ... open the file containing the old configurations of the path
       !
       INQUIRE( FILE = broy_file, EXIST = exists )
       !
       IF ( reset_broyden ) exists = .FALSE.
       !
       IF ( exists ) THEN
          !
          OPEN( UNIT = iunbroy, FILE = broy_file, STATUS = "UNKNOWN" )
          !
          READ( UNIT = iunbroy , FMT = * ) k
          READ( UNIT = iunbroy , FMT = * ) s
          !
          CLOSE( UNIT = iunbroy )
          !
          k = MIN( k + 1, broyden_ndim )
          !
       ELSE
          !
          s = 0.D0
          !
          k = 1
          !
          reset_broyden = .FALSE.
          !
       END IF
       !
       ! ... Broyden update
       !
       s(:,k) = - J0 * g(:)
       !
       DO j = 1, k - 2
          !
          s(:,k) = s(:,k) + ( s(:,j) .dot. s(:,k) ) / &
                            ( s(:,j) .dot. s(:,j) ) * s(:,j+1)
          !
       END DO
       !
       IF ( k > 1 ) THEN
          !
          s(:,k) = ( s(:,k-1) .dot. s(:,k-1) ) / &
                   ( s(:,k-1) .dot. ( s(:,k-1) - s(:,k) ) ) * s(:,k)
          !
       END IF
       !
       IF ( ( s(:,k) .dot. g(:) ) > 0.D0 ) THEN
          !
          ! ... uphill step :  reset history
          !
          k = 1
          !
          s = 0.D0
          !
          s(:,k) = - J0 * g(:)
          !
       END IF
       !
       s_norm = norm( s(:,k) )
       !
       s(:,k) = s(:,k) / s_norm * MIN( s_norm, step_max )
       !
       pos = pos + RESHAPE( SOURCE = s(:,k), &
                            SHAPE = (/ dim, num_of_images /) )
       !
       IF ( k == broyden_ndim ) THEN
          !
          DO j = k, 2, -1
             !
             s(:,j-1) = s(:,j)
             !
          END DO
          !
       END IF
       !
       ! ... save the file containing the history
       !
       OPEN( UNIT = iunbroy, FILE = broy_file )
       !
       WRITE( UNIT = iunbroy, FMT = * ) k
       WRITE( UNIT = iunbroy, FMT = * ) s
       !
       CLOSE( UNIT = iunbroy )
       !
       DEALLOCATE( g )
       DEALLOCATE( s )
       !
       RETURN
       !
     END SUBROUTINE broyden
     !
END MODULE path_opt_routines
