!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
MODULE minimization_routines
  !---------------------------------------------------------------------------
  !
  USE parameters,     ONLY :  DP
  USE constants,      ONLY :  AU, eV_to_kelvin, eps32  
  USE neb_variables,  ONLY :  pos, ds, grad, norm_grad
  USE miscellany
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
       IF ( norm_grad(index) >= eps32 ) THEN
          !
          pos(:,index) = pos(:,index) - ds * grad(:,index)
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
       USE neb_variables,  ONLY : vel, mass
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       ! 
       vel(:,index) = vel(:,index) - ds / 2.D0 * grad(:,index) / mass(:)
       pos(:,index) = pos(:,index) + ds * vel(:,index)
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
       USE neb_variables,  ONLY : vel, mass, damp, ldamped_dyn, lmol_dyn
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( ldamped_dyn ) THEN
          !
          vel(:,index) = damp * ( vel(:,index) - &
                                  ds / 2.D0 * grad(:,index)  / mass(:) )
          !
       ELSE IF ( lmol_dyn ) THEN
          !       
          vel(:,index) = vel(:,index) - ds / 2.D0 * grad(:,index) / mass(:)
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
       USE neb_variables,  ONLY : vel, mass, dim
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN)            :: index
       REAL (KIND=DP), DIMENSION(dim) :: force_versor
       REAL (KIND=DP)                 :: vel_component
       !
       !
       vel(:,index) = vel(:,index) - ds / 2.D0 * grad(:,index) / mass(:)
       !
       IF ( norm_grad(index) >= eps32 ) THEN
          !
          force_versor = - grad(:,index) / norm_grad(index)
          !
          vel_component = DOT_PRODUCT( vel(:,index) , force_versor )
          !
          IF ( vel_component > 0.D0 ) THEN
             ! 
             vel(:,index) = vel_component * force_versor
             !
          ELSE
             !
             vel(:,index) = 0.D0
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
     ! ... routine for fixed temperature dynamics
     !
     !---------------------------------------------------------------------- 
     SUBROUTINE thermalization( N_in , N_fin )
       !----------------------------------------------------------------------
       !
       USE neb_variables,  ONLY : vel, mass, temp, temp_req, deg_of_freedom
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
          local_temp = DOT_PRODUCT( mass(:) * vel(:,image) , &
                       vel(:,image) ) / REAL( deg_of_freedom ) 
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
       PRINT *, "temperature before thermalization", temp * AU * eV_to_kelvin
       !
       temp = 0.D0  
       !
       DO image = N_in, N_fin
          !
          temp = temp + DOT_PRODUCT( mass(:) * vel(:,image) , vel(:,image) )
          !
       END DO
       !
       temp = temp / REAL( ( N_fin - N_in ) * deg_of_freedom )      
       !
       PRINT *, "temperature after thermalization", temp * AU * eV_to_kelvin
       !
#endif
       !
       RETURN
       !
     END SUBROUTINE thermalization
     !
END MODULE minimization_routines
