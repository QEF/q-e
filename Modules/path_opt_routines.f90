!
! Copyright (C) 2003-2005 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
  USE path_variables, ONLY : pos, grad, frozen
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
       USE path_variables, ONLY : vel
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       !
       IF ( frozen(index) ) RETURN
       !
       vel(:,index) = vel(:,index) - 0.5D0 * ds * grad(:,index)
       !
       pos(:,index) = pos(:,index) + ds * vel(:,index)
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
       USE path_variables, ONLY : dim, vel, norm_grad
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       REAL (KIND=DP)      :: force_versor(dim)
       REAL (KIND=DP)      :: vel_component
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
       USE constants,      ONLY : eps8, au, bohr_radius_angs
       USE io_files,       ONLY : broy_file, iunbroy
       USE path_variables, ONLY : dim, reset_broyden, frozen, &
                                  n_im => num_of_images
       USE io_global,      ONLY : meta_ionode, meta_ionode_id
       USE mp,             ONLY : mp_bcast
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), ALLOCATABLE :: g(:), s(:,:), pos_old(:,:)
       INTEGER                     :: k, i, j, I_in, I_fin
       REAL (KIND=DP)              :: s_norm, coeff, &
                                      norm_g, norm_g_old, g_dot_s
       LOGICAL                     :: exists, accepted
       REAL (KIND=DP)              :: J0
       !
       REAL (KIND=DP), PARAMETER   :: step_max = 0.6D0
       !
       INTEGER,        PARAMETER   :: broyden_ndim = 5
       !
       !
       ALLOCATE( g( dim * n_im ) )
       ALLOCATE( s( dim * n_im, broyden_ndim ) )
       !
       ALLOCATE( pos_old( dim, n_im ) )
       !
       g = 0.D0
       !
       DO i = 1, n_im
          !
          IF ( frozen(i) ) CYCLE
          !
          I_in  = ( i - 1 ) * dim + 1
          I_fin = i * dim
          !
          g(I_in:I_fin) = grad(:,i)
          !
       END DO
       !
       norm_g = MAXVAL( ABS( g ) )
       !
       IF ( norm_g == 0.D0 ) RETURN
       !
       IF ( meta_ionode ) THEN
          !
          ! ... open the file containing the broyden's history
          !
          INQUIRE( FILE = broy_file, EXIST = exists )
          !
          IF ( exists ) THEN
             !
             OPEN( UNIT = iunbroy, FILE = broy_file, STATUS = "OLD" )
             !
             READ( UNIT = iunbroy , FMT = * ) i
             !
             reset_broyden = ( i /= n_im )
             !
          END IF
          !
          IF ( exists .AND. .NOT. reset_broyden ) THEN
             !
             READ( UNIT = iunbroy , FMT = * ) J0
             READ( UNIT = iunbroy , FMT = * ) accepted
             READ( UNIT = iunbroy , FMT = * ) norm_g_old
             READ( UNIT = iunbroy , FMT = * ) pos_old
             READ( UNIT = iunbroy , FMT = * ) k
             READ( UNIT = iunbroy , FMT = * ) s
             !
             IF ( accepted ) THEN
                !
                g_dot_s = ( g(:) .dot. s(:,k) ) / norm( s(:,k) )
                !
                ! ... here we check wether the previous step has to be 
                ! ... accepted or not
                !
                IF ( g_dot_s < 0.D0 ) THEN
                   !
                   PRINT '(5X,"case 1")'
                   !
                   accepted = .TRUE.
                   !
                ELSE
                   !
                   IF ( norm_g < norm_g_old ) THEN
                      !
                      PRINT '(5X,"case 2")'
                      !
                      accepted = .TRUE.
                      !
                      J0 = J0 * 1.1D0
                      !
                   ELSE
                      !
                      PRINT '(5X,"case 3")'
                      !
                      J0 = J0 * 0.8D0
                      !
                      accepted = .FALSE.
                      !
                      pos(:,1:n_im) = pos_old
                      !
                      norm_g = norm_g_old
                      !
                   END IF
                   !
                END IF
                !
             ELSE
                !
                accepted = .TRUE.
                !
             END IF
             !
             IF ( accepted ) k = k + 1
             !
          ELSE 
             !
             IF ( reset_broyden ) THEN
                !
                READ( UNIT = iunbroy , FMT = * ) J0
                !
             ELSE
                !
                J0 = ds
                !
             END IF
             !             
             accepted = .TRUE.
             !
             s = 0.D0
             !
             k = 1
             !
             reset_broyden = .FALSE.
             !
          END IF
          !
          CLOSE( UNIT = iunbroy )
          !
          PRINT '(5X,"J0        = ",F10.6)', J0
          PRINT '(5X,"norm( g ) = ",F10.6)', norm_g / bohr_radius_angs * au
          PRINT '(5X,"accepted  = ",L1,/)', accepted
          !
          IF ( accepted ) THEN
             !
             ! ... Broyden's update
             !
             IF ( k > broyden_ndim ) THEN
                !
                ! ... the Broyden's subspace is swapped
                !
                k = broyden_ndim
                !
                DO j = 1, k - 1
                   !
                   s(:,j) = s(:,j+1)
                   !
                END DO
                !
             END IF
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
                coeff = ( s(:,k-1) .dot. ( s(:,k-1) - s(:,k) ) )
                !
                IF ( coeff > eps8 ) & 
                   s(:,k) = ( s(:,k-1) .dot. s(:,k-1) ) / coeff * s(:,k)
                !
             END IF
             !
             IF ( ( s(:,k) .dot. g(:) ) > 0.D0 ) THEN
                !
                ! ... uphill step :  history reset
                !
                PRINT '(5X,"BROYDEN uphill step :  history reset",/)'
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
          ELSE
             !
             ! ... here we try a shorter step
             !
             s(:,k) = 0.5D0 * s(:,k)
             !
          END IF
          !
          ! ... save the file containing the history
          !
          OPEN( UNIT = iunbroy, FILE = broy_file )
          !
          WRITE( UNIT = iunbroy, FMT = * ) n_im
          WRITE( UNIT = iunbroy, FMT = * ) J0
          WRITE( UNIT = iunbroy, FMT = * ) accepted
          WRITE( UNIT = iunbroy, FMT = * ) norm_g
          WRITE( UNIT = iunbroy, FMT = * ) pos(:,1:n_im)
          WRITE( UNIT = iunbroy, FMT = * ) k
          WRITE( UNIT = iunbroy, FMT = * ) s
          !
          CLOSE( UNIT = iunbroy )
          !
          ! ... broyden's step
          !
          pos(:,1:n_im) = pos(:,1:n_im) + RESHAPE( s(:,k), (/ dim, n_im /) )
          !
       END IF
       !
       CALL mp_bcast( pos, meta_ionode_id )  
       !
       DEALLOCATE( g )
       DEALLOCATE( s )
       DEALLOCATE( pos_old )
       !
       RETURN
       !
     END SUBROUTINE broyden
     !
END MODULE path_opt_routines
