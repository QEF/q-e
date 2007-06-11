!
! Copyright (C) 2003-2006 Quantum-ESPRESSO group
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
  ! ... the optimisation of the reaction path (NEB and SMD calculations)
  !
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : eps8, eps16
  USE path_variables, ONLY : ds, pos, grad
  USE io_global,      ONLY : meta_ionode, meta_ionode_id
  USE mp,             ONLY : mp_bcast
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: langevin, steepest_descent, quick_min, broyden
  !
  CONTAINS
     !
     !----------------------------------------------------------------------
     SUBROUTINE langevin( idx )
       !----------------------------------------------------------------------
       !
       USE path_variables, ONLY : lang
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: idx
       !
       IF ( meta_ionode ) THEN
          !
          pos(:,idx) = pos(:,idx) - ds*grad(:,idx) + lang(:,idx)
          !
       END IF
       !
       CALL mp_bcast( pos, meta_ionode_id )
       !
       RETURN
       !
     END SUBROUTINE langevin
     !
     !----------------------------------------------------------------------
     SUBROUTINE steepest_descent( idx )
       !----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: idx
       !
       IF ( meta_ionode ) THEN
          !
          pos(:,idx) = pos(:,idx) - ds*ds*grad(:,idx)
          !
       END IF
       !
       CALL mp_bcast( pos, meta_ionode_id )
       !
       RETURN
       !
     END SUBROUTINE steepest_descent
     !
     !----------------------------------------------------------------------
     SUBROUTINE quick_min( idx, istep )
       !---------------------------------------------------------------------- 
       !
       ! ... projected Verlet algorithm
       !
       USE path_variables, ONLY : dim1, posold
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: idx, istep
       !
       REAL(DP), ALLOCATABLE :: vel(:), force_versor(:), step(:)
       REAL(DP)              :: projection, norm_grad, norm_vel, norm_step
       !
       REAL(DP), PARAMETER :: max_step = 0.6_DP  ! in bohr
       !
       !
       IF ( meta_ionode ) THEN
          !
          ALLOCATE( vel( dim1 ), force_versor( dim1 ), step( dim1 ) )
          !
          vel(:) = pos(:,idx) - posold(:,idx)
          !
          norm_grad = norm( grad(:,idx) )
          !
          norm_vel = norm( vel(:) )
          !
          IF ( norm_grad > eps16 .AND. norm_vel > eps16 ) THEN
             !
             force_versor(:) = - grad(:,idx) / norm_grad
             !
             projection = ( vel(:) .dot. force_versor(:) )
             !
             vel(:) = MAX( 0.0_DP, projection ) * force_versor(:)
             !
          ELSE
             !
             vel(:) = 0.0_DP
             !
          END IF
          !
          posold(:,idx) = pos(:,idx)
          !
          step(:) = vel(:) - ds*ds*grad(:,idx)
          !
          norm_step = norm( step(:) )
          !
          step(:) = step(:) / norm_step
          !
          pos(:,idx) = pos(:,idx) + step(:) * MIN( norm_step, max_step )
          !
          DEALLOCATE( vel, force_versor, step )
          !
       END IF
       !
       CALL mp_bcast( pos,    meta_ionode_id )
       CALL mp_bcast( posold, meta_ionode_id )
       !
       RETURN
       !
     END SUBROUTINE quick_min
     !
     ! ... Broyden (rank one) optimisation
     !
     !-----------------------------------------------------------------------
     SUBROUTINE broyden()
       !-----------------------------------------------------------------------
       !
       USE control_flags,  ONLY : lsmd
       USE io_files,       ONLY : broy_file, iunbroy, iunpath
       USE path_variables, ONLY : dim1, frozen, tangent, nim => num_of_images
       !
       IMPLICIT NONE
       !
       REAL(DP), ALLOCATABLE :: t(:), g(:), s(:,:)
       INTEGER               :: k, i, j, I_in, I_fin
       REAL(DP)              :: s_norm, coeff, norm_g
       REAL(DP)              :: J0
       LOGICAL               :: exists
       !
       REAL(DP), PARAMETER   :: step_max = 0.6_DP
       INTEGER,  PARAMETER   :: broyden_ndim = 5
       !
       !
       ! ... starting guess for the inverse Jacobian
       !
       J0 = ds*ds
       !
       ALLOCATE( g( dim1*nim ) )
       ALLOCATE( s( dim1*nim, broyden_ndim ) )
       ALLOCATE( t( dim1*nim ) )
       !
       g(:) = 0.0_DP
       t(:) = 0.0_DP
       !
       DO i = 1, nim
          !
          I_in  = ( i - 1 )*dim1 + 1
          I_fin = i * dim1
          !
          IF ( frozen(i) ) CYCLE
          !
          IF ( lsmd ) t(I_in:I_fin) = tangent(:,i)
          !
          g(I_in:I_fin) = grad(:,i)
          !
       END DO
       !
       norm_g = MAXVAL( ABS( g ) )
       !
       IF ( norm_g == 0.0_DP ) RETURN
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
             ! ... if the number of images is changed the broyden history is
             ! ... reset and the algorithm starts from scratch
             !
             exists = ( i == nim )
             !
          END IF
          !
          IF ( exists ) THEN
             !
             READ( UNIT = iunbroy , FMT = * ) k
             READ( UNIT = iunbroy , FMT = * ) s(:,:)
             !
             k = k + 1
             !
          ELSE 
             !
             s(:,:) = 0.0_DP
             !
             k = 1
             !
          END IF
          !
          CLOSE( UNIT = iunbroy )
          !
          ! ... Broyden's update
          !
          IF ( k > broyden_ndim ) THEN
             !
             ! ... the Broyden's subspace is swapped and the projection of 
             ! ... s along the current tangent is removed (this last thing 
             ! ... in the smd case only, otherwise t = 0.0_DP)
             !
             k = broyden_ndim
             !
             DO j = 1, k - 1
                !
                s(:,j) = s(:,j+1) - t(:) * ( s(:,j+1) .dot. t(:) )
                !
             END DO
             !
          END IF
          !
          s(:,k) = - J0 * g(:)
          !
          IF ( k > 1 ) THEN
             !
             DO j = 1, k - 2
                !
                s(:,k) = s(:,k) + ( s(:,j) .dot. s(:,k) ) / &
                                  ( s(:,j) .dot. s(:,j) ) * s(:,j+1)
                !
             END DO
             !
             coeff = ( s(:,k-1) .dot. ( s(:,k-1) - s(:,k) ) )
             !
             IF ( coeff > eps8 ) THEN
                !
                s(:,k) = ( s(:,k-1) .dot. s(:,k-1) ) / coeff * s(:,k)
                !
             END IF
             !
          END IF
          !
          IF ( ( s(:,k) .dot. g(:) ) > 0.0_DP ) THEN
             !
             ! ... uphill step :  history reset
             !
             WRITE( UNIT = iunpath, &
                    FMT = '(/,5X,"broyden uphill step : history is reset",/)' )
             !
             k = 1
             !
             s(:,:) = 0.0_DP
             s(:,k) = - J0 * g(:)
             !
          END IF
          !
          s_norm = norm( s(:,k) )
          !
          s(:,k) = s(:,k) / s_norm * MIN( s_norm, step_max )
          !
          ! ... save the file containing the history
          !
          OPEN( UNIT = iunbroy, FILE = broy_file )
          !
          WRITE( UNIT = iunbroy, FMT = * ) nim
          WRITE( UNIT = iunbroy, FMT = * ) k
          WRITE( UNIT = iunbroy, FMT = * ) s
          !
          CLOSE( UNIT = iunbroy )
          !
          ! ... broyden's step
          !
          pos(:,1:nim) = pos(:,1:nim) + RESHAPE( s(:,k), (/ dim1, nim /) )
          !
       END IF
       !
       CALL mp_bcast( pos, meta_ionode_id )
       !
       DEALLOCATE( t )
       DEALLOCATE( g )
       DEALLOCATE( s )
       !
       RETURN
       !
     END SUBROUTINE broyden
     !
END MODULE path_opt_routines
