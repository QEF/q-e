!
! Copyright (C) 2003-2006 Quantum ESPRESSO group
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
  USE mp_world,       ONLY : world_comm
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: langevin, steepest_descent, quick_min, broyden, broyden2
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
       CALL mp_bcast( pos, meta_ionode_id, world_comm )
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
       CALL mp_bcast( pos, meta_ionode_id, world_comm )
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
       CALL mp_bcast( pos,    meta_ionode_id, world_comm )
       CALL mp_bcast( posold, meta_ionode_id, world_comm )
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
       USE path_variables,  ONLY : lsmd
       USE path_io_units_module,  ONLY : broy_file, iunbroy, iunpath
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
       CALL mp_bcast( pos, meta_ionode_id, world_comm )
       !
       DEALLOCATE( t )
       DEALLOCATE( g )
       DEALLOCATE( s )
       !
       RETURN
       !
     END SUBROUTINE broyden
     !
     ! ... Broyden (rank one) optimisation - second attempt
     !
     !-----------------------------------------------------------------------
     SUBROUTINE broyden2()
       !-----------------------------------------------------------------------
#define DEBUG
       !
       USE path_variables,  ONLY : lsmd
       USE path_io_units_module, ONLY : broy_file, iunbroy, iunpath
       USE path_variables, ONLY : dim1, frozen, tangent, nim => num_of_images
       !
       IMPLICIT NONE
       !
       REAL(DP), PARAMETER   :: step_max = 0.6_DP
       INTEGER,  PARAMETER   :: broyden_ndim = 5
       !
       REAL(DP), ALLOCATABLE :: dx(:,:), df(:,:), x(:), f(:)
       REAL(DP), ALLOCATABLE :: x_last(:), f_last(:), mask(:)
       REAL(DP), ALLOCATABLE :: b(:,:), c(:), work(:)
       INTEGER, ALLOCATABLE  :: iwork(:)
       !
       REAL(DP)              :: x_norm, gamma0, J0, d2, d2_estimate
       LOGICAL               :: exists
       INTEGER               :: i, I_in, I_fin, info, j, niter
       !
       ! ... starting guess for the inverse Jacobian
       !
       J0 = ds*ds
       !
       ALLOCATE( dx( dim1*nim, broyden_ndim ), df( dim1*nim, broyden_ndim ) )
       ALLOCATE( x( dim1*nim ), f( dim1*nim ) )
       ALLOCATE( x_last( dim1*nim ), f_last( dim1*nim ), mask( dim1*nim ) )
       !
       ! define mask to skip frozen images 
       !
       mask(:) = 0.0_DP
       DO i = 1, nim
          I_in  = ( i - 1 )*dim1 + 1
          I_fin = i * dim1
          IF ( frozen(i) ) CYCLE
          mask(I_in:I_fin) = 1.0_DP
       END DO
       !
       ! copy current positions and gradients in local arrays
       !
       DO i = 1, nim
          I_in  = ( i - 1 )*dim1 + 1
          I_fin = i * dim1
          f(I_in:I_fin) =-grad(:,i)
          x(I_in:I_fin) = pos(:,i)
       END DO
       !
       ! only meta_ionode execute this part
       !
       IF ( meta_ionode ) THEN 
          d2 = DOT_PRODUCT( f(:),mask(:)*f(:) )
#ifdef DEBUG
          WRITE (*,*) " CURRENT ACTUAL D2 = ", d2
#endif
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
             READ( UNIT = iunbroy , FMT = * ) niter, d2_estimate
             READ( UNIT = iunbroy , FMT = * ) df(:,:), dx(:,:)
             READ( UNIT = iunbroy , FMT = * ) f_last(:), x_last(:)
             niter = min(broyden_ndim, niter + 1)
             !
             if (d2 > 2.0_DP * d2_estimate ) then
#ifdef DEBUG
                write (*,*) " bad D2 estimate ... reset history "
#endif
                niter = 1
                df(:,:) = 0.0_DP
                dx(:,:) = 0.0_DP
             end if
          ELSE 
             !
             df(:,:) = 0.0_DP
             dx(:,:) = 0.0_DP
             niter = 0
             !
          END IF
          CLOSE( UNIT = iunbroy )
          !
          ! ... Broyden's update
          !
          ! shift previous history, automatically discarding oldest iterations
          !
          DO i = broyden_ndim, 2, -1 
             df(:,i) = df(:,i-1)
             dx(:,i) = dx(:,i-1)
          END DO
          !
          ! and update it with last increment 
          !
          IF (niter > 0 ) THEN
            df(:,1) = f(:) - f_last(:)
            dx(:,1) = x(:) - x_last(:)
          END IF
          ! save for later use
          f_last(:) = f(:)
          x_last(:) = x(:)
          !
          x(:) = 0.0_DP
          IF ( niter > 0 ) THEN
             !
             ALLOCATE (b(niter,niter), c(niter), work(niter), iwork(niter))
             !
             ! create the matrix and the right-hand side of the liner system
             !
             b(:,:) = 0.0_DP
             c(:) = 0.0_DP
             DO i = 1,niter
                DO j = 1,niter
                   b(i,j) = DOT_PRODUCT(df(:,i),mask(:)*df(:,j))
                END DO
                c(i) = DOT_PRODUCT(f(:),mask(:)*df(:,i))
             END DO
             !
             ! solve the linear system
             !
             CALL DSYTRF( 'U', niter, b, niter, iwork, work, niter, info )
             CALL errore( 'broyden', 'factorization', abs(info) )
             CALL DSYTRI( 'U', niter, b, niter, iwork, work, info )
             CALL errore( 'broyden', 'DSYTRI', abs(info) ) 
             FORALL( i = 1:niter, j = 1:niter, j > i ) b(j,i) = b(i,j)
             !
             ! set the best correction vector and gradient
             !
             DO i = 1, niter
               gamma0 = DOT_PRODUCT( b(1:niter,i), c(1:niter) )
               call DAXPY(dim1*nim, -gamma0, dx(:,i),1, x,1)
               call DAXPY(dim1*nim, -gamma0, df(:,i),1, f,1)
             END DO
             !
             DEALLOCATE (b,c,work,iwork)
             !
          END IF
          d2 = DOT_PRODUCT( f(:), mask(:)*f(:) )
          x(:) =  mask(:) * ( x(:) + J0 * f(:) )
          x_norm = norm(x)
          x(:) = x_last(:) + x(:) * min ( 1.0_DP, step_max/x_norm)
#ifdef DEBUG
          WRITE (*,*) " ESTIMATED NEXT D2 = ", d2
          IF (x_norm > step_max)  &
             WRITE (*,*) " x_norm = ", x_norm, step_max
#endif
          !
          ! ... save the file containing the history
          !
          OPEN( UNIT = iunbroy, FILE = broy_file )
          !
          WRITE( UNIT = iunbroy, FMT = * ) nim
          WRITE( UNIT = iunbroy, FMT = * ) niter, d2
          WRITE( UNIT = iunbroy, FMT = * ) df(:,:), dx(:,:)
          WRITE( UNIT = iunbroy, FMT = * ) f_last(:), x_last(:)
          !
          CLOSE( UNIT = iunbroy )
          !
          ! ... copy broyden's step on the position array ...
          !
          pos(:,1:nim) =  RESHAPE( x, (/ dim1, nim /) )
          !
       END IF
       !
       ! ... and distribute it
       !
       CALL mp_bcast( pos, meta_ionode_id, world_comm )
       !
       DEALLOCATE( df, dx, f, x, f_last, x_last, mask )
       !
       RETURN
       !
     END SUBROUTINE broyden2
     !
END MODULE path_opt_routines
