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
  USE path_variables, ONLY : ds
  USE path_variables, ONLY : pos, grad, precond_grad
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
       pos(:,index) = pos(:,index) - ds * grad(:,index)
       !
       IF ( llangevin ) pos(:,index) = pos(:,index) + lang(:,index)
       !
       RETURN
       !
     END SUBROUTINE steepest_descent 
     !
     !----------------------------------------------------------------------
     SUBROUTINE quick_min( index, istep )
       !---------------------------------------------------------------------- 
       !
       USE path_variables, ONLY : dim, posold
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index, istep
       !
       REAL(DP), ALLOCATABLE :: vel(:), force_versor(:), step(:)
       REAL(DP)              :: projection, norm_pgrad, norm_vel, norm_step
       !
       REAL(DP), PARAMETER :: max_step = 0.6D0  ! in bohr
       !
       !
       ALLOCATE( vel( dim ), force_versor( dim ), step( dim ) )
       !
       CALL grad_precond( index, istep )
       !
       vel(:) = pos(:,index) - posold(:,index)
       !
       norm_pgrad = norm( precond_grad(:,index) )
       !
       norm_vel = norm( vel(:) )
       !
       IF ( norm_pgrad > eps16 .AND. norm_vel > eps16 ) THEN
          !
          force_versor(:) = - precond_grad(:,index) / norm_pgrad
          !
          projection = ( vel(:) .dot. force_versor(:) )
          !
          vel(:) = MAX( 0.D0, projection ) * force_versor(:)
          !
       ELSE
          !
          vel(:) = 0.D0
          !
       END IF
       !
       posold(:,index) = pos(:,index)
       !
       step(:) = vel(:) - ds**2 * precond_grad(:,index)
       !
       norm_step = norm( step(:) )
       !
       step(:) = step(:) / norm_step
       !
       pos(:,index) = pos(:,index) + step(:) * MIN( norm_step, max_step )
       !
       DEALLOCATE( vel, force_versor, step )
       !
       RETURN
       !
     END SUBROUTINE quick_min
     !
     !-----------------------------------------------------------------------
     SUBROUTINE grad_precond( i, istep )
       !-----------------------------------------------------------------------
       !
       ! ... this routine computes an estimate of H^-1 by using the BFGS
       ! ... algorithm and the preconditioned gradient  pg = H^-1 * g
       !
       USE path_variables, ONLY : dim, use_precond
       USE io_files,       ONLY : iunpath, iunbfgs, tmp_dir, prefix
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: i, istep
       !
       REAL(DP), ALLOCATABLE      :: pos_p(:)
       REAL(DP), ALLOCATABLE      :: grad_p(:)
       REAL(DP), ALLOCATABLE      :: inv_hess(:,:)
       REAL(DP), ALLOCATABLE      :: y(:), s(:)
       REAL(DP), ALLOCATABLE      :: Hy(:), yH(:)
       REAL(DP)                   :: sdoty, pg_norm
       CHARACTER(LEN=256)         :: bfgs_file
       CHARACTER(LEN=6), EXTERNAL :: int_to_char
       LOGICAL                    :: file_exists
       !
       INTEGER,  PARAMETER :: nrefresh    = 25
       REAL(DP), PARAMETER :: max_pg_norm = 0.6D0
       !
       !
       IF ( .NOT. use_precond ) THEN
          !
          precond_grad(:,i) = grad(:,i)
          !
          RETURN
          !
       END IF
       !
       ALLOCATE( pos_p(  dim ) )
       ALLOCATE( grad_p( dim ) )
       ALLOCATE( y( dim ), s( dim ) )
       ALLOCATE( inv_hess( dim, dim ) )
       ALLOCATE( Hy( dim ), yH( dim ) )       
       !
       bfgs_file = TRIM( tmp_dir ) // TRIM( prefix ) // "_" // &
                   TRIM( int_to_char( i ) ) // "/" // TRIM( prefix ) // '.bfgs'
       !
       INQUIRE( FILE = TRIM( bfgs_file ), EXIST = file_exists )
       !
       IF ( file_exists ) THEN
          !
          OPEN( UNIT = iunbfgs, &
                FILE = TRIM( bfgs_file ), STATUS = 'OLD', ACTION = 'READ' )
          !
          READ( iunbfgs, * ) pos_p
          READ( iunbfgs, * ) grad_p
          READ( iunbfgs, * ) inv_hess
          !
          CLOSE( UNIT = iunbfgs )
          !
          ! ... the approximate inverse hessian is reset to one every nrefresh
          ! ... iterations: this is one to clean-up the memory of the starting
          ! ... configuration
          !
          IF ( MOD( nrefresh, istep ) == 0 ) inv_hess(:,:) = identity( dim )
          !
          ! ... BFGS update
          !
          s(:) = pos(:,i)  - pos_p(:)
          y(:) = grad(:,i) - grad_p(:)
          !
          sdoty = ( s(:) .dot. y(:) )
          !
          IF ( sdoty > eps8 ) THEN
             !
             Hy(:) = ( inv_hess(:,:) .times. y(:) )
             yH(:) = ( y(:) .times. inv_hess(:,:) )
             !
             inv_hess = inv_hess + 1.D0 / sdoty * &
                        ( ( 1.D0 + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                          ( matrix( s, yH ) +  matrix( Hy, s ) ) )
             !
          END IF
          !
       ELSE
          !
          inv_hess(:,:) = identity( dim )
          !
       END IF
       !
       precond_grad(:,i) = ( inv_hess(:,:) .times. grad(:,i) )
       !
       IF ( ( precond_grad(:,i) .dot. grad(:,i) ) < 0.D0 ) THEN
          !
          WRITE( UNIT = iunpath, FMT = '(5X,/,"image ",I3)' ) i
          WRITE( UNIT = iunpath, &
                 FMT = '(5X,"uphill step: resetting bfgs history",/)' )
          !
          precond_grad(:,i) = grad(:,i)
          !
          inv_hess(:,:) = identity( dim )
          !
       END IF
       !
       pg_norm = norm( precond_grad(:,i) )
       !
       precond_grad(:,i) = precond_grad(:,i) / pg_norm
       precond_grad(:,i) = precond_grad(:,i) * MIN( pg_norm, max_pg_norm )
       !
       OPEN( UNIT = iunbfgs, &
             FILE = TRIM( bfgs_file ), STATUS = 'UNKNOWN', ACTION = 'WRITE' )
       !
       WRITE( iunbfgs, * ) pos(:,i)
       WRITE( iunbfgs, * ) grad(:,i)
       WRITE( iunbfgs, * ) inv_hess(:,:)
       !
       CLOSE( UNIT = iunbfgs )
       !
       DEALLOCATE( pos_p )
       DEALLOCATE( grad_p )
       DEALLOCATE( inv_hess )
       DEALLOCATE( y, s )
       DEALLOCATE( Hy, yH )
       !
       RETURN
       !
     END SUBROUTINE grad_precond     
     !
     ! ... Broyden (rank one) optimisation
     !
     !-----------------------------------------------------------------------
     SUBROUTINE broyden()
       !-----------------------------------------------------------------------
       !
       USE constants,      ONLY : eps8
       USE control_flags,  ONLY : lsmd
       USE io_files,       ONLY : broy_file, iunbroy, iunpath
       USE path_variables, ONLY : dim, frozen, tangent, &
                                  n_im => num_of_images
       USE io_global,      ONLY : meta_ionode, meta_ionode_id
       USE mp,             ONLY : mp_bcast
       !
       IMPLICIT NONE
       !
       REAL(DP), ALLOCATABLE :: t(:), g(:), s(:,:)
       INTEGER               :: k, i, j, I_in, I_fin
       REAL(DP)              :: s_norm, coeff, norm_g
       LOGICAL               :: exists
       !
       REAL(DP), PARAMETER   :: step_max = 0.6D0
       INTEGER,  PARAMETER   :: broyden_ndim = 5
       !
       !
       ALLOCATE( g( dim * n_im ) )
       ALLOCATE( s( dim * n_im, broyden_ndim ) )
       ALLOCATE( t( dim * n_im ) )
       !
       g = 0.D0
       t = 0.D0
       !
       DO i = 1, n_im
          !
          I_in  = ( i - 1 ) * dim + 1
          I_fin = i * dim
          !
          IF ( lsmd ) t(I_in:I_fin) = tangent(:,i)
          !
          IF ( frozen(i) ) CYCLE
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
             ! ... if the number of images is changed the broyden history is
             ! ... reset and the algorithm starts from scratch
             !
             exists = ( i == n_im )
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
             s(:,:) = 0.D0
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
             ! ... the Broyden's subspace is swapped and s is projected
             ! ... orthogonally to the current tangent (this last thing 
             ! ... in the smd case only, otherwise t = 0.D0)
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
          s(:,k) = - ds * g(:)
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
             ELSE
                !
                s(:,k) = - ds * g(:)
                !
             END IF
             !
          END IF
          !
          IF ( ( s(:,k) .dot. g(:) ) > 0.D0 ) THEN
             !
             ! ... uphill step :  history reset
             !
             WRITE( UNIT = iunpath, &
                    FMT = '(5X,"BROYDEN uphill step :  history reset",/)' )
             !
             k = 1
             !
             s(:,:) = 0.D0
             s(:,k) = - ds * g(:)
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
          WRITE( UNIT = iunbroy, FMT = * ) n_im
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
       DEALLOCATE( t )
       DEALLOCATE( g )
       DEALLOCATE( s )
       !
       RETURN
       !
     END SUBROUTINE broyden
     !
END MODULE path_opt_routines
