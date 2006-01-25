!
! Copyright (C) 2003-2005 Quantum-ESPRESSO group
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
  ! ... Written by Carlo Sbraccia ( 2003-2005 )
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
       vel(:,index) = vel(:,index) - 0.5D0 * ds * precond_grad(:,index)
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
       vel(:,index) = vel(:,index) - 0.5D0 * ds * precond_grad(:,index)
       !
       IF ( ldamped_dyn ) vel(:,index) = damp * vel(:,index)
       !
       RETURN
       !
     END SUBROUTINE velocity_Verlet_second_step
     !
     !----------------------------------------------------------------------
     SUBROUTINE quick_min_second_step( index )
       !----------------------------------------------------------------------
       !
       USE path_variables, ONLY : dim, vel
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       REAL(DP)            :: force_versor(dim)
       REAL(DP)            :: norm_pgrad
       !
       !
       vel(:,index) = vel(:,index) - 0.5D0 * ds * precond_grad(:,index)
       !
       norm_pgrad = norm( precond_grad(:,index) )
       !
       IF ( norm_pgrad > eps16 ) THEN
          !
          force_versor = - precond_grad(:,index) / norm_pgrad
          !
          vel(:,index) = force_versor * &
                         MAX( 0.D0, ( vel(:,index) .dot. force_versor ) )
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
     !-----------------------------------------------------------------------
     SUBROUTINE grad_precond( index )
       !-----------------------------------------------------------------------
       !
       ! ... this routine computes an estimate of H^-1 by using the BFGS
       ! ... algorithm and the preconditioned gradient  pg = H^-1 * g
       !
       USE path_variables, ONLY : dim, use_precond, frozen
       USE io_files,       ONLY : iunpath, iunbfgs, tmp_dir, prefix
       USE parser,         ONLY : int_to_char
       USE basic_algebra_routines
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: index
       !
       REAL(DP), ALLOCATABLE :: pos_p(:)
       REAL(DP), ALLOCATABLE :: grad_p(:)
       REAL(DP), ALLOCATABLE :: inv_hess(:,:)
       REAL(DP), ALLOCATABLE :: y(:), s(:)
       REAL(DP), ALLOCATABLE :: Hs(:), Hy(:), yH(:)
       REAL(DP)              :: sdoty, p_grad_norm
       CHARACTER(LEN=256)    :: bfgs_file
       LOGICAL               :: file_exists
       !
       REAL(DP), PARAMETER :: p_grad_norm_max = 0.6D0
       !
       !
       IF ( .NOT. use_precond ) THEN
          !
          precond_grad(:,index) = grad(:,index)
          !
          RETURN
          !
       END IF
       !
       IF ( frozen(index) ) RETURN
       !
       ALLOCATE( pos_p( dim ) )
       ALLOCATE( grad_p( dim ) )
       ALLOCATE( y( dim ), s( dim ) )
       ALLOCATE( inv_hess( dim, dim ) )
       ALLOCATE( Hs( dim ), Hy( dim ), yH( dim ) )       
       !
       bfgs_file = TRIM( tmp_dir ) // TRIM( prefix ) // "_" // &
                   TRIM( int_to_char( index ) )//"/"//TRIM( prefix )//'.bfgs'
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
          ! ... BFGS update
          !
          s(:) = pos(:,index)  - pos_p(:)
          y(:) = grad(:,index) - grad_p(:)
          !
          sdoty = ( s(:) .dot. y(:) )
          !
          IF ( sdoty > eps8 ) THEN
             !
             Hs(:) = ( inv_hess(:,:) .times. s(:) )
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
       precond_grad(:,index) = ( inv_hess(:,:) .times. grad(:,index) )
       !
       IF ( ( precond_grad(:,index) .dot. grad(:,index) ) < 0.D0 ) THEN
          !
          WRITE( UNIT = iunpath, FMT = '(5X,/,"image ",I3)' ) index
          WRITE( UNIT = iunpath, &
                 FMT = '(5X,"uphill step: resetting bfgs history",/)' )
          !
          precond_grad(:,index) = grad(:,index)
          !
          inv_hess(:,:) = identity( dim )
          !
       END IF
       !
       p_grad_norm = norm( precond_grad(:,index) )
       !
       precond_grad(:,index) = precond_grad(:,index) / &
                               p_grad_norm * MIN( p_grad_norm, p_grad_norm_max )
       !
       OPEN( UNIT = iunbfgs, &
             FILE = TRIM( bfgs_file ), STATUS = 'UNKNOWN', ACTION = 'WRITE' )
       !
       WRITE( iunbfgs, * ) pos(:,index)
       WRITE( iunbfgs, * ) grad(:,index)
       WRITE( iunbfgs, * ) inv_hess(:,:)
       !
       CLOSE( UNIT = iunbfgs )
       !
       DEALLOCATE( pos_p )
       DEALLOCATE( grad_p )
       DEALLOCATE( inv_hess )
       DEALLOCATE( y, s )
       DEALLOCATE( Hs, Hy, yH )
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
