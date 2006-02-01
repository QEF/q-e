!
! Copyright (C) 2003-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define __MIXEDBFGSMS
!
!----------------------------------------------------------------------------
MODULE bfgs_module
  !----------------------------------------------------------------------------
  !
  ! ... Ionic relaxation through the Broyden-Fletcher-Goldfarb-Shanno 
  ! ... minimization scheme and a "trust radius" line search based on the
  ! ... Wolfe conditions ( bfgs() subroutine ).
  !
  ! ... Written by Carlo Sbraccia ( 5/12/2003 )
  !
  ! ... references :  
  !
  ! ... 1) Roger Fletcher, Practical Methods of Optimization, John Wiley and 
  ! ...    Sons, Chichester, 2nd edn, 1987. 
  ! ... 2) Salomon R. Billeter, Alexander J. Turner, Walter Thiel, 
  ! ...    Phys. Chem. Chem. Phys. 2, 2177 (2000).
  ! ... 3) Salomon R. Billeter, Alessandro Curioni, Wanda Andreoni,
  ! ...    Comput. Mat. Science 27, 437, (2003).
  ! ... 4) Ren Weiqing, PhD Thesis: Numerical Methods for the Study of Energy
  ! ...    Landscapes and Rare Events. 
  !
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : iunbfgs, prefix
  USE constants, ONLY : eps16
  !
  USE basic_algebra_routines  
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: bfgs, terminate_bfgs
  !
  ! ... public variables
  !   
  PUBLIC :: bfgs_ndim,        &
            trust_radius_max, &
            trust_radius_min, &
            trust_radius_ini, &
            w_1,              &
            w_2
  !
  ! ... global variables
  !
  SAVE
  !
  REAL(DP), ALLOCATABLE :: &
      pos_p(:),               &! positions at the previous iteration
      grad_p(:),              &! gradients at the previous iteration
      inv_hess(:,:),          &! inverse of the hessian matrix ( updated
                               ! using the BFGS formula )
      step(:),                &! the last bfgs step
      step_old(:),            &! old bfgs steps
      pos_old(:,:),           &! list of m old positions
      grad_old(:,:),          &! list of m old gradients
      pos_best(:)              ! best extrapolated positions
  REAL(DP) :: &   
      trust_radius,           &! displacement along the bfgs direction
      trust_radius_old,       &! old displacement along the bfgs direction
      energy_p                 ! energy at the previous iteration
  INTEGER :: &
      scf_iter,               &! number of scf iterations
      bfgs_iter                ! number of bfgs iterations
  !
  ! ... default values for all these variables are set in 
  ! ... Modules/read_namelist.f90 (SUBROUTINE ions_defaults)
  !
  INTEGER :: &
      bfgs_ndim                ! dimension of the subspace for L-BFGS
                               ! fixed to 1 for standard BFGS algorithm
  REAL(DP)  :: &
      trust_radius_max,       &! maximum allowed displacement
      trust_radius_min,       &! minimum allowed displacement
      trust_radius_ini         ! initial displacement
  REAL(DP)  :: &
      w_1,                    &! parameters for Wolfe conditions
      w_2                      ! parameters for Wolfe conditions
  !
  ! ... Note that trust_radius_max, trust_radius_min, trust_radius_ini,
  ! ... w_1, w_2, bfgs_ndim have a default value, but can also be assigned 
  ! ... in the input.
  !
  CONTAINS
    !
    ! ... public methods :
    !
    !------------------------------------------------------------------------
    SUBROUTINE bfgs( pos, energy, grad, scratch, stdout, energy_thr, grad_thr, &
                     energy_error, grad_error, step_accepted, conv_bfgs )
      !------------------------------------------------------------------------
      !
      ! ... list of input/output arguments : 
      !  
      !  pos            : vector containing 3N coordinates of the system ( x )
      !  energy         : energy of the system ( V(x) )
      !  grad           : vector containing 3N components of ( grad( V(x) ) ) 
      !  scratch        : scratch directory
      !  stdout         : unit for standard output
      !  energy_thr     : treshold on energy difference for BFGS convergence
      !  grad_thr       : treshold on grad difference for BFGS convergence
      !                    the largest component of grad( V(x) ) is considered
      !  energy_error   : energy difference | V(x_i) - V(x_i-1) |
      !  grad_error     : the largest component of 
      !                    | grad(V(x_i)) - grad(V(x_i-1)) | 
      !  step_accepted  : .TRUE. if a new BFGS step is done
      !  conv_bfgs      : .TRUE. if BFGS convergence has been achieved
      !
      IMPLICIT NONE
      !
      ! ... input/output arguments
      !
      REAL(DP),         INTENT(INOUT) :: pos(:)
      REAL(DP),         INTENT(INOUT) :: energy       
      REAL(DP),         INTENT(INOUT) :: grad(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: stdout
      REAL(DP),         INTENT(IN)    :: energy_thr, grad_thr  
      REAL(DP),         INTENT(OUT)   :: energy_error, grad_error       
      LOGICAL,          INTENT(OUT)   :: step_accepted, conv_bfgs
      !
      ! ... local variables
      !
      INTEGER  :: dim, i, j
      LOGICAL  :: lwolfe
      REAL(DP) :: dE0s, den
      !
      REAL(DP), ALLOCATABLE :: res(:,:), overlap(:,:), work(:)
      INTEGER,  ALLOCATABLE :: iwork(:)
      INTEGER               :: k, k_m, info
      REAL(DP)              :: gamma0
      !
      !
      dim = SIZE( pos )
      !
      ! ... conditional work-space allocation
      !   
      IF ( .NOT. ALLOCATED( grad_old ) ) ALLOCATE( grad_old( dim, bfgs_ndim ) )
      IF ( .NOT. ALLOCATED( pos_old ) )  ALLOCATE( pos_old(  dim, bfgs_ndim ) )
      !
      IF ( .NOT. ALLOCATED( inv_hess ) ) ALLOCATE( inv_hess( dim, dim ) )
      !
      IF ( .NOT. ALLOCATED( pos_p ) )    ALLOCATE( pos_p(  dim ) )
      IF ( .NOT. ALLOCATED( grad_p ) )   ALLOCATE( grad_p( dim ) )
      IF ( .NOT. ALLOCATED( step ) )     ALLOCATE( step(      dim ) )
      IF ( .NOT. ALLOCATED( step_old ) ) ALLOCATE( step_old(  dim ) )
      IF ( .NOT. ALLOCATED( pos_best ) ) ALLOCATE( pos_best(  dim ) )
      !
      ! ... the BFGS file is read
      !
      CALL read_bfgs_file( pos, grad, energy, scratch, dim, stdout )
      !
      scf_iter = scf_iter + 1
      !
      IF ( scf_iter == 1 ) &
         WRITE( UNIT = stdout, FMT = '(/,5X,"BFGS Geometry Optimization")' )
      !
      energy_error = ABS( energy_p - energy )
      grad_error   = 0.D0
      !
      ! ... convergence is checked here
      !
      grad_error = MAXVAL( ABS( grad(:) ) )
      conv_bfgs  = energy_error < energy_thr      
      conv_bfgs  = conv_bfgs .AND. ( grad_error < grad_thr )
      !
      IF ( conv_bfgs ) RETURN
      !
      ! ... some output is written
      !
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"number of scf cycles",T30,"= ",I3)' ) scf_iter
      WRITE( UNIT = stdout, &
           & FMT = '(  5X,"number of bfgs steps",T30,"= ",I3,/)' ) bfgs_iter
      IF ( scf_iter > 1 ) &
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"energy old",T30,"= ",F18.10," ryd")' ) energy_p
      WRITE( UNIT = stdout, &
           & FMT = '(5X,"energy new",T30,"= ",F18.10," ryd",/)' ) energy
      !
      ! ... the bfgs algorithm starts here
      !
      IF ( ( energy > energy_p ) .AND. ( scf_iter > 1 ) ) THEN
         !
         ! ... the previous step is rejected, line search goes on
         !
         step_accepted = .FALSE.          
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"CASE: energy_new > energy_old",/)' )
         !
         ! ... the new trust radius is obtained with a quadratic interpolation
         !
         ! ... E(s) = a*s*s + b*s + c      ( we use E(0), dE(0), E(s') )
         !
         ! ... s_min = - 0.5 * ( dE(0)*s'*s' ) / ( E(s') - E(0) - dE(0)*s' )
         !
         dE0s = ( grad_p(:) .dot. step_old )
         !
         den = energy - energy_p - dE0s
         !              
         IF ( den > eps16 ) THEN
            !
            trust_radius = - 0.5D0 * dE0s * trust_radius_old / den
            !
         ELSE
            !
            ! ... no quadratic interpolation is possible
            !
            trust_radius = 0.5D0 * trust_radius_old
            !
         END IF
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius          
         !
         IF ( trust_radius < trust_radius_min ) THEN
            !
            ! ... the history is reset
            !                  
            WRITE( UNIT = stdout, FMT = '(/,5X,"resetting bfgs history",/)' )
            !
            IF ( trust_radius_old == trust_radius_min ) THEN
               !
               ! ... the history has already been reset at the previous step :
               ! ... something is going wrong
               !
               WRITE( UNIT = stdout, FMT = '(5X,"Notice: bfgs history", &
                      & " already reset at previous step",/)' )
               !
            END IF
            !
            inv_hess = identity( dim )
            !
            step = - grad
            !
            trust_radius = trust_radius_min
            !
         ELSE
            !
            ! ... values from the last succeseful bfgs step are restored
            !
            pos    = pos_p
            energy = energy_p
            grad   = grad_p
            !
            ! ... old bfgs direction ( normalized ) is recovered
            !
            step = step_old / trust_radius_old
            !
         END IF
         !
      ELSE    
         !
         ! ... a new bfgs step is done
         !
         bfgs_iter = bfgs_iter + 1
         !
         IF ( bfgs_iter > 1 ) THEN
            !
            step_accepted = .TRUE.
            !
            WRITE( UNIT = stdout, &
                 & FMT = '(5X,"CASE: energy_new < energy_old",/)' )
            !
            CALL check_wolfe_conditions( lwolfe, energy, grad )
            !
            IF ( lwolfe ) THEN
               !
               WRITE( UNIT = stdout, &
                    & FMT = '(5X,"Wolfe conditions satisfied",/)' )
               !
            ELSE
               !
               WRITE( UNIT = stdout, &
                    & FMT = '(5X,"Wolfe conditions not satisfied",/)' )
               !
            END IF
            !
            CALL update_inverse_hessian( pos, grad, dim, stdout )
            !
            IF ( bfgs_ndim > 1 ) THEN
               !
               ! ... GDIIS extrapolation
               !
               k   = MIN( bfgs_iter, bfgs_ndim )
               k_m = k + 1
               !
               ALLOCATE( res( dim, k ) )
               !
               ALLOCATE( overlap( k_m, k_m ) )
               !
               ALLOCATE( work(  k_m ) )
               ALLOCATE( iwork( k_m ) )
               !
               ! ... the new direction is added to the workspace
               !
               DO i = bfgs_ndim, 2, -1
                  !
                  pos_old(:,i)  = pos_old(:,i-1)
                  grad_old(:,i) = grad_old(:,i-1)
                  !
               END DO
               !
               pos_old(:,1)  = pos(:)
               grad_old(:,1) = grad(:)
               !
               ! ... |res_i> = H^-1 \times |g_i>
               !
               CALL DGEMM( 'N', 'N', dim, k, dim, 1.D0, inv_hess, &
                           dim, grad_old, dim, 0.D0, res, dim )
               !
               ! ... overlap_ij = <res_i|res_j>
               !
               CALL DGEMM( 'T', 'N', k, k, dim, 1.D0, res, &
                           dim, res, dim, 0.D0, overlap, k_m )
               !
               overlap( :,   k_m ) = 1.D0
               overlap( k_m, k_m ) = 0.D0
               !
               ! ... overlap is inverted
               !
               CALL DSYTRF( 'U', k_m, overlap, k_m, iwork, work, k_m, info )
               CALL DSYTRI( 'U', k_m, overlap, k_m, iwork, work, info )
               !
               FORALL( i = 1 : k_m, j = 1 : k_m, j > i )
                  !
                  ! ... overlap is symmetrised
                  ! 
                  overlap(j,i) = overlap(i,j)
                  !
               END FORALL
               !
               work = (/ ( 0.D0, i = 1, k ), 1.D0 /)
               !
               pos_best = 0.D0
               step     = 0.D0
               !
               DO i = 1, k
                  !
                  gamma0 = SUM( overlap(:,i) * work(:) )
                  !
                  pos_best = pos_best + gamma0 * pos_old(:,i)
                  !
                  step = step - gamma0 * res(:,i)
                  !
               END DO
               !
               ! ... the step must be consistent with the old positions
               !
               step = step + ( pos_best - pos )
               !
               IF ( ( grad .dot. step ) > 0.D0 ) THEN
                  !
                  ! ... if the extrapolated direction is uphill the last
                  ! ... gradient only is uded
                  !
                  step = - ( inv_hess .times. grad )
                  !
               END IF
               !
               DEALLOCATE( res, overlap, work, iwork )
               !
            ELSE
               !
               step = - ( inv_hess .times. grad )
               !
            END IF
            !
         ELSE
            !
            IF ( grad_error < 0.01D0 ) &
               trust_radius_ini = MIN( 0.2D0, trust_radius_ini )
            !
            step_accepted = .FALSE.
            !
            step = - grad
            !
         END IF
         !
         IF ( ( grad .dot. step ) > 0.D0 ) THEN
            !
            WRITE( UNIT = stdout, &
                   FMT = '(5X,"uphill step: resetting bfgs history",/)' )
            !
            step = - grad
            !
            inv_hess = identity( dim )
            !
         END IF
         !  
         ! ... the new trust radius is computed
         !
         IF ( bfgs_iter == 1 ) THEN
            !
            trust_radius =  trust_radius_ini
            !
         ELSE
            !
            trust_radius =  trust_radius_old
            !
            CALL compute_trust_radius( lwolfe, energy, grad, dim, stdout )
            !   
         END IF
         !
         IF ( conv_bfgs ) RETURN
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius
         !
      END IF
      !
      ! ... step along the bfgs direction
      !
      IF ( norm( step ) < eps16 ) THEN
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"WARNING : norm( step )",T30,"= ",F18.10)' ) &
              norm( step )
         !
         step = - grad
         !
      ELSE
         !
         step = trust_radius * step / norm( step )
         !
      END IF
      !
      ! ... informations needed for the next iteration are saved
      ! ... this must be done before positions update
      !
      CALL write_bfgs_file( pos, energy, grad, scratch )                
      !
      ! ... positions are updated
      !
      pos = pos + step
      !
    END SUBROUTINE bfgs
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_bfgs_file( pos, grad, energy, scratch, dim, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(IN)    :: pos(:)
      REAL(DP),         INTENT(IN)    :: grad(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: dim
      INTEGER,          INTENT(IN)    :: stdout
      REAL(DP),         INTENT(INOUT) :: energy
      !
      ! ... local variables
      !
      INTEGER             :: rank1, rank2
      CHARACTER (LEN=256) :: bfgs_file, hess_file
      LOGICAL             :: file_exists
      !
      !
      bfgs_file = TRIM( scratch ) // TRIM( prefix ) // '.bfgs'
      !
      INQUIRE( FILE = TRIM( bfgs_file ) , EXIST = file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... bfgs is restarted from file
         !
         OPEN( UNIT = iunbfgs, FILE = TRIM( bfgs_file ), &
               STATUS = 'UNKNOWN', ACTION = 'READ' )
         !
         READ( iunbfgs, * ) pos_p
         READ( iunbfgs, * ) grad_p
         READ( iunbfgs, * ) scf_iter
         READ( iunbfgs, * ) bfgs_iter
         READ( iunbfgs, * ) energy_p
         READ( iunbfgs, * ) pos_old
         READ( iunbfgs, * ) grad_old
         READ( iunbfgs, * ) inv_hess
         !     
         CLOSE( UNIT = iunbfgs )
         !
         trust_radius_old = norm( pos(:) - pos_p(:) )
         !
         step_old = ( pos(:) - pos_p(:) ) / trust_radius_old
         !
      ELSE
         !
         ! ... bfgs initialization
         !
         pos_p            = 0.D0
         grad_p           = 0.D0
         scf_iter         = 0
         bfgs_iter        = 0
         energy_p         = energy
         step_old         = 0.D0
         trust_radius_old = trust_radius_ini
         !
         pos_old(:,:)     = 0.D0
         grad_old(:,:)    = 0.D0
         !
         pos_old(:,1)     = pos
         grad_old(:,1)    = grad
         !
         inv_hess   = identity( dim )
         !
         hess_file = TRIM( scratch ) // TRIM( prefix ) // '.hess_in'
         !
         INQUIRE( FILE = TRIM( hess_file ) , EXIST = file_exists )
         !
         IF ( file_exists ) THEN
            !
            OPEN( UNIT = iunbfgs, FILE = TRIM( hess_file ), &
                  STATUS = 'UNKNOWN', ACTION = 'READ' )  
            !
            READ( iunbfgs, * ) rank1, rank2
            !
            IF ( ( rank1 == rank2 ) .AND. ( rank1 == dim ) ) THEN
               !
               WRITE( UNIT = stdout, &
                    & FMT = '(/,5X,"Reading the approximate inverse ", &
                    & "hessian from file",/)' )                
               !
               READ( iunbfgs, * ) inv_hess
               !
            END IF
            !
            CLOSE( UNIT = iunbfgs )                 
            !
         END IF
         !
      END IF
      !
    END SUBROUTINE read_bfgs_file
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_bfgs_file( pos, energy, grad, scratch )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(IN) :: pos(:)       
      REAL(DP),         INTENT(IN) :: energy       
      REAL(DP),         INTENT(IN) :: grad(:)       
      CHARACTER(LEN=*), INTENT(IN) :: scratch
      !
      !
      OPEN( UNIT = iunbfgs, FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs', &
            STATUS = 'UNKNOWN', ACTION = 'WRITE' )  
      !
      WRITE( iunbfgs, * ) pos
      WRITE( iunbfgs, * ) grad
      WRITE( iunbfgs, * ) scf_iter
      WRITE( iunbfgs, * ) bfgs_iter
      WRITE( iunbfgs, * ) energy
      WRITE( iunbfgs, * ) pos_old
      WRITE( iunbfgs, * ) grad_old
      WRITE( iunbfgs, * ) inv_hess
      ! 
      CLOSE( UNIT = iunbfgs )
      !
    END SUBROUTINE write_bfgs_file
    !
    !------------------------------------------------------------------------
    SUBROUTINE update_inverse_hessian( pos, grad, dim, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN)  :: pos(:)
      REAL(DP), INTENT(IN)  :: grad(:)   
      INTEGER,  INTENT(IN)  :: dim
      INTEGER,  INTENT(IN)  :: stdout
      !
      ! ... local variables
      !
      REAL(DP), ALLOCATABLE :: y(:), s(:)
      REAL(DP), ALLOCATABLE :: Hs(:), Hy(:), yH(:), aux(:)
      REAL(DP), ALLOCATABLE :: H_bfgs(:,:), H_ms(:,:)
      REAL(DP)              :: sdoty, coeff
      !
      !
      ALLOCATE( y( dim ), s( dim ) )
      ALLOCATE( Hs( dim ), Hy( dim ), yH( dim ), aux( dim ) )
      ALLOCATE( H_bfgs( dim, dim ), H_ms( dim, dim ) )
      !
      s(:) = pos(:)  - pos_p(:)
      y(:) = grad(:) - grad_p(:)
      !
      sdoty = ( s(:) .dot. y(:) )
      !
      IF ( ABS( sdoty ) < eps16 ) THEN
         !
         ! ... the history is reset
         !
         WRITE( stdout, '(/,5X,"WARINIG: unexpected behaviour in ", &
              & "update_inverse_hessian")' )
         WRITE( stdout, '(5X,"         resetting bfgs history",/)' )
         !
         inv_hess = identity( dim )
         !
         RETURN
         !
      END IF
      !
      Hs(:) = ( inv_hess .times. s(:) )
      Hy(:) = ( inv_hess .times. y(:) )
      yH(:) = ( y(:) .times. inv_hess )
      !
#if defined (__DFP)
      !
      ! ... DFP update
      !
      inv_hess = inv_hess + matrix( y(:), y(:) ) / sdoty - &
                            matrix( Hs(:), Hs(:) ) / ( s(:) .dot. Hs(:) )
      !
#endif
      !
#if defined (__BFGS)
      !
      ! ... BFGS update
      !
      inv_hess = inv_hess + 1.D0 / sdoty * &
                      ( ( 1.D0 + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                        ( matrix( s, yH ) +  matrix( Hy, s ) ) )
#endif
      !
#if defined (__MIXEDBFGSMS)
      !
      ! ... Bofill ( = BFGS + Murtag-Sargent ) update
      !
      aux(:) = y(:) - Hs(:)
      !
      coeff = ABS( s .dot. aux ) / ( norm( s ) * norm( aux ) )
      !
      H_bfgs = 1.D0 / sdoty * &
               ( ( 1.D0 + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                 ( matrix( s, yH ) +  matrix( Hy, s ) ) )
      !
      H_ms = matrix( aux, aux ) / ( s .dot. aux )
      !
      inv_hess = inv_hess + coeff * H_bfgs + ( 1.D0 - coeff ) * H_ms
      !
#endif
      !
      DEALLOCATE( y, s )
      DEALLOCATE( Hs, Hy, yH, aux )
      DEALLOCATE( H_bfgs, H_ms )
      !
      RETURN
      !
    END SUBROUTINE update_inverse_hessian
    !
    !------------------------------------------------------------------------
    SUBROUTINE check_wolfe_conditions( lwolfe, energy, grad )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN)  :: energy       
      REAL(DP), INTENT(IN)  :: grad(:)              
      LOGICAL,  INTENT(OUT) :: lwolfe
      !
      !
      lwolfe = ( energy - energy_p ) < w_1 * ( grad_p .dot. step_old )
      !
      lwolfe = lwolfe .AND. &
               ABS( grad .dot. step_old ) > - w_2 * ( grad_p .dot. step_old )
      !
    END SUBROUTINE check_wolfe_conditions
    !
    !------------------------------------------------------------------------
    SUBROUTINE compute_trust_radius( lwolfe, energy, grad, dim, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,  INTENT(IN)  :: lwolfe
      REAL(DP), INTENT(IN)  :: energy    
      REAL(DP), INTENT(IN)  :: grad(:)                         
      INTEGER,  INTENT(IN)  :: dim   
      INTEGER,  INTENT(IN)  :: stdout
      !
      ! ... local variables
      !
      REAL(DP) :: a
      LOGICAL  :: ltest
      !
      !
      ltest = ( energy - energy_p ) < w_1 * ( grad_p .dot. step_old )
      !
      ltest = ltest .AND. ( norm( step ) > trust_radius_old )
      !
      IF ( ltest ) THEN
         !
         a = 1.5D0
         !
      ELSE
         !
         a = 1.1D0
         !
      END IF
      !
      IF ( lwolfe ) THEN
         !
         trust_radius = MIN( trust_radius_max, 2.D0 * a * trust_radius_old )
         !
      ELSE
         !
         trust_radius = MIN( trust_radius_max, &
                             a * trust_radius_old, norm( step ) )
         !
      END IF
      !
      IF ( trust_radius < trust_radius_min ) THEN
         !
         ! ... the history is reset
         !
         WRITE( UNIT = stdout, FMT = '(5X,"resetting bfgs history",/)' )
         !
         inv_hess = identity( dim )
         !
         step = - grad
         !
         trust_radius = trust_radius_min
         !
      END IF
      !
    END SUBROUTINE compute_trust_radius
    !
    !------------------------------------------------------------------------
    SUBROUTINE terminate_bfgs( energy, stdout, scratch )
      !------------------------------------------------------------------------
      !
      USE io_files, ONLY : prefix, delete_if_present
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(IN) :: energy  
      INTEGER,          INTENT(IN) :: stdout         
      CHARACTER(LEN=*), INTENT(IN) :: scratch       
      !       
      !
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"bfgs converged in ",I3," scf cycles and ", &
           &         I3," bfgs steps")' ) scf_iter, bfgs_iter
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"Final energy = ",F18.10," ryd")' ) energy
      !
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"Saving the approximate inverse hessian",/)' )
      !
      OPEN( UNIT = iunbfgs, FILE = TRIM( scratch ) // TRIM( prefix ) // &
          & '.hess_out', STATUS = 'UNKNOWN', ACTION = 'WRITE' )  
      !
      WRITE( iunbfgs, * ) SHAPE( inv_hess )
      WRITE( iunbfgs, * ) inv_hess
      !
      CLOSE( UNIT = iunbfgs )
      !
      DEALLOCATE( pos_p )
      DEALLOCATE( grad_p )
      DEALLOCATE( pos_old )
      DEALLOCATE( grad_old )
      DEALLOCATE( inv_hess )
      DEALLOCATE( step )
      DEALLOCATE( step_old )
      DEALLOCATE( pos_best )
      !
      CALL delete_if_present( TRIM( scratch ) // TRIM( prefix ) // '.bfgs' )
      !
    END SUBROUTINE terminate_bfgs
    !
END MODULE bfgs_module
