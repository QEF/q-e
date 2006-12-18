!
! Copyright (C) 2003-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE bfgs_module
  !----------------------------------------------------------------------------
  !
  ! ... Ionic relaxation through the Broyden-Fletcher-Goldfarb-Shanno 
  ! ... minimization scheme and a "trust radius" line search based on the
  ! ... Wolfe conditions.
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
      bfgs_iter,              &! number of bfgs iterations
      gdiis_iter               ! number of gdiis iterations
  !
  LOGICAL :: &
      tr_min_hit               ! .TRUE. if the trust_radius has already been set
                               ! to the minimum value at the previous step
  !
  LOGICAL :: &
      conv_bfgs                ! .TRUE. when bfgs convergence has been achieved
  !
  REAL(DP), ALLOCATABLE :: &
      inv_hess0(:,:),         &! starting guess for the iverse hessian matrix
      projmat(:,:),           &! matrix that pojects out fixed ions
      dirtm(:,:),             &! transformation matrix to internal coordinates
      invtm(:,:)               ! inverse transformation matrix
  !
  ! ... default values for all these variables are set in 
  ! ... Modules/read_namelist.f90 (SUBROUTINE ions_defaults)
  !
  INTEGER :: &
      bfgs_ndim                ! dimension of the subspace for GDIIS
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
    !------------------------------------------------------------------------
    SUBROUTINE bfgs( pos, energy, grad, fixion, scratch, stdout, energy_thr, &
                     grad_thr, energy_error, grad_error, istep, nstep,       &
                     step_accepted, stop_bfgs )
      !------------------------------------------------------------------------
      !
      ! ... list of input/output arguments : 
      !  
      !  pos            : vector containing 3N coordinates of the system ( x )
      !  energy         : energy of the system ( V(x) )
      !  grad           : vector containing 3N components of grad( V(x) )
      !  fixion         : vector used to freeze a deg. of freedom
      !  scratch        : scratch directory
      !  stdout         : unit for standard output
      !  energy_thr     : treshold on energy difference for BFGS convergence
      !  grad_thr       : treshold on grad difference for BFGS convergence
      !                    the largest component of grad( V(x) ) is considered
      !  energy_error   : energy difference | V(x_i) - V(x_i-1) |
      !  grad_error     : the largest component of 
      !                    | grad(V(x_i)) - grad(V(x_i-1)) |
      !  nstep          : the maximun nuber of scf-steps
      !  step_accepted  : .TRUE. if a new BFGS step is done
      !  conv_bfgs      : .TRUE. if BFGS convergence has been achieved
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(INOUT) :: pos(:)
      REAL(DP),         INTENT(INOUT) :: energy
      REAL(DP),         INTENT(INOUT) :: grad(:)
      INTEGER,          INTENT(IN)    :: fixion(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: stdout
      REAL(DP),         INTENT(IN)    :: energy_thr, grad_thr
      INTEGER,          INTENT(OUT)   :: istep
      INTEGER,          INTENT(IN)    :: nstep
      REAL(DP),         INTENT(OUT)   :: energy_error, grad_error
      LOGICAL,          INTENT(OUT)   :: step_accepted, stop_bfgs
      !
      INTEGER  :: dim, i, j
      LOGICAL  :: lwolfe
      REAL(DP) :: dE0s, den
      !
      !
      dim = SIZE( pos )
      !
      ! ... conditional work-space allocation
      !   
      IF ( .NOT. ALLOCATED( grad_old ) )  ALLOCATE( grad_old( dim, bfgs_ndim ) )
      IF ( .NOT. ALLOCATED( pos_old ) )   ALLOCATE( pos_old(  dim, bfgs_ndim ) )
      !
      IF ( .NOT. ALLOCATED( inv_hess ) )  ALLOCATE( inv_hess(  dim, dim ) )
      IF ( .NOT. ALLOCATED( inv_hess0 ) ) ALLOCATE( inv_hess0( dim, dim ) )
      IF ( .NOT. ALLOCATED( projmat ) )   ALLOCATE( projmat(   dim, dim ) )
      !
      IF ( .NOT. ALLOCATED( pos_p ) )     ALLOCATE( pos_p(    dim ) )
      IF ( .NOT. ALLOCATED( grad_p ) )    ALLOCATE( grad_p(   dim ) )
      IF ( .NOT. ALLOCATED( step ) )      ALLOCATE( step(     dim ) )
      IF ( .NOT. ALLOCATED( step_old ) )  ALLOCATE( step_old( dim ) )
      IF ( .NOT. ALLOCATED( pos_best ) )  ALLOCATE( pos_best( dim ) )
      !
      IF ( .NOT. ALLOCATED( dirtm ) )     ALLOCATE( dirtm( dim, dim ) )
      IF ( .NOT. ALLOCATED( invtm ) )     ALLOCATE( invtm( dim, dim ) )
      !
      ! ... the BFGS file is read (pos & grad are converted here to
      ! ... internal coordinates )
      !
      CALL read_bfgs_file( pos, grad, fixion, energy, scratch, dim, stdout )
      !
      scf_iter = scf_iter + 1
      istep    = scf_iter
      !
      ! ... convergence is checked here
      !
      energy_error = ABS( energy_p - energy )
      grad_error   = MAXVAL( ABS( grad(:) ) )
      !
      conv_bfgs = energy_error < energy_thr      
      conv_bfgs = conv_bfgs .AND. ( grad_error < grad_thr )
      !
      stop_bfgs = conv_bfgs .OR. ( scf_iter >= nstep )
      !
      IF ( stop_bfgs ) THEN
         !
         ! ... quick return ( positions and forces are back-converted to
         ! ... cartesian coordinates )
         !
         pos(:)  = to_cart_coords( pos,  dim )
         grad(:) = to_cart_coords( grad, dim )
         !
         RETURN
         !
      END IF
      !
      ! ... some output is written
      !
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"number of scf cycles",T30,"= ",I3)' ) scf_iter
      WRITE( UNIT = stdout, &
           & FMT = '(5X,"number of bfgs steps",T30,"= ",I3,/)' ) bfgs_iter
      IF ( scf_iter > 1 ) &
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"energy old",T30,"= ",F18.10," Ry")' ) energy_p
      WRITE( UNIT = stdout, &
           & FMT = '(5X,"energy new",T30,"= ",F18.10," Ry",/)' ) energy
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
         ! ... the new trust radius is obtained by quadratic interpolation
         !
         ! ... E(s) = a*s*s + b*s + c      ( we use E(0), dE(0), E(s') )
         !
         ! ... s_min = - 0.5*( dE(0)*s'*s' ) / ( E(s') - E(0) - dE(0)*s' )
         !
         dE0s = ( grad_p(:) .dot. step_old(:) )
         !
         den = energy - energy_p - dE0s
         !              
         IF ( den > eps16 ) THEN
            !
            trust_radius = - 0.5D0 * dE0s * trust_radius_old / den
            !
         ELSE
            !
            ! ... no quadratic interpolation is possible: we use bisection
            !
            trust_radius = 0.5D0 * trust_radius_old
            !
         END IF
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius          
         !
         ! ... values from the last succeseful bfgs step are restored
         !
         pos(:)  = pos_p(:)
         energy  = energy_p
         grad(:) = grad_p(:)         
         !
         IF ( trust_radius < trust_radius_min ) THEN
            !
            ! ... the history is reset ( the history can be reset at most two
            ! ... consecutive times )
            !
            WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"trust_radius < trust_radius_min")' )
            WRITE( UNIT = stdout, FMT = '(/,5X,"resetting bfgs history",/)' )
            !
            IF ( tr_min_hit ) THEN
               !
               ! ... the history has already been reset at the previous step :
               ! ... something is going wrong
               !
               CALL errore( 'bfgs', &
                            'bfgs history already reset at previous step', 1 )
               !
            END IF
            !
            CALL reset_bfgs()
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            !
            trust_radius = trust_radius_min
            !
            tr_min_hit = .TRUE.
            !
         ELSE
            !
            ! ... old bfgs direction ( normalized ) is recovered
            !
            step(:) = step_old(:) / trust_radius_old
            !
            tr_min_hit = .FALSE.
            !
         END IF
         !
      ELSE    
         !
         ! ... a new bfgs step is done
         !
         bfgs_iter = bfgs_iter + 1
         !
         IF ( bfgs_iter == 1 ) THEN
            !
            ! ... first iteration
            !
            IF ( grad_error < 0.01D0 ) &
               trust_radius_ini = MIN( 0.2D0, trust_radius_ini )
            !
            step_accepted = .FALSE.
            !
         ELSE
            !
            step_accepted = .TRUE.
            !
            WRITE( UNIT = stdout, &
                 & FMT = '(5X,"CASE: energy_new < energy_old",/)' )
            !
            CALL check_wolfe_conditions( lwolfe, energy, grad )
            !
            CALL update_inverse_hessian( pos, grad, dim, stdout )
            !
         END IF
         !
         IF ( bfgs_ndim > 1 ) THEN
            !
            ! ... GDIIS extrapolation
            !
            CALL gdiis_step()
            !
         ELSE
            !
            ! ... standard Newton-Raphson step
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            !
         END IF         
         !
         IF ( ( grad(:) .dot. step(:) ) > 0.D0 ) THEN
            !
            WRITE( UNIT = stdout, &
                   FMT = '(5X,"uphill step: resetting bfgs history",/)' )
            !
            CALL reset_bfgs()
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            !
         END IF
         !  
         ! ... the new trust radius is computed
         !
         IF ( bfgs_iter == 1 ) THEN
            !
            trust_radius = trust_radius_ini
            !
            tr_min_hit = .FALSE.
            !
         ELSE
            !
            trust_radius = trust_radius_old
            !
            CALL compute_trust_radius( lwolfe, energy, grad, dim, stdout )
            !   
         END IF
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius
         !
      END IF
      !
      ! ... step along the bfgs direction
      !
      IF ( norm( step(:) ) < eps16 ) &
         CALL errore( 'bfgs', 'NR step-length unreasonably short', 1 )
      !
      step(:) = trust_radius * step(:) / norm( step(:) )
      !
      ! ... information required by next iteration is saved here ( this must
      ! ... be done before positions are updated )
      !
      CALL write_bfgs_file( pos, energy, grad, scratch )                
      !
      ! ... positions are updated and then converted to cartesian coordinates
      !
      pos(:) = pos(:) + step(:)
      !
      pos(:)  = to_cart_coords( pos,  dim )
      grad(:) = to_cart_coords( grad, dim )
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE gdiis_step()
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          REAL(DP), ALLOCATABLE :: res(:,:), overlap(:,:), work(:)
          INTEGER,  ALLOCATABLE :: iwork(:)
          INTEGER               :: k, k_m, info
          REAL(DP)              :: gamma0
          !
          !
          gdiis_iter = gdiis_iter + 1
          !
          k   = MIN( gdiis_iter, bfgs_ndim )
          k_m = k + 1
          !
          ALLOCATE( res( dim, k ) )
          ALLOCATE( overlap( k_m, k_m ) )
          ALLOCATE( work( k_m ), iwork( k_m ) )
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
          CALL DGEMM( 'N', 'N', dim, k, dim, 1.D0, &
                      inv_hess, dim, grad_old, dim, 0.D0, res, dim )
          !
          ! ... overlap_ij = <grad_i|res_j>
          !
          CALL DGEMM( 'T', 'N', k, k, dim, 1.D0, &
                      res, dim, res, dim, 0.D0, overlap, k_m )
          !
          overlap( :, k_m) = 1.D0
          overlap(k_m, : ) = 1.D0
          overlap(k_m,k_m) = 0.D0
          !
          ! ... overlap is inverted via Bunch-Kaufman diagonal pivoting method
          !
          CALL DSYTRF( 'U', k_m, overlap, k_m, iwork, work, k_m, info )
          CALL DSYTRI( 'U', k_m, overlap, k_m, iwork, work, info )
          !
          ! ... overlap is symmetrised
          !
          FORALL( i = 1:k_m, &
                  j = 1:k_m, j > i ) overlap(j,i) = overlap(i,j)
          !
          pos_best(:) = 0.D0
          step(:)     = 0.D0
          !
          DO i = 1, k
             !
             gamma0 = overlap(k_m,i)
             !
             pos_best(:) = pos_best(:) + gamma0*pos_old(:,i)
             !
             step(:) = step(:) - gamma0*res(:,i)
             !
          END DO
          !
          ! ... the step must be consistent with the old positions
          !
          step(:) = step(:) + ( pos_best(:) - pos(:) )
          !
          IF ( ( grad(:) .dot. step(:) ) > 0.D0 ) THEN
             !
             ! ... if the extrapolated direction is uphill use only the
             ! ... last gradient and reset gdiis history
             !
             step(:) = - ( inv_hess(:,:) .times. grad(:) )
             !
             gdiis_iter = 0
             !
          END IF
          !
          DEALLOCATE( res, overlap, work, iwork )
          !
        END SUBROUTINE gdiis_step
        !      
    END SUBROUTINE bfgs
    !
    !------------------------------------------------------------------------
    SUBROUTINE reset_bfgs()
      !------------------------------------------------------------------------
      !
      ! ... inv_hess0 contains the initial guess for the inverse hessian matrix
      ! ... in internal units
      !
      inv_hess(:,:) = inv_hess0(:,:)
      !
      gdiis_iter = 0
      !
    END SUBROUTINE reset_bfgs
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_bfgs_file( pos, grad, fixion, energy, scratch, dim, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(INOUT) :: pos(:)
      REAL(DP),         INTENT(INOUT) :: grad(:)
      INTEGER,          INTENT(IN)    :: fixion(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: dim
      INTEGER,          INTENT(IN)    :: stdout
      REAL(DP),         INTENT(INOUT) :: energy
      !
      CHARACTER(LEN=256) :: bfgs_file
      LOGICAL            :: file_exists
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
         READ( iunbfgs, * ) tr_min_hit
         READ( iunbfgs, * ) inv_hess0
         READ( iunbfgs, * ) dirtm
         !     
         CLOSE( UNIT = iunbfgs )
         !
         invtm(:,:) = TRANSPOSE( dirtm(:,:) )
         !
         pos(:)  = to_internal_coords( pos,  dim )
         grad(:) = to_internal_coords( grad, dim )
         !
         trust_radius_old = norm( pos(:) - pos_p(:) )
         !
         step_old = ( pos(:) - pos_p(:) ) / trust_radius_old
         !
      ELSE
         !
         ! ... bfgs initialization
         !
         WRITE( UNIT = stdout, FMT = '(/,5X,"BFGS Geometry Optimization")' )
         !
         CALL generate_hess( pos, dim )
         !
         CALL init_internal_coords( dim, fixion )
         !
         pos(:)  = to_internal_coords( pos,  dim )
         grad(:) = to_internal_coords( grad, dim )
         !
         pos_p     = 0.D0
         grad_p    = 0.D0
         scf_iter  = 0
         bfgs_iter = 0
         gdiis_iter= 0
         energy_p  = energy
         step_old  = 0.D0
         !
         trust_radius_old = trust_radius_ini
         !
         pos_old(:,:)  = 0.D0
         grad_old(:,:) = 0.D0
         !
         tr_min_hit = .FALSE.
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
      WRITE( iunbfgs, * ) tr_min_hit
      WRITE( iunbfgs, * ) inv_hess0
      WRITE( iunbfgs, * ) dirtm
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
      REAL(DP), ALLOCATABLE :: y(:), s(:)
      REAL(DP), ALLOCATABLE :: Hy(:), yH(:)
      REAL(DP), ALLOCATABLE :: H_bfgs(:,:)
      REAL(DP)              :: sdoty
      !
      !
      ALLOCATE( y( dim ), s( dim ), Hy( dim ), yH( dim ), H_bfgs( dim, dim ) )
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
         WRITE( stdout, '(/,5X,"WARINIG: unexpected ", &
                         &     "behaviour in update_inverse_hessian")' )
         WRITE( stdout, '(  5X,"         resetting bfgs history",/)' )
         !
         CALL reset_bfgs()
         !
         RETURN
         !
      END IF
      !
      Hy(:) = ( inv_hess .times. y(:) )
      yH(:) = ( y(:) .times. inv_hess )
      !
      ! ... BFGS update
      !
      inv_hess = inv_hess + 1.D0 / sdoty * &
                      ( ( 1.D0 + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                        ( matrix( s, yH ) +  matrix( Hy, s ) ) )
      !
      DEALLOCATE( y, s, Hy, yH, H_bfgs )
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
         IF ( tr_min_hit ) THEN
            !
            ! ... the history has already been reset at the previous step :
            ! ... something is going wrong
            !
            CALL errore( 'bfgs', 'history already reset at previous step', 1 )
            !
         END IF
         !
         WRITE( UNIT = stdout, &
                FMT = '(5X,"small trust_radius: resetting bfgs history",/)' )
         !
         CALL reset_bfgs()
         !
         step(:) = - ( inv_hess(:,:) .times. grad(:) )
         !
         trust_radius = trust_radius_min
         !
         tr_min_hit = .TRUE.
         !
      ELSE
         !
         tr_min_hit = .FALSE.
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
      IF ( conv_bfgs ) THEN
         !
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"bfgs converged in ",I3," scf cycles and ", &
              &         I3," bfgs steps")' ) scf_iter, bfgs_iter
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"Final energy = ",F18.10," Ry")' ) energy
         !
         CALL delete_if_present( TRIM( scratch ) // TRIM( prefix ) // '.bfgs' )
         !
      ELSE
         !
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"The maximum number of steps has been reached.")' )
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         !
      END IF
      !
      IF ( ALLOCATED( pos_p ) )      DEALLOCATE( pos_p )
      IF ( ALLOCATED( grad_p ) )     DEALLOCATE( grad_p )
      IF ( ALLOCATED( pos_old ) )    DEALLOCATE( pos_old )
      IF ( ALLOCATED( grad_old ) )   DEALLOCATE( grad_old )
      IF ( ALLOCATED( inv_hess ) )   DEALLOCATE( inv_hess )
      IF ( ALLOCATED( inv_hess0 ) )  DEALLOCATE( inv_hess0 )
      IF ( ALLOCATED( step) )        DEALLOCATE( step )
      IF ( ALLOCATED( step_old ) )   DEALLOCATE( step_old )
      IF ( ALLOCATED( pos_best ) )   DEALLOCATE( pos_best )
      IF ( ALLOCATED( dirtm ) )      DEALLOCATE( dirtm )
      IF ( ALLOCATED( invtm ) )      DEALLOCATE( invtm )
      !
    END SUBROUTINE terminate_bfgs
    !
    !------------------------------------------------------------------------
    SUBROUTINE init_internal_coords( dim, fixion )
      !------------------------------------------------------------------------
      !
      ! ... inv_hess(:,:) contains the initial estimate of the hessian matrix
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: dim
      INTEGER, INTENT(IN) :: fixion(:)
      !
      INTEGER               :: i, info, lwork
      REAL(DP), ALLOCATABLE :: w(:), work(:)
      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! ... model hessian disabled (to enable it, comment out these lines)
!!
inv_hess(:,:)  = identity( dim )
dirtm(:,:)     = identity( dim )
invtm(:,:)     = identity( dim )
inv_hess0(:,:) = inv_hess(:,:)
RETURN
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      lwork = 3*dim
      !
      ALLOCATE( w( dim ), work( lwork ) )
      !
      dirtm(:,:) = inv_hess(:,:)
      !
      CALL DSYEV( 'V', 'U', dim, dirtm, dim, w, work, lwork, info )
      !
      CALL errore( 'init_internal_coords', &
                   'error in the diagonalisation of starting hessian', info )
      !
      invtm(:,:) = TRANSPOSE( dirtm(:,:) )
      !
      inv_hess(:,:) = 0.D0
      FORALL( i = 4:dim ) inv_hess(i,i) = 1.D0 / w(i)
      !
      DEALLOCATE( w, work )
      !
      IF ( dim <= 12 ) THEN
         !
         PRINT *, "DIAGONAL INVERSE HESSIAN"
         !
         DO i = 1, dim
            !
            PRINT '(12(1X,F10.6))', inv_hess(i,:)
            !
         END DO
         !
         PRINT *, "DIRECT TRANSFORMATION MATRIX"
         !
         DO i = 1, dim
            !
            PRINT '(12(1X,F10.6))', dirtm(i,:)
            !
         END DO
         !
         PRINT *, "INVERSE TRANSFORMATION MATRIX"
         !
         DO i = 1, dim
            !
            PRINT '(12(1X,F10.6))', invtm(i,:)
            !
         END DO
         !
      END IF
      !
      projmat(:,:) = 0.D0
      !
      FORALL( i = 1:dim ) projmat(i,i) = DBLE( fixion(i) )
      !
      projmat(:,:) = MATMUL( dirtm(:,:), MATMUL( projmat(:,:), invtm(:,:) ) )
      !
      inv_hess(:,:) = MATMUL( projmat(:,:), inv_hess(:,:) )
      !
      inv_hess0(:,:) = inv_hess(:,:)

      IF ( dim <= 12 ) THEN
         !
         PRINT *, "INVERSE HESSIAN AGAIN"
         !
         DO i = 1, dim
            !
            PRINT '(12(1X,F10.6))', inv_hess(i,:)
            !
         END DO
         !
      END IF
      !
    END SUBROUTINE init_internal_coords
    !
    !------------------------------------------------------------------------
    FUNCTION to_internal_coords( vect_in, dim ) RESULT( vect_out )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN) :: dim
      REAL(DP), INTENT(IN) :: vect_in(dim)
      REAL(DP)             :: vect_out(dim)
      !
      vect_out(:) = ( dirtm(:,:) .times. vect_in(:) )
      !
    END FUNCTION to_internal_coords
    !
    !------------------------------------------------------------------------
    FUNCTION to_cart_coords( vect_in, dim ) RESULT( vect_out )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN) :: dim
      REAL(DP), INTENT(IN) :: vect_in(dim)
      REAL(DP)             :: vect_out(dim)
      !
      vect_out(:) = ( invtm(:,:) .times. vect_in(:) )
      !
    END FUNCTION to_cart_coords
    !
    !------------------------------------------------------------------------
    SUBROUTINE generate_hess( pos, dim )
      !------------------------------------------------------------------------
      !
      ! ... the hessian matrix is termporarly stored into inv_hess(:,:)
      !
      USE constants,          ONLY : rytoev
      USE constraints_module, ONLY : pbc
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN) :: dim
      REAL(DP), INTENT(IN) :: pos(:)
      !
      REAL(DP), ALLOCATABLE :: tau(:,:)
      REAL(DP)              :: k_ij, dist_ij
      INTEGER               :: i, j, in, jn
      INTEGER               :: nat
      REAL(DP)              :: r_cov_max
      !
      REAL(DP), PARAMETER :: r_cov = 2.23D0
      REAL(DP), PARAMETER :: k0    = 6174.45D0
      !
      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! ... model hessian disabled (to enable it, comment out these lines)
!!      
inv_hess(:,:) = identity( dim ); RETURN
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      nat = dim / 3
      !
      r_cov_max = r_cov
      !
      ALLOCATE( tau( 3, nat ) )
      !
      tau(:,:) = RESHAPE( pos, (/ 3, nat /) )
      !
      inv_hess(:,:) = 0.D0
      !
      DO in = 1, nat
         !
         i = 3*in - 2
         !
         DO jn = 1, nat
            !
            IF ( jn == in ) CYCLE
            !
            j = 3*jn - 2
            !
            dist_ij = norm( pbc( tau(:,in) - tau(:,jn) ) )
            !
            IF ( dist_ij > 6.D0*r_cov_max ) CYCLE
            !
            k_ij = k0*( ( r_cov + r_cov ) / dist_ij )**8
            !
            inv_hess(i:i+2,i:i+2) = inv_hess(i:i+2,i:i+2) + k_ij*identity( 3 )
            !
            inv_hess(j:j+2,i:i+2) = - k_ij*identity( 3 )
            !
         END DO
         !
      END DO
      !
      IF ( dim <= 12 ) THEN
         !
         PRINT *, "HESSIAN"
         !
         DO i = 1, dim
            !
            PRINT '(12(1X,F10.3))', inv_hess(i,:)
            !
         END DO
         !
      END IF
      !
      DEALLOCATE( tau )
      !
    END SUBROUTINE generate_hess
    !
END MODULE bfgs_module
