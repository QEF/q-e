!
! Copyright (C) 2003-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE bfgs_module
  !----------------------------------------------------------------------------
  !
  ! ... ionic relaxation through Broyden-Fletcher-Goldfarb-Shanno 
  ! ... minimization and a "trust radius" line search based on 
  ! ... Wolfe conditions ( bfgs() subroutine )
  ! ... A linear scaling BFGS is also implemented ( lin_bfgs() subroutine )
  ! ... Both subroutines are called with the same list of arguments
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
  USE kinds,       ONLY : DP
  USE io_files,    ONLY : iunbfgs
  !
  USE basic_algebra_routines  
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: bfgs,             &
            lin_bfgs
  !
  ! ... public variables
  !   
  PUBLIC :: lbfgs_ndim,       &
            trust_radius_max, &
            trust_radius_min, &
            trust_radius_ini, &
            trust_radius_end, &
            w_1,              &
            w_2
  !
  ! ... global variables
  !
  REAL(KIND=DP), ALLOCATABLE :: &
      pos_old(:,:),             &! list of m old positions ( m = 1 for 
                                 ! standard BFGS algorithm )
      inverse_hessian(:,:),     &! inverse of the hessian matrix (updated via
                                 ! BFGS formula)
      bfgs_step(:),             &! bfgs direction
      bfgs_step_old(:),         &! old bfgs direction
      gradient_old(:,:)          ! list of m old gradients ( m = 1 for 
                                 ! standard BFGS algorithm )
  INTEGER :: &
      lbfgs_ndim = 4             ! dimension of the subspace for L-BFGS
                                 ! fixed to 1 for standard BFGS algorithm
  REAL(KIND=DP) :: &   
      trust_radius,             &! displacement along the bfgs direction
      trust_radius_old,         &! old displacement along the bfgs direction
      energy_old                 ! old energy
  INTEGER :: &
      scf_iter,                 &! number of scf iterations
      bfgs_iter,                &! number of bfgs iterations
      lin_iter                   ! number of line search iterations    
  REAL(KIND=DP)  :: &
      trust_radius_max = 0.5D0, &! maximum allowed displacement
      trust_radius_min = 1.D-5, &! minimum allowed displacement
      trust_radius_ini = 0.5D0, &! initial displacement
      trust_radius_end = 1.D-7   ! bfgs stops when trust_radius is less than
                                 ! this value
  REAL(KIND=DP)  :: &
      w_1 = 0.5D-1,             &! parameters for Wolfe conditions
      w_2 = 0.5D0                ! parameters for Wolfe conditions
  !
  ! ... Note that m, trust_radius_max, trust_radius_min, trust_radius_ini,
  ! ... trust_radius_end, w_1, w_2 have a default value, but can also be 
  ! ... assigned in input
  !  
  ! Workaround for IBM xlf compiler bug
  !
  INTEGER, PUBLIC :: bfgs_xlf_bug
  !
  CONTAINS
     !
     !
     ! ... public methods :
     !
     !-----------------------------------------------------------------------
     SUBROUTINE bfgs( pos, energy, gradient, scratch, stdout, energy_thr, &
                      gradient_thr, energy_error, gradient_error,         &
                      step_accepted, conv_bfgs )
       !-----------------------------------------------------------------------
       !
       ! ... list of input/output arguments : 
       !  
       !  pos            : vector containing 3N coordinates of the system ( x )
       !  energy         : energy of the system ( V(x) )
       !  gradient       : vector containing 3N components of ( grad( V(x) ) ) 
       !  scratch        : scratch diercotry
       !  stdout         : unit for standard output
       !  energy_thr     : treshold on energy difference for BFGS convergence
       !  gradient_thr   : treshold on gradient difference for BFGS convergence
       !                    the largest component of grad( V(x) ) is considered
       !  energy_error   : energy difference | V(x_i) - V(x_i-1) |
       !  gradient_error : the largest component of 
       !                    | grad(V(x_i)) - grad(V(x_i-1)) | 
       !  step_accepted  : .TRUE. if a new BFGS step is done
       !  conv_bfgs      : .TRUE. if BFGS convergence has been achieved
       !
       USE constants,  ONLY : eps16
       !
       IMPLICIT NONE
       !
       ! ... input/output arguments
       !
       REAL(KIND=DP),     INTENT(INOUT) :: pos(:)
       REAL(KIND=DP),     INTENT(INOUT) :: energy       
       REAL(KIND=DP),     INTENT(INOUT) :: gradient(:)
       CHARACTER (LEN=*), INTENT(IN)    :: scratch
       INTEGER,           INTENT(IN)    :: stdout   
       REAL(KIND=DP),     INTENT(IN)    :: energy_thr, gradient_thr  
       REAL(KIND=DP),     INTENT(OUT)   :: energy_error, gradient_error       
       LOGICAL,           INTENT(OUT)   :: step_accepted, conv_bfgs
       !
       ! ... local variables
       !
       INTEGER  :: dim, i
       LOGICAL  :: lwolfe
       !
       !
       dim = SIZE( pos )
       !
       ! ... lbfgs_ndim is forced to be equal to 1 ( the complete inverse 
       ! ... hessian  matrix is stored )
       !
       lbfgs_ndim  = 1
       !
       ALLOCATE( pos_old( dim, lbfgs_ndim ) )
       ALLOCATE( inverse_hessian( dim, dim ) )
       ALLOCATE( bfgs_step( dim ) )              
       ALLOCATE( bfgs_step_old( dim ) )
       ALLOCATE( gradient_old( dim, lbfgs_ndim ) )
       !       
       CALL read_bfgs_file( pos, energy, gradient, scratch, dim )
       !
       scf_iter = scf_iter + 1       
       !       
       conv_bfgs = ( ( energy_old - energy ) < energy_thr )
       !
       energy_error   = ABS( energy_old - energy )
       gradient_error = 0.D0
       !
       DO i = 1, dim
          !
          conv_bfgs = ( conv_bfgs .AND. ( ABS( gradient(i) ) < gradient_thr ) )
          !
          gradient_error = MAX( gradient_error, ABS( gradient(i) ) )
          !
       END DO       
       !
       ! ... as long as the first two scf iterations have been performed the 
       ! ... error on the energy is redefined as a "large" number.
       !
       IF ( scf_iter < 2 ) energy_error = 1000.D0
       !
       IF ( conv_bfgs ) THEN
          !
          CALL terminate_bfgs( energy, stdout, scratch )
          !
          RETURN
          !
       END IF
       !
       WRITE( UNIT = stdout, &
            & FMT = '(/,5X,"number of ionic steps",T30,"= ",I3)' ) scf_iter
       WRITE( UNIT = stdout, &
            & FMT = '(  5X,"number of bfgs  steps",T30,"= ",I3,/)' ) bfgs_iter
       IF ( scf_iter > 1 ) &
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"energy old",T30,"= ",F18.10," ryd")' ) energy_old
       WRITE( UNIT = stdout, &
            & FMT = '(5X,"energy new",T30,"= ",F18.10," ryd",/)' ) energy
       !
       IF ( energy > energy_old ) THEN
          !
          ! ... the previous step is rejected, line search goes on
          !
          step_accepted = .FALSE.          
          !	  
          lin_iter = lin_iter + 1
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"CASE: energy_new > energy_old",/)' )
          !
          ! ... the old trust radius is reduced by a factor 2
          !
          trust_radius = 0.5D0 * trust_radius_old
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr",/)' ) &
              trust_radius
          !
          IF ( trust_radius < trust_radius_min ) THEN
             !
             ! ... the history is reset
             !     
             WRITE( UNIT = stdout, FMT = '(/,5X,"resetting bfgs history",/)' )
             !
             inverse_hessian = identity(dim)
             !
             bfgs_step = - ( inverse_hessian .times. gradient )
             !
             trust_radius = trust_radius_min
             !
          ELSE 
             !
             ! ... values of the last succeseful bfgs step are restored
             !
             pos      = pos_old(:,1)
             energy   = energy_old
             gradient = gradient_old(:,1)
             !
             ! ... old bfgs direction ( normalized ) is recovered
             !
             bfgs_step = bfgs_step_old / trust_radius_old
             !
          END IF   
          !
       ELSE    
          !
          ! ... a new bfgs step is done
          !
          step_accepted = .TRUE.
          !
          lin_iter  = 1
          bfgs_iter = bfgs_iter + 1
          !
          IF ( bfgs_iter > 1 ) THEN
             !
             WRITE( UNIT = stdout, &
                  & FMT = '(5X,"CASE: energy_new < energy_old",/)' )
             !
             CALL check_wolfe_conditions( lwolfe, energy, gradient )
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
             CALL update_inverse_hessian( gradient, dim, stdout )
             !
          END IF
          !
          ! ... bfgs direction ( not normalized ) 
          !
          bfgs_step = - ( inverse_hessian .times. gradient )
          !
          IF ( ( gradient .dot. bfgs_step ) > 0.D0 ) THEN
             !  
             ! ... bfgs direction is reversed if not downhill
             !
             bfgs_step = - bfgs_step
             !
             WRITE( UNIT = stdout, FMT = '(5X,"search direction reversed",/)' )
             !
             ! ... the history is reset
             !     
             WRITE( UNIT = stdout, FMT = '(5X,"resetting bfgs history",/)' )
             !
             inverse_hessian = identity(dim)
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
             CALL compute_trust_radius( lwolfe, energy, gradient, dim, &
                                        stdout, conv_bfgs )
             !
          END IF
          !
          ! ... if trust_radius < trust_radius_end convergence is achieved
          !
          IF ( conv_bfgs ) THEN
             !
             CALL terminate_bfgs( energy, stdout, scratch )
             !
             RETURN
             !
          END IF
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr",/)' ) &
              trust_radius
          !
       END IF  
       !
       ! ... step along the bfgs direction
       !
       IF ( norm( bfgs_step ) < eps16 ) THEN
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"WARNING : norm( bfgs_step )",T30,"= ",F18.10)' ) &
              norm( bfgs_step )
          !
          bfgs_step = - gradient
          !
       ELSE
          !
          bfgs_step = trust_radius * bfgs_step / norm( bfgs_step )
          !
       END IF 
       !
       ! ... informations needed for the next iteration are saved
       ! ... this must be done before positions update
       !
       CALL write_bfgs_file( pos, energy, gradient, scratch )                
       !
       ! ... positions are updated
       !
       pos = pos + bfgs_step
       !
       DEALLOCATE( pos_old )   
       DEALLOCATE( inverse_hessian )
       DEALLOCATE( bfgs_step )       
       DEALLOCATE( bfgs_step_old )
       DEALLOCATE( gradient_old )       
       !
     END SUBROUTINE bfgs
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE lin_bfgs( pos, energy, gradient, scratch, stdout, energy_thr, &
                          gradient_thr, energy_error, gradient_error,         &
                          step_accepted, conv_bfgs )
       !-----------------------------------------------------------------------
       !
       ! ... list of input/output arguments : 
       !  
       !  pos            : vector containing 3N coordinates of the system ( x )
       !  energy         : energy of the system ( V(x) )
       !  gradient       : vector containing 3N components of ( grad( V(x) ) ) 
       !  scratch        : scratch diercotry
       !  stdout         : unit for standard output
       !  energy_thr     : treshold on energy difference for BFGS convergence
       !  gradient_thr   : treshold on gradient difference for BFGS convergence
       !                    the largest component of grad( V(x) ) is considered
       !  energy_error   : energy difference | V(x_i) - V(x_i-1) |
       !  gradient_error : the largest component of 
       !                    | grad(V(x_i)) - grad(V(x_i-1)) | 
       !  step_accepted  : .TRUE. if a new BFGS step is done
       !  conv_bfgs      : .TRUE. if BFGS convergence has been achieved       
       !
       USE constants,  ONLY : eps16       
       !
       IMPLICIT NONE
       !
       !
       ! ... input/output arguments
       !
       REAL(KIND=DP),     INTENT(INOUT) :: pos(:)
       REAL(KIND=DP),     INTENT(INOUT) :: energy       
       REAL(KIND=DP),     INTENT(INOUT) :: gradient(:)
       CHARACTER (LEN=*), INTENT(IN)    :: scratch
       INTEGER,           INTENT(IN)    :: stdout   
       REAL(KIND=DP),     INTENT(IN)    :: energy_thr, gradient_thr  
       REAL(KIND=DP),     INTENT(OUT)   :: energy_error, gradient_error       
       LOGICAL,           INTENT(OUT)   :: step_accepted, conv_bfgs
       !
       ! ... local variables
       !
       INTEGER   :: dim, i
       LOGICAL   :: lwolfe
       !
       !
       dim = SIZE( pos )
       !
       ALLOCATE( pos_old( dim, lbfgs_ndim ) )
       ALLOCATE( gradient_old( dim, lbfgs_ndim ) )       
       ALLOCATE( bfgs_step( dim ) )              
       ALLOCATE( bfgs_step_old( dim ) )
       !       
       CALL read_lbfgs_file( pos, energy, gradient, scratch, dim )
       !
       scf_iter = scf_iter + 1       
       !       
       ! ... convergence is checked
       !
       conv_bfgs = ( ( energy_old - energy ) < energy_thr )
       !
       energy_error   = ABS( energy_old - energy )
       gradient_error = 0.D0
       !
       DO i = 1, dim
          !
          conv_bfgs = ( conv_bfgs .AND. ( ABS( gradient(i) ) < gradient_thr ) )
          !
          gradient_error = MAX( gradient_error, ABS( gradient(i) ) )
          !
       END DO       
       !
       ! ... as long as the first two scf iterations have been performed the 
       ! ... error on the energy is redefined as a "large" number.
       !
       IF ( scf_iter < 2 ) energy_error = 1000.D0
       !
       IF ( conv_bfgs ) THEN
          !
          ! ... convergence has been achieved
          !
          CALL terminate_bfgs( energy, stdout, scratch )
          !
          RETURN
          !
       END IF
       !
       WRITE( UNIT = stdout, &
            & FMT = '(/,5X,"number of ionic steps",T30,"= ",I3)' ) scf_iter
       WRITE( UNIT = stdout, &
            & FMT = '(  5X,"number of bfgs  steps",T30,"= ",I3,/)' ) bfgs_iter
       IF ( scf_iter > 1 ) &
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"energy old",T30,"= ",F18.10," ryd")' ) energy_old
       WRITE( UNIT = stdout, &
            & FMT = '(5X,"energy new",T30,"= ",F18.10," ryd",/)' ) energy
       !
       IF ( energy > energy_old ) THEN
          !
          ! ... the previous step is rejected, line search goes on
          !
          step_accepted = .FALSE.          
          !	  
          lin_iter = lin_iter + 1
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"CASE: energy_new > energy_old",/)' )
          !
          ! ... the old trust radius is reduced by a factor 2
          !
          trust_radius = 0.5D0 * trust_radius_old
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr",/)' ) &
              trust_radius
          !
          IF ( trust_radius < trust_radius_min ) THEN
             !
             ! ... the history is reset
             !     
             WRITE( UNIT = stdout, FMT = '(5X,"resetting bfgs history",/)' )
             !
             pos_old      = 0.D0
             gradient_old = 0.D0
             !
             bfgs_step = - gradient
             !
             trust_radius = trust_radius_min
             !
          ELSE 
             !
             ! ... values of the last succeseful bfgs step are restored
             !
             pos       = pos_old(:,lbfgs_ndim)
             energy    = energy_old
             gradient  = gradient_old(:,lbfgs_ndim)
             !
             ! ... old bfgs direction ( normalized ) is recovered
             !
             bfgs_step = bfgs_step_old / trust_radius_old
             !
          END IF   
          !
       ELSE    
          !
          ! ... a new bfgs step is done
          !
          step_accepted = .TRUE.
          !
          lin_iter  = 1
          bfgs_iter = bfgs_iter + 1
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"CASE: energy_new < energy_old",/)' )
          !
          ! ... Wolfe conditions and hessian update are needed after
          ! ... the first bfgs iteration
          !
          IF ( bfgs_iter == 1 ) THEN
             !
             bfgs_step = - gradient
             !
          ELSE
             !
             CALL check_wolfe_conditions( lwolfe, energy, gradient )
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
             CALL lbfgs_update( pos, gradient, dim )
             !
          END IF   
          !
          IF ( ( gradient .dot. bfgs_step ) > 0.D0 ) THEN
             !  
             ! ... bfgs direction is reversed if not downhill
             !
             bfgs_step = - bfgs_step
             !
             WRITE( UNIT = stdout, FMT = '(5X,"search direction reversed")' )
             !
             ! ... the history is reset
             !     
             WRITE( UNIT = stdout, FMT = '(5X,"resetting bfgs history",/)' )
             !
             pos_old      = 0.D0
             gradient_old = 0.D0
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
             CALL compute_trust_radius( lwolfe, energy, gradient, dim, &
                                        stdout, conv_bfgs )
             !
          END IF
          !
          ! ... if trust_radius < trust_radius_end convergence is achieved
          ! ... this should be a "rare event"
          !
          IF ( conv_bfgs ) THEN
             !
             CALL terminate_bfgs( energy, stdout, scratch )
             !
             RETURN
             !
          END IF
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr",/)' ) &
              trust_radius
          !
       END IF  
       !
       ! ... step along the bfgs direction
       !
       IF ( norm( bfgs_step ) < eps16 ) THEN
          !
          WRITE( UNIT = stdout, &
               & FMT = '(5X,"WARNING : norm( bfgs_step )",T30,"= ",F18.10)' ) &
              norm( bfgs_step )
          !
          bfgs_step = - gradient
          !
       ELSE
          !
          bfgs_step = trust_radius * bfgs_step / norm( bfgs_step )
          !
       END IF 
       !
       ! ... informations needed for the next iteration are saved
       ! ... this must be done before positions update
       !
       CALL write_lbfgs_file( pos, energy, gradient, scratch )                       
       !
       ! ... positions are updated for a new scf calculation
       !
       pos = pos + bfgs_step
       !
       DEALLOCATE( pos_old )
       DEALLOCATE( gradient_old )          
       DEALLOCATE( bfgs_step )              
       DEALLOCATE( bfgs_step_old )    
       !
     END SUBROUTINE lin_bfgs     
     !
     !
     ! ... private methods :
     !
     !-----------------------------------------------------------------------
     SUBROUTINE read_bfgs_file( pos, energy, gradient, scratch, dim )
       !-----------------------------------------------------------------------
       !
       USE io_files, ONLY : prefix
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP),     INTENT(INOUT) :: pos(:)
       REAL(KIND=DP),     INTENT(INOUT) :: energy       
       REAL(KIND=DP),     INTENT(INOUT) :: gradient(:)       
       CHARACTER (LEN=*), INTENT(IN)    :: scratch
       INTEGER,           INTENT(IN)    :: dim
       !
       ! ... local variables
       !
       CHARACTER (LEN=256) :: bfgs_file
       LOGICAL             :: file_exists
       !
       !
       bfgs_file = TRIM( scratch ) // TRIM( prefix ) //'.bfgs'
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
          READ( iunbfgs, * ) scf_iter
          READ( iunbfgs, * ) bfgs_iter
          READ( iunbfgs, * ) lin_iter
          READ( iunbfgs, * ) pos_old
          READ( iunbfgs, * ) energy_old
          READ( iunbfgs, * ) gradient_old
          READ( iunbfgs, * ) bfgs_step_old
          READ( iunbfgs, * ) trust_radius_old
          READ( iunbfgs, * ) inverse_hessian
          !     
          CLOSE( UNIT = iunbfgs )
          !
       ELSE
          !
          ! ... bfgs initialization
          !
          scf_iter                   = 0
          bfgs_iter                  = 0
          lin_iter                   = 0
          pos_old(:,lbfgs_ndim)      = pos
          energy_old                 = energy
          gradient_old(:,lbfgs_ndim) = gradient
          bfgs_step_old              = 0.D0
          trust_radius_old           = trust_radius_ini
          inverse_hessian            = identity(dim)
          !
       END IF    
       !
     END SUBROUTINE read_bfgs_file
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE read_lbfgs_file( pos, energy, gradient, scratch, dim )
       !-----------------------------------------------------------------------
       !
       USE io_files, ONLY : prefix
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP),     INTENT(INOUT) :: pos(:)
       REAL(KIND=DP),     INTENT(INOUT) :: energy       
       REAL(KIND=DP),     INTENT(INOUT) :: gradient(:)       
       CHARACTER (LEN=*), INTENT(IN)    :: scratch
       INTEGER,           INTENT(IN)    :: dim
       !
       ! ... local variables
       !
       CHARACTER (LEN=256) :: bfgs_file
       LOGICAL             :: file_exists
       !
       !
       bfgs_file = TRIM( scratch ) // TRIM( prefix ) //'.bfgs'
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
          READ( iunbfgs, * ) scf_iter
          READ( iunbfgs, * ) bfgs_iter
          READ( iunbfgs, * ) lin_iter
          READ( iunbfgs, * ) pos_old(:,1:lbfgs_ndim)
          READ( iunbfgs, * ) energy_old
          READ( iunbfgs, * ) gradient_old(:,1:lbfgs_ndim)
          READ( iunbfgs, * ) bfgs_step_old  
          READ( iunbfgs, * ) trust_radius_old
          !     
          CLOSE( UNIT = iunbfgs )
          !
       ELSE
          !
          ! ... bfgs initialization
          !
          scf_iter         = 0
          bfgs_iter        = 0
          lin_iter         = 0
          pos_old          = 0.D0
          energy_old       = energy
          gradient_old     = 0.D0
          trust_radius_old = trust_radius_ini
          bfgs_step_old    = 0.D0
          !
       END IF    
       !
     END SUBROUTINE read_lbfgs_file     
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_bfgs_file( pos, energy, gradient, scratch )
       !-----------------------------------------------------------------------
       !
       USE io_files, ONLY : prefix       
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP),     INTENT(IN) :: pos(:)       
       REAL(KIND=DP),     INTENT(IN) :: energy       
       REAL(KIND=DP),     INTENT(IN) :: gradient(:)       
       CHARACTER (LEN=*), INTENT(IN) :: scratch
       !
       !
       OPEN( UNIT = iunbfgs, FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs', &
             STATUS = 'UNKNOWN', ACTION = 'WRITE' )  
       !
       WRITE( iunbfgs, * ) scf_iter
       WRITE( iunbfgs, * ) bfgs_iter
       WRITE( iunbfgs, * ) lin_iter
       WRITE( iunbfgs, * ) pos
       WRITE( iunbfgs, * ) energy
       WRITE( iunbfgs, * ) gradient
       WRITE( iunbfgs, * ) bfgs_step
       WRITE( iunbfgs, * ) trust_radius
       WRITE( iunbfgs, * ) inverse_hessian
       ! 	     
       CLOSE( UNIT = iunbfgs )
       !
     END SUBROUTINE write_bfgs_file  
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_lbfgs_file( pos, energy, gradient, scratch )
       !-----------------------------------------------------------------------
       !
       USE io_files, ONLY : prefix       
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP),     INTENT(IN) :: pos(:)        
       REAL(KIND=DP),     INTENT(IN) :: energy       
       REAL(KIND=DP),     INTENT(IN) :: gradient(:)       
       CHARACTER (LEN=*), INTENT(IN) :: scratch
       !
       !
       OPEN( UNIT = iunbfgs, FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs', &
             STATUS = 'UNKNOWN', ACTION = 'WRITE' )  
       !
       WRITE( iunbfgs, * ) scf_iter
       WRITE( iunbfgs, * ) bfgs_iter
       WRITE( iunbfgs, * ) lin_iter
       WRITE( iunbfgs, * ) pos_old(:,2:lbfgs_ndim), pos
       WRITE( iunbfgs, * ) energy
       WRITE( iunbfgs, * ) gradient_old(:,2:lbfgs_ndim), gradient
       WRITE( iunbfgs, * ) bfgs_step
       WRITE( iunbfgs, * ) trust_radius
       ! 	     
       CLOSE( UNIT = iunbfgs )
       !
     END SUBROUTINE write_lbfgs_file          
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE update_inverse_hessian( gradient, dim, stdout )
       !-----------------------------------------------------------------------
       !
       USE constants,  ONLY : eps16
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(IN)  :: gradient(:)   
       INTEGER,       INTENT(IN)  :: dim
       INTEGER,       INTENT(IN)  :: stdout
       !
       ! ... local variables
       !
       REAL(KIND=DP)              :: y(dim)
       REAL(KIND=DP)              :: sdoty
       !
       !
       y = gradient - gradient_old(:,lbfgs_ndim)
       !
       sdoty = ( bfgs_step_old .dot. y )
       !
       IF ( ABS( sdoty ) < eps16 ) THEN
          !
          ! ... the history is reset
          !
          WRITE( stdout, '(/,5X,"WARINIG: unexpected behaviour in ", &
                              & "update_inverse_hessian")' )
          WRITE( stdout, '(5X,"         resetting bfgs history",/)' )
          !
          inverse_hessian = identity(dim)
          !
          RETURN
          !
       END IF 
       !
       inverse_hessian = inverse_hessian + &
              ( 1.D0 + ( y .dot. ( inverse_hessian .times. y ) ) / sdoty ) * &
              matrix( bfgs_step_old, bfgs_step_old ) / sdoty -               &
              ( matrix( bfgs_step_old, ( y .times. inverse_hessian ) ) +     &
                matrix( ( inverse_hessian .times. y ), bfgs_step_old ) ) / sdoty
       !
     END SUBROUTINE update_inverse_hessian
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE lbfgs_update( pos, gradient, dim )
       !-----------------------------------------------------------------------
       !
       USE constants,  ONLY : eps16
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(IN)  :: pos(:)
       REAL(KIND=DP), INTENT(IN)  :: gradient(:)         
       INTEGER,       INTENT(IN)  :: dim
       !
       ! ... local variables
       !
       INTEGER       :: i       
       REAL(KIND=DP) :: s(dim,lbfgs_ndim), y(dim,lbfgs_ndim)
       REAL(KIND=DP) :: alpha(lbfgs_ndim), sdoty(lbfgs_ndim)
       REAL(KIND=DP) :: preconditioning
       !
       !
       bfgs_step = gradient
       !
       s(:,lbfgs_ndim) = pos - pos_old(:,lbfgs_ndim) 
       y(:,lbfgs_ndim) = gradient - gradient_old(:,lbfgs_ndim) 
       !
       DO i = ( lbfgs_ndim - 1 ), 1, -1
          !
          s(:,i) = pos_old(:,i+1) - pos_old(:,i)
          y(:,i) = gradient_old(:,i+1) - gradient_old(:,i)
          !
       END DO
       !
       DO i = lbfgs_ndim, 1, -1
          !
          sdoty(i) = ( s(:,i) .dot. y(:,i) )
          !
          IF ( sdoty(i) > eps16 ) THEN
             !
             alpha(i) = ( s(:,i) .dot. bfgs_step ) / sdoty(i)
             !
          ELSE
             !   
             alpha(i) = 0.D0
             !
          END IF   
          !
          bfgs_step = bfgs_step - alpha(i) * y(:,i)
          !
       END DO 
       !
       preconditioning = ( y(:,lbfgs_ndim) .dot. y(:,lbfgs_ndim) )
       !
       IF ( preconditioning > eps16 ) THEN
          !
          bfgs_step =  sdoty(lbfgs_ndim) / preconditioning * bfgs_step
          !
       END IF        
       !
       DO i = 1, lbfgs_ndim
          !
          IF ( sdoty(i) > eps16 ) THEN
             !
             bfgs_step = bfgs_step + s(:,i) * ( alpha(i) - &
                        ( y(:,lbfgs_ndim) .dot. bfgs_step ) / sdoty(i) )
             !
          END IF   
          !
       END DO  
       !
       bfgs_step = - bfgs_step
       !
     END SUBROUTINE lbfgs_update  
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_wolfe_conditions( lwolfe, energy, gradient )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(IN)  :: energy       
       REAL(KIND=DP), INTENT(IN)  :: gradient(:)              
       LOGICAL,       INTENT(OUT) :: lwolfe
       !
       !
       lwolfe = ( energy - energy_old ) < & 
                w_1 * ( gradient_old(:,lbfgs_ndim) .dot. bfgs_step_old )
       !
       lwolfe = lwolfe .AND. &
                ( ABS( gradient .dot. bfgs_step_old ) > &
                  - w_2 * ( gradient_old(:,lbfgs_ndim) .dot. bfgs_step_old ) ) 
       !
     END SUBROUTINE check_wolfe_conditions
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_trust_radius( lwolfe, energy, gradient, dim, &
                                      stdout, conv_bfgs )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       LOGICAL,       INTENT(IN)  :: lwolfe
       REAL(KIND=DP), INTENT(IN)  :: energy    
       REAL(KIND=DP), INTENT(IN)  :: gradient(:)                         
       INTEGER,       INTENT(IN)  :: dim   
       INTEGER,       INTENT(IN)  :: stdout
       LOGICAL,       INTENT(OUT) :: conv_bfgs
       !
       ! ... local variables
       !
       REAL(KIND=DP) :: a
       LOGICAL       :: ltest
       !
       !
       ltest = ( energy - energy_old ) < &
               w_1 * ( gradient_old(:,lbfgs_ndim) .dot. bfgs_step_old )
       !       
       ltest = ltest .AND. ( norm( bfgs_step ) > trust_radius_old )
       !
       IF ( ltest ) THEN
          !
          a = 1.25D0
          !
       ELSE
          !
          a = 1.D0
          !
       END IF    
       !
       IF ( lwolfe ) THEN
          !
          trust_radius = MIN( trust_radius_max, 2.D0 * a * trust_radius_old )
          !
       ELSE
          !
          trust_radius = MIN( trust_radius_max, a * trust_radius_old, &
                              norm( bfgs_step ) )
          !
       END IF    
       !
       IF ( trust_radius < trust_radius_end  ) THEN
          !
          conv_bfgs = .TRUE.
          !
       ELSE IF ( trust_radius < trust_radius_min ) THEN
          !
          ! ... the history is reset
          !
          WRITE( UNIT = stdout, FMT = '(5X,"resetting bfgs history",/)' )
          !
          IF ( ALLOCATED( inverse_hessian ) ) THEN
             !
             inverse_hessian = identity(dim)
             !
             bfgs_step = - ( inverse_hessian .times. gradient )
             !
             trust_radius = trust_radius_min
             !
          ELSE
             !
             pos_old      = 0.D0
             gradient_old = 0.D0
             !     
             bfgs_step = - gradient
             !
             trust_radius = trust_radius_min
             !
          END IF      
          !
       END IF          
       !
     END SUBROUTINE compute_trust_radius 
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE terminate_bfgs( energy, stdout, scratch )
       !-----------------------------------------------------------------------
       !
       USE io_files, ONLY : prefix             
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP),     INTENT(IN) :: energy  
       INTEGER,           INTENT(IN) :: stdout         
       CHARACTER (LEN=*), INTENT(IN) :: scratch       
       !       
       !
       WRITE( UNIT = stdout, &
            & FMT = '(/,5X,"bfgs converged in ",I3," ionic steps and ", &
            &         I3," bfgs steps",/)' ) scf_iter, bfgs_iter
       WRITE( UNIT = stdout, &
            & FMT = '(5X,"Final energy",T30,"= ",F18.10," ryd")' ) energy
       !
       IF ( ALLOCATED( inverse_hessian ) ) THEN
          !
          WRITE( UNIT = stdout, &
               & FMT = '(/,5X,"Saving the approssimate hessian",/)' )
          !
          OPEN( UNIT = iunbfgs, FILE = TRIM( scratch ) // TRIM( prefix ) // &
              & '.hess', STATUS = 'UNKNOWN', ACTION = 'WRITE' )  
          !
          WRITE( iunbfgs, * ) SHAPE( inverse_hessian )
          WRITE( iunbfgs, * ) inverse_hessian
          ! 	     
          CLOSE( UNIT = iunbfgs )       
          !
          DEALLOCATE( pos_old )   
          DEALLOCATE( inverse_hessian )
          DEALLOCATE( bfgs_step )       
          DEALLOCATE( bfgs_step_old )
          DEALLOCATE( gradient_old ) 
          !
       ELSE
          !
          DEALLOCATE( pos_old )
          DEALLOCATE( gradient_old )          
          DEALLOCATE( bfgs_step )              
          DEALLOCATE( bfgs_step_old )  
          !
       END IF
       !    
       OPEN( UNIT = iunbfgs, &
             FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs' )
       CLOSE( UNIT = iunbfgs, STATUS = 'DELETE' )    
       !
     END SUBROUTINE terminate_bfgs
     !
END MODULE bfgs_module
