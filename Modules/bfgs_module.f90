!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE basic_algebra_routines
  !----------------------------------------------------------------------------
  !
  USE parameters,  ONLY : DP
  !
  IMPLICIT NONE
  !
  INTERFACE OPERATOR( .dot. )
     !
     MODULE PROCEDURE internal_dot_product
     !
  END INTERFACE
  !  
  INTERFACE OPERATOR( * )
     !
     MODULE PROCEDURE matrix_times_vector, vector_times_matrix
     !
  END INTERFACE
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     PURE FUNCTION internal_dot_product( vector1, vector2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector1(:), vector2(:)
       REAL (KIND=DP)             :: internal_dot_product
       !
       !
       internal_dot_product = DOT_PRODUCT( vector1 , vector2 )
       !
     END FUNCTION internal_dot_product
     !
     !     
     !----------------------------------------------------------------------- 
     PURE FUNCTION norm( vector )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector(:)
       REAL (KIND=DP)             :: norm
       !
       !
       norm = SQRT( vector .dot. vector )
       !
     END FUNCTION norm
     !
     !
     !-----------------------------------------------------------------------
     PURE FUNCTION matrix_times_vector( matrix , vector )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector(:)
       REAL (KIND=DP), INTENT(IN) :: matrix(:,:)
       REAL (KIND=DP)             :: matrix_times_vector(SIZE( vector ))
       INTEGER                    :: i, dim
       !
       !
       dim = SIZE( vector )
       !
       DO i = 1, dim
          !
          matrix_times_vector(i) = matrix(i,:) .dot. vector(:) 
          !
       END DO
       !
     END FUNCTION  matrix_times_vector
     !
     !
     !-----------------------------------------------------------------------
     PURE FUNCTION vector_times_matrix( vector , matrix )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector(:)
       REAL (KIND=DP), INTENT(IN) :: matrix(:,:)
       REAL (KIND=DP)             :: vector_times_matrix(SIZE( vector ))
       INTEGER                    :: i, dim
       !
       !
       dim = SIZE( vector )
       !
       DO i = 1, dim
          !
          vector_times_matrix(i) = vector(:) .dot. matrix(:,i)
          !
       END DO
       !
     END FUNCTION vector_times_matrix
     !
     !
     !-----------------------------------------------------------------------
     PURE FUNCTION matrix( vector1 , vector2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN) :: vector1(:), vector2(:)
       REAL (KIND=DP)             :: matrix(SIZE( vector1 ),SIZE( vector2 ))
       INTEGER                    :: i, j
       !
       !
       DO i = 1, SIZE( vector1 )
          !
          DO j = 1, SIZE( vector2 )
             !
             matrix(i,j) = vector1(i) * vector2(j)
             !
          END DO
          !
       END DO       
       !
     END FUNCTION matrix
     !
     !
     !-----------------------------------------------------------------------
     PURE FUNCTION identity( dim )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: dim
       REAL(KIND=DP)       :: identity(dim,dim)
       INTEGER             :: i
       !
       !
       identity = 0.D0
       !
       DO i = 1, dim
          !
          identity(i,i) = 1.D0
          !
       END DO
       !
     END FUNCTION identity
     !    
END MODULE basic_algebra_routines
!
!
!----------------------------------------------------------------------------
MODULE bfgs_module
  !----------------------------------------------------------------------------
  !
  ! ... ionic relaxation through Broyden-Fletcher-Goldfarb-Shanno 
  ! ... minimization and a "trust radius" line search based on 
  ! ... Wolfe conditions
  !
  ! ... references :  
  !
  ! ... 1) R. Fletcher, Practical Methods of Optimization, John Wiley and Sons,
  ! ...    Chichester, 2nd edn, 1987. 
  ! ... 2) Salomon R. Billeter, Alexander J. Turner, Walter Thiel, 
  ! ...    Phys. Chem. Chem. Phys. 2, 2177 (2000)
  ! ... 3) Salomon R. Billeter, Alessandro Curioni, Wanda Andreoni,
  ! ...    Comput. Mat. Science 27, 437, (2003)
  !
  !
  USE parameters,  ONLY : DP
  !
  USE basic_algebra_routines  
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: bfgs
  !
  REAL(KIND=DP), ALLOCATABLE :: &
      pos_old(:),               &! old positions
      inverse_hessian(:,:),     &! inverse of the hessian matrix (updated via
                                 ! BFGS formula)
      bfgs_step(:),             &! bfgs direction
      bfgs_step_old(:),         &! old bfgs direction
      gradient_old(:)            ! old gradient
  REAL(KIND=DP) :: &   
      trust_radius,             &! displacement along the bfgs direction
      trust_radius_old,         &! old displacement along the bfgs direction
      energy_old                 ! old energy
  INTEGER :: &
      iteration                  ! bfgs iteration    
  REAL(KIND=DP), PARAMETER  :: &
      trust_radius_max = 0.5D0, &! maximum allowed displacement
      trust_radius_min = 1.D-5, &! minimum allowed displacement
      trust_radius_ini = 0.5D0, &! initial displacement
      trust_radius_end = 1.D-7   ! bfgs stops when trust_radius is less than
                                 ! this value
  INTEGER, PARAMETER  :: &
      iunbfgs = 4                ! iunit for bfgs IO     
  !  
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE bfgs( pos, energy, gradient, scratch, stdout, energy_thr, &
                      gradient_thr, energy_error, gradient_error, &
                      step_accepted, conv_bfgs )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(INOUT)   :: pos(:)
       REAL(KIND=DP), INTENT(INOUT)   :: energy       
       REAL(KIND=DP), INTENT(INOUT)   :: gradient(:)
       CHARACTER (LEN=*), INTENT(IN)  :: scratch
       INTEGER, INTENT(IN)            :: stdout   
       REAL(KIND=DP), INTENT(IN)      :: energy_thr, gradient_thr  
       REAL(KIND=DP), INTENT(OUT)     :: energy_error, gradient_error       
       LOGICAL, INTENT(OUT)           :: step_accepted, conv_bfgs
       INTEGER                        :: dim, i
       LOGICAL                        :: lwolfe
       !
       !
       dim = SIZE( pos )
       !
       ALLOCATE( pos_old( dim ) )
       ALLOCATE( inverse_hessian( dim, dim ) )
       ALLOCATE( bfgs_step( dim ) )              
       ALLOCATE( bfgs_step_old( dim ) )
       ALLOCATE( gradient_old( dim ) )
       !
       CALL read_bfgs_file( scratch, dim )
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
       IF ( conv_bfgs ) THEN
          !
          CALL terminate_bfgs( energy, stdout, scratch )
          !
          RETURN
          !
       END IF
       !
       iteration = iteration + 1
       !
       WRITE( stdout, '(/,5X,"iteration  = ",I3)' )   iteration
       WRITE( stdout, '(5X,"energy new = ",F14.10)' ) energy
       WRITE( stdout, '(5X,"energy old = ",F14.10)' ) energy_old
       !
       IF ( energy > energy_old ) THEN
          !
          ! ... the previous step is rejected 
          !
          step_accepted = .FALSE.
          !
          WRITE( stdout, '(/,5X,"CASE: energy_new > energy_old",/)' )
          WRITE( stdout, '(/,5X,"gradient .dot. bfgs_step > 0.D0 : ",L1,/)' ) &
              ( ( gradient .dot. bfgs_step_old ) > 0.D0 )
          !
          ! ... the old trust radius is reduced by a factor 2
          !
          trust_radius = 0.5D0 * trust_radius_old
          !
          WRITE( stdout, '(5X,"trust_radius = ",F14.10)' ) trust_radius
          !
          IF ( trust_radius < trust_radius_min ) THEN
             !
             ! ... the history is resetted
             !     
             WRITE( stdout, '(/,5X,"resetting bfgs history",/)' )
             !
             inverse_hessian = identity(dim)
             !
             bfgs_step = - inverse_hessian * gradient
             !
             trust_radius = trust_radius_ini
             !
          ELSE 
             !
             ! ... values saved of the previous step are restored
             !
             pos       = pos_old
             energy    = energy_old
             gradient  = gradient_old
             !
             ! ... old bfgs direction (normalized) is recovered
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
          WRITE( stdout, '(/,5X,"CASE: energy_new < energy_old",/)' )
          !
          CALL check_wolfe_conditions( lwolfe, energy, gradient )
          !
          WRITE( stdout, '(5X,"lwolfe       = ",L1)' ) lwolfe
          !
          CALL update_inverse_hessian( gradient, dim )
          !
          ! ... bfgs direction (not normalized) 
          !
          bfgs_step = - inverse_hessian * gradient
          !
          IF ( ( gradient .dot. bfgs_step ) > 0.D0 ) THEN
             !  
             ! ... bfgs direction is reversed if not downhill
             !
             bfgs_step = - bfgs_step
             !
             WRITE( stdout, '(/,5X,"search direction reversed")' )
             !
          END IF   
          !  
          ! ... the new trust radius is computed
          !
          CALL compute_trust_radius( lwolfe, energy, gradient, dim, &
                                     stdout, conv_bfgs )
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
          WRITE( stdout, '(5X,"trust_radius = ",F14.10)' ) trust_radius
          !
       END IF  
       !
       ! ... actual positions are stored
       !
       pos_old = pos
       !
       ! ... bfgs step
       !
       bfgs_step = trust_radius * bfgs_step / norm( bfgs_step )
       !
       ! ... positions are updated
       !
       pos = pos + bfgs_step
       !
       CALL write_bfgs_file( energy, gradient, scratch )         
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
     SUBROUTINE read_bfgs_file( scratch, dim )
       !-----------------------------------------------------------------------
       !
       USE io_files, ONLY : prefix
       !
       IMPLICIT NONE
       !
       CHARACTER (LEN=*), INTENT(IN) :: scratch
       INTEGER, INTENT(IN)           :: dim
       CHARACTER (LEN=256)           :: bfgs_file
       LOGICAL                       :: file_exists
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
          READ( iunbfgs, * ) iteration
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
          iteration        = 0
          pos_old          = 0.D0
          energy_old       = 0.D0
          gradient_old     = 0.D0
          bfgs_step_old    = 0.D0
          trust_radius_old = trust_radius_ini
          inverse_hessian  = identity(dim)
          !
       END IF    
       !
     END SUBROUTINE read_bfgs_file
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE write_bfgs_file( energy, gradient, scratch )
       !-----------------------------------------------------------------------
       !
       USE io_files, ONLY : prefix       
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(IN)     :: energy       
       REAL(KIND=DP), INTENT(IN)     :: gradient(:)       
       CHARACTER (LEN=*), INTENT(IN) :: scratch
       !
       !
       OPEN( UNIT = iunbfgs, FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs', &
             STATUS = 'UNKNOWN', ACTION = 'WRITE' )  
       !
       WRITE( iunbfgs, * ) iteration
       WRITE( iunbfgs, * ) pos_old
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
     SUBROUTINE update_inverse_hessian( gradient, dim )
       !-----------------------------------------------------------------------
       !
       USE constants,  ONLY : eps16
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(IN)  :: gradient(:)   
       INTEGER, INTENT(IN)        :: dim                       
       REAL(KIND=DP)              :: gamma(dim)
       REAL(KIND=DP)              :: sdotgamma
       !
       !
       gamma = gradient - gradient_old
       !
       sdotgamma = bfgs_step_old .dot. gamma 
       !
       IF ( ABS( sdotgamma ) < eps16 ) THEN
          !
          inverse_hessian = inverse_hessian
          !
       ELSE
          !
          inverse_hessian = inverse_hessian + ( identity(dim) + &
            ( gamma .dot. ( inverse_hessian * gamma ) ) / sdotgamma ) * &
            matrix( bfgs_step_old, bfgs_step_old ) / sdotgamma - &
            ( matrix( bfgs_step_old, ( gamma * inverse_hessian ) ) + &
              matrix( ( inverse_hessian * gamma ), bfgs_step_old ) ) / sdotgamma
          !
       END IF  
       !
     END SUBROUTINE update_inverse_hessian
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_wolfe_conditions( lwolfe, energy, gradient )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(IN) :: energy       
       REAL(KIND=DP), INTENT(IN) :: gradient(:)              
       REAL(KIND=DP), PARAMETER  :: w_1 = 1.0D-4
       REAL(KIND=DP), PARAMETER  :: w_2 = 0.9D0
       !
       LOGICAL, INTENT(OUT)  :: lwolfe
       !
       lwolfe = ( energy - energy_old ) < & 
                w_1 * ( gradient_old .dot. bfgs_step_old )
       !
       lwolfe = lwolfe .AND. &
                ( ( gradient .dot. bfgs_step_old ) > &
                  w_2 * ( gradient_old .dot. bfgs_step_old ) )
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
       LOGICAL, INTENT(IN)           :: lwolfe
       REAL(KIND=DP), INTENT(IN)     :: energy    
       REAL(KIND=DP), INTENT(IN)     :: gradient(:)                         
       INTEGER, INTENT(IN)           :: dim   
       INTEGER, INTENT(IN)           :: stdout
       LOGICAL, INTENT(OUT)          :: conv_bfgs
       REAL(KIND=DP)                 :: a
       LOGICAL                       :: ltest
       !
       !
       ltest = ( energy - energy_old ) < &
               1.0D-4 * ( gradient_old .dot. bfgs_step_old )
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
          ! ... the history is resetted
          !
          WRITE( stdout, '(/,5X,"resetting bfgs history",/)' )
          !
          inverse_hessian = identity(dim)
          !
          bfgs_step = - inverse_hessian * gradient
          !
          trust_radius = trust_radius_ini
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
       REAL(KIND=DP), INTENT(IN)     :: energy  
       INTEGER, INTENT(IN)           :: stdout         
       CHARACTER (LEN=*), INTENT(IN) :: scratch       
       !       
       !
       WRITE( stdout, '(/,5X,"bfgs converged in ",I3," iterations")' ) &
           iteration
       WRITE( stdout, '(/,5X,"Final energy: ",F14.10," ryd"/)' ) energy
       !
       OPEN( UNIT = iunbfgs, &
             FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs' )
       CLOSE( UNIT = iunbfgs, STATUS = 'DELETE' )
       !
       DEALLOCATE( pos_old )   
       DEALLOCATE( inverse_hessian )
       DEALLOCATE( bfgs_step )       
       DEALLOCATE( bfgs_step_old )
       DEALLOCATE( gradient_old )     
       !
     END SUBROUTINE terminate_bfgs
     !
END MODULE bfgs_module
