!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE neb_routines
  !-----------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the NEB implementation into the PWSCF code
  ! ... Written by Carlo Sbraccia ( 04-11-2003 )
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE constants,  ONLY : AU, BOHR_RADIUS_ANGS, eV_to_kelvin
  !
  PRIVATE
  !
  PUBLIC :: initialize_neb, search_mep
  !   
  CONTAINS
    !
    !    
    !-----------------------------------------------------------------------
    SUBROUTINE initialize_neb()
      !-----------------------------------------------------------------------
      !
      USE control_flags,    ONLY : istep
      USE input_parameters, ONLY : pos, nat, restart_mode, calculation, &
                                   minimization_scheme, climbing
      USE io_files,         ONLY : prefix, iunneb, neb_file, &
                                   dat_file, int_file, xyz_file, axsf_file
      USE brilz,            ONLY : alat
      USE neb_variables,    ONLY : pos_      => pos, &
                                   climbing_ => climbing, &
                                   vel, num_of_images, dim, PES, PES_gradient, &
                                   elastic_gradient, tangent, grad, norm_grad, &
                                   error, mass, free_minimization, CI_scheme,  &
                                   optimization, k, k_min, k_max,  Emax_index, &
                                   VEC_scheme, ds, neb_thr, lquick_min ,       &
                                   ldamped_dyn, lmol_dyn
      USE neb_variables,    ONLY : neb_dyn_allocation   
      USE parser,           ONLY : int_to_char
      USE io_routines,      ONLY : read_restart
      USE formats,          ONLY : stringfmt   
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER                     :: i
      REAL (KIND=DP), ALLOCATABLE :: d_R(:)
      !
      ! ... end of local variables
      !
      !    
      ! ... NEB internal variables are set
      !
      ! ... output files are set
      !
      neb_file  = TRIM( prefix ) // ".neb"
      dat_file  = TRIM( prefix ) // ".dat"
      int_file  = TRIM( prefix ) // ".int"
      xyz_file  = TRIM( prefix ) // ".xyz"
      axsf_file = TRIM( prefix ) // ".axsf"
      !
      ! ... istep is initialized to zero
      !
      istep = 0
      !
      ! ... the dimension of all neb arrays is set 
      ! ... ( It corresponds to the dimension of the configurational space )
      !
      dim = 3 * nat
      !  
      ! ... dynamical allocation of arrays      
      !
      CALL neb_dyn_allocation()
      !
      ! ... coordinates must be in bohr ( pwscf uses alat units )
      !
      pos_ = pos(1:dim,:) * alat   
      !
      ! ... all other arrays are initialized 
      !
      PES              = 0.D0
      PES_gradient     = 0.D0
      elastic_gradient = 0.D0
      tangent          = 0.D0
      grad             = 0.D0
      norm_grad        = 0.D0
      error            = 0.D0
      mass             = 1.D0
      k                = k_min
      !
      IF ( ALLOCATED( climbing ) ) THEN
         !
         climbing_ = climbing(:)
         !
      ELSE
         !
         climbing_ = .FALSE.
         !
      END IF  
      !
      free_minimization = .FALSE.
      !
      IF ( optimization ) THEN      
         !
         free_minimization(1)             = .TRUE.
         free_minimization(num_of_images) = .TRUE.
         !
      END IF      
      !
      IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
         !     
         vel = 0.D0
         !
      END IF
      !
      ! ... initial path is read ( restart_mode == "restart" ) 
      ! ... or generated ( restart_mode = "from_scratch" )
      !
      IF ( restart_mode == "restart" ) THEN
         !
         CALL read_restart()
         !
         IF ( CI_scheme == "highest-TS" ) THEN
            !
            climbing_ = .FALSE.
            !
            climbing_(Emax_index) = .TRUE.
            !
         ELSE IF ( CI_scheme == "all-SP" ) THEN
            !
            CALL search_stationary_points()
            !
         END IF
         !
      ELSE
         !
         ALLOCATE( d_R(dim) )        
         !
         d_R = ( pos_(:,num_of_images) - pos_(:,1) )
         ! 
         d_R = d_R / REAL( ( num_of_images - 1 ), KIND = DP )
         !
         DO i = 2, ( num_of_images - 1 )
            !
            pos_(:,i) = pos_(:,( i - 1 )) + d_R(:)
            !
         END DO
         !
         DEALLOCATE( d_R )
         !
      END IF    
      !
      CALL compute_deg_of_freedom()
      !
      ! ... details of the calculation are written on output
      !
      WRITE( UNIT = iunneb, FMT = stringfmt ) &
          "calculation", TRIM( calculation )
      WRITE( UNIT = iunneb, FMT = stringfmt ) &
          "restart_mode", TRIM( restart_mode )
      WRITE( UNIT = iunneb, FMT = stringfmt ) &
          "CI_scheme", TRIM( CI_scheme )
      WRITE( UNIT = iunneb, FMT = stringfmt ) &
          "VEC_scheme", TRIM( VEC_scheme )
      WRITE( UNIT = iunneb, FMT = stringfmt ) &
          "minimization_scheme", TRIM( minimization_scheme )
      WRITE( UNIT = iunneb, &
             FMT = '(5X,"optimization",T35," = ",L1))' ) optimization
      WRITE( UNIT = iunneb, &
             FMT = '(5X,"num_of_images",T35," = ",I3)' ) num_of_images
      WRITE( UNIT = iunneb, &
             FMT = '(5X,"ds",T35," = ",F6.4)' ) ds
      WRITE( UNIT = iunneb, &
             FMT = '(5X,"k_max",T35," = ",F6.4)' ) k_max
      WRITE( UNIT = iunneb, &
             FMT = '(5X,"k_min",T35," = ",F6.4)' ) k_min
      WRITE( UNIT = iunneb, &
             FMT = '(5X,"neb_thr",T35," = ",F6.4)' ) neb_thr     
      WRITE( UNIT = iunneb, FMT = '(/)' )
      !
      RETURN
      !
      CONTAINS
         !
         SUBROUTINE compute_deg_of_freedom 
           !
           USE basis,            ONLY :  nat
           USE input_parameters, ONLY :  if_pos
           USE neb_variables,    ONLY :  deg_of_freedom
           !
           IMPLICIT NONE
           !
           INTEGER    :: ia
           !
           !
           deg_of_freedom = 0
           !
           DO ia = 1, nat
              !
              IF ( if_pos(1,ia) == 1 ) deg_of_freedom = deg_of_freedom + 1
              IF ( if_pos(2,ia) == 1 ) deg_of_freedom = deg_of_freedom + 1
              IF ( if_pos(3,ia) == 1 ) deg_of_freedom = deg_of_freedom + 1
              !
           END DO
           !
         END SUBROUTINE compute_deg_of_freedom
         !
    END SUBROUTINE initialize_neb
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE search_mep()
      !-----------------------------------------------------------------------
      !
      USE control_flags, ONLY : time_max, istep, nstep
      USE io_files,      ONLY : iunneb, iunexit, exit_file     
      USE formats,       ONLY : run_output, run_output_T_const
      USE neb_variables, ONLY : num_of_images, dim, pos, PES, error,       &
                                climbing, optimization,  CI_scheme,        &
                                Emax_index, temp, Emax, neb_thr, conv_neb, &
                                suspended_image, lsteep_des, lquick_min ,  &
                                ldamped_dyn, lmol_dyn 
      USE io_routines,   ONLY : write_restart, write_dat_files, write_output 
#if defined (__PARA)
      USE para,          ONLY : me, mypool
      USE mp,            ONLY : mp_barrier      
#endif
      USE minimization_routines
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      REAL (KIND=DP)      :: err, tcpu, langevin_action
      INTEGER             :: image
      LOGICAL             :: stat
      INTEGER             :: N_in, N_fin
      LOGICAL             :: file_exists
      !
      ! ... end of local variables
      ! ... external functions
      !
      REAL (kind=DP), EXTERNAL :: get_clock
      !
      !
      conv_neb = .FALSE.
      !
      IF ( istep == nstep ) THEN
         !
         CALL compute_tangent()
         ! 
         CALL gradient()
         !
         CALL compute_error( err )
         !
         CALL write_dat_files()
         !
         suspended_image = 0
         !
         RETURN
         !  
      END IF 
      ! 
      IF ( optimization ) THEN
         !
         N_in  = 1
         N_fin = num_of_images
         !     
      ELSE
         !
         N_in  = 2
         N_fin = num_of_images - 1
         ! 
      END IF
      !
      minimization: DO
         !
         ! ... the restart file is written 
         !
         CALL write_restart()
         !
         ! ... minimization step is done only in case of no suspended images
         ! ... when the simulation is started from scratch all gradients are 
         ! ... zero and the elastic band remains in the old position.
         !
         IF ( suspended_image == 0 ) THEN
            !
            first_minimization_loop: DO image = N_in, N_fin
               !
               IF ( lsteep_des ) THEN
                  !
                  CALL steepest_descent( image )
                  !
               ELSE IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
                  !
                  CALL velocity_Verlet_first_step( image )
                  !
               END IF
               !
            END DO first_minimization_loop 
            !
         END IF
         !
         ! ... the programs checks if the user has required a soft
         ! ... exit or if if maximum CPU time has been exceeded
         !
         INQUIRE( FILE = TRIM( exit_file ), EXIST = file_exists )
         !
         tcpu = get_clock( 'PWSCF' )
         !
         IF ( file_exists ) THEN
            !
#if defined (__PARA)
            !
            ! ... all jobs are syncronized
            !
            CALL mp_barrier()
            !
            IF ( me == 1 .AND. mypool == 1 ) THEN
#endif
               OPEN( UNIT = iunexit, FILE = TRIM( exit_file ), STATUS = "OLD" )
               CLOSE( UNIT = iunexit, STATUS = "DELETE" )
#if defined (__PARA)
            END IF
#endif       
            !  
            WRITE( UNIT = iunneb, &
                   FMT = '(/,5X,"WARNING :  soft exit required",/, &
                   & 5X,"stopping in searc_mep() ...",/)' )  
            !  
            CALL stop_pw( .FALSE. )
            !
         ELSE IF ( tcpu > time_max ) THEN
            ! 
            WRITE( UNIT = iunneb, &
                   FMT = '(5X,"Maximum CPU time exceeded",2F15.2)') &
                tcpu, time_max
            !  
            CALL stop_pw( .FALSE. )
            !
         END IF
         !
         ! ... energies and gradients acting on each image of the elastic band 
         ! ... are computed calling a driver that perfoms the scf calculations
         !
         IF ( istep == 0 ) THEN
            !
            CALL born_oppenhimer_PES( .TRUE., stat )
            !
         ELSE
            !
            CALL born_oppenhimer_PES( optimization, stat )
            !
         END IF
         !
         IF ( .NOT. stat ) THEN
            !
            EXIT minimization
            !
         END IF
         !
         IF ( CI_scheme == "highest-TS" ) THEN
            ! 
            climbing = .FALSE.
            !  
            climbing(Emax_index) = .TRUE.
            !
         ELSE IF ( CI_scheme == "all-SP" ) THEN
            !
            CALL search_stationary_points()
            !
         END IF
         !
         CALL compute_tangent()
         CALL gradient()
         ! 
         ! ... a second minimization step is needed for those algorithms 
         ! ... based on a velocity Verlet scheme
         !
         IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            second_minimization_loop: DO image = N_in, N_fin 
               !
               IF ( lquick_min ) THEN
                  !
                  CALL quick_min_second_step( image )
                  !
               ELSE IF ( ldamped_dyn .OR. lmol_dyn ) THEN
                  !    
                  CALL velocity_Verlet_second_step( image )
                  !
               END IF
               !
            END DO second_minimization_loop
            !
            IF ( lmol_dyn ) CALL thermalization( N_in , N_fin   )
            !
         END IF
         !
         CALL compute_error( err )
         !
         CALL write_dat_files()
         !
         CALL compute_action( langevin_action )
         !
         istep = istep + 1
         !
         ! ... informations are written on the standard output
         !
         IF ( lmol_dyn ) THEN
            !
            WRITE( UNIT = iunneb, FMT = '(/,"langevin action = ",F16.8,/)' ) &
                langevin_action
            !
            WRITE( UNIT = iunneb, FMT = run_output_T_const ) &
                istep,                    &
                temp * AU * eV_to_kelvin, &
                err * ( AU / BOHR_RADIUS_ANGS ) 
            !
         ELSE
            !
            WRITE( UNIT = iunneb, FMT = '(/,"langevin action = ",F16.8,/)' ) &
                langevin_action
            !            
            WRITE( UNIT = iunneb, FMT = run_output ) &
                istep,                  &
                ( Emax - PES(1) ) * AU, &
                err * ( AU / BOHR_RADIUS_ANGS )
            !   
         END IF
         !
         CALL write_output()
         !
         ! ... the program checks if the convergence has been achieved
         ! 
         IF ( ( err * AU / BOHR_RADIUS_ANGS ) <= neb_thr )  THEN
            !
            conv_neb = .TRUE.
            !
            EXIT minimization
            !  
         END IF
         !
         ! ... the programs checks if the maximum number of iterations has 
         ! ... been reached or if the user has required a soft exit
         !
         IF ( istep >= nstep ) THEN
            !
            suspended_image = 0
            !
            EXIT minimization
            !  
         END IF 
         !
         suspended_image = 0
         ! 
      END DO minimization
      !
      RETURN
      !
    END SUBROUTINE search_mep
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE compute_action( langevin_action )
      !------------------------------------------------------------------------
      !
      USE neb_variables,          ONLY : num_of_images, pos, PES_gradient
      USE basic_algebra_routines, ONLY : norm
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      REAL (KIND=DP), INTENT(OUT) :: langevin_action
      !
      ! ... local variables
      !
      INTEGER          :: image
      !
      ! ... end of local variables
      !       
      !  
      langevin_action = 0.D0
      !
      DO image = 1, ( num_of_images - 1 )  
         !  
	 langevin_action = langevin_action + &
                           norm( pos(:,(image+1)) - pos(:,image) ) * &
                           ( norm( PES_gradient(:,image+1) ) + &
                             norm( PES_gradient(:,image) ) )
         !
      END DO
      !
    END SUBROUTINE compute_action
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE compute_tangent()
      !------------------------------------------------------------------------
      !
      USE neb_variables,          ONLY : num_of_images, tangent
      USE basic_algebra_routines, ONLY : norm
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER :: image   
      !
      ! ... end of local variables
      !
      !
      tangent = 0
      !
      DO image = 2, ( num_of_images - 1 )
         !
         ! ... tangent to the path ( normalized )
         !
         !tangent(:,image) = path_tangent( image )
         CALL path_tangent_( image, tangent(:,image) )
         !
         tangent(:,image) = tangent(:,image) / norm( tangent(:,image) )
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE compute_tangent
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE elastic_constants()
      !------------------------------------------------------------------------
      ! 
      USE constants,              ONLY : pi, eps32
      USE neb_variables,          ONLY : num_of_images, Emax, Emin, &
                                         k_max, k_min, k, PES, PES_gradient, &
                                         VEC_scheme
      USE basic_algebra_routines, ONLY : norm
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER        :: image      
      REAL (KIND=DP) :: delta_E
      REAL (KIND=DP) :: norm_grad_V, norm_grad_V_min, norm_grad_V_max
      !
      ! ... end of local variables
      !
      !
      IF ( VEC_scheme == "energy-weighted" ) THEN
         !
         delta_E = Emax - Emin     
         !
         IF ( delta_E <= eps32 ) THEN
            !
            k = k_min
            !
            RETURN
            !
         END IF
         !
         DO image = 1, num_of_images 
            !
            k(image) = 0.25D0 * ( ( k_max + k_min ) -  ( k_max - k_min ) * &
                       COS( pi * ( PES(image) - Emin ) / delta_E ) )
            !
         END DO
         !
      ELSE
         !
         norm_grad_V_min = + 1.0D32
         norm_grad_V_max = - 1.0D32
         !
         DO image = 1, num_of_images 
            !  
            norm_grad_V = norm( PES_gradient(:,image) )
            !
            IF ( norm_grad_V < norm_grad_V_min ) norm_grad_V_min = norm_grad_V
            IF ( norm_grad_V > norm_grad_V_max ) norm_grad_V_max = norm_grad_V
            !
         END DO   
         !    
         DO image = 1, num_of_images 
            !
            norm_grad_V = norm( PES_gradient(:,image) )
            !
            k(image) = 0.25D0 * ( ( k_max + k_min ) - ( k_max - k_min ) * &
                       COS( pi * ( norm_grad_V - norm_grad_V_min ) / & 
                       ( norm_grad_V_max - norm_grad_V_min ) ) )
            !
         END DO      
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE elastic_constants
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE gradient()
      !------------------------------------------------------------------------
      !
      USE supercell,              ONLY : pbc
      USE neb_variables,          ONLY : pos, grad, norm_grad,              &
                                         elastic_gradient, PES_gradient, k, &
                                         num_of_images, free_minimization,  &
                                         climbing, tangent, lmol_dyn
      USE basic_algebra_routines, ONLY : norm   
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER :: i
      !
      ! ... end of local variables
      !
      CALL elastic_constants
      !
      gradient_loop: DO i = 1, num_of_images
         !
         IF ( ( i > 1 ) .AND. ( i < num_of_images ) ) THEN
            !
            IF ( lmol_dyn ) THEN
               !
               ! ... elastic gradient ( variable elastic consatnt is used )
               ! ... note that this is NOT the NEB recipe
               !
               elastic_gradient = &
                       ( ( k(i) + k(i-1) ) * pbc( pos(:,i) - pos(:,(i-1)) ) - &
                         ( k(i) + k(i+1) ) * pbc( pos(:,(i+1)) - pos(:,i) ) )
               !
            ELSE
               !
               ! ... elastic gradient only along the path ( variable elastic
               ! ... consatnt is used ) NEB recipe
               !
               elastic_gradient = tangent(:,i) * &
               ( ( k(i) + k(i-1) ) * norm( pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
                 ( k(i) + k(i+1) ) * norm( pbc( pos(:,(i+1)) - pos(:,i) ) ) )
               !
            END IF
            !
         END IF
         !
         ! ... total gradient on each image ( climbing image is used if needed )
         ! ... only the component of the PES gradient orthogonal to the path is 
         ! ... taken into account
         !
         grad(:,i) = PES_gradient(:,i)
         !
         IF ( climbing(i) ) THEN
            !
            grad(:,i) = grad(:,i) - 2.D0 * tangent(:,i) * &
                        DOT_PRODUCT( PES_gradient(:,i) , tangent(:,i) ) 
            ! 
         ELSE IF ( ( .NOT. free_minimization(i) ) .AND. &
                   ( i > 1 ) .AND. ( i < num_of_images ) ) THEN
            !
            grad(:,i) = grad(:,i) + elastic_gradient - tangent(:,i) * &
                        DOT_PRODUCT( PES_gradient(:,i) , tangent(:,i) )
            !      
         END IF
         ! 
         norm_grad(i) = norm( grad(:,i) )
         !  
      END DO gradient_loop
      !
      RETURN
      !
    END SUBROUTINE gradient
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE search_stationary_points()
      !-----------------------------------------------------------------------
      !
      USE neb_variables, ONLY : num_of_images, PES, climbing, &
                                free_minimization, optimization 
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER         :: i
      REAL (KIND=DP)  :: V_p, V_h, V_n
      !
      ! ... end of local variables
      !
      !
      climbing          = .FALSE.
      free_minimization = .FALSE.
      !
      IF ( optimization ) THEN
         !
         free_minimization(1)             = .TRUE.
         free_minimization(num_of_images) = .TRUE.
         !
      END IF    
      !
      V_p = PES(1)
      V_h = PES(2)
      !
      DO i = 2, ( num_of_images - 1 )
         !
         V_n = PES(i+1)
         !
         IF ( ( V_h > V_p ) .AND. &
              ( V_h > V_n ) ) climbing(i) = .TRUE.
         IF ( ( V_h < V_p ) .AND. &
              ( V_h < V_n ) ) free_minimization(i) = .TRUE.
         !
         V_p = V_h
         V_h = V_n
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE search_stationary_points
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_error( err )
      !-----------------------------------------------------------------------
      !
      USE neb_variables, ONLY : num_of_images, optimization, &
                                error, norm_grad
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      REAL (KIND=DP), INTENT(OUT)  :: err
      !
      ! ... local variables
      !
      INTEGER                      :: N_in, N_fin
      INTEGER                      :: i
      !
      ! ... end of local variables
      !
      !
      err = 0.D0
      !
      IF ( optimization ) THEN
         !
         N_in  = 1
         N_fin = num_of_images
         !
      ELSE
         !
         N_in  = 2
         N_fin = ( num_of_images - 1 )      
         !   
      END IF   
      !
      DO i = N_in, N_fin
         !
         ! ... the error is given by the norm of the 
         ! ... gradient ( PES + SPRINGS ).
         !
         error(i) = norm_grad(i)
         !
         IF ( error(i) > err ) err = error(i)
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE compute_error
    !
    !
    !-----------------------------------------------------------------------
    !FUNCTION path_tangent( index )
    SUBROUTINE path_tangent_( index, path_tangent )
      !-----------------------------------------------------------------------
      !
      USE supercell,     ONLY : pbc
      USE neb_variables, ONLY : pos, dim, PES
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER, INTENT(IN) :: index  
      REAL (KIND=DP)      :: path_tangent(dim)                
      !
      ! ... local variables
      !
      REAL (KIND=DP) :: V_previous, V_actual, V_next
      REAL (KIND=DP) :: abs_next, abs_previous
      REAL (KIND=DP) :: delta_V_max, delta_V_min
      !
      ! ... end of local variables
      !
      !
      V_previous = PES( index - 1 )
      V_actual   = PES( index )
      V_next     = PES( index + 1 )
      !
      IF ( ( V_next > V_actual ) .AND. ( V_actual > V_previous ) ) THEN
         !
         path_tangent = pbc( pos(:,( index + 1 )) - pos(:,index) )
         !
      ELSE IF ( ( V_next < V_actual ) .AND. ( V_actual < V_previous ) ) THEN
         !
         path_tangent = pbc( pos(:,index) - pos(:,( index - 1 )) )
         !
      ELSE
         !
         abs_next     = ABS( V_next - V_actual ) 
         abs_previous = ABS( V_previous - V_actual ) 
         !
         delta_V_max = MAX( abs_next , abs_previous ) 
         delta_V_min = MIN( abs_next , abs_previous )
         !
         IF ( V_next > V_previous ) THEN
            !
            path_tangent = &
                pbc( pos(:,( index + 1 )) - pos(:,index) ) * delta_V_max + & 
                pbc( pos(:,index) - pos(:,( index - 1 )) ) * delta_V_min
            !
         ELSE IF ( V_next < V_previous ) THEN
            !
            path_tangent = &
                pbc( pos(:,( index + 1 )) - pos(:,index) ) * delta_V_min + & 
                pbc( pos(:,index) - pos(:,( index - 1 )) ) * delta_V_max
            !
         ELSE
            !
            path_tangent = &
                pbc( pos(:,( index + 1 )) - pos(:,( index - 1 )) ) 
            !
         END IF
         !
      END IF 
      !
      RETURN
      !
    !END FUNCTION path_tangent
    END SUBROUTINE path_tangent_
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE born_oppenhimer_PES( flag, stat )
      !-----------------------------------------------------------------------
      !
      USE neb_variables, ONLY : num_of_images, Emax_index, Emin, Emax, &
                                PES, PES_gradient, suspended_image
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      LOGICAL, INTENT(IN)   :: flag
      LOGICAL, INTENT(OUT)  :: stat            
      !
      ! ... local variables
      !
      INTEGER               :: i, image
      INTEGER               :: N_in, N_fin
      ! 
      ! ... end of local variables
      !
      !
      IF ( flag ) THEN
         !
         N_in = 1
         N_fin = num_of_images
         !
      ELSE
         !
         N_in = 2
         N_fin = ( num_of_images - 1 )
         !
      END IF
      !
      IF ( suspended_image /= 0 ) N_in = suspended_image
      !
      Emin = + 1.0D16
      Emax = - 1.0D16
      !
      CALL compute_scf( N_in, N_fin, stat )  
      !
      IF ( .NOT. stat ) RETURN
      !
      DO image = 1, num_of_images
         !  
         IF ( PES(image) <= Emin ) Emin = PES(image)
         !
         IF ( PES(image) >= Emax ) THEN
            !
            Emax = PES(image)         
            !
            Emax_index = image
            !
         END IF
         !
      END DO 
      !
      RETURN
      !
    END SUBROUTINE born_oppenhimer_PES
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_scf( N_in, N_fin, stat  )
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : if_pos, sp_pos, startingwfc, startingpot
      USE constants,        ONLY : e2
      USE control_flags,    ONLY : time_max, conv_elec, ethr
      USE brilz,            ONLY : alat
      USE basis,            ONLY : tau, ityp, nat, &
                                   startingwfc_ => startingwfc, &
                                   startingpot_ => startingpot    
      USE ener,             ONLY : etot
      USE force_mod,        ONLY : force
      USE relax,            ONLY : if_pos_ => if_pos
      USE extfield,         ONLY : tefield, forcefield
      USE io_files,         ONLY : prefix, tmp_dir, &
                                   iunneb, iunexit, exit_file
      USE formats,          ONLY : scf_fmt
      USE neb_variables,    ONLY : pos, PES, PES_gradient, num_of_images, &
                                   dim, suspended_image
      USE parser,           ONLY : int_to_char
#if defined (__PARA)
      USE para,             ONLY : me, mypool
      USE mp,               ONLY : mp_barrier
#endif        
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER, INTENT(IN)  :: N_in, N_fin
      LOGICAL, INTENT(OUT) :: stat      
      !
      ! ... local variables definition
      !
      INTEGER              :: image, ia
      REAL (KIND=DP)       :: tcpu 
      CHARACTER (LEN=80)   :: tmp_dir_saved
      LOGICAL              :: file_exists, opnd 
      !
      ! ... end of local variables definition
      !
      ! ... external functions definition
      !
      REAL (KIND=DP), EXTERNAL :: get_clock
      !
      ! ... end of external functions definition
      !
      !
      stat = .TRUE.
      !
      tmp_dir_saved = tmp_dir
      ! 
      DO image = N_in, N_fin
         !
         suspended_image = image
         !
         INQUIRE( FILE = TRIM( exit_file ), EXIST = file_exists )
         !
         tcpu = get_clock( 'PWSCF' )
         !
         IF ( file_exists ) THEN
            !
#if defined (__PARA)
            !
            ! ... all jobs are syncronized
            !
            CALL mp_barrier()
            !
            IF ( me == 1 .AND. mypool == 1 ) THEN
#endif    
               OPEN( UNIT = iunexit, FILE = TRIM( exit_file ), STATUS = "OLD" )
               CLOSE( UNIT = iunexit, STATUS = "DELETE" )
#if defined (__PARA)
            END IF
#endif       
            !  
            WRITE( iunneb, '(/,5X,"WARNING :  soft exit required",/, &
                             & 5X,"stopping in compute_scf()...",/)' )
            !   
            stat = .FALSE.
            !
            RETURN
            !    
         ELSE IF ( tcpu > time_max ) THEN
            !
            WRITE( iunneb, '(5X,"Maximum CPU time exceeded",2F15.2)' ) &
                tcpu, time_max
            !
            stat = .FALSE.
            !
            RETURN
            !       
         END IF
         !
         tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                   TRIM( int_to_char( image ) ) // "/" 
         !
         WRITE( UNIT = iunneb, FMT = scf_fmt ) tcpu, image
         !
         CALL clean_pw()
         !
         CALL close_files()
         !
         ! ... unit stdout is connected to the appropriate file
         !
#if defined (__PARA)
         IF ( me == 1 .AND. mypool == 1 ) THEN
#endif  
            INQUIRE( UNIT = stdout, OPENED = opnd )
            IF ( opnd ) CLOSE( UNIT = stdout )
            OPEN( UNIT = stdout, FILE = TRIM( tmp_dir )//'PW.out', &
                  STATUS = 'UNKNOWN', POSITION = 'APPEND' )
#if defined (__PARA)
         END IF 
#endif 
         !
         IF ( .NOT. ALLOCATED( tau ) )      ALLOCATE( tau( 3, nat ) )
         IF ( .NOT. ALLOCATED( ityp ) )     ALLOCATE( ityp( nat ) )
         IF ( .NOT. ALLOCATED( force ) )    ALLOCATE( force( 3, nat ) )  
         IF ( .NOT. ALLOCATED( if_pos_ ) )  ALLOCATE( if_pos_( 3, nat ) )
         IF ( tefield  .AND. .NOT. ALLOCATED( forcefield ) ) &
                                           ALLOCATE( forcefield( 3, nat ) )
         !
         DO ia = 1, nat
            !
            ! ... tau is in alat units ( pos is in bohr )
            !
            tau(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),image) / alat
            !
         END DO
         !
         if_pos_(:,:) = if_pos(:,1:nat) 
         ityp(:)      = sp_pos(1:nat)
         !
         CALL init_run()
         !
         ! ... potential and wavefunctions are extrapolated
         !
         CALL update_pot()
         !
         CALL electrons()
         !
         IF ( .NOT. conv_elec ) THEN
            !
            WRITE( iunneb, '(/,5X,"WARNING :  scf convergence NOT achieved",/, &
                             & 5X,"stopping in compute_scf()...",/)' )
            !   
            stat = .FALSE.
            !
            RETURN
            !
         END IF   
         !
         CALL forces()
         !
         ! ... energy is converted from rydberg to hartree
         !
         PES(image) = etot / e2
         !
         ! ... gradients are converted from ( rydberg / bohr ) 
         ! ... to ( hartree / bohr )
         !
         PES_gradient(:,image) = - RESHAPE( SOURCE = force, &
                                            SHAPE = (/ dim /) ) / e2 
         !
         ! ... input values are restored at the end of each iteration
         !
         ethr = 0.D0
         startingpot_ = startingpot
         startingwfc_ = startingwfc
         !
         CALL reset_k_points()
         !
      END DO         
      !
      tmp_dir = tmp_dir_saved
      !
      suspended_image = 0
      !
      RETURN
      !
    END SUBROUTINE compute_scf  
    !
END MODULE neb_routines
