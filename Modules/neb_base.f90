!
! Copyright (C) 2003-2004 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE neb_base
  !-----------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the NEB implementation into the PWSCF-FPMD-CPV codes
  ! ... Written by Carlo Sbraccia ( 04-11-2003 )
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE constants,  ONLY : AU, BOHR_RADIUS_ANGS, eV_to_kelvin
  !
  PRIVATE
  !
  PUBLIC :: initialize_neb, compute_action, compute_tangent
  PUBLIC :: elastic_constants, gradient, search_stationary_points
  PUBLIC :: compute_error, path_tangent_
  PUBLIC :: born_oppenheimer_PES, search_mep
  !   
  CONTAINS
    !
    ! ... module procedures    
    !
    !-----------------------------------------------------------------------
    SUBROUTINE initialize_neb( prog )
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : pos, nat, restart_mode, calculation, &
                                   minimization_scheme, climbing, nstep
      USE io_files,         ONLY : prefix, iunneb, neb_file, &
                                   dat_file, int_file, xyz_file, axsf_file
      USE cell_base,        ONLY : alat
      USE neb_variables,    ONLY : pos_      => pos, &
                                   climbing_ => climbing, &
                                   vel, num_of_images, dim, PES, PES_gradient, &
                                   elastic_gradient, tangent, grad, norm_grad, &
                                   error, mass, free_minimization, CI_scheme,  &
                                   optimization, k, k_min, k_max,  Emax_index, &
                                   VEC_scheme, ds, neb_thr, lquick_min ,       &
                                   ldamped_dyn, lmol_dyn, nstep_neb, istep_neb
      USE neb_variables,    ONLY : neb_dyn_allocation   
      USE parser,           ONLY : int_to_char
      USE io_routines,      ONLY : read_restart
      USE formats,          ONLY : stringfmt   
      USE io_global,        ONLY : ionode
      !
      IMPLICIT NONE

      CHARACTER(LEN=2) :: prog   ! ... specify the calling program

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
      istep_neb = 0
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
      IF( prog == 'PW' ) THEN
        !
        pos_ = pos(1:dim,:) * alat   
        !
      ELSE
        !
        pos_ = pos(1:dim,:)
        !
      END IF
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
         ! ... consistency between the input value of nstep and the value
         ! ... of nstep_neb read from the restart_file is checked
         !
         IF ( nstep == 0 ) THEN
            !
            istep_neb = 0
            nstep_neb = nstep
            !
         END IF   
         !
         IF ( nstep > nstep_neb ) nstep_neb = nstep
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
         ! ... linear interpolation
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
      ! ... the actual number of degrees of freedom is computed
      !
      CALL compute_deg_of_freedom()
      !
      ! ... details of the calculation are written on output (only by ionode)
      !
      IF ( ionode ) THEN
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
      END IF
      !
      RETURN
      !
      CONTAINS
         !
         SUBROUTINE compute_deg_of_freedom()
           !
           USE ions_base,        ONLY :  nat
           USE input_parameters, ONLY :  if_pos
           USE neb_variables,    ONLY :  deg_of_freedom
           !
           IMPLICIT NONE
           !
           INTEGER :: ia
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
         !!! tangent(:,image) = path_tangent( image )
         !!! workaround for ifc8 compiler internal error
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
      USE neb_variables,          ONLY : pos, num_of_images, Emax, Emin, &
                                         k_max, k_min, k, PES, PES_gradient, &
                                         VEC_scheme, elastic_gradient, tangent
      USE supercell,              ONLY : pbc
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER        :: i
      REAL (KIND=DP) :: F_ortho_max, F_ortho_max_i, &
                        F_para_max_i, F_para_max
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
         DO i = 1, num_of_images 
            !
            k(i) = 0.25D0 * ( ( k_max + k_min ) -  ( k_max - k_min ) * &
                   COS( pi * ( PES(i) - Emin ) / delta_E ) )
            !
         END DO
         !
      ELSE
         !
         norm_grad_V_min = + 1.0D32
         norm_grad_V_max = - 1.0D32
         !
         DO i = 1, num_of_images 
            !  
            norm_grad_V = norm( PES_gradient(:,i) )
            !
            IF ( norm_grad_V < norm_grad_V_min ) norm_grad_V_min = norm_grad_V
            IF ( norm_grad_V > norm_grad_V_max ) norm_grad_V_max = norm_grad_V
            !
         END DO   
         !    
         DO i = 1, num_of_images 
            !
            norm_grad_V = norm( PES_gradient(:,i) )
            !
            k(i) = 0.25D0 * ( ( k_max + k_min ) - ( k_max - k_min ) * &
                   COS( pi * ( norm_grad_V - norm_grad_V_min ) / & 
                   ( norm_grad_V_max - norm_grad_V_min ) ) )
            !
         END DO      
         !
      END IF
      !
      F_ortho_max = 0.D0
      F_para_max  = 0.D0
      !
      DO i = 2, ( num_of_images - 1 )
         !
         F_ortho_max_i = MAXVAL( ABS( PES_gradient(:,i) - tangent(:,i) * &
                                    ( PES_gradient(:,i) .dot. tangent(:,i) ) ) )
         !
         elastic_gradient = tangent(:,i) * &
                ( ( k(i) + k(i-1) ) * norm( pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
                  ( k(i) + k(i+1) ) * norm( pbc( pos(:,(i+1)) - pos(:,i) ) ) )
         !
         F_para_max_i = MAXVAL( ABS( elastic_gradient(:) ) )
         !
         IF ( F_ortho_max_i > F_ortho_max ) F_ortho_max = F_ortho_max_i
         !
         IF ( F_para_max_i > F_para_max ) F_para_max = F_para_max_i
         !
      END DO
      !
      PRINT '(/5X,"F_ortho_max = ",F10.6)', F_ortho_max
      PRINT '(5X,"F_para_max  = ",F10.6)', F_para_max
      PRINT '(5X,"ALPHA       = ",F10.6)', F_ortho_max / F_para_max
      !
      k     = k * F_ortho_max / F_para_max
      k_max = k_max * F_ortho_max / F_para_max
      k_min = k_min * F_ortho_max / F_para_max
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
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER :: i
      !
      ! ... end of local variables
      !
      !
      CALL elastic_constants()
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
                                    ( PES_gradient(:,i) .dot. tangent(:,i) ) 
            ! 
         ELSE IF ( ( .NOT. free_minimization(i) ) .AND. &
                   ( i > 1 ) .AND. ( i < num_of_images ) ) THEN
            !
            grad(:,i) = elastic_gradient + PES_gradient(:,i) - &
                        tangent(:,i) * ( PES_gradient(:,i) .dot. tangent(:,i) )
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
    !!! FUNCTION path_tangent( index )
    !!! workaround for ifc8 compiler internal error
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
    !!! END FUNCTION path_tangent
    !!! workaround for ifc8 compiler internal error
    END SUBROUTINE path_tangent_
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_PES( flag, stat )
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
    END SUBROUTINE born_oppenheimer_PES
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE search_mep()
      !-----------------------------------------------------------------------
      !
      USE io_files,      ONLY : iunneb, iunexit, exit_file
      USE formats,       ONLY : run_output, run_output_T_const
      USE neb_variables, ONLY : num_of_images, dim, pos, PES, error,       &
                                climbing, optimization,  CI_scheme,        &
                                Emax_index, temp, Emax, neb_thr, conv_neb, &
                                suspended_image, lsteep_des, lquick_min ,  &
                                ldamped_dyn, lmol_dyn, istep_neb, nstep_neb
      USE io_routines,   ONLY : write_restart, write_dat_files, write_output
      USE check_stop,    ONLY : check_stop_now
      USE io_global,     ONLY : ionode
      USE minimization_routines
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      REAL (KIND=DP)      :: err
      INTEGER             :: image
      LOGICAL             :: stat
      INTEGER             :: N_in, N_fin
      LOGICAL             :: file_exists
      !
      ! ... external functions
      !
      REAL (kind=DP), EXTERNAL :: get_clock
      !
      !
      conv_neb = .FALSE.
      !
      IF ( istep_neb == nstep_neb ) THEN
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
         IF ( check_stop_now( ) ) THEN
            !
            CALL stop_pw( .FALSE. )
            !
         END IF
         !
         ! ... energies and gradients acting on each image of the elastic band
         ! ... are computed calling a driver that perfoms the scf calculations
         !
         IF ( istep_neb == 0 ) THEN
            !
            CALL born_oppenheimer_PES( .TRUE., stat )
            !
         ELSE
            !
            CALL born_oppenheimer_PES( optimization, stat )
            !
         END IF
         !
         IF ( .NOT. stat ) THEN
            !
            conv_neb = .FALSE.
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
         istep_neb = istep_neb + 1
         !
         ! ... informations are written on the standard output
         !
         IF ( ionode ) THEN
            !
            IF ( lmol_dyn ) THEN
               !
               WRITE( UNIT = iunneb, FMT = run_output_T_const ) &
                   istep_neb,                    &
                   temp * AU * eV_to_kelvin, &
                   err * ( AU / BOHR_RADIUS_ANGS )
               !
            ELSE
               !
               WRITE( UNIT = iunneb, FMT = run_output ) &
                   istep_neb,                  &
                   ( Emax - PES(1) ) * AU, &
                   err * ( AU / BOHR_RADIUS_ANGS )
               !
            END IF
            !
            CALL write_output()
            !
         END IF
         !   
         ! ... the program checks if the convergence has been achieved
         !
         IF ( ( err * AU / BOHR_RADIUS_ANGS ) <= neb_thr )  THEN
            !
            IF ( ionode ) &
               WRITE( UNIT = iunneb, &
                      FMT = '(/,5X,"NEB: convergence achieved in ",I3, &
                             &     " iterations" )' ) istep_neb
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
         IF ( istep_neb >= nstep_neb ) THEN
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
END MODULE neb_base
