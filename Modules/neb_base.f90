!
! Copyright (C) 2003-2004 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define USE_ELASTIC_CONSTANTS_RESCALING
!#define DEBUG_ELASTIC_CONSTANTS
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
  USE constants,  ONLY : au, bohr_radius_angs, eV_to_kelvin
  !
  PRIVATE
  !
  PUBLIC :: io_neb_start, io_neb_stop
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
    SUBROUTINE io_neb_start()
      !-----------------------------------------------------------------------
      !
      USE io_global,  ONLY : stdout, io_global_start
      USE mp_global,  ONLY : me_image, my_image_id, root_image
      !
      IMPLICIT NONE
      !
      !
      ! ... the I/O node is set again according to the number of parallel
      ! ... images that have been required: for each parallel image there
      ! ... is only one node that does I/O
      !
      CALL io_global_start( my_image_id, root_image )
      !
      ! ... stdout is connected to a file ( different for each image ) 
      ! ... via unit 17 ( only root_image performes I/O )
      !
      IF ( me_image == root_image ) stdout = 17
      !
      RETURN
      !
    END SUBROUTINE io_neb_start
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE io_neb_stop()
      !-----------------------------------------------------------------------
      !
      USE io_global,  ONLY : stdout, io_global_start
      USE mp_global,  ONLY : mpime, root
      !
      IMPLICIT NONE
      !
      !
      ! ... the original I/O node set 
      !
      CALL io_global_start( mpime, root )
      !
      ! ... stdout is reconnected to standard output
      !
      stdout = 6
      !
      RETURN
      !
    END SUBROUTINE io_neb_stop
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE initialize_neb( prog )
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : pos, restart_mode, calculation, &
                                   minimization_scheme, climbing, nstep, ds, &
                                   input_images
      USE control_flags,    ONLY : conv_elec
      USE ions_base,        ONLY : nat                                   
      USE io_files,         ONLY : prefix, iunneb, neb_file, &
                                   dat_file, int_file, xyz_file, axsf_file
      USE cell_base,        ONLY : alat
      USE neb_variables,    ONLY : ds_       => ds, &
                                   pos_      => pos, &
                                   climbing_ => climbing, &
                                   pos_old, grad_old, frozen, reset_vel,       &
                                   vel, num_of_images, dim, PES, PES_gradient, &
                                   elastic_gradient, tangent, grad, norm_grad, &
                                   error, mass, free_minimization, CI_scheme,  &
                                   optimization, k, k_min, k_max,  Emax_index, &
                                   neb_thr, lquick_min, lmol_dyn, ldamped_dyn, &
                                   nstep_neb, istep_neb, suspended_image,      &
                                   vel_zeroed
      USE neb_variables,    ONLY : neb_dyn_allocation   
      USE parser,           ONLY : int_to_char
      USE io_routines,      ONLY : read_restart
      USE formats,          ONLY : stringfmt   
      USE io_global,        ONLY : ionode
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      ! ... input variables
      !
      CHARACTER (LEN=2) :: prog   
        ! ... specify the calling program
      !
      ! ... local variables
      !      
      INTEGER                     :: i, j
      REAL (KIND=DP)              :: s, inter_image_dist
      REAL (KIND=DP), ALLOCATABLE :: d_R(:,:), path_length(:)
      CHARACTER (LEN=20)          :: num_of_images_char, i_char
      !
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
      conv_elec = .TRUE.
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
      ! ... all other arrays are initialized 
      !
      ds_              = ds
      PES              = 0.D0
      PES_gradient     = 0.D0
      elastic_gradient = 0.D0
      tangent          = 0.D0
      grad             = 0.D0
      grad_old         = 0.D0
      norm_grad        = 0.D0
      error            = 0.D0
      mass             = 1.D0
      k                = k_min
      frozen           = .FALSE.
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
         vel_zeroed = .FALSE.
         !
      END IF
      !
      ! ... initial path is read ( restart_mode == "restart" ) 
      ! ... or generated ( restart_mode = "from_scratch" )
      !
      IF ( restart_mode == "restart" ) THEN
         !
         ALLOCATE( path_length( num_of_images - 1 ) )
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
         ! ... path length is computed here
         !
         DO i = 1, ( num_of_images - 1 )
            !
            path_length(i) = norm( pos_(:,i+1) - pos_(:,i) )
            !
         END DO
         !
         inter_image_dist = SUM( path_length(:) ) / REAL( num_of_images - 1 )
         !
         IF ( suspended_image == 0 ) THEN
            !
            ! ... NEB forces are computed for the first iteration
            !
            CALL compute_tangent()
            CALL gradient()
            !
         END IF   
         !
      ELSE
         !
         ! ... linear interpolation
         !
         ALLOCATE( path_length( input_images - 1 ) )
         ALLOCATE( d_R( dim, ( input_images - 1 ) ) )
         !
         DO i = 1, ( input_images - 1 )
            !
            d_R(:,i) = ( pos(1:dim,i+1) - pos(1:dim,i) )
            !
            path_length(i) = norm( d_R(:,i) )
            !
         END DO   
         !
         inter_image_dist = SUM( path_length(:) ) / REAL( num_of_images - 1  )
         !
         FORALL( i = 1: ( input_images - 1 ) )
            !
            d_R(:,i) = d_R(:,i) / path_length(i)
            !
         END FORALL   
         !
         pos_(:,1) = pos(1:dim,1)
         !
         i = 1
         s = 0.D0
         !
         DO j = 2, ( num_of_images - 1 )
            !
            s = s + inter_image_dist
            !
            IF ( s > path_length(i) ) THEN
               !
               s = s - path_length(i)
               !
               i = i + 1
               !
            END IF   
            !
            IF ( i >= input_images ) &
               CALL errore( 'initialize_neb', ' i >= input_images ', i )
            !
            pos_(:,j) = pos(1:dim,i) + s * d_R(:,i)
            !
         END DO
         !
         pos_(:,num_of_images) = pos(1:dim,input_images)
         !
         ! ... coordinates must be in bohr ( pwscf uses alat units )
         !
         IF ( prog == 'PW' ) THEN
            !
            path_length(:) = path_length(:) * alat
            !
            inter_image_dist = inter_image_dist * alat
            !
            pos_(:,:) = pos_(:,:) * alat
            !
         END IF
         !
         pos_old(:,:) = pos_(:,:)
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
         num_of_images_char = int_to_char( num_of_images )
         !
         WRITE( UNIT = iunneb, FMT = stringfmt ) &
             "calculation", TRIM( calculation )
         WRITE( UNIT = iunneb, FMT = stringfmt ) &
             "restart_mode", TRIM( restart_mode )
         WRITE( UNIT = iunneb, FMT = stringfmt ) &
             "CI_scheme", TRIM( CI_scheme )
         WRITE( UNIT = iunneb, FMT = stringfmt ) &
             "minimization_scheme", TRIM( minimization_scheme )
         WRITE( UNIT = iunneb, FMT = stringfmt ) &
             "num_of_images", TRIM( num_of_images_char )
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"optimization",T35," = ",L1))' ) optimization
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"reset_vel",T35," = ",L1))' ) reset_vel
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"ds",T35," = ",F6.4," a.u.")' ) ds
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"k_max",T35," = ",F6.4," a.u.")' ) k_max
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"k_min",T35," = ",F6.4," a.u.")' ) k_min
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"neb_thr",T35," = ",F6.4," eV / A")' ) neb_thr
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"initial path length",&
                      & T35," = ",F6.3," bohr")' ) SUM( path_length(:) )   
         WRITE( UNIT = iunneb, &
                FMT = '(5X,"initial inter-image distance", &
                      & T35," = ",F6.3," bohr")' ) inter_image_dist
         !
         IF ( CI_scheme == "manual" ) THEN
            !
            DO i = 1, num_of_images
               !
               IF ( climbing(i) ) THEN
                  !
                  i_char = int_to_char( i )
                  !
                  WRITE( UNIT = iunneb, &
                         FMT = stringfmt ) "climbing image", TRIM( i_char )
                  !
               END IF       
               !
            END DO
            !
         END IF
         !
         WRITE( UNIT = iunneb, FMT = '(/)' )
         !
      END IF
      !
      DEALLOCATE( path_length )
      !
      RETURN
      !
      CONTAINS
         !
         !-------------------------------------------------------------------
         SUBROUTINE compute_deg_of_freedom()
           !-------------------------------------------------------------------
           !
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
    SUBROUTINE elastic_constants()
      !------------------------------------------------------------------------
      ! 
      USE constants,              ONLY : pi, eps32
      USE neb_variables,          ONLY : pos, num_of_images, Emax, Emin, &
                                         k_max, k_min, k, PES, PES_gradient, &
                                         elastic_gradient, tangent
      USE supercell,              ONLY : pbc
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER       :: i
      REAL(KIND=DP) :: F_ortho_max, F_ortho_max_i, &
                       F_para_max_i, F_para_max, rescale_coeff        
      REAL(KIND=DP) :: delta_E, delta_norm_grad
      REAL(KIND=DP) :: k_sum, k_diff
      REAL(KIND=DP) :: norm_grad_V, norm_grad_V_min, norm_grad_V_max
      !
      ! ... local parameters
      !
      REAL(KIND=DP), PARAMETER :: k_minimal = 0.1D0
        ! minimum allowed input elastic constant
      REAL(KIND=DP), PARAMETER :: rescale_coeff_min = 0.5D0, &
                                  rescale_coeff_max = 1.5D0
        ! minimum allowed rescaling coefficient (  50% )
        ! maximum allowed rescaling coefficient ( 150% )
      !
      !
      rescale_coeff = k_max / k_min
      !
      k_min = MAX( k_min, k_minimal )
      k_max = MAX( k_max, k_minimal * rescale_coeff )
      !
      delta_E = Emax - Emin
      !
      k_sum  = k_max + k_min
      k_diff = k_max - k_min
      !
      k(:) = k_min
      !
      IF ( delta_E > eps32 ) THEN
         !
         DO i = 1, num_of_images 
            !
            k(i) = 0.5D0 * ( k_sum - k_diff * &
                             COS( pi * ( PES(i) - Emin ) / delta_E ) )
            !
         END DO
         !
      END IF
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
      delta_norm_grad = norm_grad_V_max - norm_grad_V_min
      !
      IF ( delta_norm_grad > eps32 ) THEN
         !
         DO i = 1, num_of_images 
            !
            norm_grad_V = norm( PES_gradient(:,i) )
            !
            k(i) = k(i) + 0.5D0 * ( k_sum - k_diff * &
                          COS( pi * ( norm_grad_V - norm_grad_V_min ) / & 
                          delta_norm_grad ) )
            !
         END DO
         !
      END IF
      !
      k(:) = 0.5D0 * k(:)
      !
      F_ortho_max = 0.D0
      F_para_max  = 0.D0
      !
      DO i = 2, ( num_of_images - 1 )
         !
         F_ortho_max_i = MAXVAL( ABS( PES_gradient(:,i) - tangent(:,i) * &
                                    ( PES_gradient(:,i) .dot. tangent(:,i) ) ) )
         !
         elastic_gradient(:) = tangent(:,i) * 0.5D0 * &
                ( ( k(i) + k(i-1) ) * norm( pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
                  ( k(i) + k(i+1) ) * norm( pbc( pos(:,(i+1)) - pos(:,i) ) ) )
         !
         F_para_max_i = MAXVAL( ABS( elastic_gradient(:) ) )
         !
         IF ( F_ortho_max_i > F_ortho_max ) F_ortho_max = F_ortho_max_i
         IF ( F_para_max_i  > F_para_max  ) F_para_max  = F_para_max_i
         !
      END DO
      !
      rescale_coeff = MAX( ( F_ortho_max / F_para_max ), rescale_coeff_min )
      rescale_coeff = MIN( rescale_coeff, rescale_coeff_max )
      !
#if defined (USE_ELASTIC_CONSTANTS_RESCALING)
      !
      k     = k * rescale_coeff
      k_max = k_max * rescale_coeff
      k_min = k_min * rescale_coeff
      !
#endif
      !
#if defined (DEBUG_ELASTIC_CONSTANTS)
      !
      PRINT '(/5X,"F_ortho_max = ",F10.6 )', F_ortho_max
      PRINT '( 5X,"F_para_max  = ",F10.6 )', F_para_max
      PRINT '( 5X,"ALPHA       = ",F10.6/)', rescale_coeff
      !
      DO i = 1, num_of_images
         !
         PRINT '(F8.4)', k(i)
         !
      END DO
      !
#endif      
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
               elastic_gradient = 0.5D0 * &
                        ( ( k(i) + k(i-1) ) * pbc( pos(:,i) - pos(:,(i-1)) ) - &
                          ( k(i) + k(i+1) ) * pbc( pos(:,(i+1)) - pos(:,i) ) )
               !
            ELSE
               !
               ! ... elastic gradient only along the path ( variable elastic
               ! ... consatnt is used ) NEB recipe
               !
               elastic_gradient = tangent(:,i) * 0.5D0 * &
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
    SUBROUTINE compute_error( err_out )
      !-----------------------------------------------------------------------
      !
      USE neb_variables, ONLY : num_of_images, optimization, &
                                neb_thr, error, grad, frozen
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      REAL (KIND=DP), OPTIONAL, INTENT(OUT) :: err_out
      !
      ! ... local variables
      !
      INTEGER        :: N_in, N_fin
      INTEGER        :: i
      REAL (KIND=DP) :: err_max
      !
      !
      err_max = 0.D0
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
      DO i = 1, num_of_images
         !
         ! ... the error is given by the largest component of the gradient 
         ! ... vector ( PES + SPRINGS )
         !
         error(i) = MAXVAL( ABS( grad(:,i) ) )
         !
         IF ( ( error(i) > err_max ) .AND. &
              ( i >= N_in .AND. i <= N_fin ) ) err_max = error(i)
         !
      END DO
      !
      ! ... an image is "frozen" if the error is less than 50% of the
      ! ... largest error
      !
      frozen(:) = ( error(:) < MAX( 0.5D0 * err_max, &
                                    neb_thr * bohr_radius_angs / au ) )
      !
      IF ( PRESENT( err_out ) ) err_out = err_max
      !
      RETURN
      !
    END SUBROUTINE compute_error
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
      !
      tangent = 0
      !
      DO image = 2, ( num_of_images - 1 )
         !
         ! ... tangent to the path ( normalized )
         !
         !!!tangent(:,image) = path_tangent( image )
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
    !-----------------------------------------------------------------------
    !!!FUNCTION path_tangent( index )
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
    !!!END FUNCTION path_tangent
    !!! workaround for ifc8 compiler internal error
    END SUBROUTINE path_tangent_
    !
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_PES( flag, stat )
      !------------------------------------------------------------------------
      !
      USE neb_variables, ONLY : num_of_images, Emax_index, Emin, Emax, &
                                PES, suspended_image
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      LOGICAL, INTENT(IN)  :: flag
      LOGICAL, INTENT(OUT) :: stat
      !
      ! ... local variables
      !
      INTEGER  :: i, image
      INTEGER  :: N_in, N_fin
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
      CALL compute_scf( N_in, N_fin, stat )
      !
      IF ( .NOT. stat ) RETURN
      !
      Emin       = MINVAL( PES(:) )
      Emax       = MAXVAL( PES(:) )
      Emax_index = MAXLOC( PES(:), 1 )
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
      USE io_files,      ONLY : iunneb
      USE formats,       ONLY : run_output, run_output_T_const
      USE neb_variables, ONLY : num_of_images, dim, pos, PES, error,        &
                                climbing, optimization,  CI_scheme, frozen, &
                                Emax_index, temp, Emax, neb_thr, conv_neb,  &
                                suspended_image, lsteep_des, lquick_min ,   &
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
      REAL (KIND=DP) :: err
      INTEGER        :: image
      LOGICAL        :: stat
      INTEGER        :: N_in, N_fin
      LOGICAL        :: file_exists
      !
      ! ... external functions
      !
      REAL (KIND=DP), EXTERNAL :: get_clock
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
               IF ( frozen(image) ) CYCLE first_minimization_loop
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
         ! ... istep_neb is updated after a self-consistency step
         !
         istep_neb = istep_neb + 1         
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
         CALL compute_error( err )
         !
         ! ... a second minimization step is needed for those algorithms
         ! ... based on a velocity Verlet scheme
         !
         IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            second_minimization_loop: DO image = N_in, N_fin
               !
               IF ( frozen(image) ) CYCLE second_minimization_loop
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
         ! ... information is written on the files
         !
         CALL write_dat_files()
         !
         ! ... information is written on the standard output
         !
         IF ( ionode ) THEN
            !
            IF ( lmol_dyn ) THEN
               !
               WRITE( UNIT = iunneb, FMT = run_output_T_const ) &
                   istep_neb,                    &
                   temp * au * eV_to_kelvin, &
                   err * ( au / bohr_radius_angs )
               !
            ELSE
               !
               WRITE( UNIT = iunneb, FMT = run_output ) &
                   istep_neb,                  &
                   ( Emax - PES(1) ) * au, &
                   err * ( au / bohr_radius_angs )
               !
            END IF
            !
            CALL write_output()
            !
         END IF
         !   
         ! ... the program checks if the convergence has been achieved
         !
         IF ( ( err * au / bohr_radius_angs ) <= neb_thr )  THEN
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
         ! ... been reached
         !
         IF ( istep_neb >= nstep_neb ) THEN
            !
            IF ( ionode ) &
               WRITE( UNIT = iunneb, &
                      FMT = '(/,5X,"NEB: reached the maximum number of ", &
                             &     "steps")' )
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
