!
! Copyright (C) 2003-2004 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define USE_ELASTIC_CONSTANTS_RESCALING
!
!---------------------------------------------------------------------------
MODULE path_base
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the implementation of "NEB" and "SMD" methods into the 
  ! ... PWSCF-FPMD-CPV codes
  !
  ! ... Written by Carlo Sbraccia ( 2003-2004 )
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps32, pi, au, bohr_radius_angs, eV_to_kelvin
  !
  USE basic_algebra_routines
  !
  PRIVATE
  !
  PUBLIC :: initialize_path
  PUBLIC :: search_mep
  !
  CONTAINS
    !
    ! ... module procedures    
    !
    !-----------------------------------------------------------------------
    SUBROUTINE initialize_path( prog )
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : pos, restart_mode, calculation, &
                                   opt_scheme, climbing, nstep, input_images
      USE control_flags,    ONLY : conv_elec, lneb, lsmd
      USE ions_base,        ONLY : nat, if_pos
      USE io_files,         ONLY : prefix, iunpath, path_file, &
                                   dat_file, int_file, xyz_file, axsf_file
      USE cell_base,        ONLY : alat
      USE path_variables,   ONLY : pos_ => pos, &
                                   istep_path, nstep_path, dim, num_of_images, &
                                   pes, grad_pes, grad_proj, tangent, error,   &
                                   path_length, path_thr, deg_of_freedom, ds,  &
                                   first_last_opt, reset_vel, react_coord,     &
                                   llangevin, temp_req, use_multistep
      USE path_variables,   ONLY : climbing_ => climbing,                  &
                                   CI_scheme, vel, grad, elastic_grad,     &
                                   norm_grad, k, k_min, k_max, Emax_index, &
                                   vel_zeroed, frozen, pos_old, grad_old
      USE path_variables,   ONLY : num_of_modes, ft_vel_zeroed, ft_pos,ft_pes,&
                                   ft_vel, ft_grad, Nft, ft_coeff, ft_frozen,  &
                                   ft_error, norm_ft_grad, Nft_smooth,         &
                                   ft_pos_old, ft_grad_old
      USE path_formats,     ONLY : summary_fmt   
      USE io_global,        ONLY : meta_ionode
      USE parser,           ONLY : int_to_char
      USE path_io_routines, ONLY : read_restart
      USE path_variables,   ONLY : path_allocation
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
      INTEGER                     :: i, j, mode
      REAL (KIND=DP)              :: inter_image_dist
      REAL (KIND=DP), ALLOCATABLE :: d_R(:,:), image_spacing(:)
      CHARACTER (LEN=20)          :: num_of_images_char, nstep_path_char
      !
      !    
      ! ... output files are set
      !
      path_file = TRIM( prefix ) // ".path"
      dat_file  = TRIM( prefix ) // ".dat"
      int_file  = TRIM( prefix ) // ".int"
      xyz_file  = TRIM( prefix ) // ".xyz"
      axsf_file = TRIM( prefix ) // ".axsf"
      !
      ! ... istep is initialized to zero
      !
      istep_path = 0
      conv_elec  = .TRUE.
      !
      ! ... the dimension of all "path" arrays is set here
      ! ... ( It corresponds to the dimension of the configurational space )
      !
      dim = 3 * nat
      !
      IF ( lsmd ) THEN
         !
         ! ... some coefficients for string dynamics
         !
         Nft = ( num_of_images - 1 )
         !
         Nft_smooth = 1000
         !
         num_of_modes = ( Nft - 1 )
         !
         ft_coeff = 2.D0 / DBLE( Nft )
         !
      END IF
      !  
      ! ... dynamical allocation of arrays and initialization
      !
      IF ( lneb ) THEN
         !
         CALL path_allocation( 'neb' )
         !
         vel          = 0.D0
         pes          = 0.D0
         grad_pes     = 0.D0
         elastic_grad = 0.D0
         tangent      = 0.D0
         grad         = 0.D0
         norm_grad    = 0.D0
         error        = 0.D0
         k            = k_min
         pos_old      = 0.D0
         grad_old     = 0.D0
         frozen       = .FALSE.
         vel_zeroed   = .FALSE.
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
      ELSE IF ( lsmd ) THEN
         !
         CALL path_allocation( 'smd' )
         !
         pes           = 0.D0
         grad_pes      = 0.D0
         grad_proj     = 0.D0
         tangent       = 0.D0
         ft_vel_zeroed = .FALSE.
         !
         IF ( first_last_opt ) THEN
            !
            error       = 0.D0
            frozen      = .FALSE.
            vel         = 0.D0
            vel_zeroed  = .FALSE.
            grad        = 0.D0
            norm_grad   = 0.D0
            pos_old     = 0.D0
            grad_old    = 0.D0
            !
         END IF
         !
         ! ... fourier components
         !
         ft_pos       = 0.D0
         ft_pes       = 0.D0
         ft_vel       = 0.D0
         ft_grad      = 0.D0
         norm_ft_grad = 0.D0
         ft_error     = 0.D0
         ft_pos_old   = 0.D0
         ft_grad_old  = 0.D0
         ft_frozen    = .FALSE.
         !
      END IF
      !
      ! ... initial path is read ( restart_mode == "restart" ) 
      ! ... or generated ( restart_mode = "from_scratch" )
      !
      IF ( restart_mode == "restart" ) THEN
         !
         ALLOCATE( image_spacing( num_of_images - 1 ) )
         !
         CALL read_restart()
         !
         ! ... consistency between the input value of nstep and the value
         ! ... of nstep_path read from the restart_file is checked
         !
         IF ( nstep == 0 ) THEN
            !
            istep_path = 0
            nstep_path = nstep
            !
         END IF   
         !
         IF ( nstep > nstep_path ) nstep_path = nstep
         !
         ! ... path length is computed here
         !
         DO i = 1, ( num_of_images - 1 )
            !
            image_spacing(i) = norm( pos_(:,i+1) - pos_(:,i) )
            !
         END DO
         !
         path_length = SUM( image_spacing(:) )
         !
         inter_image_dist = SUM( image_spacing(:) ) / DBLE( num_of_images - 1 )
         !
      ELSE
         !
         CALL initial_guess()
         !
      END IF
      !
      DEALLOCATE( image_spacing )
      !
      ! ... the actual number of degrees of freedom is computed
      !
      deg_of_freedom = 0
      !
      DO i = 1, nat
         !
         IF ( if_pos(1,i) == 1 ) deg_of_freedom = deg_of_freedom + 1
         IF ( if_pos(2,i) == 1 ) deg_of_freedom = deg_of_freedom + 1
         IF ( if_pos(3,i) == 1 ) deg_of_freedom = deg_of_freedom + 1
         !
      END DO
      !
      ! ... details of the calculation are written on output (only by ionode)
      !
      IF ( meta_ionode ) THEN
         !
         nstep_path_char   = int_to_char( nstep_path )
         num_of_images_char = int_to_char( num_of_images )
         !
         WRITE( UNIT = iunpath, FMT = * )
         !
         WRITE( UNIT = iunpath, FMT = summary_fmt ) &
             "calculation", TRIM( calculation )
         !
         WRITE( UNIT = iunpath, FMT = summary_fmt ) &
             "restart_mode", TRIM( restart_mode )
         !
         IF ( lneb ) &
            WRITE( UNIT = iunpath, FMT = summary_fmt ) &
                "CI_scheme", TRIM( CI_scheme )
         !
         WRITE( UNIT = iunpath, FMT = summary_fmt ) &
             "opt_scheme", TRIM( opt_scheme )
         !
         WRITE( UNIT = iunpath, FMT = summary_fmt ) &
             "num_of_images", TRIM( num_of_images_char )
         !
         WRITE( UNIT = iunpath, FMT = summary_fmt ) &
             "nstep", TRIM( nstep_path_char )   
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"first_last_opt",T35," = ",1X,L1))' ) first_last_opt
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"reset_vel",T35," = ",1X,L1))' ) reset_vel
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"use_multistep",T35," = ",1X,L1))' ) use_multistep
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"ds",T35," = ",1X,F6.4," a.u.")' ) ds
         !
         IF ( lneb ) THEN
            !
            WRITE( UNIT = iunpath, &
                   FMT = '(5X,"k_max",T35," = ",1X,F6.4," a.u.")' ) k_max
            WRITE( UNIT = iunpath, &
                   FMT = '(5X,"k_min",T35," = ",1X,F6.4," a.u.")' ) k_min
            !
         END IF
         !
         IF ( llangevin ) &
            WRITE( UNIT = iunpath, &
                   FMT = '(5X,"required temperature",T35, &
                          &" = ",F6.1," K")' ) temp_req * eV_to_kelvin * au
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"path_thr",T35," = ",1X,F6.4," eV / A")' ) path_thr
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"initial path length",&
                      & T35," = ",F7.4," bohr")' ) path_length  
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"initial inter-image distance", &
                      & T35," = ",F7.4," bohr")' ) inter_image_dist
         !
      END IF
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE initial_guess()
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          REAL (KIND=DP) :: s
          !
          ! ... linear interpolation
          !
          ALLOCATE( image_spacing( input_images - 1 ) )
          ALLOCATE( d_R( dim,    ( input_images - 1 ) ) )
          !
          DO i = 1, ( input_images - 1 )
             !
             d_R(:,i) = ( pos(1:dim,i+1) - pos(1:dim,i) )
             !
             image_spacing(i) = norm( d_R(:,i) )
             !
          END DO   
          !
          path_length = SUM( image_spacing(:) )
          !
          inter_image_dist = path_length / DBLE( num_of_images - 1  )
          !
          FORALL( i = 1: ( input_images - 1 ) )
             !
             d_R(:,i) = d_R(:,i) / image_spacing(i)
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
             IF ( s > image_spacing(i) ) THEN
                !
                s = s - image_spacing(i)
                !
                i = i + 1
                !
             END IF   
             !
             IF ( i >= input_images ) &
                CALL errore( 'initialize_path', ' i >= input_images ', i )
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
             path_length = path_length * alat
             !
             inter_image_dist = inter_image_dist * alat
             !
             pos_(:,:) = pos_(:,:) * alat
             !
          END IF
          !
          RETURN
          !
        END SUBROUTINE initial_guess
        !
    END SUBROUTINE initialize_path
    !
    ! ... neb specific routines
    !
    !------------------------------------------------------------------------
    SUBROUTINE elastic_constants()
      !------------------------------------------------------------------------
      ! 
      USE path_variables,  ONLY : pos, num_of_images, Emax, Emin, k_max, &
                                  k_min, k, pes, grad_pes, elastic_grad, &
                                  tangent
      USE supercell,       ONLY : pbc
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
                             COS( pi * ( pes(i) - Emin ) / delta_E ) )
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
         norm_grad_V = norm( grad_pes(:,i) )
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
            norm_grad_V = norm( grad_pes(:,i) )
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
         F_ortho_max_i = MAXVAL( ABS( grad_pes(:,i) - tangent(:,i) * &
                                    ( grad_pes(:,i) .dot. tangent(:,i) ) ) )
         !
         elastic_grad(:) = tangent(:,i) * 0.5D0 * &
                ( ( k(i) + k(i-1) ) * norm( pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
                  ( k(i) + k(i+1) ) * norm( pbc( pos(:,(i+1)) - pos(:,i) ) ) )
         !
         F_para_max_i = MAXVAL( ABS( elastic_grad(:) ) )
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
    !------------------------------------------------------------------------
    SUBROUTINE neb_gradient()
      !------------------------------------------------------------------------
      !
      USE supercell,      ONLY : pbc
      USE path_variables, ONLY : pos, grad, norm_grad, elastic_grad, grad_pes, &
                                 k, lmol_dyn, num_of_images, climbing, tangent
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER :: i
      !
      !
      tangent = 0
      !
      DO i = 2, ( num_of_images - 1 )
         !
         ! ... tangent to the path ( normalized )
         !
         tangent(:,i) = path_tangent( i )
         !
         tangent(:,i) = tangent(:,i) / norm( tangent(:,i) )
         !
      END DO
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
               elastic_grad = 0.5D0 * &
                        ( ( k(i) + k(i-1) ) * pbc( pos(:,i) - pos(:,(i-1)) ) - &
                          ( k(i) + k(i+1) ) * pbc( pos(:,(i+1)) - pos(:,i) ) )
               !
            ELSE
               !
               ! ... elastic gradient only along the path ( variable elastic
               ! ... consatnt is used ) NEB recipe
               !
               elastic_grad = tangent(:,i) * 0.5D0 * &
                ( ( k(i) + k(i-1) ) * norm( pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
                  ( k(i) + k(i+1) ) * norm( pbc( pos(:,(i+1)) - pos(:,i) ) ) )
               !
            END IF
            !
         END IF
         !
         ! ... total gradient on each image ( climbing image is used if needed )
         ! ... only the component of the pes gradient orthogonal to the path is 
         ! ... taken into account
         !
         grad(:,i) = grad_pes(:,i)
         !
         IF ( climbing(i) ) THEN
            !
            grad(:,i) = grad(:,i) - 2.D0 * tangent(:,i) * &
                                            ( grad_pes(:,i) .dot. tangent(:,i) )
            ! 
         ELSE IF ( ( i > 1 ) .AND. ( i < num_of_images ) ) THEN
            !
            grad(:,i) = elastic_grad + grad_pes(:,i) - &
                             tangent(:,i) * ( grad_pes(:,i) .dot. tangent(:,i) )
            !
         END IF
         ! 
         norm_grad(i) = norm( grad(:,i) )
         !  
      END DO gradient_loop
      !
      RETURN
      !
    END SUBROUTINE neb_gradient
    !
    ! ... smd specific routines
    !
    !-----------------------------------------------------------------------
    SUBROUTINE update_num_of_images()
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : num_of_images_inp => num_of_images
      USE path_variables,   ONLY : istep_path, num_of_images, num_of_modes,  &
                                   Nft, ft_coeff, ft_vel_zeroed, pos, pes,   &
                                   use_multistep, grad_pes, err_max, error,  &
                                   frozen, vel, vel_zeroed, grad, norm_grad, &
                                   path_thr
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP) :: chamge_image_thr
      REAL (KIND=DP) :: multistep_coeff = 3.D0
      INTEGER        :: new_num_of_images
      LOGICAL        :: images_updated
      INTEGER        :: init_num_of_images = 3
      !
      !
      IF ( .NOT. use_multistep ) RETURN
      !
      images_updated = .FALSE.
      !
      chamge_image_thr = multistep_coeff * &
                         DBLE( num_of_images_inp - num_of_images ) * path_thr
      !
      IF ( istep_path == 0 ) THEN
         !
         ! ... initialization
         !
         CALL redispose_last_image( init_num_of_images )
         !
         num_of_images = init_num_of_images
         !
         images_updated = .TRUE.
         !
      ELSE IF ( err_max < chamge_image_thr ) THEN
         !
         new_num_of_images = MIN( num_of_images_inp, num_of_images + 2 )
         !
         CALL redispose_last_image( new_num_of_images )
         !
         IF ( new_num_of_images > num_of_images ) THEN
            !
            images_updated = .TRUE.
            !
            num_of_images = new_num_of_images
            !
            ft_vel_zeroed = .TRUE.
            !
         END IF
         !
      END IF
      !
      IF ( images_updated ) THEN
         !
         ! ... reciprocal space dimensions updated
         !
         Nft = ( num_of_images - 1 )
         !
         ft_coeff = 2.D0 / DBLE( Nft )
         !
      END IF
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE redispose_last_image( n )
          !--------------------------------------------------------------------
          !
          USE path_variables, ONLY : first_last_opt
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN) :: n
          !
          !
          pos(:,n)      = pos(:,num_of_images)
          pes(n)        = pes(num_of_images)
          grad_pes(:,n) = grad_pes(:,num_of_images)
          !
          IF ( first_last_opt ) THEN
             !
             error(n)      = error(num_of_images)
             frozen(n)     = frozen(num_of_images)
             vel(:,n)      = vel(:,num_of_images)
             vel_zeroed(n) = vel_zeroed(num_of_images)
             grad(:,n)     = grad(:,num_of_images)
             norm_grad(n ) = norm_grad(num_of_images)
             !
          END IF
          !
          RETURN
          !
        END SUBROUTINE redispose_last_image
        !
    END SUBROUTINE update_num_of_images
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_path_length()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim, path_length, pos, ft_pos, Nft, &
                                 Nft_smooth ,num_of_images, num_of_modes
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP), ALLOCATABLE :: r_h(:), r_n(:), delta_pos(:)
      REAL (KIND=DP)              :: x, delta_x
      INTEGER                     :: i, j, n
      !
      !
      ALLOCATE( r_h(       dim ) )
      ALLOCATE( r_n(       dim ) )
      ALLOCATE( delta_pos( dim ) )
      !
      delta_pos(:) = ( pos(:,num_of_images) - pos(:,1) )
      !
      path_length = 0.D0
      !
      r_h(:) = pos(:,1)
      !
      delta_x = 1.D0 / DBLE( Nft_smooth * Nft )
      !
      DO i = 1, Nft
         !
         DO j = 1, Nft_smooth
            !
            x = delta_x * DBLE( Nft_smooth * ( i - 1 ) + j )
            !
            r_n(:) = pos(:,1) + x * delta_pos(:)
            !
            DO n = 1, num_of_modes
               !
               r_n(:) = r_n(:) + ft_pos(:,n) * SIN( DBLE( n ) * pi * x )
               !
            END DO
            !
            path_length = path_length + norm( r_n - r_h )
            !
            r_h(:) = r_n(:)
            !
         END DO
         !
      END DO
      !
      r_n(:) = pos(:,num_of_images)
      !
      path_length = path_length + norm( r_n - r_h )
      !
      DEALLOCATE( r_h )
      DEALLOCATE( r_n )
      DEALLOCATE( delta_pos)
      !
      RETURN
      !
    END SUBROUTINE compute_path_length
    !
    !-----------------------------------------------------------------------
    SUBROUTINE to_real_space()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : num_of_modes, num_of_images, dim, &
                                 pos, ft_pos, Nft, Nft_smooth, path_length
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP), ALLOCATABLE :: r_h(:), r_n(:), delta_pos(:)
      REAL (KIND=DP)              :: x, delta_x, s, s_image
      INTEGER                     :: i, j, n, image
      !
      !
      ALLOCATE( r_h(       dim ) )
      ALLOCATE( r_n(       dim ) )
      ALLOCATE( delta_pos( dim ) )
      !
      delta_pos(:) = ( pos(:,num_of_images) - pos(:,1) )
      !
      s = 0.D0
      !
      image = 1
      !
      s_image = path_length / DBLE( Nft )
      !
      r_h(:) = pos(:,1)
      !
      delta_x = 1.D0 / DBLE( Nft_smooth * Nft )
      !
      DO i = 1, Nft
         !
         DO j = 1, Nft_smooth
            !
            x = delta_x * DBLE( Nft_smooth * ( i - 1 ) + j )
            !
            r_n(:) = pos(:,1) + x * delta_pos(:)
            !
            DO n = 1, num_of_modes
               !
               r_n(:) = r_n(:) + ft_pos(:,n) * SIN( DBLE( n ) * pi * x )
               !
            END DO
            !
            s = s + norm( r_n - r_h )
            !
            IF ( s >= s_image ) THEN
               !
               image = image + 1
               !
               pos(:,image) = r_n(:)
               !
               s_image = DBLE( image ) * path_length / DBLE( Nft )
               !
            END IF
            !
            r_h(:) = r_n(:)
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( r_h )
      DEALLOCATE( r_n )
      DEALLOCATE( delta_pos)
      !
      RETURN
      !
    END SUBROUTINE to_real_space
    !
    !------------------------------------------------------------------------
    SUBROUTINE to_reciprocal_space()
      !------------------------------------------------------------------------
      !
      ! ... functions that are fourier transformed are defined here :
      !
      ! ... f(j), j = 1,...,N   is the function in real space (it starts from 1)
      !
      ! ... where :  f(1) /= f(N)
      !
      ! ... f_star(j) = f(j+1) - f(0) + j/(N-1)*( f(N) - f(0) )
      !
      ! ... with index running in  j = 0,..., N-1  and :
      !
      ! ... f_star(0) = f_star(N-1) = 0   so that f_star(j) is stroed between
      ! ...                               0  and  Nft - 1 ( = N - 2 )
      !
      USE path_variables, ONLY : dim, num_of_images, num_of_modes, Nft, &
                                 Nft_smooth, pos, ft_pos, ft_coeff, pos_star
      !
      IMPLICIT NONE
      !
      INTEGER        :: j, n
      REAL (KIND=DP) :: x, coeff
      !
      !
      ft_pos  = 0.D0
      !
      DO j = 0, ( Nft - 1 )
         !
         x = DBLE( j ) / DBLE( Nft )
         !
         pos_star(:,j) = pos(:,j+1) - pos(:,1) - &
                         x * ( pos(:,num_of_images) - pos(:,1) )
         !
      END DO
      !
      ! ... fourier components of pos_star are computed
      !
      DO n = 1, num_of_modes
         !
         coeff = DBLE( n ) * pi / DBLE( Nft )
         !
         DO j = 0, ( Nft - 1 )
            !
            x = DBLE( j )
            !
            ft_pos(:,n) = ft_pos(:,n) + &
                          pos_star(:,j) * SIN( coeff * x )
            !
         END DO
         !
      END DO
      !
      ! ... normalization
      !
      ft_pos = ft_pos * ft_coeff
      !
      RETURN
      !
    END SUBROUTINE to_reciprocal_space
    !
    !-----------------------------------------------------------------------
    SUBROUTINE smd_gradient()
      !-----------------------------------------------------------------------
      !
      ! ... functions that are fourier transformed are defined here :
      !
      ! ... f(j), j = 1,...,N   is the function in real space (it starts from 1)
      !
      ! ... where :  f(1) /= f(N)
      !
      ! ... f_star(j) = f(j+1) - f(0) + j/(N-1)*( f(N) - f(0) )
      !
      ! ... with index running in  j = 0,..., N-1  and :
      !
      ! ... f_star(0) = f_star(N-1) = 0   so that f_star(j) is stroed between
      ! ...                               0  and  Nft - 1 ( = N - 2 )
      !
      USE ions_base,      ONLY : if_pos
      USE path_variables, ONLY : dim, num_of_images, num_of_modes, Nft,  &
                                 Nft_smooth, pes, grad_pes, grad_proj,   &
                                 tangent, ft_pes, ft_grad, norm_ft_grad, &
                                 ft_coeff, pes_star, grad_proj_star,     &
                                 llangevin, lang_proj, lang_proj_star,   &
                                 ft_lang, grad, first_last_opt, ds
      !
      IMPLICIT NONE
      !
      INTEGER                     :: j, n
      REAL (KIND=DP)              :: x, coeff
      REAL (KIND=DP), ALLOCATABLE :: phi(:)
      !
      !
      IF ( first_last_opt ) THEN
         !s
         grad(:,:) = grad_pes(:,:)
         !
      END IF
      !
      ft_pes  = 0.D0
      ft_grad = 0.D0
      !
      IF ( llangevin ) THEN
         !
         ALLOCATE( phi( dim ) )
         !
         ft_lang = 0.D0
         !
      END IF
      !
      ! ... we compute the tangent to the path and project pes gradients 
      ! ... (in real space)
      !
      DO j = 1, num_of_images
         !
         ! ... tangent to the path ( normalized )
         !
         tangent(:,j) = path_tangent( j )
         !
         tangent(:,j) = tangent(:,j) / norm( tangent(:,j) )
         !
         ! ... projection of the pes gradients
         !
         grad_proj(:,j) = grad_pes(:,j) - &
                          tangent(:,j) * ( tangent(:,j) .dot. grad_pes(:,j) )
         !
         IF ( llangevin ) THEN
            !
            ! ... the random term used in langevin dynamics is generated here
            !
            phi(:) = gaussian_vect() * DBLE( RESHAPE( if_pos, (/ dim /) ) )
            !
            IF ( first_last_opt ) &
               grad(:,j) = grad(:,j) - phi(:) / SQRT( ds )
            !
            lang_proj(:,j) = phi(:) - &
                             tangent(:,j) * ( tangent(:,j) .dot. phi(:) )
            !
         END IF
         !
      END DO
      !
      ! ... here we compute pes_star, grad_proj_star and lang_proj_star
      !
      DO j = 0, ( Nft - 1 )
         !
         x = DBLE( j ) / DBLE( Nft )
         !
         pes_star(j) = pes(j+1) - pes(1) - &
                       x * ( pes(num_of_images) - pes(1) )
         !
         grad_proj_star(:,j) = grad_proj(:,j+1) - grad_proj(:,1) - x * &
                               ( grad_proj(:,num_of_images) - grad_proj(:,1) )
         !
         IF ( llangevin ) THEN
            !            
            lang_proj_star(:,j) = lang_proj(:,j+1) - lang_proj(:,1) - x * &
                                 ( lang_proj(:,num_of_images) - lang_proj(:,1) )
            !
         END IF
         !
      END DO
      !
      ! ... here we compute fourier components for pes_star and grad_proj_star
      !
      DO n = 1, num_of_modes
         !
         coeff = DBLE( n ) * pi / DBLE( Nft )
         !
         DO j = 0, ( Nft - 1 )
            !
            x = DBLE( j )
            !
            ft_pes(n) = ft_pes(n) + pes_star(j) * SIN( coeff * x )
            !
            ft_grad(:,n) = ft_grad(:,n) + &
                           grad_proj_star(:,j) * SIN( coeff * x )
            !
            IF ( llangevin ) &
               ft_lang(:,n) = ft_lang(:,n) + &
                              lang_proj_star(:,j) * SIN( coeff * x )
            !
         END DO
         !
         ! ... and the norm of the fourier component of the gradient
         !
         norm_ft_grad(n) = norm( ft_grad(:,n) )
         !
      END DO
      !
      ! ... normalizations
      !
      ft_pes  = ft_pes  * ft_coeff
      ft_grad = ft_grad * ft_coeff      
      !
      norm_ft_grad = norm_ft_grad * ft_coeff
      !
      IF ( llangevin ) THEN
         !
         DEALLOCATE( phi )
         !
         ft_lang = ft_lang * ft_coeff
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE smd_gradient
    !
    ! ... shared routines
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_error( err_out )
      !-----------------------------------------------------------------------
      !
      USE control_flags,  ONLY : lneb, lsmd
      USE path_variables, ONLY : num_of_images, num_of_modes, grad, &
                                 first_last_opt, path_thr, error, ft_error, &
                                 ft_pes, ft_grad, frozen, ft_frozen
      USE mp_global, ONLY : mpime
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      REAL (KIND=DP), OPTIONAL, INTENT(OUT) :: err_out
      !
      ! ... local variables
      !
      INTEGER        :: i, n
      INTEGER        :: N_in, N_fin
      REAL (KIND=DP) :: err_max
      !
      !
      err_max = 0.D0
      !
      IF ( lneb ) THEN
         !
         IF ( first_last_opt ) THEN
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
            error(i) = MAXVAL( ABS( grad(:,i) ) ) / bohr_radius_angs * au
            !
         END DO
         !
         err_max = MAXVAL( error(N_in:N_fin), 1 )
         !
         frozen(:) = ( error(:) < MAX( 0.5D0 * err_max, path_thr ) )
         !
      ELSE IF ( lsmd ) THEN
         !
         DO n = 1, num_of_modes
            !
            ! ... the error ( in eV / A ) is given by the largest component 
            ! ... of the gradient vector (in reciprocal space)
            !
            ft_error(n) = MAXVAL( ABS( ft_grad(:,n) ) ) / bohr_radius_angs * au
            !
            IF ( ft_error(n) > err_max ) err_max = ft_error(n)
            !
         END DO
         !
         err_max = MAXVAL( ft_error, 1 )
         !
         IF ( first_last_opt ) THEN
            !
            error(1) = MAXVAL( ABS( grad(:,1) ) )
            !
            error(num_of_images) = MAXVAL( ABS( grad(:,num_of_images) ) )
            !
            error = error / bohr_radius_angs * au
            !
            frozen(:) = ( error(:) < path_thr )
            !
            err_max = MAX( error(1), error(num_of_images), err_max )
            !
         END IF
         !
      END IF
      !
      IF ( PRESENT( err_out ) ) err_out = err_max
      !
      RETURN
      !
    END SUBROUTINE compute_error
    !
    !-----------------------------------------------------------------------
    FUNCTION path_tangent( index )
      !-----------------------------------------------------------------------
      !
      USE supercell,      ONLY : pbc
      USE control_flags,  ONLY : lneb, lsmd
      USE path_variables, ONLY : pos, dim, num_of_modes, num_of_images, &
                                 pes, path_length, ft_pos, Nft
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
      INTEGER        :: n
      REAL (KIND=DP) :: x, pi_n
      REAL (KIND=DP) :: V_previous, V_actual, V_next
      REAL (KIND=DP) :: abs_next, abs_previous
      REAL (KIND=DP) :: delta_V_max, delta_V_min
      !
      !
      IF ( lneb ) THEN
         !
         ! ... NEB definition of the tangent
         !
         V_previous = pes( index - 1 )
         V_actual   = pes( index )
         V_next     = pes( index + 1 )
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
      ELSE IF ( lsmd ) THEN
         !
         ! ... tangent from fourier interpolation
         !
         x = DBLE( index - 1 ) / DBLE( Nft )
         !
         path_tangent(:) = ( pos(:,num_of_images) - pos(:,1) )
         !
         DO n = 1, num_of_modes
            !
            pi_n = pi * DBLE( n )
            !
            path_tangent(:) = path_tangent(:) + &
                              ft_pos(:,n) * pi_n * COS( pi_n * x )
            !
         END DO
         !
         path_tangent(:) = path_tangent(:) / path_length
         !
      END IF
      !
      RETURN
      !
    END FUNCTION path_tangent
    !
    !-----------------------------------------------------------------------
    FUNCTION gaussian_vect()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim, temp_req
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP) :: gaussian_vect(dim)
      REAL (KIND=DP) :: x1, x2, w, coeff
      INTEGER        :: i
      !
      REAL (KIND=DP), EXTERNAL :: rndm
      !
      !
      coeff = SQRT( 2.D0 * temp_req )
      !
      DO i = 1, dim, 2
         !
         gaussian_loop: DO
            !
            x1 = 2.D0 * rndm() - 1.D0
            x2 = 2.D0 * rndm() - 1.D0
            !
            w = x1 * x1 + x2 * x2
            !
            IF ( w < 1.D0 ) EXIT gaussian_loop
            !
         END DO gaussian_loop
         !
         w = SQRT( ( - 2.D0 * LOG( w ) ) / w )
         !
         gaussian_vect(i) = x1 * w * coeff
         !
         IF ( i >= dim ) EXIT
         !
         gaussian_vect(i+1) = x2 * w * coeff
         !
      END DO
      !
      RETURN
      !
    END FUNCTION gaussian_vect
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_saddle()
      !-----------------------------------------------------------------------
      !
      ! ... the transition state configuration and the forward activation 
      ! ... energy are computed here
      !
      USE control_flags,    ONLY : lneb, lsmd
      USE path_variables,   ONLY : dim, num_of_images, num_of_modes, Nft, &
                                   Nft_smooth, pes, ft_pes, pos, ft_pos,  &
                                   activation_energy, Emax_index
      USE path_io_routines, ONLY : write_ts_config
      !
      IMPLICIT NONE
      !
      INTEGER                     :: i, j, n
      REAL (KIND=DP)              :: ener, delta_e, x, x_ts, delta_x
      REAL (KIND=DP), ALLOCATABLE :: pos_ts(:)
      !
      !
      IF ( lneb ) THEN
         !
         activation_energy = ( pes(Emax_index) - pes(1) ) * au
         !
      ELSE IF ( lsmd ) THEN
         !
         ALLOCATE( pos_ts( dim ) )
         !
         delta_x = 1.D0 / DBLE( Nft_smooth * Nft )
         !
         delta_e = pes(num_of_images) - pes(1)
         !
         activation_energy = 0.D0
         !
         x_ts = 0.D0
         !
         DO i = 1, Nft
            !
            DO j = 0, ( Nft_smooth - 1 )
               !
               x = delta_x * DBLE( Nft_smooth * ( i - 1 ) + j ) 
               !
               ener = x * delta_e
               !
               DO n = 1, num_of_modes
                  !
                  ener = ener + ft_pes(n) * SIN( DBLE( n ) * pi * x )
                  !
               END DO
               !
               ener = ener * au
               !
               IF ( ener > activation_energy ) THEN
                  !
                  activation_energy = ener
                  !
                  x_ts = x
                  !
               END IF
               !
            END DO
            !
         END DO
         !
         pos_ts(:) = pos(:,1) + x_ts * ( pos(:,num_of_images) - pos(:,1) )
         !
         DO n = 1, num_of_modes
            !
            pos_ts(:) = pos_ts(:) + ft_pos(:,n) * SIN( DBLE( n ) * pi * x_ts )
            !
         END DO
         !
         CALL write_ts_config( pos_ts )
         !
         DEALLOCATE( pos_ts )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE find_saddle
    !
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_pes( stat )
      !------------------------------------------------------------------------
      !
      USE control_flags,  ONLY : lneb
      USE path_variables, ONLY : num_of_images, suspended_image, istep_path, &
                                 pes, first_last_opt, Emin , Emax, Emax_index, &
                                 frozen
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      LOGICAL, INTENT(OUT) :: stat
      !
      ! ... local variables
      !
      INTEGER        :: N_in, N_fin, i
      REAL (KIND=DP) :: val
      !
      !
      IF ( istep_path == 0 ) THEN
         !
         N_in  = 1
         N_fin = num_of_images
         !
      ELSE IF ( first_last_opt ) THEN
         !
         N_in  = 1
         N_fin = num_of_images
         !
         IF ( frozen(1) ) N_in = 2
         !
         IF ( frozen(num_of_images) ) N_fin = ( num_of_images - 1 )
         !
      ELSE
         !
         N_in  = 2
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
      IF ( lneb ) THEN
         !
#if defined (__PGI)
         !
         Emax_index = 1
         !
         Emax = pes(1)
         Emin = pes(1)
         !   
         DO i = 2, num_of_images
            !
            val = pes(i)
            !
            IF ( val < Emin ) Emin = val
            !
            IF ( val > Emax ) THEN
               !
               Emax = val
               !
               Emax_index = i
               !
            END IF
            !
         END DO
         !
#else
         !
         Emin       = MINVAL( pes(:) )
         Emax       = MAXVAL( pes(:) )
         Emax_index = MAXLOC( pes(:), 1 )
         !
#endif
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE born_oppenheimer_pes
    !
    !-----------------------------------------------------------------------
    SUBROUTINE search_mep()
      !-----------------------------------------------------------------------
      !
      USE control_flags,    ONLY : lneb, lsmd
      USE path_variables,   ONLY : conv_path, istep_path, nstep_path, &
                                   lquick_min, ldamped_dyn, lmol_dyn, &
                                   suspended_image, err_max
      USE path_variables,   ONLY : climbing, CI_scheme, Emax_index
      USE path_io_routines, ONLY : write_restart, write_dat_files, write_output
      USE check_stop,       ONLY : check_stop_now
      USE io_files,         ONLY : iunpath
      USE io_global,        ONLY : meta_ionode
      USE path_formats,     ONLY : scf_iter_fmt
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER :: mode, image
      LOGICAL :: stat
      !
      ! ... external functions
      !
      REAL (KIND=DP), EXTERNAL :: get_clock
      !
      !
      conv_path = .FALSE.
      !
      CALL search_mep_init()
      !
      IF ( istep_path == nstep_path ) THEN
         !
         CALL write_dat_files()
         !
         CALL write_output()
         !
         suspended_image = 0
         !
         CALL write_restart()
         !
         RETURN
         !
      END IF
      !
      ! ... path optimization loop
      !
      optimization: DO
         !
         IF ( meta_ionode ) &
            WRITE( UNIT = iunpath, FMT = scf_iter_fmt ) istep_path + 1
         !
         ! ... the restart file is written (in real space)
         !
         CALL write_restart()
         !
         IF ( suspended_image == 0 ) THEN
            !
            ! ... minimization step is done only in case of no suspended images
            ! ... when the simulation is started from scratch all gradients are
            ! ... zero and fourier components are not optimized.
            !
            CALL first_opt_step()
            !
         END IF
         !
         IF ( check_stop_now() ) THEN
            !
            ! ... the programs checks if the user has required a soft
            ! ... exit or if if maximum CPU time has been exceeded
            !
            CALL write_restart()
            !
            CALL stop_pw( .FALSE. )
            !
         END IF
         !
         ! ... energies and gradients acting on each image of the path (in real
         ! ... space) are computed calling a driver for the scf calculations
         !
         CALL born_oppenheimer_PES( stat )
         !
         IF ( .NOT. stat ) THEN
            !
            conv_path = .FALSE.
            !
            EXIT optimization
            !
         END IF         
         !
         ! ... istep_path is updated after a self-consistency step
         !
         istep_path = istep_path + 1
         !
         IF ( lneb ) THEN
            !
            IF ( CI_scheme == "highest-TS" ) THEN
               !
               climbing = .FALSE.
               !
               climbing(Emax_index) = .TRUE.
               !
            END IF
            !
            CALL neb_gradient()
            !
         ELSE IF ( lsmd ) THEN
            !
            ! ... reparametrized path in reciprocal space
            !
            CALL to_reciprocal_space()
            !
            ! ... the fourier components of the pes and of the projected 
            ! ... gradient are computed here
            !
            CALL smd_gradient()
            !
         END IF
         !
         ! ... the transition state configuration and the forward activation 
         ! ... energy are computed here
         !
         CALL find_saddle()
         !
         ! ... the error is computed here
         !
         CALL compute_error( err_max )
         !
         IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            ! ... a second minimization step is needed for those algorithms
            ! ... based on a velocity Verlet scheme 
            !
            CALL second_opt_step()
            !
         END IF
         !
         ! ... information is written on the files
         !
         CALL write_dat_files()
         !
         ! ... information is written on the standard output
         !
         CALL write_output()
         !
         ! ... exit conditions
         !
         IF ( check_exit( err_max ) ) EXIT optimization
         !
         suspended_image = 0
         !
      END DO optimization
      !
      ! ... the restart file is written before exit (again in real space)
      !
      CALL write_restart()
      !
      RETURN
      !
    END SUBROUTINE search_mep
    !
    !------------------------------------------------------------------------
    SUBROUTINE search_mep_init()
      !------------------------------------------------------------------------
      !
      USE control_flags,    ONLY : lneb, lsmd
      USE path_variables,   ONLY : istep_path, nstep_path, &
                                   suspended_image, reset_vel
      USE path_variables,   ONLY : ft_pos, ft_grad, ft_pos_old, ft_grad_old
      USE path_io_routines, ONLY : write_dat_files, write_output
      !
      IMPLICIT NONE
      !
      !
      IF ( suspended_image == 0 ) THEN
         !
         IF ( lneb ) THEN
            !
            ! ... neb forces
            !
            CALL neb_gradient()
            !
            ! ... the error is computed only when the run is not starting 
            ! ... from scratch or when reset_vel == .FALSE. 
            !
            IF ( ( istep_path > 0 ) .AND. &
                 ( .NOT. reset_vel ) ) CALL compute_error()
            !
         ELSE IF ( lsmd ) THEN
            !
            ! ... the fourier components of pos are computed here
            !
            CALL to_reciprocal_space()
            !
            ! ... the path-length is computed here
            !
            CALL compute_path_length()
            !
            ! ... the fourier components of the projected gradient
            ! ... are computed here
            !
            CALL smd_gradient()
            !
            ! ... the error is computed only when the run is not starting 
            ! ... from scratch
            !
            IF ( istep_path > 0 ) CALL compute_error()
            !
            ! ... ft_pos_old  and  ft_grad_old  are initialized here
            !
            ft_pos_old(:,:)  = ft_pos(:,:)
            ft_grad_old(:,:) = ft_grad(:,:)
            !
         END IF
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE search_mep_init
    !
    !------------------------------------------------------------------------
    FUNCTION check_exit( err_max )
      !------------------------------------------------------------------------
      !
      USE input_parameters, ONLY : num_of_images_inp => num_of_images
      USE control_flags,    ONLY : lneb, lsmd
      USE io_files,         ONLY : iunpath
      USE io_global,        ONLY : ionode
      USE path_variables,   ONLY : path_thr, istep_path, nstep_path, &
                                   conv_path, suspended_image, &
                                   num_of_images, llangevin, lmol_dyn
      USE path_formats,     ONLY : final_fmt
      !
      IMPLICIT NONE
      !
      LOGICAL                    :: check_exit
      REAL (KIND=DP), INTENT(IN) :: err_max
      LOGICAL                    :: exit_condition
      !
      !
      check_exit = .FALSE.
      !
      ! ... the program checks if the convergence has been achieved
      !
      exit_condition = ( .NOT. ( llangevin .OR. lmol_dyn )  ) .AND. & 
                       ( num_of_images == num_of_images_inp ) .AND. &
                       ( err_max <= path_thr )
                       
      !
      IF ( exit_condition )  THEN
         !
         WRITE( UNIT = iunpath, FMT = final_fmt )
         !
         IF ( ionode .AND. lneb ) &
            WRITE( UNIT = iunpath, &
                   FMT = '(/,5X,"neb: convergence achieved in ",I3, &
                          &     " iterations" )' ) istep_path
         IF ( ionode .AND. lsmd ) &
            WRITE( UNIT = iunpath, &
                   FMT = '(/,5X,"smd: convergence achieved in ",I3, &
                          &     " iterations" )' ) istep_path
         !
         suspended_image = 0
         !
         conv_path = .TRUE.
         !
         check_exit = .TRUE.
         !
         RETURN
         !
      END IF
      !
      ! ... the program checks if the maximum number of iterations has
      ! ... been reached
      !
      IF ( istep_path >= nstep_path ) THEN
         !
         WRITE( UNIT = iunpath, FMT = final_fmt )
         !         
         IF ( ionode .AND. lneb ) &
            WRITE( UNIT = iunpath, &
                   FMT = '(/,5X,"neb: reached the maximum number of ", &
                          &     "steps")' )
         IF ( ionode .AND. lsmd ) &
            WRITE( UNIT = iunpath, &
                   FMT = '(/,5X,"smd: reached the maximum number of ", &
                          &     "steps")' )
         !
         suspended_image = 0
         !
         check_exit = .TRUE.
         !
         RETURN
         !
      END IF
      !
      RETURN
      !
    END FUNCTION check_exit
    !
    !------------------------------------------------------------------------
    SUBROUTINE first_opt_step()
      !------------------------------------------------------------------------
      !
      USE control_flags,  ONLY : lneb, lsmd
      USE path_variables, ONLY : istep_path, first_last_opt, &
                                 num_of_images, num_of_modes, Nft, &
                                 lsteep_des, lquick_min, ldamped_dyn, &
                                 lmol_dyn, llangevin
      USE path_opt_routines
      !
      IMPLICIT NONE
      !
      INTEGER :: image, mode
      !
      !
      IF ( istep_path == 0 ) THEN
         !
         IF ( lsmd ) THEN
            !
            ! ... the number of images is initialized
            !
            CALL update_num_of_images()
            !
            ! ... real space representation of the path :
            !
            ! ... the path-length is computed here
            !
            CALL compute_path_length()
            !
            ! ... the path in real space with the new number of images is 
            ! ... obtained interpolating with the "old" number of modes
            !
            CALL to_real_space()
            !
            ! ... the number of modes is updated (if necessary)
            !
            num_of_modes = ( Nft - 1 )
            !
         END IF
         !
         RETURN
         !
      END IF
      !
      IF ( first_last_opt ) THEN
         !
         IF ( lsteep_des .OR. llangevin ) THEN
            !
            CALL r_steepest_descent( 1 )
            CALL r_steepest_descent( num_of_images )
            !
         ELSE IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            CALL r_velocity_Verlet_first_step( 1 )
            CALL r_velocity_Verlet_first_step( num_of_images )
            !
         END IF
         !
      END IF
      !
      IF ( lneb ) THEN
         !
         first_r_opt_loop: DO image = 2, ( num_of_images - 1 )
            !
            IF ( lsteep_des ) THEN
               !
               CALL r_steepest_descent( image )
               !
            ELSE IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
               !
               CALL r_velocity_Verlet_first_step( image )
               !
            END IF
            !
         END DO first_r_opt_loop
         !
      ELSE IF ( lsmd ) THEN
         !
         ! ... the number of images is updated
         !
         CALL update_num_of_images()
         !
         first_ft_opt_loop: DO mode = 1, num_of_modes
            !
            IF ( lsteep_des .OR. llangevin ) THEN
               !
               CALL ft_steepest_descent( mode )
               !
            ELSE IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
               !
               CALL ft_velocity_Verlet_first_step( mode )
               !
            END IF
            !
         END DO first_ft_opt_loop
         !
         ! ... real space representation of the path :
         !
         ! ... the path-length is computed here
         !
         CALL compute_path_length()
         !
         ! ... the path in real space with the new number of images is 
         ! ... obtained interpolating with the "old" number of modes
         !
         CALL to_real_space()
         !
         ! ... the number of modes is updated (if necessary)
         !
         num_of_modes = ( Nft - 1 )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE first_opt_step
    !
    !------------------------------------------------------------------------
    SUBROUTINE second_opt_step()
      !------------------------------------------------------------------------
      !
      USE control_flags,  ONLY : lneb, lsmd
      USE path_variables, ONLY : first_last_opt, num_of_images, num_of_modes, &
                                 lquick_min, ldamped_dyn, lmol_dyn, llangevin
      USE path_opt_routines
      !
      IMPLICIT NONE
      !
      INTEGER :: image, mode
      !
      !
      IF ( first_last_opt ) THEN
         !
         IF ( lquick_min ) THEN
            !
            CALL r_quick_min_second_step( 1 )
            CALL r_quick_min_second_step( num_of_images )
            !
         ELSE IF ( ldamped_dyn .OR. lmol_dyn ) THEN
            !
            CALL r_velocity_Verlet_second_step( 1 )
            CALL r_velocity_Verlet_second_step( num_of_images )
            !
         END IF
         !
      END IF
      !
      IF ( lneb ) THEN
         !
         second_r_opt_loop: DO image = 2, ( num_of_images - 1 )
            !
            IF ( lquick_min ) THEN
               !
               CALL r_quick_min_second_step( image )
               !
            ELSE IF ( ldamped_dyn .OR. lmol_dyn ) THEN
               !
               CALL r_velocity_Verlet_second_step( image )
               !
            END IF
            !
         END DO second_r_opt_loop
         !
      ELSE IF ( lsmd ) THEN
         !
         second_ft_opt_loop: DO mode = 1, num_of_modes
            !
            IF ( lquick_min ) THEN
               !
               CALL ft_quick_min_second_step( mode )
               !
            ELSE IF ( ldamped_dyn .OR. lmol_dyn ) THEN
               !
               CALL ft_velocity_Verlet_second_step( mode )
               !
            END IF
            !
         END DO second_ft_opt_loop
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE second_opt_step    
    !
END MODULE path_base
