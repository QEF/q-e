!
! Copyright (C) 2003-2004 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define  USE_ELASTIC_CONSTANTS_RESCALING
!#define  DEBUG_ELASTIC_CONSTANTS
!#define  DEBUG_H3COLL
!
!-----------------------------------------------------------------------
MODULE path_base
  !-----------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the implementation of "NEB" and "STRING" methods into the 
  ! ... PWSCF-FPMD-CPV codes
  !
  ! ... Written by Carlo Sbraccia ( 2003-2004 )
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps32, pi, au, bohr_radius_angs, eV_to_kelvin
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
                                   minimization_scheme, climbing, nstep, ds, &
                                   input_images
      USE control_flags,    ONLY : conv_elec, lneb, lsmd
      USE ions_base,        ONLY : nat                                   
      USE io_files,         ONLY : prefix, iunpath, path_file, &
                                   dat_file, int_file, xyz_file, axsf_file
      USE cell_base,        ONLY : alat
      USE path_variables,   ONLY : pos_ => pos, &
                                   ds_  => ds,  &
                                   istep_path, nstep_path, dim, num_of_images, &
                                   pes, grad_pes, grad_proj, tangent, error,   &
                                   frozen, path_length, first_last_opt,        &
                                   reset_vel, path_thr
      USE path_variables,   ONLY : climbing_ => climbing,                  &
                                   CI_scheme, vel, grad, elastic_grad,     &
                                   norm_grad, k, k_min, k_max, Emax_index, &
                                   vel_zeroed
      USE path_variables,   ONLY : num_of_modes, ft_vel_zeroed, ft_pos, ft_pes,&
                                   ft_vel, ft_grad, Nft, ft_coeff, ft_frozen,  &
                                   ft_error, norm_ft_grad, Nft_smooth
      USE path_formats,     ONLY : summary_fmt   
      USE io_global,        ONLY : ionode
      USE parser,           ONLY : int_to_char
      USE path_io_routines, ONLY : read_restart
      USE path_variables,   ONLY : path_allocation
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
      INTEGER                     :: i, j, mode
      REAL (KIND=DP)              :: s, inter_image_dist
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
      ! ... the dimension of all neb arrays is set 
      ! ... ( It corresponds to the dimension of the configurational space )
      !
      dim = 3 * nat
      !      
      ! ... some coefficients for string dynamics
      !
      IF ( lsmd ) THEN
         !
         Nft = ( num_of_images - 1 )
         !
         Nft_smooth = 100
         !
         num_of_modes = ( Nft - 1 )
         !
         ft_coeff = 2.D0 / REAL( Nft )
         !
      END IF
      !  
      ! ... dynamical allocation of arrays and initialization
      !
      IF ( lneb ) THEN
         !
         CALL path_allocation( 'neb' )
         !
         ds_          = ds
         vel          = 0.D0
         pes          = 0.D0
         grad_pes     = 0.D0
         elastic_grad = 0.D0
         tangent      = 0.D0
         grad         = 0.D0
         norm_grad    = 0.D0
         error        = 0.D0
         k            = k_min
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
         ! ... mode dependent time step
         !
         FORALL( mode = 1 : ( Nft - 1 ) )
            !
            ds_(mode) = ds / SQRT( REAL( mode ) )
            !
         END FORALL
         !
         IF ( first_last_opt ) THEN
            !
            error       = 0.D0
            frozen      = .FALSE.
            vel         = 0.D0
            vel_zeroed  = .FALSE.
            grad        = 0.D0
            norm_grad   = 0.D0
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
         ! ... of nstep_neb read from the restart_file is checked
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
         inter_image_dist = SUM( image_spacing(:) ) / REAL( num_of_images - 1 )
         !
      ELSE
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
         inter_image_dist = SUM( image_spacing(:) ) / REAL( num_of_images - 1  )
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
      END IF
      !
      DEALLOCATE( image_spacing )
      !
      nstep_path_char = int_to_char( nstep_path )
      !
      ! ... the actual number of degrees of freedom is computed
      !
      CALL compute_deg_of_freedom()
      !
      IF ( lsmd ) THEN
         !
         ! ... fourier terms ( sin(n*pi*j/N) and n*pi*cos(n*pi*j/N) ) are
         ! ... computed once here and stored in a table
         !
         CALL compute_fourier_sin_and_cos()
         !
      END IF
      !
      ! ... details of the calculation are written on output (only by ionode)
      !
      IF ( ionode ) THEN
         !
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
             "minimization_scheme", TRIM( minimization_scheme )
         !
         WRITE( UNIT = iunpath, FMT = summary_fmt ) &
             "num_of_images", TRIM( num_of_images_char )
         !
         WRITE( UNIT = iunpath, FMT = summary_fmt ) &
             "nstep", TRIM( nstep_path_char )   
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"first_last_opt",T35," = ",L1))' ) first_last_opt
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"reset_vel",T35," = ",L1))' ) reset_vel
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"ds",T35," = ",F6.4," a.u.")' ) ds
         !
         IF ( lneb ) THEN
            !
            WRITE( UNIT = iunpath, &
                   FMT = '(5X,"k_max",T35," = ",F6.4," a.u.")' ) k_max
            WRITE( UNIT = iunpath, &
                   FMT = '(5X,"k_min",T35," = ",F6.4," a.u.")' ) k_min
            !
         END IF
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"path_thr",T35," = ",F6.4," eV / A")' ) path_thr
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"initial path length",&
                      & T35," = ",F6.4," bohr")' ) path_length  
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"initial inter-image distance", &
                      & T35," = ",F6.4," bohr")' ) inter_image_dist
         !
      END IF
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
           USE path_variables,   ONLY :  deg_of_freedom
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
      USE basic_algebra_routines
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
    !------------------------------------------------------------------------
    SUBROUTINE compute_fourier_sin_and_cos()
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : Nft, ft_sin, ft_cos
      !
      IMPLICIT NONE
      !
      INTEGER        :: j, n
      REAL (KIND=DP) :: x, pi_over_Nft
      !
      !
      pi_over_Nft = pi / REAL( Nft )
      !
      DO j = 0, ( Nft - 1 )
         !
         DO n = 1, ( Nft - 1 )
            !
            x = pi_over_Nft * REAL( n * j )
            !
            ft_sin(n,j) = SIN( x )
            ft_cos(n,j) = pi * REAL( n ) * COS( x )
            !
         END DO
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE compute_fourier_sin_and_cos
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_path_length()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim, path_length, pos, ft_pos, &
                                 num_of_images, Nft, Nft_smooth
      USE basic_algebra_routines
      USE path_variables, ONLY : istep_path
      USE parser
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP), ALLOCATABLE :: r_h(:), r_n(:), delta_pos(:)
      REAL (KIND=DP)              :: x, y
      INTEGER                     :: i, j, n
      !
      !
      ALLOCATE( r_h(       dim ) )
      ALLOCATE( r_n(       dim ) )
      ALLOCATE( delta_pos( dim ) )
      !
#if defined (DEBUG_H3COLL)
      !
      OPEN( UNIT = 666, &
            FILE = '2Dplot_' // TRIM( int_to_char( istep_path ) ) // '.dat' )
      !
#endif
      !
      delta_pos(:) = ( pos(:,num_of_images) - pos(:,1) )
      !
      path_length = 0.D0
      !
      r_h(:) = pos(:,1)
      !
      DO i = 1, ( num_of_images - 1 )
         !
         DO j = 0, ( Nft_smooth - 1 )
            !
            x = REAL( Nft_smooth * ( i - 1 ) + j ) / &
                REAL( Nft_smooth * ( num_of_images - 1 ) )
            !
            r_n(:) = pos(:,1) + x * delta_pos(:)
            !
            DO n = 1, ( Nft - 1 )
               !
               r_n(:) = r_n(:) + ft_pos(:,n) * SIN( REAL( n ) * pi * x )
               !
            END DO
            !
            path_length = path_length + norm( r_n - r_h )
            !
            r_h(:) = r_n(:)
            !
#if defined (DEBUG_H3COLL)
            !
            x = r_h(7) - r_h(4)
            y = r_h(4) - r_h(1)
            !
            WRITE( 666, '(2(2X,F10.8))' ) x, y
           ! IF ( j == 0 ) THEN
           !    !
           !    WRITE( 666, '(4(2X,F12.8),2X,I2)' ) &
           !        REAL( Nft_smooth * ( i - 1 ) + j ) / REAL( Nft_smooth * ( num_of_images - 1 ) ), path_length, x, y, i
           !    !
           ! ELSE
           !    !
           !    WRITE( 666, '(4(2X,F12.8))' ) &
           !        REAL( Nft_smooth * ( i - 1 ) + j ) / REAL( Nft_smooth * ( num_of_images - 1 ) ), path_length, x, y
           !    !
           ! END IF
#endif
            !
         END DO
         !
      END DO
      !
      r_n(:) = pos(:,num_of_images)
      !
      path_length = path_length + norm( r_n - r_h )
      !
#if defined (DEBUG_H3COLL)
      !
      r_h(:) = r_n(:)
      !
      x = r_h(7) - r_h(4)
      y = r_h(4) - r_h(1)
      !
      WRITE( 666, '(2(2X,F12.8))' ) x, y
     ! WRITE( 666, '(4(2X,F12.8),2X,I2)' ) 1.D0, path_length, x, y, i
      !
      CLOSE( 666 )
      !
#endif      
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
      USE path_variables, ONLY : num_of_images, Nft, path_length, &
                                 dim, pos, ft_pos, pos_star, Nft_smooth
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      INTEGER                     :: i, j, n, image
      REAL (KIND=DP)              :: x, ds      
      REAL (KIND=DP), ALLOCATABLE :: r_h(:), r_n(:), delta_pos(:)
      !
      !
      ALLOCATE( r_h(       dim ) )
      ALLOCATE( r_n(       dim ) )
      ALLOCATE( delta_pos( dim ) )
      !
      image = 1 
      !
      delta_pos(:) = ( pos(:,num_of_images) - pos(:,1) )
      !
      ds = 0.D0
      !
      r_h(:) = pos(:,1)
      !
      DO i = 1, ( num_of_images - 1 )
         !
         DO j = 0, ( Nft_smooth - 1 )
            !
            x = REAL( Nft_smooth * ( i - 1 ) + j ) / &
                REAL( Nft_smooth * ( num_of_images - 1 ) )
            !
            r_n(:) = pos(:,1) + x * delta_pos(:)
            !
            DO n = 1, ( Nft - 1 )
               !
               r_n(:) = r_n(:) + ft_pos(:,n) * SIN( REAL( n ) * pi * x )
               !
            END DO
            !
            ds = ds + norm( r_n - r_h )
            !
            IF ( ds > image * path_length / ( num_of_images - 1 ) ) THEN
               !
               image = image + 1
               !
               pos(:,image) = r_h(:)
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
      DEALLOCATE( delta_pos )
      !
      RETURN
      !
    END SUBROUTINE to_real_space
    !
    !------------------------------------------------------------------------
    SUBROUTINE to_reciprocal_space()
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim, num_of_images, Nft, pos, pes, &
                                 ft_pos, ft_pes, ft_sin, ft_coeff, &
                                 pos_star, pes_star
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER        :: j, n
      REAL (KIND=DP) :: ratio, delta_pes
      !
      ! ... initialization
      !
      ft_pos  = 0.D0
      ft_pes  = 0.D0      
      !      
      !
      ! ... functions that are fourier transformed are defined here :
      !
      ! ... f(j), j = 1,...,N   is the function in real space (it starts from 1)
      !
      ! ... f(1) /= f(N)
      !
      ! ... f_star(j) = f(j+1) - f(0) + j/(N-1) * ( f(N) - f(0) )
      !
      ! ... j = 0,...,N-1 
      !
      ! ... f_star(0) = f_star(N-1) = 0   so that f_star(j) is stroed between
      ! ...                               0  and  Nft = N - 1
      !
      DO j = 0, ( Nft - 1 )
         !
         ratio = REAL( j ) / REAL( Nft )
         !
         pos_star(:,j) = pos(:,j+1) - pos(:,1) - &
                         ratio * ( pos(:,num_of_images) - pos(:,1) )
         !
         pes_star(j) = pes(j+1) - pes(1) - &
                       ratio * ( pes(num_of_images) - pes(1) )
         !
      END DO
      !
      ! ... fourier components of pos_star and pes_star are computed
      !
      DO n = 1, ( Nft - 1 )
         !
         DO j = 0, ( Nft - 1 )
            !
            ft_pos(:,n) = ft_pos(:,n) + pos_star(:,j) * ft_sin(n,j)
            !
            ft_pes(n) = ft_pes(n) + pes_star(j) * ft_sin(n,j)
            !
         END DO
         !
      END DO
      !
      ! ... normalization
      !
      ft_pos = ft_pos * ft_coeff
      ft_pes = ft_pes * ft_coeff
      !
      RETURN
      !
    END SUBROUTINE to_reciprocal_space
    !
    !-----------------------------------------------------------------------
    SUBROUTINE smd_gradient()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim, num_of_images, Nft, grad_pes, &
                                 grad_proj, tangent, ft_grad, norm_ft_grad, &
                                 ft_sin, ft_coeff, grad_proj_star
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      INTEGER        :: j, n
      REAL (KIND=DP) :: ratio
      !
      !
      ft_grad = 0.D0
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
      END DO
      !
      ! ... here we compute grad_proj_star
      !
      DO j = 0, ( Nft - 1 )
         !
         ratio = REAL( j ) / REAL( Nft )
         !
         grad_proj_star(:,j) = grad_proj(:,j+1) - &
                               grad_proj(:,1) - ratio * &
                               ( grad_proj(:,num_of_images) - grad_proj(:,1) )
         !
      END DO
      !
      ! ... here we compute fourier components for grad_proj_star
      !
      DO n = 1, ( Nft - 1 )
         !
         DO j = 0, ( Nft - 1 )
            !
            ft_grad(:,n) = ft_grad(:,n) + grad_proj_star(:,j) * ft_sin(n,j)
            !
         END DO
         !
         ! ... and the norm of the fourier component
         !
         norm_ft_grad(n) = norm( ft_grad(:,n) )
         !
      END DO
      !
      ! ... normalization
      !
      ft_grad = ft_grad * ft_coeff
      !
      norm_ft_grad = norm_ft_grad * ft_coeff
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
        ! ft_frozen(:) = ( ft_error(:) < MAX( 0.5D0 * err_max, path_thr ) )
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
      USE path_variables, ONLY : pos, dim, Nft, num_of_images, &
                                 pes, path_length, ft_pos, ft_cos
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
      REAL (KIND=DP) :: path_length_inv
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
         path_length_inv = 1.D0 / path_length
         !
         path_tangent(:) = path_length_inv * &
                           ( pos(:,num_of_images) - pos(:,1) )
         !
         DO n = 1, ( Nft - 1 )
            !
            path_tangent(:) = path_tangent(:) + &
                              path_length_inv * ft_pos(:,n) * ft_cos(n,index)
            !
         END DO
         !
      END IF
      !
      RETURN
      !
    END FUNCTION path_tangent
    !
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_pes( stat )
      !------------------------------------------------------------------------
      !
      USE control_flags,  ONLY : lneb
      USE path_variables, ONLY : num_of_images, suspended_image, istep_path, &
                                 pes, first_last_opt, Emin , Emax, Emax_index
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      LOGICAL, INTENT(OUT) :: stat
      !
      ! ... local variables
      !
      INTEGER  :: N_in, N_fin
      !
      !
      IF ( first_last_opt .OR. ( istep_path == 0 ) ) THEN
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
      IF ( lneb ) THEN
         !
         Emin       = MINVAL( pes(:) )
         Emax       = MAXVAL( pes(:) )
         Emax_index = MAXLOC( pes(:), 1 )
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
      USE io_files,         ONLY : iunpath
      USE control_flags,    ONLY : lneb, lsmd
      USE path_variables,   ONLY : conv_path, istep_path, nstep_path, temp, &
                                   activation_energy, first_last_opt, pos,  &
                                   lsteep_des, lquick_min, ldamped_dyn,     &
                                   lmol_dyn, suspended_image, path_thr,     &
                                   num_of_images, grad_pes
      USE path_variables,   ONLY : climbing, CI_scheme, Emax_index, Emax
      USE path_variables,   ONLY : ft_pos, ft_grad, ft_pos_old, ft_grad_old, &
                                   num_of_modes                   
      USE path_io_routines, ONLY : write_restart, write_dat_files, write_output
      USE check_stop,       ONLY : check_stop_now
      USE io_global,        ONLY : ionode
      USE path_opt_routines
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      REAL (KIND=DP) :: err_max
      INTEGER        :: mode, image
      LOGICAL        :: stat
      !
      ! ... external functions
      !
      REAL (KIND=DP), EXTERNAL :: get_clock
      !
      !
      conv_path = .FALSE.
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
            ! ... from scratch
            !
            IF ( istep_path > 0 ) CALL compute_error()
            !
         ELSE IF ( lsmd ) THEN
            !
            ! ... the fourier components of pos, pes are computed here
            !
            CALL to_reciprocal_space()
            !
            ! ... the fourier components of the projected gradient
            ! ... are computed here
            !
            CALL compute_path_length()
            !
            CALL smd_gradient()
            !
            IF ( first_last_opt ) grad(:,:) = grad_pes(:,:)
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
      IF ( istep_path == nstep_path ) THEN
         !
         CALL write_dat_files()
         !
         CALL write_output()
         !
         suspended_image = 0
         !
         RETURN
         !
      END IF
      !
      minimization: DO
         !
         ! ... the restart file is written (in real space)
         !
         CALL write_restart()
         !
         ! ... minimization step is done only in case of no suspended images
         ! ... when the simulation is started from scratch all gradients are
         ! ... zero and fourier components are not optimized.
         !
         IF ( suspended_image == 0 ) THEN
            !
            IF ( first_last_opt ) THEN
               !
               IF ( lsteep_des ) THEN
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
               first_ft_opt_loop: DO mode = 1, num_of_modes
                  !
                  IF ( lsteep_des ) THEN
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
               ! ... the length of the new path is computed here
               !
               CALL compute_path_length()
               !
               ! ... back to real space
               !
               CALL to_real_space()
               !
            END IF
            !
         END IF
         !
         ! ... the programs checks if the user has required a soft
         ! ... exit or if if maximum CPU time has been exceeded
         !
         IF ( check_stop_now() ) THEN
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
            EXIT minimization
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
            ! ... the fourier components of pos, pes are computed here
            !
            CALL to_reciprocal_space()
            !
            ! ... the fourier components of the projected gradient
            ! ... are computed here
            !
            CALL compute_path_length()
            !
            CALL smd_gradient()
            !
            IF ( first_last_opt ) grad(:,:) = grad_pes(:,:)
            !
         END IF
         !
         ! ... the error is computed here
         !
         CALL compute_error( err_max )
         !
         ! ... a second minimization step is needed for those algorithms
         ! ... based on a velocity Verlet scheme
         !
         IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
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
         END IF
         !
         ! ... information is written on the files ( activation energy is
         ! ... computed in this routine )
         !
         CALL write_dat_files()
         !
         ! ... information is written on the standard output
         !
         CALL write_output()
         !   
         ! ... the program checks if the convergence has been achieved
         !
         IF ( err_max <= path_thr )  THEN
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
            conv_path = .TRUE.
            !
            EXIT minimization
            !
         END IF
         !
         ! ... the programs checks if the maximum number of iterations has
         ! ... been reached
         !
         IF ( istep_path >= nstep_path ) THEN
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
            EXIT minimization
            !
         END IF
         !
         suspended_image = 0
         !
      END DO minimization
      !
      ! ... the restart file is written before exit (again in real space)
      !
      CALL write_restart()
      !
      RETURN
      !
    END SUBROUTINE search_mep
    !
END MODULE path_base
