!
! Copyright (C) 2003-2005 PWSCF-FPMD-CPV group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!#define NEW_TANGENT
!#define NEW_VEC
#define IPRINT 1
!
!---------------------------------------------------------------------------
MODULE path_base
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the implementation of "NEB" and "SMD" methods into the 
  ! ... PWSCF-FPMD-CPV codes
  !
  ! ... Written by Carlo Sbraccia ( 2003-2005 )
  !
  USE io_files,  ONLY : iunpath
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
      USE io_files,         ONLY : prefix, tmp_dir, path_file, dat_file, &
                                   int_file, xyz_file, axsf_file, broy_file
      USE cell_base,        ONLY : alat
      USE path_variables,   ONLY : pos_ => pos, &
                                   istep_path, nstep_path, dim, num_of_images, &
                                   pes, grad_pes, tangent, error, path_length, &
                                   path_thr, deg_of_freedom, ds, react_coord,  &
                                   first_last_opt, reset_vel, llangevin,       &
                                   temp_req, use_freezing, tune_load_balance,  &
                                   lbroyden
      USE path_variables,   ONLY : climbing_ => climbing,                  &
                                   CI_scheme, vel, grad, elastic_grad,     &
                                   norm_grad, k, k_min, k_max, Emax_index, &
                                   vel_zeroed, frozen, pos_old, grad_old,  &
                                   lquick_min, ldamped_dyn, lmol_dyn
      USE path_variables,   ONLY : num_of_modes, pos_in_av, pos_fin_av, &
                                   ft_pos, ft_pos_av, Nft, fixed_tan,   &
                                   ft_coeff, Nft_smooth, use_multistep, &
                                   free_energy, pos_in_h, pos_fin_h,    &
                                   ft_pos_h, av_counter
      USE path_formats,     ONLY : summary_fmt
      USE mp_global,        ONLY : nimage
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
      INTEGER                     :: i
      REAL (KIND=DP)              :: inter_image_dist, k_ratio
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
      broy_file = TRIM( tmp_dir ) // TRIM( prefix ) // ".broyden"
      !
      ! ... istep is initialised to zero
      !
      istep_path = 0
      conv_elec  = .TRUE.
      !
      ! ... the dimension of all "path" arrays is set here
      ! ... ( It corresponds to the dimension of the configurational space )
      !
      dim = 3 * nat
      !
      IF ( nimage > 1 ) THEN
         !
         ! ... the automatic tuning of the load balance in 
         ! ... image-parallelisation is switched off
         !
         tune_load_balance = .FALSE.
         !
         ! ... freezing allowed only with the automatic tuning of 
         ! ... the load balance
         !
         use_freezing = tune_load_balance
         !
      END IF
      !
      IF ( lneb ) THEN
         !
         ! ... elastic constants are rescaled here on
         ! ... the base of the input time step ds ( m = 4 ) :
         !
         k_ratio = k_min / k_max
         !
         IF ( lbroyden ) THEN
            !
            k_max = ( pi / 2.D0 )**2 / 16.D0
            !
         ELSE
            !
            k_max = ( pi / ds )**2 / 16.D0
            !
         END IF
         !
         k_min = k_max * k_ratio
         !
      ELSE IF ( lsmd ) THEN
         !
         av_counter = 0
         !
         IF ( free_energy) THEN
            !
            fixed_tan      = .TRUE.
            first_last_opt = .TRUE.
            !
            use_multistep  = .FALSE.
            use_freezing   = .FALSE.
            !
         END IF
         !
         ! ... some coefficients for string dynamics
         !
         Nft = ( num_of_images - 1 )
         !
         Nft_smooth = 50
         !
         num_of_modes = ( Nft - 1 )
         !
         ft_coeff = 2.D0 / DBLE( Nft )
         !
      END IF
      !  
      ! ... dynamical allocation of arrays and initialisation
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
         pes        = 0.D0
         grad_pes   = 0.D0
         tangent    = 0.D0
         error      = 0.D0
         vel        = 0.D0
         grad       = 0.D0
         norm_grad  = 0.D0
         pos_old    = 0.D0
         grad_old   = 0.D0
         !
         frozen     = .FALSE.
         vel_zeroed = .FALSE.
         !
         ! ... fourier components
         !
         ft_pos     = 0.D0
         ft_pos_av  = 0.D0
         ft_pos_h   = 0.D0
         !
         pos_in_av  = 0.D0
         pos_fin_av = 0.D0
         pos_in_h   = 0.D0
         pos_fin_h  = 0.D0
         !
      END IF
      !
      ! ... initial path is read from file ( restart_mode == "restart" ) 
      ! ... or generated from the input images ( restart_mode = "from_scratch" )
      ! ... It is alway read from file in the case of "free-energy" calculations
      !
      IF ( free_energy .OR. restart_mode == "restart" ) THEN
         !
         ALLOCATE( image_spacing( num_of_images - 1 ) )
         !
         CALL read_restart()
         !
         IF ( free_energy ) THEN
            !
            pos_in_av  = pos_(:,1)
            pos_fin_av = pos_(:,num_of_images)
            !
         END IF
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
                FMT = '(5X,"use_freezing",T35," = ",1X,L1))' ) use_freezing
                
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"fixed_tan",T35," = ",1X,L1))' ) fixed_tan
                
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
          INTEGER        :: i, j
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
          DO i = 1, input_images - 1
             !
             d_R(:,i) = d_R(:,i) / image_spacing(i)
             !
          END DO
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
    ! ... neb-specific routines
    !
    !-----------------------------------------------------------------------
    FUNCTION neb_tangent( index )
      !-----------------------------------------------------------------------
      !
      USE supercell,      ONLY : pbc
      USE path_variables, ONLY : pos, dim, num_of_images, pes
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER, INTENT(IN) :: index
      REAL (KIND=DP)      :: neb_tangent(dim)
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
      ! ... NEB definition of the tangent
      !
      IF ( index == 1 ) THEN
         !
         neb_tangent = pbc( pos(:,( index + 1 )) - pos(:,index) )
         !
         RETURN
         !
      ELSE IF ( index == num_of_images ) THEN
         !
         neb_tangent = pbc( pos(:,index ) - pos(:,( index - 1 )) )
         !
         RETURN
         !
      END IF
      !
      V_previous = pes( index - 1 )
      V_actual   = pes( index )
      V_next     = pes( index + 1 )
      !
      IF ( ( V_next > V_actual ) .AND. ( V_actual > V_previous ) ) THEN
         !
         neb_tangent = pbc( pos(:,( index + 1 )) - pos(:,index) )
         !
      ELSE IF ( ( V_next < V_actual ) .AND. ( V_actual < V_previous ) ) THEN
         !
         neb_tangent = pbc( pos(:,index) - pos(:,( index - 1 )) )
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
            neb_tangent = &
                   pbc( pos(:,( index + 1 )) - pos(:,index) ) * delta_V_max + & 
                   pbc( pos(:,index) - pos(:,( index - 1 )) ) * delta_V_min
            !
         ELSE IF ( V_next < V_previous ) THEN
            !
            neb_tangent = &
                   pbc( pos(:,( index + 1 )) - pos(:,index) ) * delta_V_min + &
                   pbc( pos(:,index) - pos(:,( index - 1 )) ) * delta_V_max
            !
         ELSE
            !
            neb_tangent = pbc( pos(:,( index + 1 )) - pos(:,( index - 1 )) ) 
            !
         END IF
         !
      END IF
      !
      RETURN
      !
    END FUNCTION neb_tangent
    !
    !------------------------------------------------------------------------
    SUBROUTINE elastic_constants()
      !------------------------------------------------------------------------
      ! 
      USE path_variables,  ONLY : pos, num_of_images, Emax, Emin, &
                                  k_max, k_min, k, pes, dim
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      INTEGER       :: i
      REAL(KIND=DP) :: delta_E
      REAL(KIND=DP) :: k_sum, k_diff
      !
      REAL(KIND=DP) :: omega( num_of_images ), v1( dim ), v2( dim )
      REAL(KIND=DP) :: omega_max
      !
      !
      k_sum  = k_max + k_min
      k_diff = k_max - k_min
      !
      k(:) = k_min
      !
#if defined (NEW_VEC)
      !
      ! ... here the curvature is computed
      !
      omega = 0.D0
      !
      DO i = 2, num_of_images - 1
         !
         v1 = ( pos(:,i+1) - pos(:,i) ) / norm( pos(:,i+1) - pos(:,i) )
         v2 = ( pos(:,i) - pos(:,i-1) ) / norm( pos(:,i) - pos(:,i-1) )
         !
         omega(i) = norm( v1 - v2 ) / norm( pos(:,i+1) - pos(:,i-1) )
         !
      END DO
      !
      omega_max = MAXVAL( omega )
      !
      IF ( omega_max > eps32 ) THEN
         !
         DO i = 1, num_of_images 
            !
            k(i) = 0.5D0 * ( k_sum - k_diff * COS( pi * omega(i) / omega_max ) )
            !
         END DO
         !
      END IF
      !
#else
      !
      delta_E = Emax - Emin
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
#endif
      !
      k(:) = 0.5D0 * k(:)
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
      USE path_variables, ONLY : pos, grad, norm_grad, elastic_grad,   &
                                 grad_pes, k, lmol_dyn, num_of_images, &
                                 climbing, tangent
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
      DO i = 1, num_of_images
         !
         ! ... tangent to the path ( normalised )
         !
         tangent(:,i) = neb_tangent( i )
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
            ! ... elastic gradient only along the path ( variable elastic
            ! ... consatnt is used ) NEB recipe
            !
            elastic_grad = tangent(:,i) * 0.5D0 * &
                ( ( k(i) + k(i-1) ) * norm( pbc( pos(:,i) - pos(:,(i-1)) ) ) - &
                  ( k(i) + k(i+1) ) * norm( pbc( pos(:,(i+1)) - pos(:,i) ) ) )
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
    ! ... smd-specific routines
    !
    !-----------------------------------------------------------------------
    SUBROUTINE update_num_of_images()
      !-----------------------------------------------------------------------
      !
      USE input_parameters, ONLY : num_of_images_inp => num_of_images
      USE path_variables,   ONLY : istep_path, num_of_images, &
                                   path_thr, Nft, ft_coeff, pos, pes, &
                                   use_multistep, grad_pes, err_max,  &
                                   frozen, vel, vel_zeroed, reset_broyden
      USE io_global,        ONLY : meta_ionode
      !
      IMPLICIT NONE
      !
      REAL (KIND=DP) :: change_image_thr
      REAL (KIND=DP) :: multistep_coeff = 2.D0
      INTEGER        :: new_num_of_images
      LOGICAL        :: images_updated
      INTEGER        :: N_in, N_fin
      INTEGER        :: init_num_of_images = 3
      !
      !
      IF ( .NOT. use_multistep ) RETURN
      !
      images_updated = .FALSE.
      !
      change_image_thr = multistep_coeff * &
                         DBLE( num_of_images_inp - num_of_images ) * path_thr
      !
      IF ( istep_path == 0 ) THEN
         !
         ! ... initialisation
         !
         IF ( meta_ionode ) &
            WRITE( UNIT = iunpath, &
                 & FMT = '(5X,"initial number of images = ",I3,/)' ) &
                init_num_of_images
         !
         CALL redispose_last_image( init_num_of_images )
         !
         num_of_images = init_num_of_images
         !
         images_updated = .TRUE.
         !
      ELSE IF ( err_max < change_image_thr ) THEN
         !
         new_num_of_images = MIN( num_of_images_inp, num_of_images + 2 )
         !
         IF ( new_num_of_images > num_of_images ) THEN
            !
            IF ( meta_ionode ) &
               WRITE( UNIT = iunpath, &
                    & FMT = '(5X,"new number of images = ",I3,/)' ) &
                   new_num_of_images
            !
            CALL redispose_last_image( new_num_of_images )
            !
            N_in  = 2
            N_fin = new_num_of_images - 1
            !
            vel(:,N_in:N_fin) = 0.D0
            !
            frozen(N_in:N_fin) = .FALSE.
            !
            vel_zeroed(N_in:N_fin) = .FALSE.
            !
            images_updated = .TRUE.
            !
            num_of_images = new_num_of_images
            !
            reset_broyden = .TRUE.
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
          USE path_variables, ONLY : error, grad, norm_grad, pos_old, grad_old
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
          error(n)      = error(num_of_images)
          vel(:,n)      = vel(:,num_of_images)
          grad(:,n)     = grad(:,num_of_images)
          norm_grad(n)  = norm_grad(num_of_images)
          !
          pos_old(:,n)  = pos_old(:,num_of_images)
          grad_old(:,n) = grad_old(:,num_of_images)
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
      REAL (KIND=DP)              :: x, delta_x, s, s_image, pi_n
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
               s_image = s_image + path_length / DBLE( Nft )
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
      REAL (KIND=DP) :: x, coeff, inv_Nft
      !
      !
      inv_Nft = 1.D0 / DBLE( Nft )
      !
      DO j = 0, ( Nft - 1 )
         !
         x = DBLE( j ) * inv_Nft
         !
         pos_star(:,j) = pos(:,j+1) - pos(:,1) - &
                         x * ( pos(:,num_of_images) - pos(:,1) )
         !
      END DO
      !
      ! ... fourier components of pos_star are computed
      !
      ft_pos  = 0.D0
      !
      DO n = 1, num_of_modes
         !
         coeff = DBLE( n ) * pi * inv_Nft
         !
         DO j = 0, ( Nft - 1 )
            !
            x = DBLE( j )
            !
            ft_pos(:,n) = ft_pos(:,n) + pos_star(:,j) * SIN( coeff * x )
            !
         END DO
         !
      END DO
      !
      ! ... normalisation
      !
      ft_pos = ft_pos * ft_coeff
      !
      RETURN
      !
    END SUBROUTINE to_reciprocal_space
    !
    !-----------------------------------------------------------------------
    SUBROUTINE smd_tangent( image )
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : pos, ft_pos, pos_in_av, pos_fin_av, &
                                 ft_pos_av,  num_of_modes, num_of_images, &
                                 tangent
      !
      IMPLICIT NONE
      !
      INTEGER,        INTENT(IN) :: image
      INTEGER                    :: n
      REAL (KIND=DP)             :: s, pi_n
      !
      !
      s = DBLE( image - 1 ) / DBLE( num_of_images - 1 )
      !
#if defined (NEW_TANGENT)
      !
      tangent(:,image) = ( pos_fin_av - pos_in_av )
      !
      DO n = 1, num_of_modes
         !
         pi_n = pi * DBLE( n )
         !
         tangent(:,image) = tangent(:,image) + &
                            ft_pos_av(:,n) * pi_n * COS( pi_n * s )
         !
      END DO
      !
#else
      !
      tangent(:,image) = ( pos(:,num_of_images) - pos(:,1) )
      !
      DO n = 1, num_of_modes
         !
         pi_n = pi * DBLE( n )
         !
         tangent(:,image) = tangent(:,image) + &
                            ft_pos(:,n) * pi_n * COS( pi_n * s )
         !
      END DO
      !
#endif
      !
      tangent(:,image) = tangent(:,image) / norm( tangent(:,image) )
      !
      RETURN
      !
    END SUBROUTINE smd_tangent
    !
    !-----------------------------------------------------------------------
    SUBROUTINE update_history()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : istep_path, av_counter, num_of_images,    &
                                 pos, ft_pos, pos_in_av, pos_fin_av,       &
                                 ft_pos_av, pos_in_h, pos_fin_h, ft_pos_h, &
                                 history_ndim
      !
      IMPLICIT NONE
      !
      INTEGER, SAVE :: ipos
      !
      !
      av_counter = MIN( av_counter + 1, history_ndim )
      !
      ipos = MOD( istep_path, history_ndim ) + 1
      !
      ! ... history update
      !
      pos_in_h(:,ipos)  = pos(:,1)
      pos_fin_h(:,ipos) = pos(:,num_of_images)
      !
      ft_pos_h(:,:,ipos) = ft_pos(:,:)
      !
      ! ... average path
      !
      pos_in_av  = SUM( ARRAY = pos_in_h,  DIM = 2 ) / av_counter
      pos_fin_av = SUM( ARRAY = pos_fin_h, DIM = 2 ) / av_counter
      !
      ft_pos_av  = SUM( ARRAY = ft_pos_h, DIM = 3 ) / av_counter
      !
    END SUBROUTINE update_history
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
      USE path_variables, ONLY : dim, num_of_images, grad_pes,   &
                                 tangent, llangevin, lang, grad, &
                                 norm_grad, first_last_opt, fixed_tan
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      !
      !
      IF ( .NOT. fixed_tan ) THEN
         !
         CALL to_reciprocal_space()
         !
         CALL update_history()
         !
      END IF
      !
      ! ... we project pes gradients and gaussian noise
      !
      DO i = 1, num_of_images
         !
         CALL smd_tangent( i )
         !
         IF ( llangevin ) THEN
            !
            ! ... the random term used in langevin dynamics is generated here
            !
            lang(:,i) = gaussian_vect() * DBLE( RESHAPE( if_pos, (/ dim /) ) )
            !
         END IF
         !
         IF ( fixed_tan .OR. &
              ( i > 1 ) .AND. ( i < num_of_images ) ) THEN
            !
            ! ... projection of the pes gradients 
            !
            grad(:,i) = grad_pes(:,i) - &
                        tangent(:,i) * ( tangent(:,i) .dot. grad_pes(:,i) )
            !
            IF ( llangevin ) THEN
               !
               lang(:,i) = lang(:,i) - &
                           tangent(:,i) * ( tangent(:,i) .dot. lang(:,i) )
               !
            END IF
            !
         ELSE
            !
            grad(:,i) = grad_pes(:,i)
            !
         END IF
         !
         norm_grad(i) = norm( grad(:,i) )
         !
      END DO
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
      USE path_variables, ONLY : num_of_images, grad, llangevin, &
                                 use_freezing, first_last_opt,   &
                                 path_thr, error, frozen
      USE mp_global,      ONLY : nimage, inter_image_comm
      USE mp,             ONLY : mp_bcast
      USE io_global,      ONLY : meta_ionode, meta_ionode_id
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
      INTEGER        :: N_in, N_fin, free_me, num_of_scf_images
      REAL (KIND=DP) :: err_max, val      
      !
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
         ! ... vector ( PES + SPRINGS in the neb case )
         !
         error(i) = MAXVAL( ABS( grad(:,i) ) ) / bohr_radius_angs * au
         !
      END DO
      !
      err_max = MAXVAL( error(N_in:N_fin), 1 )
      !
      IF ( use_freezing ) THEN
         !
         frozen(:) = ( error(:) < MAX( 0.5D0 * err_max, path_thr ) )
         !
      ELSE
         !
         frozen = .FALSE.
         !
      END IF
      !
      IF ( nimage > 1 .AND. use_freezing ) THEN
         !
         IF ( meta_ionode ) THEN
            !
            ! ... in the case of image-parallelisation the number of images
            ! ... to be optimised must be larger than nimage
            !
            IF ( nimage > ( N_fin - N_in + 1 ) ) &
               CALL errore( 'search_MEP', &
                          & 'nimage is larger than the number of images ', 1 )
            !
            find_scf_images: DO
               !
               num_of_scf_images = COUNT( .NOT. frozen )
               !
               IF ( num_of_scf_images >= nimage ) EXIT find_scf_images
               !
               free_me = MAXLOC( error, 1, frozen(N_in:N_fin) )
               !
               frozen(free_me) = .FALSE.
               !
            END DO find_scf_images
            !
         END IF
         !
         CALL mp_bcast( frozen, meta_ionode_id, inter_image_comm )
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
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_pes( stat )
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : num_of_images, suspended_image,  &
                                 istep_path, pes, first_last_opt, &
                                 Emin , Emax, Emax_index, frozen, &
                                 error
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
      !
      !
      IF ( istep_path == 0 .OR. first_last_opt ) THEN
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
      IF ( suspended_image /= 0 ) N_in = suspended_image
      !
      CALL compute_scf( N_in, N_fin, stat )
      !
      IF ( .NOT. stat ) RETURN
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
      Emin       = MINVAL( pes(1:num_of_images) )
      Emax       = MAXVAL( pes(1:num_of_images) )
      Emax_index = MAXLOC( pes(1:num_of_images), 1 )
      !
#endif
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
      USE path_variables,   ONLY : conv_path, istep_path, nstep_path,  &
                                   lquick_min, ldamped_dyn, lmol_dyn,  &
                                   suspended_image, activation_energy, &
                                   err_max, num_of_modes, Nft, pes,    &
                                   climbing, CI_scheme, Emax_index, fixed_tan
      USE path_io_routines, ONLY : write_restart, write_dat_files, write_output
      USE check_stop,       ONLY : check_stop_now
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
      ! ... path optimisation loop
      !
      optimisation: DO
         !
         IF ( meta_ionode ) &
            WRITE( UNIT = iunpath, FMT = scf_iter_fmt ) istep_path + 1
         !
         ! ... the restart file is written (in real space)
         !
         IF ( MOD( istep_path, IPRINT ) == 0 ) CALL write_restart()
         !
         IF ( suspended_image == 0 ) THEN
            !
            ! ... minimisation step is done only in case of no suspended images
            ! ... when the simulation is started from scratch all gradients are
            ! ... zero.
            !
            CALL first_opt_step()
            !
            IF ( lsmd .AND. .NOT. fixed_tan ) THEN
               !
               ! ... fourier components of the path
               !
               CALL to_reciprocal_space()
               !
               ! ... the path-length is computed here
               !
               CALL compute_path_length()
               !
               ! ... real space representation of the path :
               !
               CALL update_num_of_images()
               !               
               ! ... the path in real space with the new number of images is 
               ! ... obtained interpolating with the "old" number of modes
               !
               ! ... here the parametrisation of the path is enforced
               !
               CALL to_real_space()
               !
               ! ... the number of modes is updated (if necessary)
               !
               num_of_modes = ( Nft - 1 )
               !
            END IF
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
            EXIT optimisation
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
            ! ... the projected gradients are computed here
            !
            CALL smd_gradient()
            !
         END IF
         !
         ! ... the forward activation energy is computed here
         !
         activation_energy = ( pes(Emax_index) - pes(1) ) * au
         !
         ! ... the error is computed here
         !
         CALL compute_error( err_max )
         !
         IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            ! ... a second minimisation step is needed for those algorithms
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
         IF ( check_exit( err_max ) ) EXIT optimisation
         !
         suspended_image = 0
         !
      END DO optimisation
      !
      ! ... the restart file is written before exit
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
      USE control_flags,  ONLY : lneb, lsmd
      USE path_variables, ONLY : istep_path, suspended_image, reset_vel, &
                                 fixed_tan, ft_pos_av, ft_pos, path_length, &
                                 path_length_av
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
            ! ... back to real space
            !
            CALL to_real_space()
            !
            ! ... projected gradients are computed here
            !
            CALL smd_gradient()
            !
            ! ... the error is computed only when the run is not starting 
            ! ... from scratch
            !
            IF ( istep_path > 0 ) CALL compute_error()
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
      USE io_global,        ONLY : meta_ionode
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
         IF ( meta_ionode ) THEN
            !
            WRITE( UNIT = iunpath, FMT = final_fmt )
            !
            IF ( lneb ) &
               WRITE( UNIT = iunpath, &
                      FMT = '(/,5X,"neb: convergence achieved in ",I3, &
                             &     " iterations" )' ) istep_path
            IF ( lsmd ) &
               WRITE( UNIT = iunpath, &
                      FMT = '(/,5X,"smd: convergence achieved in ",I3, &
                             &     " iterations" )' ) istep_path
            !
         END IF
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
         IF ( meta_ionode ) THEN
            !
            WRITE( UNIT = iunpath, FMT = final_fmt )
            !         
            IF ( lneb ) &
               WRITE( UNIT = iunpath, &
                      FMT = '(/,5X,"neb: reached the maximum number of ", &
                             &     "steps")' )
            IF ( lsmd ) &
               WRITE( UNIT = iunpath, &
                      FMT = '(/,5X,"smd: reached the maximum number of ", &
                             &     "steps")' )
            !
         END IF
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
      USE path_variables, ONLY : first_last_opt, num_of_images, &
                                 lsteep_des, lquick_min, ldamped_dyn, &
                                 lmol_dyn, lbroyden, llangevin, istep_path
      USE path_opt_routines
      !
      IMPLICIT NONE
      !
      INTEGER :: image
      !
      IF ( lbroyden ) THEN
         !
         CALL broyden()
         !
      END IF
      !
      IF ( first_last_opt ) THEN
         !
         IF ( lsteep_des .OR. llangevin ) THEN
            !
            CALL steepest_descent( 1 )
            CALL steepest_descent( num_of_images )
            !
         ELSE IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            CALL velocity_Verlet_first_step( 1 )
            CALL velocity_Verlet_first_step( num_of_images )
            !
         END IF
         !
      END IF
      !
      DO image = 2, ( num_of_images - 1 )
         !
         IF ( lsteep_des .OR. llangevin ) THEN
            !
            CALL steepest_descent( image )
            !
         ELSE IF ( lquick_min .OR. ldamped_dyn .OR. lmol_dyn ) THEN
            !
            CALL velocity_Verlet_first_step( image )
            !
         END IF
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE first_opt_step
    !
    !------------------------------------------------------------------------
    SUBROUTINE second_opt_step()
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : first_last_opt, num_of_images, &
                                 lquick_min, ldamped_dyn, lmol_dyn
      USE path_opt_routines
      !
      IMPLICIT NONE
      !
      INTEGER :: image
      !
      !
      IF ( first_last_opt ) THEN
         !
         IF ( lquick_min ) THEN
            !
            CALL quick_min_second_step( 1 )
            CALL quick_min_second_step( num_of_images )
            !
         ELSE IF ( ldamped_dyn .OR. lmol_dyn ) THEN
            !
            CALL velocity_Verlet_second_step( 1 )
            CALL velocity_Verlet_second_step( num_of_images )
            !
         END IF
         !
      END IF
      !
      DO image = 2, ( num_of_images - 1 )
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
      END DO
      !
      RETURN
      !
    END SUBROUTINE second_opt_step    
    !
END MODULE path_base
