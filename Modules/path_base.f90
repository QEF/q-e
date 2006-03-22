!
! Copyright (C) 2003-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!---------------------------------------------------------------------------
MODULE path_base
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all subroutines and functions needed for
  ! ... the implementation of "NEB" and "SMD" methods into the 
  ! ... PWSCF-FPMD-CPV codes
  !
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
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
      USE input_parameters,   ONLY : pos_ => pos, &
                                     restart_mode, calculation, opt_scheme, &
                                     climbing, nstep, input_images
      USE control_flags,      ONLY : conv_elec, lneb, lsmd, lcoarsegrained
      USE ions_base,          ONLY : nat, amass, ityp, if_pos
      USE constraints_module, ONLY : nconstr
      USE io_files,           ONLY : prefix, tmp_dir, path_file, dat_file, &
                                     int_file, xyz_file, axsf_file, broy_file
      USE cell_base,          ONLY : alat
      USE path_variables,     ONLY : climbing_ => climbing,                 &
                                     pos, istep_path, nstep_path, dim,      &
                                     num_of_images, pes, grad_pes, tangent, &
                                     error, path_length, path_thr, mass,    &
                                     deg_of_freedom, ds, use_masses,        &
                                     first_last_opt, frozen, llangevin,     &
                                     temp_req, use_freezing, k,             &
                                     tune_load_balance, grad, CI_scheme,    &
                                     posold, elastic_grad, k_min, k_max,    &
                                     num_of_modes, ft_pos, Nft, fixed_tan,  &
                                     ft_coeff, Nft_smooth, use_multistep,   &
                                     use_fourier, use_precond
      USE path_formats,       ONLY : summary_fmt
      USE mp_global,          ONLY : nimage
      USE io_global,          ONLY : meta_ionode
      USE path_io_routines,   ONLY : read_restart
      USE path_variables,     ONLY : path_allocation
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=2), INTENT(IN) :: prog
      !
      INTEGER                    :: i
      REAL(DP)                   :: inter_image_dist
      REAL(DP), ALLOCATABLE      :: d_R(:,:), image_spacing(:)
      CHARACTER(LEN=20)          :: nim_char, nstep_path_char
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      LOGICAL                    :: file_exists
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
      IF ( lcoarsegrained ) THEN
         !
         dim = nconstr
         !
         use_masses = .FALSE.
         !
      ELSE
         !
         dim = 3 * nat
         !
      END IF
      !
      IF ( nimage > 1 ) THEN
         !
         ! ... the automatic tuning of the load balance in 
         ! ... image-parallelisation is switched off by default
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
      IF ( lsmd ) THEN
         !
         IF ( use_fourier ) THEN
            !
            ! ... some coefficients for Fourier string dynamics
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
      END IF
      !  
      ! ... dynamical allocation of arrays
      !
      IF ( lneb ) THEN
         !
         CALL path_allocation( 'neb' )
         !
      ELSE IF ( lsmd ) THEN
         !
         CALL path_allocation( 'smd' )
         !
      END IF
      !
      IF ( use_masses ) THEN
         !
         ! ... mass weighted coordinates are used
         !
         DO i = 1, nat
            !
            mass(3*i-2) = amass(ityp(i))
            mass(3*i-1) = amass(ityp(i))
            mass(3*i-0) = amass(ityp(i))
            !
         END DO
         !
      ELSE
         !
         mass = 1.D0
         !
      END IF
      !
      ! ... initial path is read from file ( restart_mode == "restart" ) or
      ! ... generated from the input images ( restart_mode = "from_scratch" )
      ! ... It is always read from file in the case of "free-energy" 
      ! ... calculations
      !
      IF ( restart_mode == "restart" ) THEN
         !
         INQUIRE( FILE = path_file, EXIST = file_exists )
         !
         IF ( .NOT. file_exists ) restart_mode = "from_scratch"
         !
      END IF
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
            image_spacing(i) = norm( pos(:,i+1) - pos(:,i) )
            !
         END DO
         !
         path_length = SUM( image_spacing(:) )
         !
         inter_image_dist = SUM( image_spacing(:) ) / DBLE( num_of_images - 1 )
         !
         DEALLOCATE( image_spacing )
         !
      ELSE
         !
         CALL initial_guess()
         !
      END IF
      !
      ! ... initialization of the allocatable arrays
      !
      posold       = pos
      pes          = 0.D0
      grad_pes     = 0.D0
      elastic_grad = 0.D0
      tangent      = 0.D0
      grad         = 0.D0
      error        = 0.D0
      frozen       = .FALSE.
      !
      IF ( lneb ) THEN
         !
         k = k_min
         !
         climbing_ = climbing(1:num_of_images)
         !
      ELSE IF ( lsmd ) THEN
         !
         IF ( use_fourier ) THEN
            !
            ! ... fourier components of the path
            !
            ft_pos = 0.D0
            !
         END IF
         !
      END IF
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
         nstep_path_char = int_to_char( nstep_path )
         nim_char        = int_to_char( num_of_images )
         !
         WRITE( iunpath, * )
         WRITE( iunpath, summary_fmt ) "calculation",   TRIM( calculation )
         WRITE( iunpath, summary_fmt ) "restart_mode",  TRIM( restart_mode )
         WRITE( iunpath, summary_fmt ) "opt_scheme",    TRIM( opt_scheme )
         WRITE( iunpath, summary_fmt ) "num_of_images", TRIM( nim_char )
         WRITE( iunpath, summary_fmt ) "nstep",         TRIM( nstep_path_char )
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"first_last_opt",T35," = ",1X,L1))' ) first_last_opt
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"coarse-grained phase-space",T35, " = ",1X,L1))' ) &
               lcoarsegrained
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"use_freezing",T35," = ",1X,L1))' ) use_freezing
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"ds",T35," = ",1X,F6.4," a.u.")' ) ds
         !
         IF ( lneb ) THEN
            !
            WRITE( iunpath, summary_fmt ) "CI_scheme", TRIM( CI_scheme )
            !
            WRITE( UNIT = iunpath, &
                   FMT = '(5X,"k_max",T35," = ",1X,F6.4," a.u.")' ) k_max
            WRITE( UNIT = iunpath, &
                   FMT = '(5X,"k_min",T35," = ",1X,F6.4," a.u.")' ) k_min
            !
         END IF
         !
         IF ( lsmd ) THEN
            !
            WRITE( UNIT = iunpath, &
                FMT = '(5X,"use_multistep",T35," = ",1X,L1))' ) use_multistep
            !
            WRITE( UNIT = iunpath, &
                FMT = '(5X,"fixed_tan",T35," = ",1X,L1))' ) fixed_tan
            !
            IF ( llangevin ) &
               WRITE( UNIT = iunpath, &
                      FMT = '(5X,"required temperature",T35, &
                             &" = ",F6.1," K")' ) temp_req * eV_to_kelvin * au
            !
         END IF
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
          REAL(DP) :: s
          INTEGER  :: i, j
          !
          ! ... linear interpolation ( pos_ indicate the input coordinates)
          !
          ALLOCATE( image_spacing( input_images - 1 ) )
          ALLOCATE( d_R( dim,    ( input_images - 1 ) ) )
          !
          DO i = 1, ( input_images - 1 )
             !
             d_R(1:dim,i) = ( pos_(1:dim,i+1) - pos_(1:dim,i) )
             !
             image_spacing(i) = norm( d_R(:,i) )
             !
          END DO   
          !
          path_length = SUM( image_spacing(:) )
          !
          inter_image_dist = path_length / DBLE( num_of_images - 1  )
          !
          DO i = 1, ( input_images - 1 )
             !
             d_R(:,i) = d_R(:,i) / image_spacing(i)
             !
          END DO
          !
          pos(1:dim,1) = pos_(1:dim,1)
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
             pos(:,j) = pos_(1:dim,i) + s * d_R(:,i)
             !
          END DO
          !
          pos(:,num_of_images) = pos_(1:dim,input_images)
          !
          IF ( prog == 'PW' .AND. .NOT. lcoarsegrained ) THEN
             !
             ! ... coordinates must be in bohr ( pwscf uses alat units )
             !
             path_length = path_length * alat
             !
             inter_image_dist = inter_image_dist * alat
             !
             pos(:,:) = pos(:,:) * alat
             !
          END IF
          !
          DEALLOCATE( image_spacing, d_R )
          !
          RETURN
          !
        END SUBROUTINE initial_guess
        !
    END SUBROUTINE initialize_path
    !
    !-----------------------------------------------------------------------
    SUBROUTINE real_space_tangent( i )
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : pos, num_of_images, pes, tangent
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: i
      !
      REAL(DP) :: V_previous, V_actual, V_next
      REAL(DP) :: abs_next, abs_previous
      REAL(DP) :: delta_V_max, delta_V_min
      !
      !
      ! ... NEB definition of the tangent
      !
      IF ( i == 1 ) THEN
         !
         tangent(:,i) = pos(:,i+1) - pos(:,i)
         !
         RETURN
         !
      ELSE IF ( i == num_of_images ) THEN
         !
         tangent(:,i) = pos(:,i) - pos(:,i-1)
         !
         RETURN
         !
      END IF
      !
      V_previous = pes( i - 1 )
      V_actual   = pes( i )
      V_next     = pes( i + 1 )
      !
      IF ( ( V_next > V_actual ) .AND. ( V_actual > V_previous ) ) THEN
         !
         tangent(:,i) = pos(:,i+1) - pos(:,i)
         !
      ELSE IF ( ( V_next < V_actual ) .AND. ( V_actual < V_previous ) ) THEN
         !
         tangent(:,i) = pos(:,i) - pos(:,i-1)
         !
      ELSE
         !
         abs_next     = ABS( V_next     - V_actual ) 
         abs_previous = ABS( V_previous - V_actual ) 
         !
         delta_V_max = MAX( abs_next, abs_previous ) 
         delta_V_min = MIN( abs_next, abs_previous )
         !
         IF ( V_next > V_previous ) THEN
            !
            tangent(:,i) = ( pos(:,i+1) - pos(:,i) ) * delta_V_max + & 
                           ( pos(:,i) - pos(:,i-1) ) * delta_V_min
            !
         ELSE IF ( V_next < V_previous ) THEN
            !
            tangent(:,i) = ( pos(:,i+1) - pos(:,i) ) * delta_V_min + &
                           ( pos(:,i) - pos(:,i-1) ) * delta_V_max
            !
         ELSE
            !
            tangent(:,i) = pos(:,i+1) - pos(:,i-1)
            !
         END IF
         !
      END IF
      !
      tangent(:,i) = tangent(:,i) / norm( tangent(:,i) )
      !
      RETURN
      !
    END SUBROUTINE real_space_tangent
    !
    !------------------------------------------------------------------------
    SUBROUTINE elastic_constants()
      !------------------------------------------------------------------------
      ! 
      USE path_variables, ONLY : num_of_images, Emax, Emin, &
                                 k_max, k_min, k, pes
      !
      IMPLICIT NONE
      !
      INTEGER  :: i
      REAL(DP) :: delta_E
      REAL(DP) :: k_sum, k_diff
      !
      !
      ! ... standard neb ( with springs )
      !
      k_sum  = k_max + k_min
      k_diff = k_max - k_min
      !
      k(:) = k_min
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
      USE path_variables,    ONLY : pos, grad, elastic_grad, grad_pes, k, &
                                    num_of_images, climbing, mass, tangent
      !
      IMPLICIT NONE
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
            ! ... elastic gradient only along the path ( variable elastic
            ! ... consatnt is used ) NEB recipe
            !
            elastic_grad = tangent(:,i) * 0.5D0 * &
                       ( ( k(i) + k(i-1) ) * norm( pos(:,i) - pos(:,(i-1)) ) - &
                         ( k(i) + k(i+1) ) * norm( pos(:,(i+1)) - pos(:,i) ) )
            !
         END IF
         !
         ! ... total gradient on each image ( climbing image is used if needed )
         ! ... only the component of the pes gradient orthogonal to the path is 
         ! ... taken into account
         !
         grad(:,i) = grad_pes(:,i) / SQRT( mass(:) )
         !
         IF ( climbing(i) ) THEN
            !
            grad(:,i) = grad(:,i) - 2.D0 * tangent(:,i) * &
                                    ( grad(:,i) .dot. tangent(:,i) )
            ! 
         ELSE IF ( ( i > 1 ) .AND. ( i < num_of_images ) ) THEN
            !
            grad(:,i) = elastic_grad + grad(:,i) - &
                        tangent(:,i) * ( grad(:,i) .dot. tangent(:,i) )
            !
         END IF
         !
      END DO gradient_loop
      !
      RETURN
      !
    END SUBROUTINE neb_gradient
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fourier_tangent( i )
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : pos, ft_pos, num_of_modes, num_of_images, &
                                 tangent
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: i
      INTEGER             :: n
      REAL(DP)            :: s, pi_n
      !
      !
      s = DBLE( i - 1 ) / DBLE( num_of_images - 1 )
      !
      tangent(:,i) = ( pos(:,num_of_images) - pos(:,1) )
      !
      DO n = 1, num_of_modes
         !
         pi_n = pi * DBLE( n )
         !
         tangent(:,i) = tangent(:,i) + ft_pos(:,n) * pi_n * COS( pi_n * s )
         !
      END DO
      !
      tangent(:,i) = tangent(:,i) / norm( tangent(:,i) )
      !
      RETURN
      !
    END SUBROUTINE fourier_tangent
    !
    !-----------------------------------------------------------------------
    SUBROUTINE smd_gradient()
      !-----------------------------------------------------------------------
      !
      USE ions_base,         ONLY : if_pos
      USE path_variables,    ONLY : dim, mass, num_of_images, grad_pes, &
                                    tangent, llangevin, lang, grad, ds, &
                                    temp_req
      USE random_numbers,    ONLY : gauss_dist
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      !
      !
      ! ... we project pes gradients and gaussian noise
      !
      DO i = 1, num_of_images
         !
         IF ( llangevin ) THEN
            !
            ! ... the random term used in langevin dynamics is generated here
            !
            lang(:,i) = gauss_dist( 0.D0, SQRT( 2.D0*temp_req*ds ), dim )
            !
            lang(:,i) = lang(:,i) * DBLE( RESHAPE( if_pos, (/ dim /) ) )
            !
         END IF
         !
         grad(:,i) = grad_pes(:,i) / SQRT( mass(:) )
         !
         IF ( ( i > 1 ) .AND. ( i < num_of_images ) ) THEN
            !
            ! ... projection of the pes gradients 
            !
            grad(:,i) = grad(:,i) - &
                        tangent(:,i) * ( tangent(:,i) .dot. grad(:,i) )
            !
            IF ( llangevin ) THEN
               !
               lang(:,i) = lang(:,i) - &
                           tangent(:,i) * ( tangent(:,i) .dot. lang(:,i) )
               !
            END IF
            !
         END IF
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
    SUBROUTINE compute_tangent()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : fixed_tan, use_fourier, num_of_images, &
                                 istep_path
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      !
      !
      IF ( fixed_tan .AND. istep_path > 0 ) RETURN
      !
      DO i = 1, num_of_images
         !
         IF ( use_fourier ) THEN
            !
            CALL fourier_tangent( i )
            !
         ELSE
            !
            CALL real_space_tangent( i )
            !
         END IF
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE compute_tangent
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_error( err_out )
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : pos, posold, num_of_images, grad, &
                                 use_freezing, first_last_opt, path_thr, &
                                 error, frozen, lquick_min
      USE mp_global,      ONLY : nimage, inter_image_comm
      USE mp,             ONLY : mp_bcast
      USE io_global,      ONLY : meta_ionode, meta_ionode_id
      !
      IMPLICIT NONE
      !
      REAL(DP), OPTIONAL, INTENT(OUT) :: err_out
      !
      INTEGER  :: i
      INTEGER  :: N_in, N_fin, free_me, num_of_scf_images
      REAL(DP) :: err_max
      !
      !
      IF ( first_last_opt ) THEN
         !
         N_in  = 1
         N_fin = num_of_images
         !
         frozen = .FALSE.
         !
      ELSE
         !
         N_in  = 2
         N_fin = ( num_of_images - 1 )      
         !
         frozen = .FALSE.
         !
         ! ... the first and the last images are always frozen
         !
         frozen( N_in  - 1 ) = .TRUE.
         frozen( N_fin + 1 ) = .TRUE.
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
         frozen(N_in:N_fin) = ( error(N_in:N_fin) < &
                                MAX( 0.5D0 * err_max, path_thr ) )
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
               num_of_scf_images = COUNT( .NOT. frozen(N_in:N_fin) )
               !
               IF ( num_of_scf_images >= nimage ) EXIT find_scf_images
               !
               free_me = MAXLOC( error(N_in:N_fin), 1, frozen(N_in:N_fin) )
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
      IF ( use_freezing .AND. lquick_min ) THEN
         !
         ! ... the old positions of the frozen images are set to the present
         ! ... position (equivalent to resetting the velocity)
         !
         FORALL( i = N_in:N_fin, frozen(i) ) posold(:,i) = pos(:,i)
         !
      END IF
      !
      IF ( PRESENT( err_out ) ) err_out = err_max
      !
      RETURN
      !
    END SUBROUTINE compute_error
    !
    !------------------------------------------------------------------------
    SUBROUTINE fe_profile()
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : ni => num_of_images
      USE path_variables, ONLY : pos, pes, grad_pes, &
                                 Emin, Emax, Emax_index
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      !
      !
      pes(:) = 0.D0
      !
      DO i = 2, ni
         !
         pes(i) = pes(i-1) + 0.5D0 * ( ( pos(:,i) - pos(:,i-1) ) .dot. &
                                       ( grad_pes(:,i) + grad_pes(:,i-1) ) )
         !
      END DO
      !
      Emin       = MINVAL( pes(1:ni) )
      Emax       = MAXVAL( pes(1:ni) )
      Emax_index = MAXLOC( pes(1:ni), 1 )
      !
      RETURN
      !
    END SUBROUTINE fe_profile
    !
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_pes( stat )
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : num_of_images, suspended_image,  &
                                 istep_path, pes, first_last_opt, &
                                 Emin, Emax, Emax_index
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(OUT) :: stat
      INTEGER              :: N_in, N_fin
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
      Emin       = MINVAL( pes(1:num_of_images) )
      Emax       = MAXVAL( pes(1:num_of_images) )
      Emax_index = MAXLOC( pes(1:num_of_images), 1 )
      !
      RETURN
      !
    END SUBROUTINE born_oppenheimer_pes
    !
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_fes( stat )
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : num_of_images, suspended_image, &
                                 istep_path, first_last_opt
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(OUT) :: stat
      INTEGER              :: N_in, N_fin
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
      CALL compute_fes_grads( N_in, N_fin, stat )
      !
      RETURN
      !
    END SUBROUTINE born_oppenheimer_fes
    !
    !-----------------------------------------------------------------------
    SUBROUTINE check_domain()
      !-----------------------------------------------------------------------
      !
      USE path_variables,     ONLY : pos, num_of_images, &
                                     istep_path, first_last_opt
      USE constraints_module, ONLY : target
      USE metadyn_base,       ONLY : impose_domain_constraints
      !
      IMPLICIT NONE
      !
      INTEGER :: N_in, N_fin, i
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
      DO i = N_in, N_fin
         !
         target(:) = pos(:,i)
         !
         CALL impose_domain_constraints()
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE check_domain
    !
    !-----------------------------------------------------------------------
    SUBROUTINE search_mep()
      !-----------------------------------------------------------------------
      !
      USE control_flags,    ONLY : lneb, lsmd, lcoarsegrained
      USE path_variables,   ONLY : conv_path, istep_path, nstep_path,  &
                                   suspended_image, activation_energy, &
                                   err_max, num_of_modes, Nft, pes,    &
                                   climbing, CI_scheme, Emax_index,    &
                                   fixed_tan, use_fourier, pos, ft_pos
      USE path_io_routines, ONLY : write_restart, write_dat_files, write_output
      USE check_stop,       ONLY : check_stop_now
      USE io_global,        ONLY : meta_ionode
      USE path_formats,     ONLY : scf_iter_fmt
      USE path_reparametrisation
      !
      IMPLICIT NONE
      !
      INTEGER :: mode
      LOGICAL :: stat
      !
      REAL(DP), EXTERNAL :: get_clock
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
         IF ( check_stop_now() ) THEN
            !
            ! ... the programs checks if the user has required a soft
            ! ... exit or if if maximum CPU time has been exceeded
            !
            CALL write_restart()
            !
            conv_path = .FALSE.
            !
            RETURN
            !
         END IF
         !
         ! ... normalised tangent of the path
         !
         CALL compute_tangent()
         !
         ! ... energies and gradients acting on each image of the path (in real
         ! ... space) are computed calling a driver for the scf calculations
         !
         IF ( lcoarsegrained ) THEN
            !
            CALL born_oppenheimer_fes( stat )
            !
         ELSE
            !
            CALL born_oppenheimer_pes( stat )
            !
         END IF
         !
         IF ( .NOT. stat ) THEN
            !
            conv_path = .FALSE.
            !
            EXIT optimisation
            !
         END IF         
         !
         ! ... istep_path is updated after a self-consistency step has been
         ! ... completed
         !
         istep_path = istep_path + 1
         !
         IF ( lcoarsegrained ) CALL fe_profile()
         !
         IF ( lneb ) THEN
            !
            IF ( CI_scheme == "auto" ) THEN
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
         ! ... the error is computed here (frozen images are also set here)
         !
         CALL compute_error( err_max )
         !
         ! ... information is written on the files
         !
         CALL write_dat_files()
         !
         ! ... information is written on the standard output
         !
         CALL write_output()
         !
         ! ... the restart file is written
         !
         CALL write_restart()
         !
         ! ... exit conditions
         !
         IF ( check_exit( err_max ) ) EXIT optimisation
         !
         suspended_image = 0
         !
         ! ... if convergence is not yet achieved, the path is optimised
         !
         CALL optimisation_step()
         !
         IF ( lcoarsegrained ) CALL check_domain()
         !
         IF ( lsmd .AND. .NOT. fixed_tan ) THEN
            !
            IF ( use_fourier ) THEN
               !
               ! ... fourier components of the path
               !
               CALL to_reciprocal_space( pos, ft_pos )
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
               ! ... the new fourier components
               !
               CALL to_reciprocal_space( pos, ft_pos )
               !
            ELSE
               !
               CALL spline_reparametrisation()
               !
            END IF
            !
         END IF
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
      USE control_flags,  ONLY : lsmd
      USE path_variables, ONLY : istep_path, suspended_image, &
                                 use_fourier, pos, ft_pos
      USE path_reparametrisation
      !
      IMPLICIT NONE
      !
      !
      IF ( istep_path == 0 .OR. suspended_image /= 0 ) RETURN
      !
      IF ( lsmd ) THEN
         !
         IF ( use_fourier ) THEN
            !
            ! ... the fourier components of the path are computed here
            !
            CALL to_reciprocal_space( pos, ft_pos )
            !
            ! ... the path-length is computed here
            !
            CALL compute_path_length()
            !
            ! ... back to real space
            !
            CALL to_real_space()
            !
            ! ... the new fourier components
            !
            CALL to_reciprocal_space( pos, ft_pos )
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
                                   num_of_images, llangevin
      USE path_formats,     ONLY : final_fmt
      !
      IMPLICIT NONE
      !
      LOGICAL              :: check_exit
      REAL(DP), INTENT(IN) :: err_max
      LOGICAL              :: exit_condition
      !
      !
      check_exit = .FALSE.
      !
      ! ... the program checks if the convergence has been achieved
      !
      exit_condition = ( .NOT. llangevin .AND. & 
                         ( num_of_images == num_of_images_inp ) .AND. &
                         ( err_max <= path_thr ) )
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
    SUBROUTINE optimisation_step()
      !------------------------------------------------------------------------
      !
      USE path_variables,    ONLY : num_of_images, frozen, lsteep_des, &
                                    lquick_min, lbroyden, llangevin, istep_path
      USE path_opt_routines, ONLY : broyden, steepest_descent, quick_min
      !
      IMPLICIT NONE
      !
      INTEGER :: image
      !
      !
      IF ( lbroyden ) THEN
         !
         CALL broyden()
         !
      ELSE
         !
         DO image = 1, num_of_images
            !
            IF ( frozen(image) ) CYCLE
            !
            IF ( lsteep_des .OR. llangevin ) THEN
               !
               CALL steepest_descent( image )
               !
            ELSE IF ( lquick_min ) THEN
               !
               CALL quick_min( image, istep_path )
               !
            END IF
            !
         END DO
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE optimisation_step
    !
END MODULE path_base
