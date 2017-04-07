!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------------
MODULE path_base
  !---------------------------------------------------------------------------
  !
  ! ... This module contains most of the subroutines and functions needed by
  ! ... the implementation of "NEB" and "SMD" methods into Quantum ESPRESSO
  !
  ! ... Other relevant files are:
  !
  ! ... path_variables.f90
  ! ... path_io_routines.f90
  ! ... path_opt_routines.f90
  ! ... path_reparametrisation.f90
  ! ... path_formats.f90
  ! ... compute_scf.f90
  !
  ! ... The code is based on the NEB algorithm described in :
  !
  ! ...  1) G. Henkelman, B.P. Uberuaga, and H. Jonsson;
  ! ...     J.Chem.Phys., 113, 9901, (2000)
  ! ...  2) G. Henkelman, and H. Jonsson; J.Chem.Phys., 113, 9978, (2000)
  !
  ! ... More details about the implementation can be found at
  !
  ! ...    http://www.sissa.it/cm/thesis/2005/sbraccia.pdf
  !
  ! ... Code written and maintained by Carlo Sbraccia ( 2003-2007 )
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps32, pi, autoev, bohr_radius_angs, eV_to_kelvin
  USE path_io_units_module,  ONLY : iunpath
  USE io_global, ONLY : meta_ionode, meta_ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
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
    SUBROUTINE initialize_path()
      !-----------------------------------------------------------------------
      !
      USE control_flags,    ONLY : conv_elec
      USE ions_base,        ONLY : amass, ityp
      USE io_files,         ONLY : prefix, tmp_dir
      USE mp_global,        ONLY : nimage
      USE path_input_parameters_module, ONLY : pos_      => pos, &
                                   climbing_ => climbing, &
                                   input_images, nstep_path_ => nstep_path
      USE path_input_parameters_module, ONLY : restart_mode
      USE path_input_parameters_module, ONLY : nat
      USE path_variables, ONLY : fix_atom_pos
      USE path_variables,   ONLY : climbing, pos, istep_path, nstep_path,    &
                                   dim1, num_of_images, pes, grad_pes, mass, &
                                   use_masses, tangent, error, path_length,  &
                                   deg_of_freedom, frozen, use_freezing, k,  &
                                   k_min, tune_load_balance, grad, posold,   &
                                   elastic_grad, pending_image, first_last_opt
      USE path_variables,   ONLY : path_allocation
      USE path_io_routines, ONLY : read_restart
      USE path_io_units_module, ONLY : path_file, dat_file, crd_file, &
                                   int_file, xyz_file, axsf_file, broy_file
      USE fcp_variables,        ONLY : lfcpopt
      USE fcp_opt_routines,     ONLY : fcp_opt_allocation
      !
      IMPLICIT NONE
      !
      INTEGER :: i, fii, lii
      LOGICAL :: file_exists
      !
      ! ... output files are set
      !
      path_file = TRIM( prefix ) // ".path"
      dat_file  = TRIM( prefix ) // ".dat"
      int_file  = TRIM( prefix ) // ".int"
      crd_file  = TRIM( prefix ) // ".crd"
      xyz_file  = TRIM( prefix ) // ".xyz"
      axsf_file = TRIM( prefix ) // ".axsf"
      !
      broy_file = TRIM( tmp_dir ) // TRIM( prefix ) // ".broyden"
      !
      ! ... istep_path is initialised to zero
      !
      istep_path    = 0
      pending_image = 0
      conv_elec     = .TRUE.
      !
      ! ... the dimension of all "path" arrays (dim1) is set here
      ! ... ( it corresponds to the dimension of the configurational space )
      !
      !
      dim1 = 3*nat
      !
      !
      IF ( nimage > 1 ) THEN
         !
         ! ... the automatic tuning of the load balance in
         ! ... image-parallelisation is switched off by default
         !
         tune_load_balance = .FALSE.
         !
         ! ... in the case of image-parallelisation the number of images
         ! ... to be optimised must be larger than nimage
         !
         IF ( first_last_opt ) THEN
            !
            fii = 1
            lii = num_of_images
            !
         ELSE
            !
            fii = 2
            lii = num_of_images - 1
            !
         END IF
         !
         IF ( nimage > ( lii - fii + 1 ) ) &
            CALL errore( 'initialize_path', 'nimage is ' // &
                       & 'larger than the available number of images', 1 )
         !
      END IF
      !
      ! ... dynamical allocation of arrays
      !
      CALL path_allocation()
      if ( lfcpopt ) CALL fcp_opt_allocation()
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
         mass = 1.0_DP
         !
      END IF
      !
      ! ... initialization of the allocatable arrays
      !
      pos(:,1:input_images) = pos_(1:dim1,1:input_images)
      !
      pes          = 0.0_DP
      grad_pes     = 0.0_DP
      elastic_grad = 0.0_DP
      tangent      = 0.0_DP
      grad         = 0.0_DP
      error        = 0.0_DP
      frozen       = .FALSE.
      !
      k = k_min
      !
      IF ( ALLOCATED( climbing_ ) ) THEN
         !
         climbing = climbing_
         !
      ELSE
         !
         climbing = .FALSE.
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
         IF ( meta_ionode ) THEN
            !
            INQUIRE( FILE = path_file, EXIST = file_exists )
            !
            IF ( .NOT. file_exists ) THEN
               !
               WRITE( iunpath, &
                      & '(/,5X,"restart file ''",A,"'' not found: ", &
                      &   /,5X,"starting from scratch")' )  TRIM( path_file )
               !
               restart_mode = "from_scratch"
               !
            END IF
            !
         END IF
         !
         CALL mp_bcast( restart_mode, meta_ionode_id, world_comm )
         !
      END IF
      !
      IF ( restart_mode == "restart" ) THEN
         !
         CALL read_restart()
         !
         ! ... consistency between the input value of nstep and the value
         ! ... of nstep_path read from the restart_file is checked
         !
         IF ( nstep_path_ == 0 ) THEN
            !
            istep_path = 0
            nstep_path = nstep_path_
            !
         END IF
         !
         IF ( nstep_path_ > nstep_path ) nstep_path = nstep_path_
         !
         ! ... in case first_last_opt has been set to true, reset the frozen
         ! ... array to false (all the images have to be optimized, at least
         ! ... on the first iteration)
         !
         IF ( first_last_opt ) frozen = .FALSE.
         !
         ! ... path length is computed here
         !
         path_length = 0.0_DP
         !
         DO i = 1, ( num_of_images - 1 )
            !
            path_length = path_length + norm( pos(:,i+1) - pos(:,i) )
            !
         END DO
         !
      ELSE
         !
         CALL initial_guess()
         !
         posold(:,:) = pos(:,:)
         !
      END IF
      !
      ! ... the actual number of degrees of freedom is computed
      !
      deg_of_freedom = 0
      !
      DO i = 1, nat
         !
         IF ( fix_atom_pos(1,i) == 1 ) deg_of_freedom = deg_of_freedom + 1
         IF ( fix_atom_pos(2,i) == 1 ) deg_of_freedom = deg_of_freedom + 1
         IF ( fix_atom_pos(3,i) == 1 ) deg_of_freedom = deg_of_freedom + 1
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE initialize_path
    !
    !--------------------------------------------------------------------
    SUBROUTINE initial_guess()
      !--------------------------------------------------------------------
      !
      ! ... linear interpolation
      !
      USE path_input_parameters_module, ONLY : input_images
      USE path_variables,   ONLY : pos, dim1, num_of_images, path_length
      USE path_io_units_module,         ONLY : iunpath
      !
      IMPLICIT NONE
      !
      REAL(DP) :: s
      INTEGER  :: i, j
      LOGICAL  :: tooclose
      REAL(DP), ALLOCATABLE :: pos_n(:,:), dr(:,:), image_spacing(:)
      !
      !
      IF ( meta_ionode ) THEN
         !
         ALLOCATE( pos_n( dim1, num_of_images ) )
         ALLOCATE( dr( dim1, input_images - 1 ) )
         ALLOCATE( image_spacing( input_images - 1 ) )
         !
         tooclose = .false.
         DO i = 1, input_images - 1
            !
            dr(:,i) = ( pos(:,i+1) - pos(:,i) )
            !
            image_spacing(i) = norm( dr(:,i) )
            tooclose = tooclose .OR. ( image_spacing(i) < 0.01 )
            !
         END DO
         !
         path_length = SUM( image_spacing(:) )
         !
         DO i = 1, input_images - 1
            !
            dr(:,i) = dr(:,i) / image_spacing(i)
            !
         END DO
         !
         pos_n(:,1) = pos(:,1)
         !
         i = 1
         s = 0.0_DP
         !
         DO j = 2, num_of_images - 1
            !
            s = s + path_length / DBLE( num_of_images - 1 )
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
               CALL errore( 'initialize_path', 'i >= input_images', i )
            !
            pos_n(:,j) = pos(:,i) + s * dr(:,i)
            !
         END DO
         !
         pos_n(:,num_of_images) = pos(:,input_images)
         !
         pos(:,:) = pos_n(:,:)
         !
         path_length = 0.0_DP
         !
         DO i = 1, num_of_images - 1
            !
            path_length = path_length + norm( pos(:,i+1) - pos(:,i) )
            !
         END DO
         !
         WRITE( UNIT = iunpath, &
                 FMT = '(/,5X,"initial path length",&
                        & T35," = ",F7.4," bohr")' ) path_length
         !
         WRITE( UNIT = iunpath, &
                FMT = '(5X,"initial inter-image distance",T35," = ",F7.4, &
                       &" bohr")' ) path_length / DBLE( num_of_images - 1 )
         !
         DEALLOCATE( image_spacing, dr, pos_n )
         !
      END IF
      !
      CALL mp_bcast( tooclose,    meta_ionode_id, world_comm )
      IF ( tooclose) CALL errore ('initial_guess', &
           ' something wrong: images are too close',1) 
      CALL mp_bcast( pos,         meta_ionode_id, world_comm )
      CALL mp_bcast( path_length, meta_ionode_id, world_comm )
      !
      RETURN
      !
    END SUBROUTINE initial_guess
    !
    !-----------------------------------------------------------------------
    FUNCTION real_space_tangent( i ) RESULT( rtan )
      !-----------------------------------------------------------------------
      !
      ! ... improved definition of the tangent (see JCP 113, 9978)
      !
      USE path_variables, ONLY : dim1, pos, num_of_images, pes
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: i
      REAL(DP)            :: rtan( dim1 )
      !
      REAL(DP) :: V_previous, V_actual, V_next
      REAL(DP) :: abs_next, abs_previous
      REAL(DP) :: delta_V_max, delta_V_min
      !
      !
      IF ( i == 1 ) THEN
         !
         rtan(:) = pos(:,i+1) - pos(:,i)
         !
         RETURN
         !
      ELSE IF ( i == num_of_images ) THEN
         !
         rtan(:) = pos(:,i) - pos(:,i-1)
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
         rtan(:) = pos(:,i+1) - pos(:,i)
         !
      ELSE IF ( ( V_next < V_actual ) .AND. ( V_actual < V_previous ) ) THEN
         !
         rtan(:) = pos(:,i) - pos(:,i-1)
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
            rtan(:) = ( pos(:,i+1) - pos(:,i) ) * delta_V_max + &
                      ( pos(:,i) - pos(:,i-1) ) * delta_V_min
            !
         ELSE IF ( V_next < V_previous ) THEN
            !
            rtan(:) = ( pos(:,i+1) - pos(:,i) ) * delta_V_min + &
                      ( pos(:,i) - pos(:,i-1) ) * delta_V_max
            !
         ELSE
            !
            rtan(:) = pos(:,i+1) - pos(:,i-1)
            !
         END IF
         !
      END IF
      !
      rtan(:) = rtan(:) / norm( rtan(:) )
      !
      RETURN
      !
    END FUNCTION real_space_tangent
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
            k(i) = 0.5_DP*( k_sum - k_diff * &
                           COS( pi * ( pes(i) - Emin ) / delta_E ) )
            !
         END DO
         !
      END IF
      !
      k(:) = 0.5_DP*k(:)
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
      IF ( meta_ionode ) THEN
         !
         CALL elastic_constants()
         !
         gradient_loop: DO i = 1, num_of_images
            !
            IF ( i > 1 .AND. i < num_of_images ) THEN
               !
               ! ... elastic gradient only along the path ( variable elastic
               ! ... consatnt is used ) NEB recipe
               !
               elastic_grad = tangent(:,i) * 0.5_DP * &
                       ( ( k(i) + k(i-1) ) * norm( pos(:,i) - pos(:,(i-1)) ) - &
                         ( k(i) + k(i+1) ) * norm( pos(:,(i+1)) - pos(:,i) ) )
               !
            END IF
            !
            ! ... total gradient on each image ( climbing image is used if
            ! ... required ) only the component of the pes gradient orthogonal
            ! ... to the path is used
            !
            grad(:,i) = grad_pes(:,i) / SQRT( mass(:) )
            !
            IF ( climbing(i) ) THEN
               !
               grad(:,i) = grad(:,i) - &
                           2.0_DP*tangent(:,i)*( grad(:,i) .dot. tangent(:,i) )
               !
            ELSE IF ( i > 1 .AND. i < num_of_images ) THEN
               !
               grad(:,i) = elastic_grad + grad(:,i) - &
                           tangent(:,i)*( grad(:,i) .dot. tangent(:,i) )
               !
            END IF
            !
         END DO gradient_loop
         !
      END IF
      !
      CALL mp_bcast( grad, meta_ionode_id, world_comm )
      !
      RETURN
      !
    END SUBROUTINE neb_gradient
    !
    !-----------------------------------------------------------------------
    SUBROUTINE smd_gradient()
      !-----------------------------------------------------------------------
      !
      USE ions_base,         ONLY : if_pos
      USE path_variables,    ONLY : dim1, mass, num_of_images, grad_pes, &
                                    tangent, llangevin, lang, grad, ds, &
                                    temp_req
      USE path_variables,    ONLY : climbing
      USE random_numbers,    ONLY : gauss_dist
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      !
      !
      IF ( meta_ionode ) THEN
         !
         grad(:,:) = 0.0_DP
         lang(:,:) = 0.0_DP
         !
         ! ... we project pes gradients and gaussian noise
         !
         DO i = 1, num_of_images
            !
            IF ( llangevin ) THEN
               !
               ! ... the random term used in langevin dynamics is generated here
               !
               lang(:,i) = gauss_dist( 0.0_DP, SQRT( 2.0_DP*temp_req*ds ), dim1 )
               !
               lang(:,i) = lang(:,i)*DBLE( RESHAPE( if_pos, (/ dim1 /) ) )
               !
            END IF
            !
            grad(:,i) = grad_pes(:,i) / SQRT( mass(:) )
            !
            IF ( climbing(i) ) THEN
               !
               grad(:,i) = grad(:,i) - &
                           2.0_DP*tangent(:,i)*( grad(:,i) .dot. tangent(:,i) )
               !
            ELSE IF ( i > 1 .AND. i < num_of_images ) THEN
               !
               ! ... projection of the pes gradients
               !
               grad(:,i) = grad(:,i) - &
                           tangent(:,i)*( grad(:,i) .dot. tangent(:,i) )
               !
               IF ( llangevin ) THEN
                  !
                  lang(:,i) = lang(:,i) - &
                              tangent(:,i)*( lang(:,i) .dot. tangent(:,i) )
                  !
               END IF
               !
            END IF
            !
         END DO
         !
      END IF
      !
      CALL mp_bcast( grad, meta_ionode_id, world_comm )
      CALL mp_bcast( lang, meta_ionode_id, world_comm )
      !
      RETURN
      !
    END SUBROUTINE smd_gradient
    !
    ! ... shared routines
    !
    !-----------------------------------------------------------------------
    FUNCTION new_tangent() RESULT( ntan )
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : dim1, num_of_images
      !
      IMPLICIT NONE
      !
      REAL(DP) :: ntan( dim1, num_of_images )
      !
      INTEGER :: i
      !
      !
      IF ( meta_ionode ) THEN
         !
         DO i = 1, num_of_images
            !
            ntan(:,i) = real_space_tangent( i )
            !
         END DO
         !
      END IF
      !
      CALL mp_bcast( ntan, meta_ionode_id, world_comm )
      !
      RETURN
      !
    END FUNCTION new_tangent
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_error( err_out )
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : pos, posold, num_of_images, grad, &
                                 use_freezing, first_last_opt, path_thr, &
                                 error, frozen, lquick_min
      USE mp_global,      ONLY : nimage
      !
      IMPLICIT NONE
      !
      REAL(DP), OPTIONAL, INTENT(OUT) :: err_out
      !
      INTEGER  :: i
      INTEGER  :: fii, lii, freed, num_of_scf_images
      REAL(DP) :: err_max
      LOGICAL  :: first
      !
      !
      IF ( first_last_opt ) THEN
         !
         fii  = 1
         lii = num_of_images
         !
         frozen = .FALSE.
         !
      ELSE
         !
         fii  = 2
         lii = num_of_images - 1
         !
         frozen = .FALSE.
         !
         ! ... the first and the last images are always frozen
         !
         frozen(fii-1) = .TRUE.
         frozen(lii+1) = .TRUE.
         !
      END IF
      !
      IF ( meta_ionode ) THEN
         !
         DO i = 1, num_of_images
            !
            ! ... the error is given by the largest component of the gradient
            ! ... vector ( PES + SPRINGS in the neb case )
            !
            error(i) = MAXVAL( ABS( grad(:,i) ) ) / bohr_radius_angs * autoev
            !
         END DO
         !
         err_max = MAXVAL( error(fii:lii), 1 )
         !
         IF ( use_freezing ) THEN
            !
            frozen(fii:lii) = ( error(fii:lii) < &
                                MAX( 0.5_DP*err_max, path_thr ) )
            !
         END IF
         !
         IF ( nimage > 1 .AND. use_freezing ) THEN
            !
            find_scf_images: DO
               !
               num_of_scf_images = COUNT( .NOT.frozen(fii:lii) )
               !
               IF ( num_of_scf_images >= nimage ) EXIT find_scf_images
               !
               first = .TRUE.
               !
               search: DO i = fii, lii
                  !
                  IF ( .NOT.frozen(i) ) CYCLE search
                  !
                  IF ( first ) THEN
                     !
                     first = .FALSE.
                     freed = i
                     !
                     CYCLE search
                     !
                  END IF
                  !
                  IF ( error(i) > error(freed) ) freed = i
                  !
               END DO search
               !
               frozen(freed) = .FALSE.
               !
            END DO find_scf_images
            !
         END IF
         !
         IF ( use_freezing .AND. lquick_min ) THEN
            !
            ! ... the old positions of the frozen images are set to the
            ! ... present position (equivalent to resetting the velocity)
            !
            FORALL( i = fii:lii, frozen(i) ) posold(:,i) = pos(:,i)
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( error,   meta_ionode_id, world_comm )
      CALL mp_bcast( err_max, meta_ionode_id, world_comm )
      CALL mp_bcast( frozen,  meta_ionode_id, world_comm )
      CALL mp_bcast( posold,  meta_ionode_id, world_comm )
      !
      IF ( PRESENT( err_out ) ) err_out = err_max
      !
      RETURN
      !
    END SUBROUTINE compute_error
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fcp_compute_error( err_out )
      !-----------------------------------------------------------------------
      !
      USE path_variables,   ONLY : num_of_images, first_last_opt
      USE fcp_variables,    ONLY : fcp_mu
      USE fcp_opt_routines, ONLY : fcp_neb_ef
      !
      IMPLICIT NONE
      !
      REAL(DP), OPTIONAL, INTENT(OUT) :: err_out
      !
      INTEGER  :: i
      INTEGER  :: fii, lii
      REAL(DP) :: err_max
      !
      !
      IF ( first_last_opt ) THEN
         !
         fii  = 1
         lii = num_of_images
         !
      ELSE
         !
         fii  = 2
         lii = num_of_images - 1
         !
      END IF
      !
      IF ( meta_ionode ) THEN
         !
         err_max = MAXVAL( ABS( fcp_mu - fcp_neb_ef(fii:lii) ), 1 )
         !
      END IF
      !
      CALL mp_bcast( err_max, meta_ionode_id, world_comm )
      !
      IF ( PRESENT( err_out ) ) err_out = err_max
      !
      RETURN
      !
    END SUBROUTINE fcp_compute_error
    !
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_pes( stat )
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : num_of_images, &
                                 pending_image, istep_path, pes, &
                                 first_last_opt, Emin, Emax, Emax_index
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(OUT) :: stat
      !
      INTEGER  :: fii, lii
      !
      !
      IF ( istep_path == 0 .OR. first_last_opt ) THEN
         !
         fii = 1
         lii = num_of_images
         !
      ELSE
         !
         fii = 2
         lii = num_of_images - 1
         !
      END IF
      !
      IF ( pending_image /= 0 ) fii = pending_image
      !
      CALL compute_scf( fii, lii, stat )
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
    SUBROUTINE fe_profile()
      !------------------------------------------------------------------------
      !
      USE path_variables, ONLY : num_of_images
      USE path_variables, ONLY : pos, pes, grad_pes, &
                                 Emin, Emax, Emax_index
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      !
      !
      pes(:) = 0.0_DP
      !
      DO i = 2, num_of_images
         !
         pes(i) = pes(i-1) + 0.5_DP*( ( pos(:,i) - pos(:,i-1) ) .dot. &
                                     ( grad_pes(:,i) + grad_pes(:,i-1) ) )
         !
      END DO
      !
      Emin       = MINVAL( pes(1:num_of_images) )
      Emax       = MAXVAL( pes(1:num_of_images) )
      Emax_index = MAXLOC( pes(1:num_of_images), 1 )
      !
      RETURN
      !
    END SUBROUTINE fe_profile
    !
    !-----------------------------------------------------------------------
    SUBROUTINE search_mep()
      !-----------------------------------------------------------------------
      !
      USE path_variables,    ONLY : lneb, lsmd
      USE path_variables,   ONLY : conv_path, istep_path, nstep_path,  &
                                   pending_image, activation_energy, &
                                   err_max, pes, climbing, CI_scheme,  &
                                   Emax_index, fixed_tan, tangent
      USE path_io_routines, ONLY : write_restart, write_dat_files, write_output
      USE path_formats,     ONLY : scf_iter_fmt
      USE fcp_variables,    ONLY : lfcpopt
      !
      USE path_reparametrisation
      !
      IMPLICIT NONE
      !
      LOGICAL :: stat
      REAL(DP) :: fcp_err_max = 0.0_DP
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
         pending_image = 0
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
         ! ... new positions are saved on file:  it has to be done here
         ! ... because, in the event of an unexpected crash the new positions
         ! ... would be lost. At this stage the forces and the energies are
         ! ... not yet known (but are not necessary for restarting); the
         ! ... restart file is written again as soon as the energies and
         ! ... forces have been computed.
         !
         CALL write_restart()
         !
         IF ( meta_ionode ) &
            WRITE( UNIT = iunpath, FMT = scf_iter_fmt ) istep_path + 1
         !
         ! ... energies and gradients acting on each image of the path (in real
         ! ... space) are computed calling a driver for the scf calculations
         !
         !
         CALL born_oppenheimer_pes( stat )
         !
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
         ! ... the new tangent is computed here :
         ! ... the improved definition of tangent requires the energies
         !
         IF ( .NOT. fixed_tan ) tangent(:,:) = new_tangent()
         !
         IF ( CI_scheme == "auto" ) THEN
            !
            climbing = .FALSE.
            !
            climbing(Emax_index) = .TRUE.
            !
         END IF
         !
         IF ( lneb ) CALL neb_gradient()
         IF ( lsmd ) CALL smd_gradient()
         !
         ! ... the forward activation energy is computed here
         !
         activation_energy = ( pes(Emax_index) - pes(1) )*autoev
         !
         ! ... the error is computed here (frozen images are also set here)
         !
         CALL compute_error( err_max )
         IF ( lfcpopt ) CALL fcp_compute_error( fcp_err_max )
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
         IF ( check_exit( err_max, fcp_err_max ) ) EXIT optimisation
         !
         ! ... if convergence is not yet achieved, the path is optimised
         !
         CALL optimisation_step()
         !
         IF ( lsmd ) CALL reparametrise()
         !
      END DO optimisation
      !
      ! ... the restart file is written before the exit
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
      USE path_variables,  ONLY : lsmd
      USE path_variables, ONLY : pending_image, tangent
      !
      USE path_reparametrisation
      !
      IMPLICIT NONE
      !
      !
      IF ( pending_image /= 0 ) RETURN
      !
      IF ( lsmd ) CALL reparametrise()
      !
      tangent(:,:) = new_tangent()
      !
      RETURN
      !
    END SUBROUTINE search_mep_init
    !
    !------------------------------------------------------------------------
    FUNCTION check_exit( err_max, fcp_err_max )
      !------------------------------------------------------------------------
      !
      USE path_input_parameters_module, ONLY : num_of_images_inp => num_of_images
      USE path_variables,    ONLY : lneb, lsmd
      USE path_variables,   ONLY : path_thr, istep_path, nstep_path, &
                                   conv_path, pending_image, &
                                   num_of_images, llangevin
      USE path_formats,     ONLY : final_fmt
      USE fcp_variables,    ONLY : lfcpopt, fcp_relax_crit
      !
      IMPLICIT NONE
      !
      LOGICAL              :: check_exit
      REAL(DP), INTENT(IN) :: err_max
      REAL(DP), INTENT(IN) :: fcp_err_max
      LOGICAL              :: exit_condition
      !
      !
      check_exit = .FALSE.
      !
      ! ... the program checks if the convergence has been achieved
      !
      exit_condition = ( .NOT.llangevin .AND. &
                         ( num_of_images == num_of_images_inp ) .AND. &
                         ( err_max <= path_thr ) )
      !
      IF ( lfcpopt .AND. fcp_err_max > fcp_relax_crit ) THEN
         exit_condition = .FALSE.
      END IF
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
         pending_image = 0
         !
         conv_path  = .TRUE.
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
         pending_image = 0
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
                                    lquick_min, lbroyden, lbroyden2, &
                                    llangevin, istep_path
      USE path_opt_routines, ONLY : quick_min, broyden, broyden2, &
                                    steepest_descent, langevin
      USE fcp_variables,     ONLY : lfcpopt
      USE fcp_opt_routines,  ONLY : fcp_line_minimisation
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
      ELSE IF (lbroyden2 ) THEN
         !
         CALL broyden2()
         !
      ELSE
         !
         DO image = 1, num_of_images
            !
            IF ( frozen(image) ) CYCLE
            !
            IF ( lsteep_des ) THEN
               !
               CALL steepest_descent( image )
               !
            ELSE IF ( llangevin ) THEN
               !
               CALL langevin( image )
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
      IF ( lfcpopt ) CALL fcp_line_minimisation()
      !
      RETURN
      !
    END SUBROUTINE optimisation_step
    !
END MODULE path_base
