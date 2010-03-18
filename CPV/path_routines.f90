!
! Copyright (C) 2002-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE path_routines
  !-----------------------------------------------------------------------
  !
  ! ... This module contains some interface subroutines needed for
  ! ... the NEB implementation into the FPMD code
  !
  ! ... Written by Carlo Sbraccia ( 2003-2006 )
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE constants,  ONLY : bohr_radius_angs, eV_to_kelvin
  !
  PRIVATE
  !
  PUBLIC :: iosys_path
  !   
  CONTAINS
    !    
    !------------------------------------------------------------------------
    SUBROUTINE iosys_path()
      !------------------------------------------------------------------------
      !
      USE input_parameters, ONLY : full_phs_path_flag, atomic_positions
      USE input_parameters, ONLY : pos, CI_scheme, opt_scheme, num_of_images, &
                                   first_last_opt, temp_req, ds, k_max,       &
                                   k_min, path_thr, restart_mode, nstep,      &
                                   calculation, use_freezing,                 &
                                   phase_space, ion_dynamics, etot_conv_thr,  &
                                   forc_conv_thr
      !
      USE path_variables, ONLY : lsteep_des, lquick_min, lbroyden, lbroyden2, &
                                 nstep_path, &
                                 num_of_images_  => num_of_images, &
                                 CI_scheme_      => CI_scheme, &
                                 first_last_opt_ => first_last_opt, &
                                 temp_req_       => temp_req, &
                                 ds_             => ds, &
                                 k_max_          => k_max, &
                                 k_min_          => k_min, &
                                 path_thr_       => path_thr, &
                                 use_freezing_   => use_freezing
      !
      USE io_files,      ONLY : prefix, tmp_dir, &
                                check_writable, delete_if_present
      USE io_global,     ONLY : ionode, ionode_id
      USE ions_base,     ONLY : nat
      USE cell_base,     ONLY : alat, a1, a2, a3
      USE mp_global,     ONLY : mpime
      USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
      USE control_flags, ONLY : lpath, lneb, lcoarsegrained, lconstrain, &
                                lmd, tprnfor
      USE metadyn_vars,  ONLY : init_metadyn_vars
      USE wrappers,      ONLY : f_mkdir
      !
      IMPLICIT NONE
      !
      INTEGER                     :: image, i, ia
      INTEGER                     :: ios
      REAL(DP), ALLOCATABLE       :: tau(:,:) 
      CHARACTER(LEN=256)          :: tmp_dir_neb
      CHARACTER(LEN=256)          :: filename
      CHARACTER (LEN=6), EXTERNAL :: int_to_char
      !
      tprnfor    = .TRUE.
      nstep_path = nstep
      !
      SELECT CASE( TRIM( phase_space ) )
      CASE( 'full' )
        !
        lcoarsegrained  = .FALSE.
        !
      CASE ( 'coarse-grained' )
        !
        lcoarsegrained  = .TRUE.
        !
      END SELECT
      !
      IF ( lcoarsegrained ) THEN
        !
        CALL init_metadyn_vars()
        !
        lmd        = .TRUE.
        lconstrain = .TRUE.
        !
        SELECT CASE( TRIM( ion_dynamics ) )
        CASE( 'verlet', 'damp' )
           !
           CONTINUE
           !
        CASE DEFAULT
           !
           CALL errore( 'iosys_path ', 'calculation=' // TRIM( calculation ) // &
                      & ': ion_dynamics=' // TRIM( ion_dynamics ) // &
                      & ' not supported', 1 )
           !
        END SELECT
        !
      END IF
      !
      IF ( num_of_images < 2 ) &
         CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                    & ': num_of_images must be at least 2', 1 )
      !
      IF ( ( CI_scheme /= "no-CI"  ) .AND. &
           ( CI_scheme /= "auto"   ) .AND. &
           ( CI_scheme /= "manual" ) ) THEN
         !
         CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                    & ': unknown CI_scheme', 1 )
         !
      END IF
      !
      ! ... initialization of logical variables
      !
      lsteep_des  = .FALSE.
      lquick_min  = .FALSE.
      lbroyden    = .FALSE.
      lbroyden2   = .FALSE.
      !
      SELECT CASE ( opt_scheme )
      CASE ( "sd" )
         !
         lsteep_des = .TRUE.
         !
      CASE ( "quick-min" )
         !
         lquick_min = .TRUE.
         !
      CASE( "broyden" )
         !
         lbroyden = .TRUE.
         !
      CASE( "broyden2" )
         !
         lbroyden2 = .TRUE.
         !
      CASE default
         !
         CALL errore( 'iosys', 'calculation = ' // TRIM( calculation ) // &
                    & ': unknown opt_scheme', 1 )
         !
      END SELECT
      !
      num_of_images_  = num_of_images
      CI_scheme_      = CI_scheme
      first_last_opt_ = first_last_opt
      temp_req_       = temp_req
      ds_             = ds
      k_max_          = k_max
      k_min_          = k_min
      path_thr_       = path_thr
      use_freezing_   = use_freezing
      !
      lpath      = .TRUE.
      lneb       = .TRUE.
      nstep_path = nstep
      !
      IF ( full_phs_path_flag ) THEN
         !
         ALLOCATE( tau( 3, nat ) )
         !
         DO image = 1, num_of_images
            !
            tau = RESHAPE( pos(1:3*nat,image), SHAPE( tau ) )
            !
            ! ... convert input atomic positions to internally used format:     
            !
            SELECT CASE ( TRIM( atomic_positions ) )
            CASE( 'alat' )
               !
               ! ... input atomic positions are divided by a0
               !
               tau(:,1:nat) = tau(:,1:nat) * alat
               !
            CASE( 'bohr' )
               !
               ! ... input atomic positions are in a.u.: do nothing
               !
               tau(:,1:nat) = tau(:,1:nat)
               !
            CASE( 'crystal' )
               !
               ! ... input atomic positions are in crystal axis ("scaled")
               !
               DO ia = 1, nat
                  !
                  DO i = 1, 3
                     !
                     tau(i,ia) = a1(i) * tau(1,ia) + &
                                 a2(i) * tau(2,ia) + &
                                 a3(i) * tau(3,ia)
                     !
                  END DO
                  !
               END DO
               !
            CASE( 'angstrom' )
               !
               ! ... atomic positions in Angstrom
               !
               tau(:,1:nat) = tau(:,1:nat) / bohr_radius_angs
               !
            CASE DEFAULT
               !
               CALL errore( 'iosys_path', ' tau_format = ' // &
                          & TRIM( atomic_positions ) // ' not implemented ', 1 )
               !
            END SELECT
            !
            pos(1:3*nat,image) = RESHAPE( tau, (/ 3 * nat /) )
            !
         END DO
         !
         DEALLOCATE( tau )
         !
      END IF
      !
      DO image = 1, num_of_images
        !
        ios = 0
        !
        tmp_dir_neb = TRIM( tmp_dir ) // "/" // TRIM( prefix ) // "_" // &
                      TRIM( int_to_char( image ) ) // '/'
        !
        ! ... a scratch directory for this image of the elastic band is
        ! ... created ( only by the master node )
        !
        IF ( ionode ) ios = f_mkdir( TRIM( tmp_dir_neb ) )
        !
        ! ... all jobs are syncronized
        !
        CALL mp_barrier()
        !
        ! ... each job checks whether the scratch directory is writable
        !
        ios = check_writable ( tmp_dir_neb, mpime )
        CALL mp_sum( ios )
        !
        IF ( ios /= 0 ) &
           CALL errore( 'outdir:', TRIM( tmp_dir_neb ) // &
                      & ' non existent or non writable', 1 )
        !
        ! ... if starting from scratch all temporary files are removed
        ! ... from outdir ( only by the master node )
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
           ! ... standard output of the self consistency is removed
           !
           IF ( ionode ) &
              CALL delete_if_present ( TRIM( tmp_dir_neb ) // 'CP.out' )
           !
        END IF
        !
      END DO
      !
      RETURN
      !
    END SUBROUTINE iosys_path
    !
END MODULE path_routines 
