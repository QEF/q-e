!
! Copyright (C) 2002-2004 PWSCF-FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE path_routines
  !-----------------------------------------------------------------------
  !
  ! ... This module contains some interface subroutines needed for
  ! ... the NEB implementation into the FPMD code
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : dbl
  USE constants,  ONLY : au, bohr_radius_angs, eV_to_kelvin
  USE path_base,  ONLY : initialize_path, search_mep
  !
  PRIVATE
  !
  PUBLIC :: iosys_path
  !   
  CONTAINS
    !    
    !-----------------------------------------------------------------------
    SUBROUTINE iosys_path()
      !-----------------------------------------------------------------------
      !
      USE path_variables, ONLY : lsteep_des, lquick_min , ldamped_dyn, &
                                 lmol_dyn, nstep_path, &
                                 num_of_images_  => num_of_images, &
                                 CI_scheme_      => CI_scheme, &
                                 first_last_opt_ => first_last_opt, &
                                 damp_           => damp, &
                                 temp_req_       => temp_req, &
                                 ds_             => ds, &
                                 k_max_          => k_max, &
                                 k_min_          => k_min, &
                                 path_thr_       => path_thr
      !
      USE input_parameters, ONLY : CI_scheme, minimization_scheme, &
                                   num_of_images, first_last_opt, damp, &
                                   temp_req, ds, k_max, k_min, path_thr
      USE input_parameters, ONLY : outdir, prefix, restart_mode, calculation
      USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc
      USE input_parameters, ONLY : rd_ht, trd_ht, cell_symmetry
      USE input_parameters, ONLY : nstep
      USE input_parameters, ONLY : max_seconds
      USE input_parameters, ONLY : ntyp, nat, na_inp, sp_pos, rd_pos, &
                                   atom_mass, atom_label, if_pos, rd_vel, &
                                   atomic_positions
      !
      USE io_global,     ONLY : ionode, ionode_id
      USE mp_global,     ONLY : mpime
      USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
      USE control_flags, ONLY : lpath, lneb
      USE parser,        ONLY : int_to_char
      USE cell_base,     ONLY : cell_base_init
      USE ions_base,     ONLY : ions_base_init
      USE check_stop,    ONLY : check_stop_init
      !
      IMPLICIT NONE
      !
      INTEGER            :: image
      INTEGER            :: ios
      CHARACTER(LEN=256) :: outdir_saved
      CHARACTER(LEN=256) :: filename
      REAL(dbl)          :: alat_
      !
      INTEGER, EXTERNAL :: C_MKDIR
      !
      IF ( num_of_images < 2 ) THEN
        CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': num_of_images must be at least 2', 1 )
      END IF
      !
      IF ( ( CI_scheme /= "no-CI"      ) .AND. &
           ( CI_scheme /= "highest-TS" ) .AND. &
           ( CI_scheme /= "all-SP"     ) .AND. &
           ( CI_scheme /= "manual"     ) ) THEN
         !
         CALL errore( ' iosys ', 'calculation=' // TRIM( calculation ) // &
                   & ': unknown CI_scheme', 1 )
         !
      END IF
      !
      lsteep_des  = .FALSE.
      lquick_min  = .FALSE.
      ldamped_dyn = .FALSE.
      lmol_dyn    = .FALSE.      
      !
      SELECT CASE ( minimization_scheme )
      CASE ( "sd" )
         !
         lsteep_des  = .TRUE.
         !
      CASE ( "quick-min" )
         !
         lquick_min  = .TRUE.
         !
      CASE ( "damped-dyn" )
         !
         ldamped_dyn = .TRUE.
         !
      CASE ( "mol-dyn" )
         !
         lmol_dyn    = .TRUE.
         !
         IF ( temp_req == 0 ) &
            WRITE( stdout,'(/,T2,"WARNING: tepm_req has not been set" )')
         !
         temp_req = temp_req / ( eV_to_kelvin * au )
         !
      CASE default
         !
         CALL errore( ' iosys ','calculation=' // TRIM( calculation ) // &
                   & ': unknown minimization_scheme', 1 )
         !
      END SELECT
      !
      CALL cell_base_init( ibrav, celldm, trd_ht, cell_symmetry, &
                           rd_ht, a, b, c, cosab, cosac, cosbc, alat_ )

      CALL ions_base_init( ntyp, nat, na_inp, sp_pos, rd_pos, rd_vel, &
                           atom_mass, atom_label, if_pos, atomic_positions  )
      !
      CALL check_stop_init( max_seconds )
      !
      num_of_images_   = num_of_images
      CI_scheme_       = CI_scheme
      first_last_opt_  = first_last_opt
      damp_            = damp
      temp_req_        = temp_req
      ds_              = ds
      k_max_           = k_max
      k_min_           = k_min
      path_thr_        = path_thr
      !
      lpath      = .TRUE.
      lneb       = .TRUE.
      nstep_path = nstep
      nstep      = 1000
      !
      outdir_saved = outdir
      !
      DO image = 1, num_of_images
        !
        ios = 0
        outdir = TRIM( outdir_saved ) // "/" // TRIM( prefix ) // "_" // &
                 TRIM( int_to_char( image ) ) // '/'
        !
        IF ( ionode ) THEN
           !
           ! ... a scratch directory for this image of the elastic band is
           ! ... created ( only by the master node )
           !
           WRITE( stdout, * ) 'Creating dir : ',TRIM( outdir )
           !
           ios = C_MKDIR( TRIM( outdir ), LEN_TRIM( outdir ) )
           !
        END IF
        !
        ! ... all jobs are syncronized
        !
        CALL mp_barrier()
        !
        ! ... each job checks whether the scratch directory is accessible
        ! ... or not
        !
        filename = TRIM( outdir ) // 'cp' // TRIM( int_to_char( mpime ) )
        !
        OPEN( UNIT = 4, FILE = TRIM( filename ) , &
              STATUS = 'UNKNOWN', FORM = 'UNFORMATTED', IOSTAT = ios )
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
        CALL mp_sum( ios )
        !
        IF ( ios /= 0 ) &
           CALL errore( 'outdir: ', TRIM( outdir ) // &
                      & ' non existent or non writable', 1 )
        !
        ! ... if starting from scratch all temporary files are removed
        ! ... from outdir ( only by the master node )
        !
        IF ( restart_mode == 'from_scratch' ) THEN
           !
           IF ( ionode ) THEN
              !
              ! ... standard output of the self consistency is removed
              !
              OPEN( UNIT = 4, FILE = TRIM( outdir ) // 'CP.out', &
                    STATUS = 'UNKNOWN' )
              CLOSE( UNIT = 4, STATUS = 'DELETE' )
              !
           END IF
           !
        END IF
        !
      END DO
      !
      outdir = outdir_saved
      !
      RETURN
      !
    END SUBROUTINE iosys_path
    !
END MODULE path_routines 
