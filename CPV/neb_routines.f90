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
  USE kinds,      ONLY : DP, dbl
  USE constants,  ONLY : au, bohr_radius_angs, eV_to_kelvin
  USE neb_base,   ONLY : initialize_neb, compute_action, compute_tangent, &
                         elastic_constants, gradient, search_stationary_points, &
                         compute_error, path_tangent_, born_oppenheimer_PES, &
                         search_mep

  !
  PRIVATE

  !
  PUBLIC :: iosys_neb
  PUBLIC :: initialize_neb
  PUBLIC :: search_mep
  !   
  CONTAINS
    !
    !    
    !-----------------------------------------------------------------------
    SUBROUTINE iosys_neb()
      !-----------------------------------------------------------------------
      USE neb_variables, ONLY : lsteep_des, lquick_min , ldamped_dyn, lmol_dyn, &
                            num_of_images_  => num_of_images, &
                            CI_scheme_      => CI_scheme, &
                            optimization_   => optimization, &
                            damp_           => damp, &
                            temp_req_       => temp_req, &
                            ds_             => ds, &
                            k_max_          => k_max, &
                            k_min_          => k_min, &
                            neb_thr_        => neb_thr, &
                            nstep_neb 
      !
      USE input_parameters, ONLY : CI_scheme, minimization_scheme, &
                            num_of_images, optimization, damp, temp_req, &
                            ds, k_max, k_min, neb_thr
      !
      USE input_parameters, ONLY : outdir, prefix, restart_mode, calculation
      USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc
      USE input_parameters, ONLY : rd_ht, trd_ht, cell_symmetry
      USE input_parameters, ONLY : nstep, max_seconds
      USE input_parameters, ONLY : ntyp , nat , na_inp , sp_pos , rd_pos , atom_mass, &
             atom_label, if_pos, rd_vel, atomic_positions
      !
      USE io_global, ONLY: ionode, ionode_id
      USE mp_global, ONLY: mpime
      USE mp,        ONLY: mp_bcast, mp_barrier, mp_sum
      USE parser, ONLY: int_to_char
      USE cell_base, ONLY: cell_base_init
      USE ions_base, ONLY: ions_base_init
      USE check_stop, ONLY: check_stop_init


      IMPLICIT NONE

      INTEGER :: image
      INTEGER :: ios
      CHARACTER(LEN=256) :: outdir_saved
      CHARACTER(LEN=256) :: filename

      REAL(dbl) :: alat_

      INTEGER, EXTERNAL :: C_MKDIR

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
      SELECT CASE ( minimization_scheme )
      !
      CASE ( "sd" )
         lsteep_des  = .TRUE.
         lquick_min  = .FALSE.
         ldamped_dyn = .FALSE.
         lmol_dyn    = .FALSE.
      CASE ( "quick-min" )
         lsteep_des  = .FALSE.
         lquick_min  = .TRUE.
         ldamped_dyn = .FALSE.
         lmol_dyn    = .FALSE.
      CASE ( "damped-dyn" )
         lsteep_des  = .FALSE.
         lquick_min  = .FALSE.
         ldamped_dyn = .TRUE.
         lmol_dyn    = .FALSE.
      CASE ( "mol-dyn" )
         lsteep_des  = .FALSE.
         lquick_min  = .FALSE.
         ldamped_dyn = .FALSE.
         lmol_dyn    = .TRUE.
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
 
      CALL cell_base_init( ibrav, celldm, trd_ht, cell_symmetry, rd_ht, &
             a, b, c, cosab, cosac, cosbc, alat_ )

      CALL ions_base_init( ntyp , nat , na_inp , sp_pos , rd_pos , rd_vel, atom_mass, &
             atom_label, if_pos, atomic_positions  )

      CALL check_stop_init( max_seconds )

      num_of_images_ = num_of_images
      CI_scheme_     = CI_scheme
      optimization_  = optimization
      damp_          = damp
      temp_req_      = temp_req
      ds_            = ds
      k_max_         = k_max
      k_min_         = k_min
      neb_thr_       = neb_thr


      nstep_neb = nstep
      nstep     = 100

      outdir_saved = outdir
      !
      DO image = 1, num_of_images
        !
        ios = 0
        outdir = TRIM( outdir_saved ) // "/" // TRIM( prefix ) //"_" // &
                 TRIM( int_to_char( image ) ) // '/'
        WRITE(*,*) 'CREATING:', outdir
        !
        IF ( ionode ) THEN
           !
           ! ... a scratch directory for this image of the elastic band is
           ! ... created ( only by the master node )
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

      RETURN
    END SUBROUTINE
    !    
    !
END MODULE neb_routines 
