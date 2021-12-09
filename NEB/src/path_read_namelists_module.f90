!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE path_read_namelists_module
  !----------------------------------------------------------------------------
  !
  !  ... this module handles the reading of input namelists
  !  ... written by: Carlo Cavazzoni
  !  --------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE path_input_parameters_module
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  REAL(DP), PARAMETER :: fcp_not_set = 1.0E+99_DP
  !
  PUBLIC :: path_read_namelist
  !
  ! ... modules needed by read_xml.f90
  !
  !  ----------------------------------------------
  !
  CONTAINS
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist PATH
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE path_defaults( )
       !-----------------------------------------------------------------------
       !
       USE path_input_parameters_module
       !
       IMPLICIT NONE
       !
       !
       ! ... ( 'full' | 'coarse-grained' )
       !
       ! ... defaults for "path" optimisations variables
       !
       restart_mode  = 'from_scratch'
       string_method  = 'neb'
       num_of_images  = 0
       first_last_opt = .FALSE.
       use_masses     = .FALSE.
       use_freezing   = .FALSE.
       opt_scheme     = 'quick-min'
       temp_req       = 0.0_DP
       ds             = 1.0_DP
       path_thr       = 0.05_DP
       CI_scheme      = 'no-CI'
       k_max          = 0.1_DP
       k_min          = 0.1_DP
       fixed_tan      = .FALSE.
       nstep_path    = 1
       !
       ! ... defaults for "FCP" optimisations variables
       !
       lfcp         = .FALSE.
       fcp_mu       = fcp_not_set
       fcp_thr      = 0.01_DP
       fcp_scheme   = 'lm'
       fcp_ndiis    = 4
       fcp_rdiis    = 1.0_DP
       fcp_max_volt = 1.0_DP
       !
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist NEB
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE path_bcast()
       !-----------------------------------------------------------------------
       !
       USE io_global, ONLY: meta_ionode_id
       USE mp,        ONLY: mp_bcast
       USE mp_world,  ONLY: world_comm
       USE path_input_parameters_module
       !
       IMPLICIT NONE
       !
       ! ... "path" variables broadcast
       !
       CALL mp_bcast( restart_mode,    meta_ionode_id, world_comm )
       CALL mp_bcast( string_method,   meta_ionode_id, world_comm ) 
       CALL mp_bcast( num_of_images,   meta_ionode_id, world_comm )
       CALL mp_bcast( first_last_opt,  meta_ionode_id, world_comm )
       CALL mp_bcast( use_masses,      meta_ionode_id, world_comm )
       CALL mp_bcast( use_freezing,    meta_ionode_id, world_comm )
       CALL mp_bcast( fixed_tan,       meta_ionode_id, world_comm )
       CALL mp_bcast( CI_scheme,       meta_ionode_id, world_comm )
       CALL mp_bcast( opt_scheme,      meta_ionode_id, world_comm )
       CALL mp_bcast( temp_req,        meta_ionode_id, world_comm )
       CALL mp_bcast( ds,              meta_ionode_id, world_comm )
       CALL mp_bcast( k_max,           meta_ionode_id, world_comm )
       CALL mp_bcast( k_min,           meta_ionode_id, world_comm )
       CALL mp_bcast( path_thr,        meta_ionode_id, world_comm )
       CALL mp_bcast( nstep_path,      meta_ionode_id, world_comm )
       !
       ! ... "FCP" variables broadcast
       !
       CALL mp_bcast( lfcp,            meta_ionode_id, world_comm )
       CALL mp_bcast( fcp_mu,          meta_ionode_id, world_comm )
       CALL mp_bcast( fcp_thr,         meta_ionode_id, world_comm )
       CALL mp_bcast( fcp_scheme,      meta_ionode_id, world_comm )
       CALL mp_bcast( fcp_ndiis,       meta_ionode_id, world_comm )
       CALL mp_bcast( fcp_rdiis,       meta_ionode_id, world_comm )
       CALL mp_bcast( fcp_max_volt,    meta_ionode_id, world_comm )
       !
       RETURN
       !
     END SUBROUTINE
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE path_checkin( )
       !-----------------------------------------------------------------------
       !
       USE path_input_parameters_module
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=20) :: sub_name = ' path_checkin '
       INTEGER           :: i
       LOGICAL           :: allowed = .FALSE.
       !
       !
       ! ... general "path" variables checkin
              IF ( ds < 0.0_DP ) &
          CALL errore( sub_name,' ds out of range ',1)
       IF ( temp_req < 0.0_DP ) &
          CALL errore( sub_name,' temp_req out of range ',1)
       !
       allowed = .FALSE.
       DO i = 1, SIZE( opt_scheme_allowed )
          IF ( TRIM( opt_scheme ) == &
               opt_scheme_allowed(i) ) allowed = .TRUE.
       END DO
       IF ( .NOT. allowed ) &
          CALL errore( sub_name, ' opt_scheme '''// &
                     & TRIM( opt_scheme )//''' not allowed ', 1 )
       !
       !
       ! ... NEB(SM) specific checkin
       !
       IF ( k_max < 0.0_DP )  CALL errore( sub_name, 'k_max out of range', 1 )
       IF ( k_min < 0.0_DP )  CALL errore( sub_name, 'k_min out of range', 1 )
       IF ( k_max < k_min ) CALL errore( sub_name, 'k_max < k_min', 1 )
       !
!       IF ( nstep_path < 1 ) CALL errore ( sub_name, 'step_path out of range', 1 )
       !
       allowed = .FALSE.
       DO i = 1, SIZE( CI_scheme_allowed )
          IF ( TRIM( CI_scheme ) == CI_scheme_allowed(i) ) allowed = .TRUE.
       END DO
       !
       IF ( .NOT. allowed ) &
          CALL errore( sub_name, ' CI_scheme ''' // &
                      & TRIM( CI_scheme ) //''' not allowed ', 1 )
       !
       !
       ! ... FCP algorithm
       !
       IF ( lfcp ) THEN
          !
          IF( fcp_mu == fcp_not_set ) &
             CALL errore( sub_name,' fcp_mu is not set ', 1 )
          !
          IF ( fcp_thr <= 0.0_DP ) &
             CALL errore( sub_name, 'fcp_thr out of range', 1 )
          !
          allowed = .FALSE.
          DO i = 1, SIZE( fcp_scheme_allowed )
             IF ( TRIM( fcp_scheme ) == fcp_scheme_allowed(i) ) allowed = .TRUE.
          END DO
          !
          IF ( .NOT. allowed ) &
             CALL errore( sub_name, ' fcp_scheme ''' // &
                         & TRIM( fcp_scheme ) //''' not allowed ', 1 )
          !
          IF ( fcp_ndiis < 1 ) &
             CALL errore( sub_name, 'fcp_ndiis out of range', 1 )
          !
          IF ( fcp_rdiis <= 0.0_DP ) &
             CALL errore( sub_name, 'fcp_rdiis out of range', 1 )
          !
          IF ( fcp_max_volt <= 0.0_DP ) &
             CALL errore( sub_name, 'fcp_max_volt out of range', 1 )
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Namelist parsing main routine
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE path_read_namelist(unit)
       !-----------------------------------------------------------------------
       !
       !  this routine reads data from standard input and puts them into
       !  module-scope variables (accessible from other routines by including
       !  this module, or the one that contains them)
       !  ----------------------------------------------
       !
       ! ... declare modules
       !
       USE io_global, ONLY : ionode, ionode_id
       USE mp,        ONLY : mp_bcast
       USE mp_world,  ONLY : world_comm
       !
       IMPLICIT NONE
       !
       ! ... declare variables
       !
       INTEGER, intent(in) :: unit
       !
       !
       ! ... declare other variables
       !
       INTEGER :: ios
       !
       ! ... end of declarations
       !
       !  ----------------------------------------------
       !
       !
       ! ... default settings for all namelists
       !
       CALL path_defaults( )
       !
       ! ... Here start reading standard input file
       !
       ! ... PATH namelist
       !
       ios = 0
       IF ( ionode ) THEN
          !
          READ( unit, path, iostat = ios )
          !
       END IF
       CALL mp_bcast( ios, ionode_id, world_comm )
       IF( ios /= 0 ) THEN
          CALL errore( ' path_read_namelists ', &
                     & ' reading namelist path ', ABS(ios) )
       END IF
       !
       CALL path_bcast( )
       CALL path_checkin( )
       !
       RETURN
       !
     END SUBROUTINE path_read_namelist
     !
END MODULE path_read_namelists_module
