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
         lfcpopt              = .FALSE.
         fcp_mu               = 0.0_DP
         fcp_relax_step       = 0.1_DP
         fcp_relax_crit       = 0.001_DP
         fcp_tot_charge_first = 0.0_DP
         fcp_tot_charge_last  = 0.0_DP
       !
       ! for reading ions namelist we need to set calculation=relax
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
       USE io_global, ONLY: ionode_id
       USE mp,        ONLY: mp_bcast
       USE mp_world,  ONLY: world_comm
       USE path_input_parameters_module
       !
       IMPLICIT NONE
       !
       ! ... "path" variables broadcast
       !
       CALL mp_bcast( restart_mode,         ionode_id, world_comm )
       CALL mp_bcast( string_method,        ionode_id, world_comm ) 
       CALL mp_bcast( num_of_images,        ionode_id, world_comm )
       CALL mp_bcast( first_last_opt,       ionode_id, world_comm )
       CALL mp_bcast( use_masses,           ionode_id, world_comm )
       CALL mp_bcast( use_freezing,         ionode_id, world_comm )
       CALL mp_bcast( fixed_tan,            ionode_id, world_comm )
       CALL mp_bcast( CI_scheme,            ionode_id, world_comm )
       CALL mp_bcast( opt_scheme,           ionode_id, world_comm )
       CALL mp_bcast( temp_req,             ionode_id, world_comm )
       CALL mp_bcast( ds,                   ionode_id, world_comm )
       CALL mp_bcast( k_max,                ionode_id, world_comm )
       CALL mp_bcast( k_min,                ionode_id, world_comm )
       CALL mp_bcast( path_thr,             ionode_id, world_comm )
       CALL mp_bcast( nstep_path,           ionode_id, world_comm )
       CALL mp_bcast( lfcpopt,              ionode_id, world_comm )
       CALL mp_bcast( fcp_mu,               ionode_id, world_comm )
       CALL mp_bcast( fcp_relax_step,       ionode_id, world_comm )
       CALL mp_bcast( fcp_relax_crit,       ionode_id, world_comm )
       CALL mp_bcast( fcp_tot_charge_first, ionode_id, world_comm )
       CALL mp_bcast( fcp_tot_charge_last,  ionode_id, world_comm )
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
