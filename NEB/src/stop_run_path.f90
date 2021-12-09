!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run_path( lflag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  ! ... Called at the end of the run with flag = .TRUE. (removes 'restart')
  ! ... or during execution with flag = .FALSE. (does not remove 'restart')
  !
  USE io_global,          ONLY : ionode, stdout
  USE mp_global,          ONLY : mp_global_end
  USE environment,        ONLY : environment_end
  USE path_variables,     ONLY : path_deallocation
  USE path_io_units_module,      ONLY : iunpath
  USE fcp_variables,             ONLY : lfcp, fcp_deallocation
  USE gcscf_variables,           ONLY : lgcscf, gcscf_deallocation
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
  LOGICAL             :: exst, opnd, flag2
  !
  !
  CALL close_files(lflag)
  !
  ! as environment_end writes final info on stdout
  ! stdout has to be redirected to iunpath (neb output)
  !
  stdout=iunpath
  !
  CALL environment_end( 'NEB' )
  !
  CALL clean_pw( .TRUE. )
  !
  CALL path_deallocation()
  IF ( lfcp ) CALL fcp_deallocation()
  IF ( lgcscf ) CALL gcscf_deallocation()
  !
  CALL mp_global_end()
  !
  IF ( .not. lflag ) THEN
     !
     STOP 1
     !
  END IF
  !
END SUBROUTINE stop_run_path
!
!----------------------------------------------------------------------------
