!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stop_cp_run()
  !----------------------------------------------------------------------------
  !! Close all files and synchronize processes before stopping.
  !
  USE environment,        ONLY : environment_end
  USE control_flags,      ONLY : lconstrain
  USE constraints_module, ONLY : deallocate_constraint
  USE mp_global,          ONLY : mp_global_end
  !
#if defined (__ENVIRON)
  USE plugin_flags,        ONLY : use_environ
  USE environ_base_module, ONLY : clean_environ
#endif
  !
  IMPLICIT NONE
  !
  CALL environment_end( 'CP' )
  !
  CALL deallocate_modules_var()
  !
  IF ( lconstrain ) CALL deallocate_constraint()
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_clean()
#endif 
#if defined (__ENVIRON)
  IF (use_environ) CALL clean_environ()
#endif
  !
  CALL mp_global_end()
  !
  STOP 0
  !
END SUBROUTINE stop_cp_run
