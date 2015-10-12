!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run()
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  !
  USE environment,        ONLY : environment_end
  USE control_flags,      ONLY : lconstrain
  USE constraints_module, ONLY : deallocate_constraint
  USE mp_global,          ONLY : mp_global_end
  !
  IMPLICIT NONE
  !
  !
  CALL environment_end( 'CP' )
  !
  CALL deallocate_modules_var()
  !
  IF ( lconstrain ) CALL deallocate_constraint()
  !
  CALL plugin_clean()
  !
  CALL mp_global_end()
  !
END SUBROUTINE stop_run

SUBROUTINE do_stop( flag )
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  !
  IF ( flag ) THEN
     STOP
  ELSE
     STOP 1
  END IF
  !
END SUBROUTINE do_stop
