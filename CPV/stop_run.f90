!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  !
  USE io_global,          ONLY : stdout, ionode
  USE control_flags,      ONLY : lpath, lconstrain, program_name
  USE io_files,           ONLY : prefix
  USE environment,        ONLY : environment_end
  USE image_io_routines,   ONLY : io_image_stop
  USE constraints_module, ONLY : deallocate_constraint
  USE mp_global,          ONLY : mp_global_end
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  LOGICAL             :: exst
  !
  !
  CALL environment_end( program_name )
  !
  IF ( lpath ) CALL io_image_stop()
  !
  CALL deallocate_modules_var()
  !
  IF ( lconstrain ) CALL deallocate_constraint()
  !
  CALL mp_global_end()
  !
  IF ( flag ) THEN
     !
     STOP
     !
  ELSE
     !
     STOP 1
     !
  END IF
  !
END SUBROUTINE stop_run
