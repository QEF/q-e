!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_ph( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Synchronize processes before stopping.
  !
  USE kinds, ONLY : DP
  USE mp_global, ONLY : mp_global_end
  USE mp_images, ONLY : nimage
  USE ph_restart,      ONLY : destroy_status_run
  USE save_ph,         ONLY : clean_input_variables
  USE environment,        ONLY : environment_end
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  !
  CALL clean_input_variables()
  !
  CALL destroy_status_run()
  !
  CALL deallocate_part()
  !
  CALL print_clock_ph()
  !
  CALL environment_end('PHONON')
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
  ENDIF
  !
END SUBROUTINE stop_ph

SUBROUTINE stop_smoothly_ph(flag)
IMPLICIT NONE
LOGICAL, INTENT(IN) :: flag

CALL collect_grid_files()

CALL close_phq(.FALSE.)

CALL stop_ph(flag)

END SUBROUTINE stop_smoothly_ph
