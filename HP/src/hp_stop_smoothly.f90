!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE hp_stop_smoothly(flag)
  !----------------------------------------------------------------------
  !
  ! Deallocate dynamical arrays, close files, and stop.
  !
  USE environment,    ONLY : environment_end
  USE mp_global,      ONLY : mp_global_end 
  USE ldaU_hp,        ONLY : code
  !
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: flag
  !
  CALL hp_clean_q (.FALSE.)
  !
  ! Deallocate some arrays
  !
  CALL hp_dealloc_1()
  CALL hp_dealloc_2()
  !
  ! Print clocks
  !
  IF (flag) THEN
     CALL print_clock_pw()
     CALL hp_print_clock()
  ENDIF
  !
  CALL environment_end(code)
  !
#if defined (__MPI)
  CALL mp_global_end()
#endif
  !
  STOP 1 
  !
END SUBROUTINE hp_stop_smoothly
