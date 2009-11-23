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
  USE mp,    ONLY : mp_end, mp_barrier
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  !
  !
  CALL print_clock_ph()
  !
  CALL mp_barrier()
  !
  CALL mp_end()
  !
  CALL deallocate_part()
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

CALL close_phq(.FALSE.)

CALL stop_ph(flag)

END SUBROUTINE stop_smoothly_ph
