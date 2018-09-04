!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE hp_close_q ( flag )
  !----------------------------------------------------------------------------
  !
  ! This subroutine closes all files.
  ! Called at the end of the run with flag=.TRUE. (removes 'recover')
  ! or during execution with flag=.FALSE. (does not remove 'recover')
  !
  USE buffers,        ONLY : close_buffer
  USE control_flags,  ONLY : twfcollect
  USE units_lr,       ONLY : iuwfc, iuatwfc
  USE ldaU_hp,        ONLY : iudwfc, iudvwfc
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  LOGICAL :: opnd
  !
  IF ( twfcollect ) THEN
     CALL close_buffer(iuwfc,'delete')
  ELSE
     CALL close_buffer(iuwfc,'keep')
  ENDIF
  !
  IF (flag) THEN
     CALL close_buffer(iudwfc,'delete')
     CALL close_buffer(iudvwfc,'delete')
  ELSE
     CALL close_buffer(iudwfc,'keep')
     CALL close_buffer(iudvwfc,'keep')
  ENDIF
  !
  CALL close_buffer(iuatwfc,'delete')
  !
  RETURN
  !
END SUBROUTINE hp_close_q
