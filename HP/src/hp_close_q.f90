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
  USE io_files,       ONLY : iunhub
  USE units_lr,       ONLY : iuwfc, iuatswfc, iudwf
  USE ldaU_hp,        ONLY : iudvwfc
  USE control_lr,     ONLY : lgamma
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  LOGICAL :: opnd
  !
  CALL close_buffer(iuwfc,'delete')
  !
  IF (flag) THEN
     CALL close_buffer(iudwf,'delete')
     CALL close_buffer(iudvwfc,'delete')
  ELSE
     CALL close_buffer(iudwf,'keep')
     CALL close_buffer(iudvwfc,'keep')
  ENDIF
  !
  CALL close_buffer(iuatswfc,'delete')
  IF (lgamma) CALL close_buffer(iunhub,'delete')
  !
  RETURN
  !
END SUBROUTINE hp_close_q
