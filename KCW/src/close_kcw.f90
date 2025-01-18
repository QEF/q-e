!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
SUBROUTINE close_kcw ()
  !-----------------------------------------------------------------
  !
  !! This routine close KC buffers.
  !
  USE buffers,               ONLY : close_buffer
  USE control_kcw,           ONLY : iuwfc_wann, iurho_wann, calculation, iuwfc_wann_allk
  USE units_lr,              ONLY : iuwfc
  !
  IMPLICIT NONE
  !
  CALL close_buffer  ( iuwfc, 'delete' )
  !
  IF (calculation /= 'wann2kcw') CALL close_buffer  ( iurho_wann,'delete')
  IF (calculation == 'wann2kcw') CALL close_buffer  ( iuwfc_wann_allk,'delete')
  IF (calculation /= 'screen' ) CALL close_buffer  ( iuwfc_wann,'delete')
  !
END SUBROUTINE close_kcw
