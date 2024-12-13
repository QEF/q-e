!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE clean_pw_kcw( )
  !-----------------------------------------------------------------------
  !
  !! This routine deallocate all the variables of pwscf and of the
  !! screen code.
  !
  USE control_kcw,     ONLY : iudvwfc, setup_pw
  USE units_lr,        ONLY : iudwf, iuwfc
  USE buffers,         ONLY : close_buffer

  USE lr_symm_base,    ONLY : nsymq
  !
  IMPLICIT NONE
  !
  IF( setup_pw ) CALL clean_pw( .FALSE. )
  CALL kcw_deallocate_q ()
  nsymq=0
  !
  ! ... Close the files
  !
  CALL close_buffer(iuwfc, 'delete')
  CALL close_buffer(iudwf, 'delete')
  CALL close_buffer(iudvwfc,'delete')
  !
RETURN
END SUBROUTINE clean_pw_kcw
