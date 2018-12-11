!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_clean_q (flag)
  !-----------------------------------------------------------------------
  !
  ! This routine deallocates the variables of PWscf and of the
  ! HP code, and resets the same variables as after reading input in
  ! hp_readin, so that it is possible to start a calculation at a new q.
  !
  USE control_flags,   ONLY : twfcollect
  USE lr_symm_base,    ONLY : nsymq
  !
  IMPLICIT NONE
  LOGICAL :: flag
  !
  twfcollect = .FALSE.
  !
  CALL clean_pw(.FALSE.)
  !
  ! Deallocate the arrays
  !
  CALL hp_dealloc_q()
  !
  nsymq = 0
  !
  ! Close the files
  !
  CALL hp_close_q (flag)
  !
  RETURN
  !
END SUBROUTINE hp_clean_q
