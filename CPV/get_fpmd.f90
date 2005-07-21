!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This function is used to distinguish whether the
! executable we are running is cp.x or fpmd.x

  FUNCTION get_program()
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=4) :: get_program
    !
    get_program = 'FPMD'
    !
    RETURN
  END FUNCTION

