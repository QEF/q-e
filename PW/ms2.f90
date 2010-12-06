MODULE ms2
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

#if defined(__MS2)

  USE iso_c_binding
  USE mp_global,        ONLY : root, world_comm, mpime
  USE mp,               ONLY : mp_bcast, mp_barrier

  IMPLICIT NONE

  PUBLIC

  LOGICAL :: MS2_enabled = .FALSE.       ! Enable MS2

  CHARACTER(LEN=256) :: MS2_handler = ' ' ! Arguments to be passed to the MS2 transport


CONTAINS 

SUBROUTINE ms2_initialization()
  PRINT *, "Placeholder for the MS2 project. This code shouldn't be used at all!"
END SUBROUTINE ms2_initialization

SUBROUTINE set_positions()
  PRINT *, "Placeholder for the MS2 project. This code shouldn't be used at all!"
END SUBROUTINE set_positions

SUBROUTINE return_forces()
  PRINT *, "Placeholder for the MS2 project. This code shouldn't be used at all!"
END SUBROUTINE return_forces

#endif

END MODULE ms2
