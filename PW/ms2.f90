MODULE ms2
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
#if defined(__MS2)
  ! This code is a safety net for people who accidentaly or due to
  ! curiosity try to compile pw.x with the __MS2 preprocessing symbol
  ! defined.
  !
  ! Since the use of error pragmas is not secure, this module will
  ! allow the user to reach the end of the compile process safely, and
  ! then present her/him a descriptive error at runtime.
  ! 
  ! Note that this measure doesn't prevent users from running ph.x (or
  ! other tools) with -D__MS2, which should produce no harm, since the
  ! net result will just be that MS2 related options in the
  ! configuration will be ignored.

  USE iso_c_binding
  USE mp_global,        ONLY : root, world_comm, mpime
  USE mp,               ONLY : mp_bcast, mp_barrier

  IMPLICIT NONE

  PUBLIC

  LOGICAL :: MS2_enabled = .FALSE.       ! Enable MS2

  CHARACTER(LEN=256) :: MS2_handler = ' ' ! Arguments to be passed to the MS2 transport

CONTAINS 

SUBROUTINE ms2_initialization()
  PRINT *, "*******************************************************************"
  PRINT *, "* This code was compiled with the __MS2 symbol defined manually   *"
  PRINT *, "* which is not the correct way to configure the MS2 system.       *"
  PRINT *, "*                                                                 *"
  PRINT *, "* In order to do so, you need to properly install the MS2 package *"
  PRINT *, "* and library, patch Quantum ESPRESSO and only then, compile      *"
  PRINT *, "*******************************************************************"
  STOP
END SUBROUTINE ms2_initialization

SUBROUTINE set_positions()
  PRINT *, "Placeholder for the MS2 project. This code shouldn't be used at all!"
END SUBROUTINE set_positions

SUBROUTINE ms2ec_add_esf()
  PRINT *, "Placeholder for the MS2 project. This code shouldn't be used at all!"
END SUBROUTINE ms2ec_add_esf

SUBROUTINE return_forces()
  PRINT *, "Placeholder for the MS2 project. This code shouldn't be used at all!"
END SUBROUTINE return_forces

#endif

END MODULE ms2
