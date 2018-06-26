!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plugin_check(calling_subroutine)
!
! This routine is used to raise an error if
! a plugin is activated
! DO NOT REMOVE THE TAGS ! ***ADDSON_NAME KIND_OF_PATCH***
!
USE plugin_flags
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: calling_subroutine
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
!Environ patch
IF (use_environ) CALL errore( calling_subroutine, 'Embedding environment cannot be used.', 1)
!Environ patch
! ***Environ CALLS END***
!
END SUBROUTINE plugin_check
