! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE plugin_clean( prog, lflag )
!! This routine is used for cleaning calls from plugins.
!
! DO NOT REMOVE THE TAGS ! ***ADDSON_NAME KIND_OF_PATCH***
!
USE plugin_flags
!
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
LOGICAL, INTENT(IN) :: lflag
!
!
!
END SUBROUTINE plugin_clean

