!
! Copyright (C) 2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_print_energies()
!----------------------------------------------------------------------------
!
! This routine is used for printing energy contrib from plugins
! DO NOT REMOVE THE TAGS ! ***ADDSON_NAME KIND_OF_PATCH***
!
USE io_global,        ONLY : stdout, ionode
USE kinds,            ONLY : DP
USE io_files,         ONLY : tmp_dir
!
USE plugin_flags
!
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
IMPLICIT NONE
!
! ***Environ VARIABLES BEGIN***
! ***Environ VARIABLES END***
!
! ***Environ CALLS BEGIN***
! ***Environ CALLS END***
!
END SUBROUTINE plugin_print_energies
