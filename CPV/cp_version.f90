!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE cp_version
        USE global_version, only : version_number
        IMPLICIT NONE
        SAVE
#if ! defined __G95
        INCLUDE 'cpver.h'
#else
        CHARACTER(LEN=70), PARAMETER :: version_date = 'Sat Jan 15 19:44:57 CET 2005'
#endif
      END MODULE cp_version
