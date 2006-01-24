!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE cp_version
        USE global_version, only : version_number
        IMPLICIT NONE
        SAVE
        INCLUDE 'cpver.h'
      END MODULE cp_version
