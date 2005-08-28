!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!=----------------------------------------------------------------------------=!
  MODULE green_functions
!=----------------------------------------------------------------------------=!

    USE kinds

    IMPLICIT NONE
    SAVE

    PRIVATE

    PUBLIC :: greenf

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

    FUNCTION greenf( x, rc )
      REAL(DP) :: greenf
      REAL(DP), INTENT(IN) :: x, rc
! ... declare external function
      REAL(DP) :: erf, erfc
      EXTERNAL erf, erfc

      greenf = erf( x / rc ) / x + erfc( x / rc ) / x
     
    END FUNCTION greenf

!=----------------------------------------------------------------------------=!
  END MODULE green_functions
!=----------------------------------------------------------------------------=!
