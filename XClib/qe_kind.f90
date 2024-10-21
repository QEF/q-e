!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
MODULE kind_l
!------------------------------------------------------------------------
!! Double precision in xc_lib.
!
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
  PUBLIC :: DP
  !
END MODULE      
