!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
MODULE util_param
!------------------------------------------------------------------------------!
!! 
!! This module is a duplication of the Modules/kind.f90 one, placed here fore 
!! convenience. 
!!   
USE parallel_include
!
CHARACTER(LEN = 5), PARAMETER :: crash_file = 'CRASH'
INTEGER, PARAMETER :: DP = selected_real_kind(14, 200)
INTEGER, PARAMETER :: i8b = selected_int_kind(18)
INTEGER, PARAMETER :: stdout = 6    ! unit connected to standard output
!
REAL(DP), PARAMETER :: eps14 = 1.0E-14_DP
REAL(DP), PARAMETER :: eps16 = 1.0E-16_DP
! 
!------------------------------------------------------------------------------!
END MODULE util_param
!------------------------------------------------------------------------------!
