!
! Copyright (C) 2002-2004 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
MODULE upf_kinds
!------------------------------------------------------------------------------!
  ! 
  IMPLICIT NONE
  SAVE
  ! ... kind definitions
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)
  INTEGER, PARAMETER :: i4b = selected_int_kind(9)
  INTEGER, PARAMETER :: i8b = selected_int_kind(18)
  PRIVATE
  PUBLIC :: i4b, i8b, sgl, DP
  !
END MODULE
