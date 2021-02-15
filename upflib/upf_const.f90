!
! Copyright (C) 2002-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE upf_const
  !----------------------------------------------------------------------------
  !
  USE upf_kinds, ONLY : DP
  !
  ! ... The constants needed everywhere
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... Mathematical constants
  ! 
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP 
  REAL(DP), PARAMETER :: tpi    = 2.0_DP * pi
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  REAL(DP), PARAMETER :: sqrt2  = 1.41421356237309504880_DP
  !
  ! ... zero up to a given accuracy
  !
  REAL(DP), PARAMETER :: eps4  = 1.0E-4_DP
  REAL(DP), PARAMETER :: eps6  = 1.0E-6_DP
  REAL(DP), PARAMETER :: eps8  = 1.0E-8_DP
  REAL(DP), PARAMETER :: eps12 = 1.0E-12_DP
  REAL(DP), PARAMETER :: eps14 = 1.0E-14_DP
  REAL(DP), PARAMETER :: eps16 = 1.0E-16_DP
  REAL(DP), PARAMETER :: eps24 = 1.0E-24_DP
  REAL(DP), PARAMETER :: eps32 = 1.0E-32_DP
  !
  ! physical constants
  !
  REAL(DP), PARAMETER :: e2 = 2.0_DP      ! the square of the electron charge
  !
END MODULE upf_const

