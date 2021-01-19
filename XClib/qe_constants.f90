!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE constants_l
  !----------------------------------------------------------------------------
  !! A subset of the constants in 'Modules' folder of QE. Here the ones
  !! needed internally in xc_lib only.
  !
  USE kind_l, ONLY : DP
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
  REAL(DP), PARAMETER :: sqrtpi = 1.77245385090551602729_DP 
  REAL(DP), PARAMETER :: sqrtpm1= 1.0_DP / sqrtpi
  REAL(DP), PARAMETER :: sqrt2  = 1.41421356237309504880_DP
  !
  !
  ! ... Physical constants, atomic units:
  !
  REAL(DP), PARAMETER :: H_PLANCK_SI    = 6.62607015E-34_DP      ! J s
  REAL(DP), PARAMETER :: HARTREE_SI     = 4.3597447222071E-18_DP ! J
  REAL(DP), PARAMETER :: BOHR_RADIUS_SI = 0.529177210903E-10_DP  ! m
  REAL(DP), PARAMETER :: C_SI           = 2.99792458E+8_DP       ! m sec^-1
  !
  REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  !
  !  Speed of light in atomic units
  REAL(DP), PARAMETER :: C_AU = C_SI / BOHR_RADIUS_SI * AU_SEC
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
  REAL(DP), PARAMETER :: gsmall = 1.0E-12_DP
  !
  REAL(DP), PARAMETER :: e2 = 2.0_DP      ! the square of the electron charge
  !
END MODULE constants_l
