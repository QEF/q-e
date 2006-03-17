!
! Copyright (C) 2002-2003 FPMD & PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE constants
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  ! ... The constants needed everywhere
  !
  IMPLICIT NONE
  !
  SAVE
  !
  REAL(DP) :: PI, tpi, FPI
  REAL(DP) :: SQRTPI, SQRTPM1, sqrt2 
  !
  ! ... Physical constants
  !
  REAL(DP), PARAMETER :: K_BOLTZMAN_SI    = 1.38066D-23  ! J K^-1 
  REAL(DP), PARAMETER :: K_BOLTZMAN_AU    = 3.1667D-6    ! Hartree K^-1 
  REAL(DP), PARAMETER :: K_BOLTZMAN_M1_AU = 315795.26D0  ! Hartree^-1 K 
  REAL(DP), PARAMETER :: FACTEM           = 315795.26D0  ! 27.212d0*11605.d0 Hartree^-1 K 
  !
  ! ... Physical constants defining the Atomic Units System
  !
  REAL(DP), PARAMETER :: BOHR_RADIUS_SI   = 0.529177D-10 ! m
  REAL(DP), PARAMETER :: BOHR_RADIUS_CM   = 0.529177D-8  ! cm
  REAL(DP), PARAMETER :: BOHR_RADIUS_ANGS = 0.529177D0   ! angstrom
  REAL(DP), PARAMETER :: ELECTRONMASS_SI  = 9.10953D-31  ! Kg
  REAL(DP), PARAMETER :: ELECTRONMASS_AMU = 5.4858D-4    ! amu
  !
  ! ... Units conversion factors
  !
  REAL(DP), PARAMETER :: ELECTRONVOLT_SI  = 1.6021892D-19  ! J  
  REAL(DP), PARAMETER :: AMU_SI           = 1.66057D-27    ! Kg
  REAL(DP), PARAMETER :: ANGSTROM_AU      = 1.889727D0     ! au
  REAL(DP), PARAMETER :: AU_TO_OHMCMM1    = 46000.0D0      ! (ohm cm)^-1
  REAL(DP), PARAMETER :: AU_KB            = 294210.0D0     ! Kbar
  REAL(DP), PARAMETER :: KB_AU            = 1.0D0/294210.0D0 ! au
  REAL(DP), PARAMETER :: AU_GPA           = 29421.0D0      ! GPa
  REAL(DP), PARAMETER :: GPA_AU           = 1.0D0/29421.0D0  ! au
  REAL(DP), PARAMETER :: AU               = 27.211652D0    ! eV
  REAL(DP), PARAMETER :: RY               = 13.605826D0    ! eV
  REAL(DP), PARAMETER :: SCMASS           = 1822.89D0   ! amu to au ( mass of a proton )
  REAL(DP), PARAMETER :: AMU_AU           = 1822.89D0   ! au
  REAL(DP), PARAMETER :: AU_TERAHERTZ     = 2.418D-5    ! THz
  REAL(DP), PARAMETER :: TERAHERTZ        = 2.418D-5    ! from au to THz
  REAL(DP), PARAMETER :: AU_SEC           = 2.4189D-17  ! sec
  REAL(DP), PARAMETER :: AU_PS            = 2.4189D-5   ! sec
  REAL(DP), PARAMETER :: DIP_DEBYE        = 2.54168     ! hartree atomic units dipole to Debye
  !     
  !
  PARAMETER( pi        = 3.14159265358979323846_DP )
  PARAMETER( tpi       = 2.0_DP * 3.14159265358979323846_DP )
  PARAMETER( fpi       = 4.0_DP * 3.14159265358979323846_DP )
  PARAMETER( sqrtpi    = 1.77245385090551602729_DP )
  PARAMETER( sqrtpm1   = 1.0_DP / sqrtpi )
  PARAMETER( sqrt2     = 1.41421356237309504880 )
  !
  !
  REAL(DP), PARAMETER :: rhothr = 1.0e-5_DP ! tolerance
  REAL(DP), PARAMETER :: gsmall = 1.0d-12
  !
  REAL(DP), PARAMETER :: e2 = 2.D0      ! the square of the electron charge
  REAL(DP), PARAMETER :: degspin = 2.D0 ! the number of spins per level
  REAL(DP), PARAMETER :: rytoev=13.6058d0      ! conversion from Ry to eV
  !
  ! ... mass conversion: a.m.u to a.u. (Ry)
  !
  REAL(DP), PARAMETER :: amconv= 1.66042D-24 / 9.1095D-28 * 0.5D0 
  !
  ! ... pressure conversion from Ry/(a.u)^3 to K
  !
  REAL(DP), PARAMETER :: uakbar= 147105.d0
  !
  ! ... zero up to a given accuracy
  !
  REAL(DP), PARAMETER :: eps4  = 1.0D-4
  REAL(DP), PARAMETER :: eps8  = 1.0D-8
  REAL(DP), PARAMETER :: eps14 = 1.0D-14
  REAL(DP), PARAMETER :: eps16 = 1.0D-16
  REAL(DP), PARAMETER :: eps32 = 1.0D-32
  !
  REAL(DP), PARAMETER :: eV_to_kelvin = 1.1604D4            ! from eV to Kelvin
  REAL(DP), PARAMETER :: ry_to_kelvin = 315642.28D0 * 0.5D0 ! from Ry to Kelvin
  !
END MODULE constants
