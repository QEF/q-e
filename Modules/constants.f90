!
! Copyright (C) 2002-2006 Quantum-Espresso group
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
  ! ... Mathematical constants
  ! 
  REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: tpi= 2.0_DP * pi
  REAL(DP), PARAMETER :: fpi= 4.0_DP * pi
  REAL(DP), PARAMETER :: sqrtpi = 1.77245385090551602729_DP 
  REAL(DP), PARAMETER :: sqrtpm1= 1.0_DP / sqrtpi
  REAL(DP), PARAMETER :: sqrt2  = 1.41421356237309504880_DP
  !
  ! ... Physical constants, SI (NIST CODATA 2002)
  !
  REAL(DP), PARAMETER :: H_PLANCK_SI      = 6.6260693D-34    ! J s
  REAL(DP), PARAMETER :: K_BOLTZMANN_SI   = 1.3806505D-23    ! J K^-1 
  REAL(DP), PARAMETER :: ELECTRON_SI      = 1.60217653D-19   ! C
  REAL(DP), PARAMETER :: ELECTRONVOLT_SI  = 1.60217653D-19   ! J  
  REAL(DP), PARAMETER :: ELECTRONMASS_SI  = 9.1093826D-31    ! Kg
  REAL(DP), PARAMETER :: HARTREE_SI       = 4.35974417D-18   ! J
  REAL(DP), PARAMETER :: RYDBERG_SI       = HARTREE_SI/2.0_DP! J
  REAL(DP), PARAMETER :: BOHR_RADIUS_SI   = 0.5291772108D-10 ! m
  REAL(DP), PARAMETER :: AMU_SI           = 1.66053886D-27   ! Kg
  !
  ! ... Physical constants, atomic units:
  ! ... AU for "Hartree" atomic units (e = m = hbar = 1)
  ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
  !
  REAL(DP), PARAMETER :: K_BOLTZMANN_AU   = K_BOLTZMANN_SI / HARTREE_SI
  REAL(DP), PARAMETER :: K_BOLTZMANN_RY   = K_BOLTZMANN_SI / RYDBERG_SI
  !
  ! ... Unit conversion factors: energy and masses
  !
  REAL(DP), PARAMETER :: AUTOEV           = HARTREE_SI / ELECTRONVOLT_SI
  REAL(DP), PARAMETER :: RYTOEV           = AUTOEV / 2.0_DP
  REAL(DP), PARAMETER :: AMU_AU           = AMU_SI / ELECTRONMASS_SI
  REAL(DP), PARAMETER :: AMU_RY           = AMU_AU / 2.0_DP
  !
  ! ... Unit conversion factors: atomic unit of time, in s and ps
  !
  REAL(DP), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(DP), PARAMETER :: AU_PS            = AU_SEC * 1.0D+12
  !
  ! ... Unit conversion factors: pressure (1 Pa = 1 J/m^3, 1GPa = 10 Kbar )
  !
  REAL(DP), PARAMETER :: AU_GPA           = HARTREE_SI / BOHR_RADIUS_SI ** 3 &
                                            / 1.0D+9 
  REAL(DP), PARAMETER :: RY_KBAR          = 10.0_dp * AU_GPA / 2.0_dp
  !
  ! ... Unit conversion factors: 1 debye = 10^-18 esu*cm 
  ! ...                                  = 3.3356409519*10^-30 C*m 
  ! ...                                  = 0.208194346 e*A
  ! ... ( 1 esu = (0.1/c) Am, c=299792458 m/s)
  !
  REAL(DP), PARAMETER :: DEBYE_SI         = 3.3356409519 * 1.0D-30 ! C*m 
  REAL(DP), PARAMETER :: AU_DEBYE         = ELECTRON_SI * BOHR_RADIUS_SI / &
	                                       DEBYE_SI
  !
  REAL(DP), PARAMETER :: eV_to_kelvin = ELECTRONVOLT_SI / K_BOLTZMANN_SI
  REAL(DP), PARAMETER :: ry_to_kelvin = RYDBERG_SI / K_BOLTZMANN_SI
  !
  ! ... zero up to a given accuracy
  !
  REAL(DP), PARAMETER :: eps4  = 1.0D-4
  REAL(DP), PARAMETER :: eps6  = 1.0D-6
  REAL(DP), PARAMETER :: eps8  = 1.0D-8
  REAL(DP), PARAMETER :: eps14 = 1.0D-14
  REAL(DP), PARAMETER :: eps16 = 1.0D-16
  REAL(DP), PARAMETER :: eps32 = 1.0D-32
  !
  REAL(DP), PARAMETER :: gsmall = 1.0d-12
  !
  REAL(DP), PARAMETER :: e2 = 2.D0      ! the square of the electron charge
  REAL(DP), PARAMETER :: degspin = 2.D0 ! the number of spins per level
  !
  !!!!!! COMPATIBIILITY
  !
  REAL(DP), PARAMETER :: amconv = AMU_RY
  REAL(DP), PARAMETER :: uakbar = RY_KBAR
  REAL(DP), PARAMETER :: bohr_radius_cm = bohr_radius_si * 100.0
  REAL(DP), PARAMETER :: BOHR_RADIUS_ANGS = bohr_radius_cm * 1.0D8
  REAL(DP), PARAMETER :: ANGSTROM_AU = 1.0/BOHR_RADIUS_ANGS
  REAL(DP), PARAMETER :: DIP_DEBYE = AU_DEBYE
  REAL(DP), PARAMETER :: AU_TERAHERTZ  = AU_PS
  REAL(DP), PARAMETER :: AU_TO_OHMCMM1 = 46000.0D0 ! (ohm cm)^-1
  !

END MODULE constants
