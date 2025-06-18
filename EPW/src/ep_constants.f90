  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-----------------------------------------------------------------------
  MODULE ep_constants
  !-----------------------------------------------------------------------
  !!
  !! Useful constants used in EPW
  !! Some constants are renamed with respect to QE.
  !!
  !! Jan 2020 SP - Use QE Modules/constants.f90 whenever possible.
  !!
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : bohr_radius_si, ry_to_cmm1, rytoev, autoev, ev_to_kelvin, ry_to_ghz, &
                         ry_to_kelvin, h_planck_si, electronvolt_si, k_boltzmann_si, pi,      &
                         EPSNOUGHT_SI, ELECTRON_SI, AMU_SI, AMU_RY
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! Mathematical constants
  !
  ! pi = 3.141592653589793238462643383279502884197169399375105820974944E0_DP
  REAL(KIND = DP), PARAMETER :: twopi    = 2.0E0_DP * pi
  REAL(KIND = DP), PARAMETER :: fpi      = 4.0E0_DP * pi
  REAL(KIND = DP), PARAMETER :: pibytwo  = pi / 2.0E0_DP
  REAL(KIND = DP), PARAMETER :: one      = 1.0E0_DP
  REAL(KIND = DP), PARAMETER :: two      = 2.0E0_DP
  REAL(KIND = DP), PARAMETER :: thre     = 3.0E0_DP
  REAL(KIND = DP), PARAMETER :: four     = 4.0E0_DP
  REAL(KIND = DP), PARAMETER :: half     = 0.5E0_DP
  REAL(KIND = DP), PARAMETER :: zero     = 0.0E0_DP
  REAL(KIND = DP), PARAMETER :: e2       = 2.0E0_DP             ! the square of the electron charge
  COMPLEX(KIND = DP), PARAMETER :: ci    = (0.0E0_DP, 1.0E0_DP)
  COMPLEX(KIND = DP), PARAMETER :: cone  = (1.0E0_DP, 0.0E0_DP)
  COMPLEX(KIND = DP), PARAMETER :: czero = (0.0E0_DP, 0.0E0_DP)
  !
  ! Unit conversion factors
  !
  REAL(KIND = DP), PARAMETER :: ang2cm   = 1.0E-8_DP
  REAL(KIND = DP), PARAMETER :: ang2m    = 1.0E-10_DP
  REAL(KIND = DP), PARAMETER :: cm2m     = 1.0E-2_DP
  REAL(KIND = DP), PARAMETER :: bohr     = bohr_radius_si * 1.0E10_DP ! 0.52917720859E-10_DP
  REAL(KIND = DP), PARAMETER :: ryd2ev   = rytoev                     ! 13.6056981
  REAL(KIND = DP), PARAMETER :: ryd2mev  = ryd2ev * 1.0E3_DP          ! 13605.6981
  REAL(KIND = DP), PARAMETER :: ha2ev    = autoev                     ! 27.2113962
  REAL(KIND = DP), PARAMETER :: rydcm1   = ry_to_cmm1                 ! ryd2ev * 8065.541 = 109737.315859
  REAL(KIND = DP), PARAMETER :: bohr2ang = bohr                       ! 0.52917720859
  REAL(KIND = DP), PARAMETER :: bohr2nm  = bohr / 10.0_DP             ! 0.052917720859
  REAL(KIND = DP), PARAMETER :: ev2cmm1  = rydcm1 / ryd2ev            ! 8065.541
  REAL(KIND = DP), PARAMETER :: mev2cmm1 = ev2cmm1 * 1E-3_DP          ! 8.065541
  REAL(KIND = DP), PARAMETER :: cmm12meV = 1.0E0_DP / mev2cmm1        ! 1.0 / 8.065541 = 0.12398424358
  REAL(KIND = DP), PARAMETER :: kelvin2eV= 1.0E0_DP / ev_to_kelvin    ! 8.6173427909d-05
  REAL(KIND = DP), PARAMETER :: kelvin2Ry= 1.0E0_DP / ry_to_kelvin    ! 1/157887.326147 = 6.33363E-6
  REAL(KIND = DP), PARAMETER :: ryd2ghz  = ry_to_ghz                  ! 3.289828d6
  REAL(KIND = DP), PARAMETER :: hbarJ    = h_planck_si / twopi        ! 1.054571800E-34 J*s
  REAL(KIND = DP), PARAMETER :: hbar     = hbarJ / electronvolt_si    ! 6.582119514E-16 eV*s
  REAL(KIND = DP), PARAMETER :: mev2ps   = 1.0E3_DP * hbar * 1.0E12   ! 1000/((1/hbar)*1e-12) = 0.6582119514
  REAL(KIND = DP), PARAMETER :: mev2invps= 1.0 / mev2ps               ! 1.51926746069
  REAL(KIND = DP), PARAMETER :: byte2Mb  = 7.62939453125E-6_DP        ! 8 / (1024 * 1024) because 8 bytes per number, value in Mb
  REAL(KIND = DP), PARAMETER :: kb       = k_boltzmann_si / electronvolt_si ! 8.6173324d-05 eV/K
  REAL(KIND = DP), PARAMETER :: cc2cb    = 6.74822779181357d24        ! cubic cm to cubic bohr
  REAL(KIND = DP), PARAMETER :: eps_vac  = EPSNOUGHT_SI               ! 8.854187817E-12  F/m (C^2/Jm)
  REAL(KIND = DP), PARAMETER :: echg     = ELECTRON_SI                ! 1.602176634E-19_DP  C
  REAL(KIND = DP), PARAMETER :: amu      = AMU_SI                     ! 1.66053906660E-27_DP  kg
  REAL(KIND = DP), PARAMETER :: sm1toryd = h_planck_si / echg / rytoev    ! convert frequencies in s-1 to Ry
  REAL(KIND = DP), PARAMETER :: ry2thz_sr  = ryd2ghz / 1.0E3_DP * twopi  ! 20670.598952748 Convert rydberg to THz for scattering rates.
  !
  ! ... zero up to a given accuracy
  !
  REAL(KIND = DP), PARAMETER :: eps2  = 1.0E-2_DP
  REAL(KIND = DP), PARAMETER :: eps4  = 1.0E-4_DP
  REAL(KIND = DP), PARAMETER :: eps5  = 1.0E-5_DP
  REAL(KIND = DP), PARAMETER :: eps6  = 1.0E-6_DP
  REAL(KIND = DP), PARAMETER :: eps8  = 1.0E-8_DP
  REAL(KIND = DP), PARAMETER :: eps10 = 1.0E-10_DP
  REAL(KIND = DP), PARAMETER :: eps12 = 1.0E-12_DP
  REAL(KIND = DP), PARAMETER :: eps14 = 1.0E-14_DP
  REAL(KIND = DP), PARAMETER :: eps16 = 1.0E-16_DP
  REAL(KIND = DP), PARAMETER :: eps20 = 1.0E-20_DP
  REAL(KIND = DP), PARAMETER :: eps24 = 1.0E-24_DP
  REAL(KIND = DP), PARAMETER :: eps32 = 1.0E-32_DP
  REAL(KIND = DP), PARAMETER :: eps40 = 1.0E-40_DP
  REAL(KIND = DP), PARAMETER :: eps80 = 1.0E-80_DP
  REAL(KIND = DP), PARAMETER :: eps160 = 1.0E-160_DP
  !
  !-----------------------------------------------------------------------
  END MODULE ep_constants
  !-----------------------------------------------------------------------
