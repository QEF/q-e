  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt . 
  !
  !-----------------------------------------------------------------------
  MODULE constants_epw
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! Mathematical constants
  ! 
  REAL(DP), PARAMETER :: pi     = 3.141592653589793238462643383279502884197169399375105820974944d0
  REAL(DP), PARAMETER :: twopi  = 2.d0 * pi
  REAL(DP), PARAMETER :: fpi    = 4.d0 * pi
  REAL(DP), PARAMETER :: pibytwo=3.141592653589793238462643383279502884197169399375105820974944d0 / 2.d0
  REAL(DP), PARAMETER :: one    = 1.d0
  REAL(DP), PARAMETER :: two    = 2.d0
  REAL(DP), PARAMETER :: zero   = 0.d0
  REAL(DP), PARAMETER :: e2     = 2.0_DP      ! the square of the electron charge
  COMPLEX(DP), PARAMETER :: ci   = (0.d0, 1.d0)
  COMPLEX(DP), PARAMETER :: cone = (1.d0, 0.d0)
  COMPLEX(DP), PARAMETER :: czero = (0.d0, 0.d0)
  !
  ! Unit conversion factors
  !
  REAL(DP), PARAMETER :: ang2cm   = 1.0d-8
  REAL(DP), PARAMETER :: ang2m    = 1.0d-10  
  REAL(DP), PARAMETER :: bohr     = 0.52917721092d0
  REAL(DP), PARAMETER :: ryd2mev  = 13605.6981d0
  REAL(DP), PARAMETER :: ryd2ev   = 13.6056981d0
  REAL(DP), PARAMETER :: rydcm1   = 13.6056981d0 * 8065.541d0
  REAL(DP), PARAMETER :: bohr2ang = 0.52917721092d0
  REAL(DP), PARAMETER :: ev2cmm1  = 8065.541d0
  REAL(DP), PARAMETER :: kelvin2eV= 8.6173427909d-05
  REAL(DP), PARAMETER :: ryd2ghz  = 3.289828d6
  REAL(DP), PARAMETER :: mev2ps   = 0.6582119514  ! 1000/((1/hbar)*1e-12)
  REAL(DP), PARAMETER :: mev2invps = 1.0 / meV2ps  
  REAL(DP), PARAMETER :: kb       = 8.6173324d-05 ! eV/K
  REAL(DP), PARAMETER :: electron_SI = 1.602176487d-19
  REAL(DP), PARAMETER :: hbar     = 6.582119514E-16 ! eV*s
  REAL(DP), PARAMETER :: hbarJ    = 1.054571800E-34 ! J*s  
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
  END MODULE constants_epw

