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
  !! 
  !! Useful constants used in EPW
  !! 
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! Mathematical constants
  ! 
  REAL(KIND = DP), PARAMETER :: pi     = 3.141592653589793238462643383279502884197169399375105820974944d0
  REAL(KIND = DP), PARAMETER :: twopi  = 2.d0 * pi
  REAL(KIND = DP), PARAMETER :: fpi    = 4.d0 * pi
  REAL(KIND = DP), PARAMETER :: pibytwo=3.141592653589793238462643383279502884197169399375105820974944d0 / 2.d0
  REAL(KIND = DP), PARAMETER :: one    = 1.d0
  REAL(KIND = DP), PARAMETER :: two    = 2.d0
  REAL(KIND = DP), PARAMETER :: zero   = 0.d0
  REAL(KIND = DP), PARAMETER :: e2     = 2.0_DP      ! the square of the electron charge
  COMPLEX(KIND = DP), PARAMETER :: ci    = (0.d0, 1.d0)
  COMPLEX(KIND = DP), PARAMETER :: cone  = (1.d0, 0.d0)
  COMPLEX(KIND = DP), PARAMETER :: czero = (0.d0, 0.d0)
  !
  ! Unit conversion factors
  !
  REAL(KIND = DP), PARAMETER :: ang2cm   = 1.0d-8
  REAL(KIND = DP), PARAMETER :: ang2m    = 1.0d-10  
  REAL(KIND = DP), PARAMETER :: cm2m     = 1.0d-2
  REAL(KIND = DP), PARAMETER :: bohr     = 0.52917721092d0
  REAL(KIND = DP), PARAMETER :: ryd2mev  = 13605.6981d0
  REAL(KIND = DP), PARAMETER :: ryd2ev   = 13.6056981d0
  REAL(KIND = DP), PARAMETER :: ha2ev    = 2.d0 * ryd2ev
  REAL(KIND = DP), PARAMETER :: rydcm1   = ryd2ev * 8065.541d0
  REAL(KIND = DP), PARAMETER :: bohr2ang = 0.52917721092d0
  REAL(KIND = DP), PARAMETER :: ev2cmm1  = 8065.541d0
  REAL(KIND = DP), PARAMETER :: mev2cmm1 = 8.065541d0
  REAL(KIND = DP), PARAMETER :: cmm12meV = 1.0d0 / 8.065541d0
  REAL(KIND = DP), PARAMETER :: kelvin2eV= 8.6173427909d-05
  REAL(KIND = DP), PARAMETER :: kelvin2Ry= 6.333627859634130e-06
  REAL(KIND = DP), PARAMETER :: ryd2ghz  = 3.289828d6
  REAL(KIND = DP), PARAMETER :: mev2ps   = 0.6582119514  ! 1000/((1/hbar)*1e-12)
  REAL(KIND = DP), PARAMETER :: mev2invps = 1.0 / meV2ps  
  REAL(KIND = DP), PARAMETER :: kb       = 8.6173324d-05 ! eV/K
  REAL(KIND = DP), PARAMETER :: electron_SI = 1.602176487d-19
  REAL(KIND = DP), PARAMETER :: hbar     = 6.582119514E-16 ! eV*s
  REAL(KIND = DP), PARAMETER :: hbarJ    = 1.054571800E-34 ! J*s 
  REAL(KIND = DP), PARAMETER :: byte2Mb  = 7.62939453125E-6 ! 8 / (1024 * 1024) because 8 bytes per number, value in Mb; FM : true if DP is REAL(8) in qe
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
  END MODULE constants_epw
  !-----------------------------------------------------------------------
