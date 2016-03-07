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
  REAL(DP), PARAMETER :: bohr     = 0.52917721092d0
  REAL(DP), PARAMETER :: ryd2mev  = 13605.6981d0
  REAL(DP), PARAMETER :: ryd2ev   = 13.6056981d0
  REAL(DP), PARAMETER :: rydcm1   = 13.6056981d0 * 8065.541d0
  REAL(DP), PARAMETER :: bohr2ang = 0.52917721092d0
  REAL(DP), PARAMETER :: ev2cmm1  = 8065.541d0
  REAL(DP), PARAMETER :: kelvin2eV= 8.6173427909d-05
  REAL(DP), PARAMETER :: ryd2ghz  = 3.289828d6
  !
  END MODULE constants_epw

