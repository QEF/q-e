!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE molecule_const
  !--------------------------------------------------------------------------
  !
  ! ... this module defines constatns used for molecule_type
  !
  USE constants, ONLY : RYDBERG_SI, BOHR_RADIUS_CM
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... avogadro number
  REAL(DP), PARAMETER :: AVOGADRO = 6.022142E+23_DP
  ! ... 1 cal = (CAL_TO_J) * J
  REAL(DP), PARAMETER :: CAL_TO_J = 4.1868_DP
  ! ... 1 bohr^3 = (BOHR3_CM3) * cm^3
  REAL(DP), PARAMETER :: BOHR3_CM3 = BOHR_RADIUS_CM * BOHR_RADIUS_CM * BOHR_RADIUS_CM
  ! ... 1 bohr^3 = (BOHR3_L) * L
  REAL(DP), PARAMETER :: BOHR3_L = BOHR3_CM3 / 1000.0_DP
  !
  ! ... units of energy
  ! ... 1 Ry = (RY_TO_KJMOL) * kJ mol^-1
  REAL(DP), PARAMETER :: RY_TO_KJMOLm1 = RYDBERG_SI * AVOGADRO / 1000.0_DP
  ! ... 1 Ry = (RY_TO_KCALMOL) * kcal mol^-1
  REAL(DP), PARAMETER :: RY_TO_KCALMOLm1 = RY_TO_KJMOLm1 / CAL_TO_J
  !
  ! ... units of density
  ! ... 1 bohr^-3 = (BOHRm3_TO_MOLCMm3) * mol cm^-3
  REAL(DP), PARAMETER :: BOHRm3_TO_MOLCMm3 = 1.0_DP / AVOGADRO / BOHR3_CM3
  ! ... 1 bohr^-3 = (BOHRm3_TO_MOLLm1) * mol L^-1
  REAL(DP), PARAMETER :: BOHRm3_TO_MOLLm1 = 1.0_DP / AVOGADRO / BOHR3_L
  !
  ! ... public components
  PUBLIC :: RY_TO_KJMOLm1
  PUBLIC :: RY_TO_KCALMOLm1
  PUBLIC :: BOHRm3_TO_MOLCMm3
  PUBLIC :: BOHRm3_TO_MOLLm1
  !
END MODULE molecule_const
