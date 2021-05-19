!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE fcp_capacitance(capacitance)
  !----------------------------------------------------------------------------
  !
  ! ... evaluate capacitance for FCP
  !
  USE cell_base,     ONLY : alat, at
  USE constants,     ONLY : fpi, e2
  USE esm,           ONLY : esm_bc, esm_w
  USE kinds,         ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: capacitance
  !
  REAL(DP) :: fac
  REAL(DP) :: epsr
  REAL(DP) :: z0
  REAL(DP) :: area_xy
  !
  ! ... set permittivity and length of z-axis
  !
  IF (TRIM(esm_bc) == 'bc2' .OR. TRIM(esm_bc) == 'bc3' .OR. TRIM(esm_bc) == 'bc4') THEN
     !
     ! ... Parallel plate capacitor
     !
     IF (TRIM(esm_bc) == 'bc2') THEN
        fac = 2.0_DP
     ELSE
        fac = 1.0_DP
     END IF
     !
     z0 = 0.5_DP * alat * at(3, 3) + esm_w
     !
     epsr = 1.0_DP
     !
  ELSE
     !
     CALL errore('fcp_capacitance', 'cannot evaluate capacitance', 1)
     !
  END IF
  !
  ! ... calculate capacitance
  !
  area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
  !
  capacitance = fac * (epsr / fpi / e2) * area_xy / z0
  !
  !
END SUBROUTINE fcp_capacitance
