!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_stress(rismt, sigma, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate stress tensor of 3D-RISM.
  ! ...
  ! ... S_sol = S_ele + S_ion + S_lj
  ! ...
  ! ... (1) S_ele: contribution from electrons, which is included
  ! ...            in Hellmann-Feynman force, is set to zero now.
  ! ...
  ! ... (2) S_ion: contribution from ions (local part)
  ! ...
  ! ... (3) S_lj:  contribution from Lennard-Jones potential
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE solute,    ONLY : get_solU_LJ_stress
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  INTEGER,         INTENT(OUT) :: ierr
  !
  REAL(DP) :: sigmaion(3, 3)
  REAL(DP) :: sigmalj (3, 3)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... ions term (local part)
  sigmaion = 0.0_DP
  CALL solvation_stress_ion(rismt, sigmaion, ierr)
  !
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... Lennad-Jones term
  sigmalj = 0.0_DP
  CALL get_solU_LJ_stress(rismt, sigmalj, ierr)
  !
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... total solvation stress
  sigma = sigmaion + sigmalj
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_stress
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_stress_ion(rismt, sigma, ierr)
  !---------------------------------------------------------------------------
  !
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE cell_base,     ONLY : alat
  USE control_flags, ONLY : gamma_only
  USE kinds,         ONLY : DP
  USE rism,          ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  INTEGER,         INTENT(OUT) :: ierr
  !
  REAL(DP)             :: fact
  REAL(DP)             :: sigmesm(3, 3)
  COMPLEX(DP), POINTER :: rhogt(:)
  LOGICAL              :: laue
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ng < rismt%gvec%ngm) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... factor of Gamma-trick
  IF (gamma_only) THEN
    fact = 2.0_DP
  ELSE
    fact = 1.0_DP
  END IF
  !
  ! ... set solvent density
  IF (laue) THEN
    rhogt => rismt%rhog_pbc
  ELSE
    rhogt => rismt%rhog
  END IF
  !
  ! ... calculate stress
  sigma = 0.0_DP
  ! TODO
  ! TODO implement short-range local potential term
  ! TODO
  !
  ! ... clear solvent density
  NULLIFY(rhogt)
  !
  ! ... add long-range from ESM(BC1)
  IF (laue) THEN
    sigmesm = 0.0_DP
    CALL solvation_esm_stress(rismt, 1.0_DP / alat, sigmesm, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
    sigma = sigma + sigmesm
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_stress_ion
