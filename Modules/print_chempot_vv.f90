!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE print_chempot_vv(rismt, lhand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... print 1D-RISM's chemical potentials
  !
  USE constants,      ONLY : RYTOEV
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE molecule_const, ONLY : RY_TO_KJMOLm1, RY_TO_KCALMOLm1
  USE rism,           ONLY : rism_type, ITYPE_1DRISM
  USE solvmol,        ONLY : nsolV, solVs, get_nsite_in_solVs, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  LOGICAL,         INTENT(IN)  :: lhand  ! if true, right-hand. if false, left-hand
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER, PARAMETER       :: LEN_LABEL = 16
  !
  INTEGER                  :: nv
  INTEGER                  :: iv1, iv2
  INTEGER                  :: ivv
  INTEGER                  :: isolV1, isolV2
  INTEGER                  :: iatom1, iatom2
  REAL(DP)                 :: rho1, rho2
  REAL(DP), ALLOCATABLE    :: uscl(:,:)
  REAL(DP), ALLOCATABLE    :: usgf(:,:)
  REAL(DP)                 :: usol_eV
  REAL(DP)                 :: usol_kJ
  REAL(DP)                 :: usol_kcal
  CHARACTER(LEN=LEN_LABEL) :: label1
  CHARACTER(LEN=LEN_LABEL) :: label2
  CHARACTER(LEN=LEN_LABEL) :: label3
  !
  ! ... number of sites in solvents
  nv = get_nsite_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_1DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nsite < (nv * (nv + 1) / 2)) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  ALLOCATE(uscl(nsolV + 1, nsolV))
  ALLOCATE(usgf(nsolV + 1, nsolV))
  uscl = 0.0_DP
  usgf = 0.0_DP
  !
  ! ... sum chemical potentials
  DO iv1 = 1, nv
    ! ... properties of site1
    isolV1 = isite_to_isolV(iv1)
    iatom1 = isite_to_iatom(iv1)
    IF (lhand) THEN
      rho1 = solVs(isolV1)%density
    ELSE
      rho1 = solVs(isolV1)%subdensity
    END IF
    !
    DO iv2 = 1, iv1
      ! ... properties of site2
      isolV2 = isite_to_isolV(iv2)
      iatom2 = isite_to_iatom(iv2)
      IF (lhand) THEN
        rho2 = solVs(isolV2)%density
      ELSE
        rho2 = solVs(isolV2)%subdensity
      END IF
      !
      ivv = iv1 * (iv1 - 1) / 2 + iv2
      !
      uscl(isolV2, isolV1) = uscl(isolV2, isolV1) + rho2 * rismt%usol(ivv)
      usgf(isolV2, isolV1) = usgf(isolV2, isolV1) + rho2 * rismt%usol_GF(ivv)
      IF (iv2 < iv1) THEN
        uscl(isolV1, isolV2) = uscl(isolV1, isolV2) + rho1 * rismt%usol(ivv)
        usgf(isolV1, isolV2) = usgf(isolV1, isolV2) + rho1 * rismt%usol_GF(ivv)
      END IF
      !
    END DO
  END DO
  !
  DO isolV1 = 1, nsolV
    uscl(nsolV + 1, isolV1) = 0.0_DP
    usgf(nsolV + 1, isolV1) = 0.0_DP
    DO isolV2 = 1, nsolV
      uscl(nsolV + 1, isolV1) = uscl(nsolV + 1, isolV1) + uscl(isolV2, isolV1)
      usgf(nsolV + 1, isolV1) = usgf(nsolV + 1, isolV1) + usgf(isolV2, isolV1)
    END DO
  END DO
  !
  ! ... write chemical potentials
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"Chemical potential of solvation")')
  !
  DO isolV1 = 1, nsolV
    DO isolV2 = 1, (nsolV + 1)
      WRITE(stdout, '()')
      !
      IF (isolV2 == 1) THEN
        label1 = solVs(isolV1)%name
      ELSE
        label1 = '          '
      END IF
      IF (isolV2 <= nsolV) THEN
        label2 = solVs(isolV2)%name
      ELSE
        label2 = 'Total     '
      END IF
      label3 = 'Closure   '
      usol_eV   = uscl(isolV2, isolV1) * RYTOEV
      usol_kJ   = uscl(isolV2, isolV1) * RY_TO_KJMOLm1
      usol_kcal = uscl(isolV2, isolV1) * RY_TO_KCALMOLm1
#if defined (__DEBUG_RISM)
      WRITE(stdout, 1000) label1, label2, label3, usol_eV, usol_kJ, usol_kcal
#else
      WRITE(stdout, 1000) label1, label2, label3, usol_kcal
#endif
      !
      label1 = '          '
      label2 = '          '
      label3 = 'GaussFluct'
      usol_eV   = usgf(isolV2, isolV1) * RYTOEV
      usol_kJ   = usgf(isolV2, isolV1) * RY_TO_KJMOLm1
      usol_kcal = usgf(isolV2, isolV1) * RY_TO_KCALMOLm1
#if defined (__DEBUG_RISM)
      WRITE(stdout, 1000) label1, label2, label3, usol_eV, usol_kJ, usol_kcal
#else
      WRITE(stdout, 1000) label1, label2, label3, usol_kcal
#endif
      !
#if defined (__DEBUG_RISM)
1000  FORMAT(5X,A10,X,A10,X,A10,X,E14.6,' eV',E14.6,' kJ/mol',E14.6,' kcal/mol')
#else
1000  FORMAT(5X,A10,X,A10,X,A10,X,E14.6,' kcal/mol')
#endif
      !
    END DO
  END DO
  !
  WRITE(stdout, '()')
  !
  ! ... deallocate memory
  DEALLOCATE(uscl)
  DEALLOCATE(usgf)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE print_chempot_vv
