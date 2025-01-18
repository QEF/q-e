!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE print_chempot_lauerism(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... print Laue-RISM's charges and chemical potentials
  !
  USE constants,      ONLY : RYTOEV
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE molecule_const, ONLY : RY_TO_KJMOLm1, RY_TO_KCALMOLm1
  USE mp,             ONLY : mp_sum
  USE rism,           ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,        ONLY : nsolV, solVs, get_nuniq_in_solVs, &
                           & iuniq_to_isite, isite_to_isolV
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER, PARAMETER       :: LEN_LABEL = 16
  !
  INTEGER                  :: nq
  INTEGER                  :: iq
  INTEGER                  :: iiq
  INTEGER                  :: iv
  INTEGER                  :: isolV
  INTEGER                  :: natom
  REAL(DP)                 :: rho
  REAL(DP), ALLOCATABLE    :: nsol(:)
  REAL(DP), ALLOCATABLE    :: qsol(:)
  REAL(DP), ALLOCATABLE    :: usgf(:)
  REAL(DP)                 :: usol_eV
  REAL(DP)                 :: usol_kJ
  REAL(DP)                 :: usol_kcal
  CHARACTER(LEN=LEN_LABEL) :: label1
  CHARACTER(LEN=LEN_LABEL) :: label2
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  ALLOCATE(nsol(nsolV + 1))
  ALLOCATE(qsol(nsolV + 1))
  ALLOCATE(usgf(nsolV + 1))
  nsol = 0.0_DP
  qsol = 0.0_DP
  usgf = 0.0_DP
  !
  ! ... sum charges and chemical potentials
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    natom = solVs(isolV)%natom
    !
    nsol(isolV) = nsol(isolV) + rismt%nsol(iiq) / DBLE(natom)
    qsol(isolV) = qsol(isolV) + rismt%qsol(iiq)
    usgf(isolV) = usgf(isolV) + rismt%usol_GF(iiq)
  END DO
  !
  nsol(nsolV + 1) = 0.0_DP
  qsol(nsolV + 1) = 0.0_DP
  usgf(nsolV + 1) = 0.0_DP
  DO isolV = 1, nsolV
    nsol(nsolV + 1) = nsol(nsolV + 1) + nsol(isolV)
    qsol(nsolV + 1) = qsol(nsolV + 1) + qsol(isolV)
    usgf(nsolV + 1) = usgf(nsolV + 1) + usgf(isolV)
  END DO
  !
  CALL mp_sum(nsol, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(qsol, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(usgf, rismt%mp_site%inter_sitg_comm)
  !
  ! ... write numbers
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"Total number of solvent (in expand-cell)")')
  !
  DO isolV = 1, nsolV
    !
    label1 = solVs(isolV)%name
    WRITE(stdout, 100) label1, nsol(isolV)
    !
100 FORMAT(5X,A10,X,F12.6)
    !
  END DO
  !
  ! ... write charges
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"Total charge of solvent (in expand-cell)")')
  !
  DO isolV = 1, (nsolV + 1)
    !
    IF (isolV <= nsolV) THEN
      label1 = solVs(isolV)%name
    ELSE
      label1 = 'Total     '
    END IF
    WRITE(stdout, 200) label1, qsol(isolV)
    !
200 FORMAT(5X,A10,X,F12.6,' e')
    !
  END DO
  !
  ! ... write chemical potentials
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"Chemical potential of solvation (in expand-cell)")')
  !
  DO isolV = 1, (nsolV + 1)
    !
    IF (isolV <= nsolV) THEN
      label1 = solVs(isolV)%name
    ELSE
      label1 = 'Total     '
    END IF
    label2 = 'GaussFluct'
    usol_eV   = usgf(isolV) * RYTOEV
    usol_kJ   = usgf(isolV) * RY_TO_KJMOLm1
    usol_kcal = usgf(isolV) * RY_TO_KCALMOLm1
#if defined (__DEBUG_RISM)
    WRITE(stdout, 300) label1, label2, usol_eV, usol_kJ, usol_kcal
#else
    WRITE(stdout, 300) label1, label2, usol_kcal
#endif
    !
#if defined (__DEBUG_RISM)
300 FORMAT(5X,A10,X,A10,X,E14.6,' eV',E14.6,' kJ/mol',E14.6,' kcal/mol')
#else
300 FORMAT(5X,A10,X,A10,X,E14.6,' kcal/mol')
#endif
    !
  END DO
  !
  WRITE(stdout, '()')
  !
  ! ... deallocate memory
  DEALLOCATE(nsol)
  DEALLOCATE(qsol)
  DEALLOCATE(usgf)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE print_chempot_lauerism
