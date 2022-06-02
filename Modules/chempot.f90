!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE chempot(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate chemical potential of solvation
  !
  USE cell_base, ONLY : omega
  USE constants, ONLY : K_BOLTZMANN_RY, fpi
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE rism,      ONLY : rism_type, get_chempot_type, ITYPE_1DRISM, ITYPE_3DRISM, CHEMPOT_GF
  USE solvmol,   ONLY : get_nuniq_in_solVs, solVs, &
                      & iuniq_to_nsite, iuniq_to_isite, isite_to_isolV
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: ichempot
  INTEGER               :: ir
  INTEGER               :: isite
  INTEGER               :: nq
  INTEGER               :: iq
  INTEGER               :: iiq
  INTEGER               :: iv
  INTEGER               :: nv
  INTEGER               :: isolV
  REAL(DP)              :: beta
  REAL(DP)              :: r
  REAL(DP)              :: dr
  REAL(DP)              :: fac
  REAL(DP)              :: rhov
  REAL(DP), ALLOCATABLE :: weight(:)
  LOGICAL               :: lweight
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_1DRISM .AND. rismt%itype /= ITYPE_3DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%itype == ITYPE_1DRISM .AND. rismt%nr /= rismt%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%itype == ITYPE_3DRISM) THEN
    nq = get_nuniq_in_solVs()
    !
    IF (rismt%mp_site%nsite < nq) THEN
      ierr = IERR_RISM_INCORRECT_DATA_TYPE
      RETURN
    END IF
  END IF
  !
  ! ... if no data, return as normally done
  IF (rismt%nsite < 1) THEN
    GOTO 100
  END IF
  !
  ! ... set variables
  ichempot = get_chempot_type(rismt)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... calculate chemical potentials
  IF (rismt%nr > 0) THEN
    !
    ! ... set integral weight
    IF (rismt%itype == ITYPE_1DRISM) THEN
      ALLOCATE(weight(rismt%nr))
      dr = rismt%rfft%rgrid(2) - rismt%rfft%rgrid(1)
!$omp parallel do default(shared) private(ir, r)
      DO ir = 1, rismt%nr
        r = rismt%rfft%rgrid(rismt%mp_task%ivec_start + ir - 1)
        weight(ir) = fpi * r * r * dr
      END DO
!$omp end parallel do
      lweight = .TRUE.
    ELSE
      ALLOCATE(weight(1))
      weight(1) = 1.0_DP
      lweight = .FALSE.
    END IF
    !
    ! ... chemical potential for each site
    DO isite = 1, rismt%nsite
      CALL chempot_for_a_site(rismt%nr, ichempot,   beta, rismt%hr(1, isite), rismt%csr(1, isite), &
                            & rismt%ulr(1, isite), weight, lweight, rismt%usol(isite))
      CALL chempot_for_a_site(rismt%nr, CHEMPOT_GF, beta, rismt%hr(1, isite), rismt%csr(1, isite), &
                            & rismt%ulr(1, isite), weight, lweight, rismt%usol_GF(isite))
    END DO
    !
    IF (rismt%itype == ITYPE_3DRISM) THEN
      ! ... weight of FFT mesh (for 3D-RISM)
      fac = omega / DBLE(rismt%dfft%nr1 * rismt%dfft%nr2 * rismt%dfft%nr3)
      rismt%usol    = fac * rismt%usol
      rismt%usol_GF = fac * rismt%usol_GF
      !
      ! ... weight of solvent density (for 3D-RISM)
      DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        iiq   = iq - rismt%mp_site%isite_start + 1
        iv    = iuniq_to_isite(1, iq)
        nv    = iuniq_to_nsite(iq)
        isolV = isite_to_isolV(iv)
        rhov  = DBLE(nv) * solVs(isolV)%density
        rismt%usol(   iiq) = rhov * rismt%usol(   iiq)
        rismt%usol_GF(iiq) = rhov * rismt%usol_GF(iiq)
      END DO
    END IF
    !
    ! ... delete integral weight
    DEALLOCATE(weight)
    !
  ELSE
    !
    rismt%usol    = 0.0_DP
    rismt%usol_GF = 0.0_DP
    !
  END IF
  !
  CALL mp_sum(rismt%usol,    rismt%mp_task%itask_comm)
  CALL mp_sum(rismt%usol_GF, rismt%mp_task%itask_comm)
  !
  ! ... normally done
100 CONTINUE
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE chempot
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_for_a_site(nr, ichempot, beta, hr, csr, ulr, wr, lweight, usol)
  !---------------------------------------------------------------------------
  !
  ! ... calculate chemical potential for a site
  !
  USE kinds, ONLY : DP
  USE rism,  ONLY : CHEMPOT_HNC, CHEMPOT_KH, CHEMPOT_GF
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  INTEGER,  INTENT(IN)  :: ichempot
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: csr(1:*)
  REAL(DP), INTENT(IN)  :: ulr(1:*)
  REAL(DP), INTENT(IN)  :: wr(1:*)
  LOGICAL,  INTENT(IN)  :: lweight
  REAL(DP), INTENT(OUT) :: usol
  !
  INTEGER               :: ir
  REAL(DP), ALLOCATABLE :: tr(:)
  !
  REAL(DP), EXTERNAL :: ddot
  !
  ALLOCATE(tr(nr))
  !
  IF (ichempot == CHEMPOT_HNC) THEN
    CALL chempot_HNC_x(nr, beta, hr, csr, ulr, tr)
    !
  ELSE IF (ichempot == CHEMPOT_KH) THEN
    CALL chempot_KH_x(nr, beta, hr, csr, ulr, tr)
    !
  ELSE IF (ichempot == CHEMPOT_GF) THEN
    CALL chempot_GF_x(nr, beta, hr, csr, ulr, tr)
    !
  ELSE
    usol = 0.0_DP
    DEALLOCATE(tr)
    RETURN
  END IF
  !
  IF (lweight) THEN
    usol = ddot(nr, tr, 1, wr, 1)
  ELSE
    usol = 0.0_DP
!$omp parallel do default(shared) private(ir) reduction(+:usol)
    DO ir = 1, nr
      usol = usol + tr(ir)
    END DO
!$omp end parallel do
  END IF
  !
  usol = usol / beta
  !
  DEALLOCATE(tr)
  !
END SUBROUTINE chempot_for_a_site
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_HNC_x(nr, beta, hr, csr, ulr, tr)
  !---------------------------------------------------------------------------
  !
  ! ... HyperNetted-Chain model
  ! ... (J.P.Hansen et al., Theory of simple liquids. Academic Press, London, 1990.)
  ! ...
  ! ...            /    [  1                   1              ]
  ! ...   kB * T * | dr [ --- h(r)^2 - c(r) - --- h(r) * c(r) ]
  ! ...            /    [  2                   2              ]
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: csr(1:*)
  REAL(DP), INTENT(IN)  :: ulr(1:*)
  REAL(DP), INTENT(OUT) :: tr(1:*)
  !
  INTEGER  :: ir
  REAL(DP) :: cr0
  REAL(DP) :: hr0
  !
!$omp parallel do default(shared) private(ir, hr0, cr0)
  DO ir = 1, nr
    hr0 = hr(ir)
    cr0 = csr(ir) - beta * ulr(ir)
    !
    tr(ir) = 0.5_DP * hr0 * hr0 - cr0 - 0.5_DP * hr0 * cr0
  END DO
!$omp end parallel do
  !
END SUBROUTINE chempot_HNC_x
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_KH_x(nr, beta, hr, csr, ulr, tr)
  !---------------------------------------------------------------------------
  !
  ! ... Kovalenko and Hirata's model
  ! ... (A.Kovalenko, F.Hirata, J. Chem. Phys. 1999, 110, 10095-10112)
  ! ...
  ! ...            /    [  1                                  1              ]
  ! ...   kB * T * | dr [ --- h(r)^2 * Theta(-h(r)) - c(r) - --- h(r) * c(r) ]
  ! ...            /    [  2                                  2              ]
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: csr(1:*)
  REAL(DP), INTENT(IN)  :: ulr(1:*)
  REAL(DP), INTENT(OUT) :: tr(1:*)
  !
  INTEGER  :: ir
  REAL(DP) :: cr0
  REAL(DP) :: hr0
  !
!$omp parallel do default(shared) private(ir, hr0, cr0)
  DO ir = 1, nr
    hr0 = hr(ir)
    cr0 = csr(ir) - beta * ulr(ir)
    !
    IF (hr0 < 0.0_DP) THEN
      tr(ir) = 0.5_DP * hr0 * hr0 - cr0 - 0.5_DP * hr0 * cr0
    ELSE
      tr(ir) = -cr0 - 0.5_DP * hr0 * cr0
    END IF
  END DO
!$omp end parallel do
  !
END SUBROUTINE chempot_KH_x
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_GF_x(nr, beta, hr, csr, ulr, tr)
  !---------------------------------------------------------------------------
  !
  ! ... Gaussian Fluctuation model
  ! ... (D.Chandler et al., J. Chem. Phys. 1984, 81, 1975-1982)
  ! ...
  ! ...            /    [           1              ]
  ! ...   kB * T * | dr [ - c(r) - --- h(r) * c(r) ]
  ! ...            /    [           2              ]
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: csr(1:*)
  REAL(DP), INTENT(IN)  :: ulr(1:*)
  REAL(DP), INTENT(OUT) :: tr(1:*)
  !
  INTEGER  :: ir
  REAL(DP) :: cr0
  REAL(DP) :: hr0
  !
!$omp parallel do default(shared) private(ir, hr0, cr0)
  DO ir = 1, nr
    hr0 = hr(ir)
    cr0 = csr(ir) - beta * ulr(ir)
    !
    tr(ir) = -cr0 - 0.5_DP * hr0 * cr0
  END DO
!$omp end parallel do
  !
END SUBROUTINE chempot_GF_x
