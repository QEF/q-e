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
SUBROUTINE normalize_lauerism(rismt, charge, lboth, expand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... normalize total correlations to remove noise of solvent charge density,
  ! ... which is defined as
  !                   ----
  !        rho(gxy,z) >    q(v) * rho(v) * h(v; gxy,z)
  !                   ----
  !                     v
  !
  ! ... NOTE: h(gxy,z) is used, not g(gxy,z),
  ! ...       to obtain electrostatic consistent chemical potential of solvation.
  !
  ! ... Variables:
  ! ...   charge: total charge of solvent system
  ! ...   expand: use expand-cell(.TRUE.) or unit-cell(.FALSE.)
  !
  USE cell_base,      ONLY : at, alat
  USE constants,      ONLY : eps8, eps12, eps24
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,          ONLY : DP
  USE mp,             ONLY : mp_sum
  USE rism,           ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, nsolV, solVs, iuniq_to_nsite, &
                           & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: charge
  LOGICAL,         INTENT(IN)    :: lboth  ! both-hands calculation, or not
  LOGICAL,         INTENT(IN)    :: expand ! expand-cell(.TRUE.) or unit-cell(.FALSE.)
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: nq
  INTEGER               :: iq
  INTEGER               :: iiq
  INTEGER               :: iv
  INTEGER               :: nv
  INTEGER               :: isolV
  INTEGER               :: iatom
  INTEGER               :: irz
  INTEGER               :: iirz
  INTEGER               :: igxy
  INTEGER               :: jgxy
  INTEGER               :: ihand
  INTEGER               :: nhand
  INTEGER,  ALLOCATABLE :: izright_tail(:)
  INTEGER,  ALLOCATABLE :: izleft_tail(:)
  REAL(DP)              :: qv
  REAL(DP)              :: rhov1
  REAL(DP)              :: rhov2
  REAL(DP)              :: qrho1
  REAL(DP)              :: qrho2
  REAL(DP)              :: vrho
  REAL(DP)              :: vqrho
  REAL(DP)              :: dz
  REAL(DP)              :: area_xy
  REAL(DP)              :: dvol
  REAL(DP)              :: vtmp
  REAL(DP)              :: ntmp
  REAL(DP)              :: nsol_
  REAL(DP)              :: msol_
  REAL(DP)              :: charge0
  REAL(DP)              :: chgwei1
  REAL(DP)              :: chgwei2
  REAL(DP)              :: hz
  REAL(DP)              :: hr0
  REAL(DP)              :: hr1
  REAL(DP)              :: hr2
  REAL(DP), ALLOCATABLE :: hwei(:,:)
  REAL(DP), ALLOCATABLE :: vsol(:,:)
  REAL(DP), ALLOCATABLE :: nsol(:,:)
  REAL(DP), ALLOCATABLE :: msol(:,:)
  REAL(DP), ALLOCATABLE :: qsol(:)
  !
  INTEGER,    PARAMETER :: HZ_EDGE  = 3  ! to avoid noise at edges of unit-cell
  REAL(DP),   PARAMETER :: HZ_THR   = 1.0E-3_DP
  REAL(DP),   PARAMETER :: HZ_SMEAR = 2.0_DP  ! in bohr
  REAL(DP),   PARAMETER :: HW_THR   = 1.0E-8_DP
  !
#if defined (__DEBUG_RISM)
  !
  CALL start_clock('3DRISM_norm')
#endif
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
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzs < rismt%dfft%nr3) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%nsite > 0) THEN
    ALLOCATE(izright_tail(rismt%nsite))
    ALLOCATE(izleft_tail( rismt%nsite))
  END IF
  IF (rismt%lfft%nrz * rismt%nsite > 0) THEN
    ALLOCATE(hwei(rismt%lfft%nrz, rismt%nsite))
  END IF
  IF (rismt%nsite > 0) THEN
    ALLOCATE(vsol(rismt%nsite, 2))
    ALLOCATE(nsol(rismt%nsite, 2))
  END IF
  ALLOCATE(msol(nsolV, 2))
  ALLOCATE(qsol(nsolV))
  !
  ! ... set variables
  dz      = rismt%lfft%zstep * alat
  area_xy = ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)) * alat * alat
  dvol    = area_xy * dz
  !
  ! ...
  ! ... Define domain of total correlations
  ! ......................................................................................
  !
  ! ... detect truncating positions for each solvent site
  IF (rismt%nsite > 0) THEN
    izright_tail = 0
    izleft_tail  = 0
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      !
      izleft_tail(iiq) = 1
      DO irz = 1, rismt%lfft%izleft_gedge
        hz   = DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq))
        IF (ABS(hz) < HZ_THR .OR. ABS(irz - rismt%lfft%izcell_start) <= HZ_EDGE) THEN
          ! NOP
        ELSE
          izleft_tail(iiq) = irz
          EXIT
        END IF
      END DO
      !
      izright_tail(iiq) = rismt%lfft%nrz
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        iirz = rismt%lfft%nrz + rismt%lfft%izright_gedge - irz
        hz   = DBLE(rismt%hsgz(iirz, iiq) + rismt%hlgz(iirz, iiq))
        IF (ABS(hz) < HZ_THR .OR. ABS(iirz - rismt%lfft%izcell_end) <= HZ_EDGE) THEN
          ! NOP
        ELSE
          izright_tail(iiq) = iirz
          EXIT
        END IF
      END DO
      !
    END IF
    !
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(izright_tail, rismt%mp_site%intra_sitg_comm)
    CALL mp_sum(izleft_tail,  rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... calculate weight for total correlations
  IF (rismt%lfft%nrz * rismt%nsite > 0) THEN
    hwei = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
!$omp parallel do default(shared) private(irz)
    DO irz = 1, rismt%lfft%izleft_gedge
      hwei(irz, iiq) = 0.5_DP * erfc(DBLE(izleft_tail(iiq) - irz ) * dz / HZ_SMEAR)
      IF (hwei(irz, iiq) < HW_THR) THEN
        hwei(irz, iiq) = 0.0_DP
      END IF
    END DO
!$omp end parallel do
    !
!$omp parallel do default(shared) private(irz)
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      hwei(irz, iiq) = 0.5_DP * erfc(DBLE(irz - izright_tail(iiq)) * dz / HZ_SMEAR)
      IF (hwei(irz, iiq) < HW_THR) THEN
        hwei(irz, iiq) = 0.0_DP
      END IF
    END DO
!$omp end parallel do
    !
  END DO
  !
  ! ... truncate hgz, hg0, hsgz, hlgz
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
    IF (expand) THEN
      ! ... expand-cell
      DO igxy = 1, rismt%lfft%ngxy
        jgxy = (igxy - 1) * rismt%nrzl
        !
!$omp parallel do default(shared) private(irz)
        DO irz = 1, rismt%lfft%izleft_gedge
          rismt%hsgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hsgz(irz + jgxy, iiq)
          rismt%hlgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hlgz(irz + jgxy, iiq)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(irz)
        DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
          rismt%hsgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hsgz(irz + jgxy, iiq)
          rismt%hlgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hlgz(irz + jgxy, iiq)
        END DO
!$omp end parallel do
        !
      END DO
      !
    ELSE
      ! ... unit-cell
      DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
        jgxy = (igxy - 1) * rismt%nrzs
        !
!$omp parallel do default(shared) private(irz, iirz)
        DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz + jgxy, iiq) = hwei(irz, iiq) * rismt%hgz(iirz + jgxy, iiq)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(irz, iirz)
        DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz + jgxy, iiq) = hwei(irz, iiq) * rismt%hgz(iirz + jgxy, iiq)
        END DO
!$omp end parallel do
        !
      END DO
      !
!$omp parallel do default(shared) private(irz)
      DO irz = 1, rismt%lfft%izleft_gedge
        rismt%hg0(irz, iiq) = hwei(irz, iiq) * rismt%hg0(irz, iiq)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        rismt%hg0(irz, iiq) = hwei(irz, iiq) * rismt%hg0(irz, iiq)
      END DO
!$omp end parallel do
      !
    END IF
    !
  END DO
  !
  ! ... volume of solvent domain
  ! ... NOTE: vol1 and vol2 are defined only if gxystart > 1
  IF (rismt%nsite > 0) THEN
    vsol(:, :) = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      !
      vtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:vtmp)
      DO irz = 1, rismt%lfft%izleft_gedge
        vtmp = vtmp + dvol * hwei(irz, iiq)
      END DO
!$omp end parallel do
      vsol(iiq, 2) = vtmp
      !
      vtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:vtmp)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        vtmp = vtmp + dvol * hwei(irz, iiq)
      END DO
!$omp end parallel do
      vsol(iiq, 1) = vtmp
      !
    END IF
    !
  END DO
  !
  ! ...
  ! ... Correct stoichiometry of solvent molecules
  ! ......................................................................................
  !
  ! ... calculate total numbers of solvent atoms
  ! ... NOTE: nsol is defined only if gxystart > 1
  IF (rismt%nsite > 0) THEN
    nsol(:, :) = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    rhov1 = solVs(isolV)%density
    rhov2 = solVs(isolV)%subdensity
    !
    IF (rismt%lfft%gxystart < 2) THEN
      CYCLE
    END IF
    !
    IF (expand) THEN
      ! ... expand-cell
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = 1, rismt%lfft%izleft_gedge
        ntmp = ntmp + (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      nsol(iiq, 2) = nsol(iiq, 2) + rhov2 * dvol * ntmp
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        ntmp = ntmp + (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      nsol(iiq, 1) = nsol(iiq, 1) + rhov1 * dvol * ntmp
      !
    ELSE
      ! ... unit-cell
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = 1, rismt%lfft%izleft_gedge
        ntmp = ntmp + &
             & (hwei(irz, iiq) * DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      nsol(iiq, 2) = nsol(iiq, 2) + rhov2 * dvol * ntmp
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        ntmp = ntmp + &
             & (hwei(irz, iiq) * DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      nsol(iiq, 1) = nsol(iiq, 1) + rhov1 * dvol * ntmp
      !
    END IF
    !
  END DO
  !
  ! ... sum numbers and charges of solvent atoms in a molecule
  ! ... NOTE: msol and qsol are defined only if gxystart > 1
  msol(:, :) = 0.0_DP
  qsol(:)    = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    qv    = solVs(isolV)%charge(iatom)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      msol(isolV, 1) = msol(isolV, 1) + DBLE(nv) * nsol(iiq, 1)
      msol(isolV, 2) = msol(isolV, 2) + DBLE(nv) * nsol(iiq, 2)
      qsol(isolV)    = qsol(isolV)    + DBLE(nv) * qv
    END IF
  END DO
  !
  CALL mp_sum(msol, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(qsol, rismt%mp_site%inter_sitg_comm)
  !
  DO isolV = 1, nsolV
    IF (solVs(isolV)%natom > 0) THEN
      msol(isolV, 1) = msol(isolV, 1) / DBLE(solVs(isolV)%natom)
      msol(isolV, 2) = msol(isolV, 2) / DBLE(solVs(isolV)%natom)
      qsol(isolV)    = qsol(isolV)    / DBLE(solVs(isolV)%natom)
    ELSE
      msol(isolV, 1) = 0.0_DP
      msol(isolV, 2) = 0.0_DP
      qsol(isolV)    = 0.0_DP
    END IF
  END DO
  !
  ! ... renormalize hgz, hg0, hsgz, hlgz
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    !
    IF (rismt%lfft%gxystart < 2) THEN
      CYCLE
    END IF
    !
    IF (lboth) THEN
      nhand = 2
    ELSE
      nhand = 1
    END IF
    !
    DO ihand = 1, nhand
      !
      nsol_ = nsol(iiq,   ihand)
      msol_ = msol(isolV, ihand)
      IF (nhand == 1) THEN
        nsol_ = nsol_ + nsol(iiq,   2)
        msol_ = msol_ + msol(isolV, 2)
      END IF
      !
      IF (ABS(msol_ - nsol_) > eps8) THEN
        rhov1 = solVs(isolV)%density
        rhov2 = solVs(isolV)%subdensity
        !
        IF (ihand == 1) THEN
          vrho = vsol(iiq, 1) * rhov1
        ELSE !IF (ihand == 2) THEN
          vrho = vsol(iiq, 2) * rhov2
        END IF
        !
        IF (nhand == 1) THEN
          vrho = vrho + vsol(iiq, 2) * rhov2
        END IF
        !
        IF (ABS(vrho) <= eps12) THEN  ! will not be occurred
          CALL errore('normalize_lauerism', 'vrho is zero', 1)
        END IF
        hr0 = (msol_ - nsol_) / vrho
        !
        IF (expand) THEN
          ! ... expand-cell
          IF (nhand == 1 .OR. ihand == 2) THEN
!$omp parallel do default(shared) private(irz)
            DO irz = 1, rismt%lfft%izleft_gedge
              ! correct only short-range (for convenience)
              rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                                 & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
            END DO
!$omp end parallel do
          END IF
          !
          IF (nhand == 1 .OR. ihand == 1) THEN
!$omp parallel do default(shared) private(irz)
            DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
              ! correct only short-range (for convenience)
              rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                                 & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
            END DO
!$omp end parallel do
          END IF
          !
        ELSE
          ! ... unit-cell
          IF (nhand == 1 .OR. ihand == 2) THEN
!$omp parallel do default(shared) private(irz, iirz)
            DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
              iirz = irz - rismt%lfft%izcell_start + 1
              rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                                 & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
            END DO
!$omp end parallel do
          END IF
          !
          IF (nhand == 1 .OR. ihand == 1) THEN
!$omp parallel do default(shared) private(irz, iirz)
            DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
              iirz = irz - rismt%lfft%izcell_start + 1
              rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                                 & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
            END DO
!$omp end parallel do
          END IF
          !
          IF (nhand == 1 .OR. ihand == 2) THEN
!$omp parallel do default(shared) private(irz)
            DO irz = 1, rismt%lfft%izleft_gedge
              rismt%hg0(irz, iiq) = rismt%hg0(irz, iiq) + hwei(irz, iiq) * hr0
            END DO
!$omp end parallel do
          END IF
          !
          IF (nhand == 1 .OR. ihand == 1) THEN
!$omp parallel do default(shared) private(irz)
            DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
              rismt%hg0(irz, iiq) = rismt%hg0(irz, iiq) + hwei(irz, iiq) * hr0
            END DO
!$omp end parallel do
          END IF
          !
        END IF
        !
      END IF
      !
    END DO
    !
  END DO
  !
  ! ...
  ! ... Correct total charge of solvent system
  ! ......................................................................................
  !
  ! ... calculate total charge
  ! ... NOTE: charge0 is defined only if gxystart > 1
  charge0 = 0.0_DP
  !
  IF (lboth) THEN
    chgwei1 = 0.0_DP
    chgwei2 = 0.0_DP
  ELSE
    chgwei1 = 1.0_DP
    chgwei2 = 1.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      charge0 = charge0 &
            & + DBLE(nv) * qsol(isolV) * (msol(isolV, 1) + msol(isolV, 2))
      !
      IF (lboth) THEN
        chgwei1 = chgwei1 + DBLE(nv) * qsol(isolV) * msol(isolV, 1)
        chgwei2 = chgwei2 + DBLE(nv) * qsol(isolV) * msol(isolV, 2)
      END IF
    END IF
  END DO
  !
  CALL mp_sum(charge0, rismt%mp_site%inter_sitg_comm)
  !
  IF (lboth) THEN
    CALL mp_sum(chgwei1, rismt%mp_site%inter_sitg_comm)
    CALL mp_sum(chgwei2, rismt%mp_site%inter_sitg_comm)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      !
      chgwei1 = ABS(chgwei1)
      chgwei2 = ABS(chgwei2)
      !
      IF ((chgwei1 + chgwei2) <= eps8) THEN
        chgwei1 = 1.0_DP
        chgwei2 = 1.0_DP
      END IF
    END IF
  END IF
  !
  ! ... vqrho = sum(vol * qv^2 * rhov^2)
  ! ... NOTE: vqrho is defined only if gxystart > 1
  vqrho = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    rhov1 = solVs(isolV)%density
    rhov2 = solVs(isolV)%subdensity
    !
    IF (rismt%lfft%gxystart > 1) THEN
      qv    = qsol(isolV)
      qrho1 = chgwei1 * qv * qv * rhov1 * rhov1
      qrho2 = chgwei2 * qv * qv * rhov2 * rhov2
      vqrho = vqrho + DBLE(nv) * (vsol(iiq, 1) * qrho1 + vsol(iiq, 2) * qrho2)
    END IF
  END DO
  !
  CALL mp_sum(vqrho, rismt%mp_site%inter_sitg_comm)
  !
  ! ... renormalize hgz, hg0, hsgz, hlgz
  IF (rismt%lfft%gxystart > 1) THEN
    !
    IF (ABS(charge - charge0) > eps8) THEN
      !
      IF (ABS(vqrho) <= eps24) THEN  ! will not be occurred
        CALL errore('normalize_lauerism', 'vqrho is zero', 1)
      END IF
      hr0 = (charge - charge0) / vqrho
      !
      DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        iiq   = iq - rismt%mp_site%isite_start + 1
        iv    = iuniq_to_isite(1, iq)
        isolV = isite_to_isolV(iv)
        rhov1 = solVs(isolV)%density
        rhov2 = solVs(isolV)%subdensity
        qv    = qsol(isolV)
        hr1   = hr0 * chgwei1 * qv * rhov1
        hr2   = hr0 * chgwei2 * qv * rhov2
        !
        IF (expand) THEN
          ! ... expand-cell
!$omp parallel do default(shared) private(irz)
          DO irz = 1, rismt%lfft%izleft_gedge
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr2, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz)
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr1, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
        ELSE
          ! ... unit-cell
!$omp parallel do default(shared) private(irz, iirz)
          DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr2, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz, iirz)
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr1, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz)
          DO irz = 1, rismt%lfft%izleft_gedge
            rismt%hg0(irz, iiq) = rismt%hg0(irz, iiq) + hwei(irz, iiq) * hr2
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz)
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
            rismt%hg0(irz, iiq) = rismt%hg0(irz, iiq) + hwei(irz, iiq) * hr1
          END DO
!$omp end parallel do
          !
        END IF
        !
      END DO
      !
    END IF
    !
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
  ! ... deallocate memory
100 CONTINUE
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(izright_tail)
    DEALLOCATE(izleft_tail)
  END IF
  IF (rismt%lfft%nrz * rismt%nsite > 0) THEN
    DEALLOCATE(hwei)
  END IF
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(vsol)
    DEALLOCATE(nsol)
  END IF
  DEALLOCATE(msol)
  DEALLOCATE(qsol)
#if defined (__DEBUG_RISM)
  !
  CALL stop_clock('3DRISM_norm')
#endif
  !
END SUBROUTINE normalize_lauerism
