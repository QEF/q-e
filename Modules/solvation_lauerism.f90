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
SUBROUTINE solvation_lauerism(rismt, charge, ireference, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... setup charge density, solvation potential and solvation energy for DFT,
  ! ... which are derived from Laue-RISM.
  ! ...
  ! ... outside of the unit-cell, charge density is approximated as
  !
  !                   ----
  !        rho(gxy,z) >    q(v) * rho(v) * h(v; gxy,z)
  !                   ----
  !                     v
  !
  ! ... Variables:
  ! ...   charge:     total charge of solvent system
  ! ...   ireference: reference position of solvation potential
  !
  USE cell_base,      ONLY : at, alat
  USE constants,      ONLY : eps8
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE lauefft,        ONLY : fw_lauefft_2xy
  USE mp,             ONLY : mp_sum
  USE rism,           ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, solVs, iuniq_to_nsite, &
                           & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: charge
  INTEGER,         INTENT(IN)    :: ireference
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iq
  INTEGER                  :: iiq
  INTEGER                  :: iv
  INTEGER                  :: nv
  INTEGER                  :: isolV
  INTEGER                  :: iatom
  INTEGER                  :: irz
  INTEGER                  :: iirz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  INTEGER                  :: kgxy
  INTEGER                  :: izright_tail
  INTEGER                  :: izleft_tail
  REAL(DP)                 :: rhov1
  REAL(DP)                 :: rhov2
  REAL(DP)                 :: qv
  REAL(DP)                 :: ntmp
  REAL(DP)                 :: dz
  REAL(DP)                 :: area_xy
  REAL(DP)                 :: vol, dvol
  REAL(DP)                 :: voltmp
  REAL(DP)                 :: fac1
  REAL(DP)                 :: fac2
  REAL(DP)                 :: charge0
  REAL(DP)                 :: chgtmp
  REAL(DP)                 :: rhog0
  REAL(DP)                 :: vsol0
  REAL(DP),    ALLOCATABLE :: wei(:)
  COMPLEX(DP), ALLOCATABLE :: ggz(:,:)
  !
  REAL(DP),    PARAMETER   :: RHO_THR   = 1.0E-16_DP
  REAL(DP),    PARAMETER   :: RHO_SMEAR = 1.0_DP  ! in bohr
  REAL(DP),    PARAMETER   :: WEI_THR   = 1.0E-32_DP
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
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
  IF(rismt%lfft%nrz > 0) THEN
    ALLOCATE(wei(rismt%lfft%nrz))
  END IF
  IF (rismt%nrzs * rismt%ngxy * rismt%nsite > 0) THEN
    ALLOCATE(ggz(rismt%nrzs * rismt%ngxy, rismt%nsite))
  END IF
  !
  ! ... set variables
  dz      = rismt%lfft%zstep * alat
  area_xy = ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)) * alat * alat
  dvol    = area_xy * dz
  !
  ! ... gr -> ggz
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    IF (rismt%nrzs * rismt%ngxy > 0) THEN
      ggz(:, iiq) = C_ZERO
    END IF
    !
    IF (rismt%nr > 0 .AND. (rismt%nrzs * rismt%ngxy) > 0) THEN
      CALL fw_lauefft_2xy(rismt%lfft, rismt%gr(:, iiq), ggz(:, iiq), rismt%nrzs, 1)
    END IF
  END DO
  !
  ! ... make nsol, qsol
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    rhov1 = DBLE(nv) * solVs(isolV)%density
    rhov2 = DBLE(nv) * solVs(isolV)%subdensity
    qv    = solVs(isolV)%charge(iatom)
    !
    rismt%nsol(iiq) = 0.0_DP
    rismt%qsol(iiq) = 0.0_DP
    !
    IF (rismt%lfft%gxystart > 1) THEN
      fac1 = rhov1 * dvol
      fac2 = rhov2 * dvol
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = 1, (rismt%lfft%izleft_start - 1)
        ntmp = ntmp + fac2 * (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      rismt%nsol(iiq) = rismt%nsol(iiq) + ntmp
      rismt%qsol(iiq) = rismt%qsol(iiq) + qv * ntmp
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz, iirz) reduction(+:ntmp)
      DO irz = rismt%lfft%izleft_start, rismt%lfft%izleft_gedge
        iirz = irz - rismt%lfft%izcell_start + 1
        ntmp = ntmp + fac2 * DBLE(ggz(iirz, iiq))
      END DO
!$omp end parallel do
      rismt%nsol(iiq) = rismt%nsol(iiq) + ntmp
      rismt%qsol(iiq) = rismt%qsol(iiq) + qv * ntmp
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz, iirz) reduction(+:ntmp)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%izright_end
        iirz = irz - rismt%lfft%izcell_start + 1
        ntmp = ntmp + fac1 * DBLE(ggz(iirz, iiq))
      END DO
!$omp end parallel do
      rismt%nsol(iiq) = rismt%nsol(iiq) + ntmp
      rismt%qsol(iiq) = rismt%qsol(iiq) + qv * ntmp
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = (rismt%lfft%izright_end + 1), rismt%lfft%nrz
        ntmp = ntmp + fac1 * (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      rismt%nsol(iiq) = rismt%nsol(iiq) + ntmp
      rismt%qsol(iiq) = rismt%qsol(iiq) + qv * ntmp
      !
    END IF
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(rismt%nsol, rismt%mp_site%intra_sitg_comm)
    CALL mp_sum(rismt%qsol, rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... make qtot
  rismt%qtot = 0.0_DP
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    rismt%qtot = rismt%qtot + rismt%qsol(iiq)
  END DO
  !
  CALL mp_sum(rismt%qtot, rismt%mp_site%inter_sitg_comm)
  !
  ! ... make rhog (which is Laue-rep.)
  IF (rismt%nrzl * rismt%ngxy > 0) THEN
    rismt%rhog(:) = C_ZERO
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    rhov1 = DBLE(nv) * solVs(isolV)%density
    rhov2 = DBLE(nv) * solVs(isolV)%subdensity
    qv    = solVs(isolV)%charge(iatom)
    !
    DO igxy = 1, rismt%ngxy
      jgxy = (igxy - 1) * rismt%nrzl
      kgxy = (igxy - 1) * rismt%nrzs
      !
!$omp parallel do default(shared) private(irz)
      DO irz = 1, (rismt%lfft%izleft_start - 1)
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov2 * (rismt%hsgz(irz + jgxy, iiq) + rismt%hlgz(irz + jgxy, iiq))
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz, iirz)
      DO irz = rismt%lfft%izleft_start, rismt%lfft%izleft_gedge
        iirz = irz - rismt%lfft%izcell_start + 1
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov2 * ggz(iirz + kgxy, iiq)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz, iirz)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%izright_end
        iirz = irz - rismt%lfft%izcell_start + 1
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov1 * ggz(iirz + kgxy, iiq)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz)
      DO irz = (rismt%lfft%izright_end + 1), rismt%lfft%nrz
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov1 * (rismt%hsgz(irz + jgxy, iiq) + rismt%hlgz(irz + jgxy, iiq))
      END DO
!$omp end parallel do
      !
    END DO
  END DO
  !
  IF (rismt%nrzl * rismt%ngxy > 0) THEN
    CALL mp_sum(rismt%rhog, rismt%mp_site%inter_sitg_comm)
  END IF
  !
  ! ... detect truncating positions
  izright_tail = 0
  izleft_tail  = 0
  !
  IF (rismt%lfft%gxystart > 1) THEN
    izleft_tail = 1
    DO irz = 1, rismt%lfft%izleft_gedge
      IF (ABS(rismt%rhog(irz)) > RHO_THR) THEN
        izleft_tail = irz
        EXIT
      END IF
    END DO
    !
    izright_tail = rismt%lfft%nrz
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      iirz = rismt%lfft%nrz + rismt%lfft%izright_gedge - irz
      IF (ABS(rismt%rhog(iirz)) > RHO_THR) THEN
        izright_tail = iirz
        EXIT
      END IF
    END DO
    !
  END IF
  !
  CALL mp_sum(izright_tail, rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(izleft_tail,  rismt%mp_site%intra_sitg_comm)
  !
  ! ... calculate weight
  IF(rismt%lfft%nrz > 0) THEN
    wei = 0.0_DP
  END IF
  !
!$omp parallel do default(shared) private(irz)
  DO irz = 1, rismt%lfft%izleft_gedge
    wei(irz) = 0.5_DP * erfc(DBLE(izleft_tail - irz ) * dz / RHO_SMEAR)
    IF (wei(irz) < WEI_THR) THEN
      wei(irz) = 0.0_DP
    END IF
  END DO
!$omp end parallel do
  !
!$omp parallel do default(shared) private(irz)
  DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
    wei(irz) = 0.5_DP * erfc(DBLE(irz - izright_tail) * dz / RHO_SMEAR)
    IF (wei(irz) < WEI_THR) THEN
      wei(irz) = 0.0_DP
    END IF
  END DO
!$omp end parallel do
  !
  ! ... evaluate volume
  vol = 0.0_DP
  !
  IF (rismt%lfft%gxystart > 1) THEN
    !
    voltmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:voltmp)
    DO irz = 1, rismt%lfft%izleft_gedge
      voltmp = voltmp + dvol * wei(irz)
    END DO
!$omp end parallel do
    vol = vol + voltmp
    !
    voltmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:voltmp)
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      voltmp = voltmp + dvol * wei(irz)
    END DO
!$omp end parallel do
    vol = vol + voltmp
    !
  END IF
  !
  CALL mp_sum(vol, rismt%mp_site%intra_sitg_comm)
  !
  ! ... evaluate total charge
  charge0 = 0.0_DP
  !
  IF (rismt%lfft%gxystart > 1) THEN
    !
    chgtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:chgtmp)
    DO irz = 1, rismt%lfft%izleft_gedge
      chgtmp = chgtmp + dvol * wei(irz) * rismt%rhog(irz)
    END DO
!$omp end parallel do
    charge0 = charge0 + chgtmp
    !
    chgtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:chgtmp)
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      chgtmp = chgtmp + dvol * wei(irz) * rismt%rhog(irz)
    END DO
!$omp end parallel do
    charge0 = charge0 + chgtmp
    !
  END IF
  !
  CALL mp_sum(charge0, rismt%mp_site%intra_sitg_comm)
  !
  ! ... renormalize rhog
  IF (rismt%lfft%gxystart > 1) THEN
    IF (ABS(vol) <= eps8) THEN  ! will not be occurred
      CALL errore('solvation_lauerism', 'vol is zero', 1)
    END IF
    rhog0 = (charge - charge0) / vol
    !
!$omp parallel do default(shared) private(irz)
    DO irz = 1, rismt%lfft%izleft_gedge
      rismt%rhog(irz) = (rismt%rhog(irz) + CMPLX(rhog0, 0.0_DP, kind=DP)) * wei(irz)
    END DO
!$omp end parallel do
    !
!$omp parallel do default(shared) private(irz)
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      rismt%rhog(irz) = (rismt%rhog(irz) + CMPLX(rhog0, 0.0_DP, kind=DP)) * wei(irz)
    END DO
!$omp end parallel do
    !
  END IF
  !
  WRITE(stdout, '(/,5X,"solvent charge ",F10.5, &
                  & ", renormalised to ",F10.5)') charge0, charge
  !
  ! ... make vpot
  CALL solvation_esm_potential(rismt, ireference, vsol0, ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... make rhog_pbc and vpot_pbc
  CALL solvation_pbc(rismt, ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... make esol
  rismt%esol = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    ! this version only supports G.F.
    !rismt%esol = rismt%esol + rismt%usol(iiq)
    rismt%esol = rismt%esol + rismt%usol_GF(iiq)
  END DO
  !
  CALL mp_sum(rismt%esol, rismt%mp_site%inter_sitg_comm)
  !
  ! ... make vsol (reference level shifting)
  rismt%vsol = vsol0
  !
  ! ... deallocate memory
  IF(rismt%lfft%nrz > 0) THEN
    DEALLOCATE(wei)
  END IF
  IF (rismt%nrzs * rismt%ngxy * rismt%nsite > 0) THEN
    DEALLOCATE(ggz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_lauerism
