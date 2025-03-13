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
SUBROUTINE suscept_laue(rism1t, rismlt, alpha, lhand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... create inter-site susceptibility for Laue-RISM (from 1D-RISM)
  ! ...
  ! ...                     1  / inf
  ! ...   x21(gxy,z'-z) = ---- | dgz [ w21(g) + rho2 * h21(g) ] * cos(gz*(z'-z))
  ! ...                    pi  / 0
  ! ...
  ! ... x21 depends on norm of gxy, and is even along z'-z.
  !
  USE constants, ONLY : pi, sqrtpi, eps12
  USE cell_base, ONLY : alat, tpiba2
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE, &
                      & IERR_RISM_1DRISM_IS_NOT_AVAIL, IERR_RISM_LARGE_LAUE_BOX
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_size, mp_rank, mp_max, mp_sum, mp_get, mp_gather, mp_bcast, mp_barrier
  USE rism,      ONLY : rism_type, ITYPE_1DRISM, ITYPE_LAUERISM
  USE solvavg,   ONLY : solvavg_init, solvavg_clear, solvavg_print, solvavg_put
  USE solvmol,   ONLY : nsolV, solVs, get_nsite_in_solVs, get_nuniq_in_solVs, &
                      & iuniq_to_nsite, iuniq_to_isite, isite_to_isolV, isite_to_iatom
  USE splinelib, ONLY : spline, splint
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)    :: rism1t
  TYPE(rism_type), INTENT(INOUT) :: rismlt
  REAL(DP),        INTENT(IN)    :: alpha  ! in bohr
  LOGICAL,         INTENT(IN)    :: lhand  ! if true, right-hand. if false, left-hand.
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: nv
  INTEGER               :: nq
  INTEGER               :: iq1, iq2
  INTEGER               :: iv1, iv2
  INTEGER               :: iw1, iw2
  INTEGER               :: isolV1, isolV2
  INTEGER               :: iatom1, iatom2
  CHARACTER(LEN=6)      :: satom1, satom2
  INTEGER               :: iiq2
  INTEGER               :: nv2, iiv2
  REAL(DP)              :: qv2
  INTEGER               :: ivv
  INTEGER               :: irz
  INTEGER               :: igz
  INTEGER               :: igxy
  INTEGER               :: jgxy
  INTEGER               :: irank
  INTEGER               :: nproc
  INTEGER,  ALLOCATABLE :: rank_map(:,:)
  INTEGER,  ALLOCATABLE :: root_spline(:)
  REAL(DP)              :: rfft_1d
  REAL(DP)              :: rfft_laue
  REAL(DP)              :: rho2
  REAL(DP)              :: rz
  REAL(DP)              :: gz, ggz
  REAL(DP)              :: ggxy
  REAL(DP)              :: gs
  REAL(DP)              :: gsmax
  REAL(DP)              :: dg
  REAL(DP)              :: pidg
  REAL(DP)              :: xgs
  REAL(DP)              :: dxg0
  REAL(DP)              :: ddxg0
  REAL(DP), ALLOCATABLE :: xg_1d(:)
  REAL(DP), ALLOCATABLE :: xg_spl(:)
  REAL(DP), ALLOCATABLE :: xg_d2y(:)
  REAL(DP), ALLOCATABLE :: gs_t(:,:)
  REAL(DP), ALLOCATABLE :: xg_t(:,:)
  REAL(DP), ALLOCATABLE :: xg_0(:,:)
  REAL(DP), ALLOCATABLE :: xgs21(:)
  REAL(DP), ALLOCATABLE :: cosgz(:,:)
  !
  REAL(DP), PARAMETER   :: LAUE_BOX_SCALE = 1.2_DP
  !
  EXTERNAL :: dgemm
  !
  ! ... number of sites in solvents
  nv = get_nsite_in_solVs()
  nq = get_nuniq_in_solVs()
  !
  IF (rism1t%itype /= ITYPE_1DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rism1t%nr /= rism1t%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rism1t%nsite < (nv * (nv + 1) / 2)) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rism1t%rfft%ngrid < rism1t%mp_task%nvec) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (.NOT. rism1t%avail) THEN
    ierr = IERR_RISM_1DRISM_IS_NOT_AVAIL
    RETURN
  END IF
  !
  ! ... check data type of Laue-RISM
  IF (rismlt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismlt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismlt%ngs < rismlt%lfft%nglxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismlt%nrzl < rismlt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... check alpha
  IF (alpha <= 0.0_DP) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... check FFT-range
  rfft_1d   = rism1t%rfft%rgrid(rism1t%rfft%ngrid)
  rfft_laue = (rismlt%lfft%zright - rismlt%lfft%zleft) * alat
  IF ((LAUE_BOX_SCALE * rfft_laue) > rfft_1d) THEN
    ierr = IERR_RISM_LARGE_LAUE_BOX
    RETURN
  END IF
  !
  ! ... allocate working memory
  ALLOCATE(rank_map(3, nq))
  ALLOCATE(root_spline(nq))
  ALLOCATE(xg_1d(rism1t%ng))
  ALLOCATE(xg_spl(rism1t%mp_task%nvec))
  ALLOCATE(xg_d2y(rism1t%mp_task%nvec))
  IF ((rism1t%rfft%ngrid * rismlt%lfft%nglxy) > 0) THEN
    ALLOCATE(gs_t(rism1t%rfft%ngrid, rismlt%lfft%nglxy))
  END IF
  IF ((rism1t%rfft%ngrid * rismlt%lfft%nglxy) > 0) THEN
    ALLOCATE(xg_t(rism1t%rfft%ngrid, rismlt%lfft%nglxy))
  END IF
  IF (rismlt%nsite > 0) THEN
    ALLOCATE(xg_0(rismlt%nsite, nq))
  END IF
  IF ((rismlt%nrzl * rismlt%ngs) > 0) THEN
    ALLOCATE(xgs21(rismlt%nrzl * rismlt%ngs))
  END IF
  IF ((rism1t%rfft%ngrid * rismlt%lfft%nrz) > 0) THEN
    ALLOCATE(cosgz(rism1t%rfft%ngrid, rismlt%lfft%nrz))
  END IF
  !
  ! ... setup roots to prepare spline
  DO iq1 = 1, nq
    rank_map(1, iq1) = mp_rank(rismlt%intra_comm)
    rank_map(2, iq1) = 0
    rank_map(3, iq1) = 0
    !
    root_spline(iq1) = 0
    IF (rismlt%mp_site%isite_start <= iq1 .AND. iq1 <= rismlt%mp_site%isite_end) THEN
      IF (rismlt%mp_site%me_sitg == rismlt%mp_site%root_sitg) THEN
        rank_map(2, iq1) = rank_map(1, iq1) + 1
        IF (rism1t%is_intra) THEN
          root_spline(iq1) = rism1t%mp_task%me_task + 1
        END IF
      END IF
    END IF
    !
    CALL mp_max(rank_map(2, iq1), rismlt%intra_comm)
    rank_map(2, iq1) = rank_map(2, iq1) - 1
    !
    CALL mp_max(root_spline(iq1), rismlt%intra_comm)
    root_spline(iq1) = root_spline(iq1) - 1
    IF (root_spline(iq1) < 0) THEN
      root_spline(iq1) = rism1t%mp_task%root_task
    END IF
    !
    IF (rism1t%is_intra .AND. rism1t%mp_task%me_task == root_spline(iq1)) THEN
      rank_map(3, iq1) = rank_map(1, iq1) + 1
    END IF
    !
    CALL mp_max(rank_map(3, iq1), rismlt%intra_comm)
    rank_map(3, iq1) = rank_map(3, iq1) - 1
  END DO
  !
  ! ... set variables
  gsmax = rism1t%rfft%ggrid(rism1t%rfft%ngrid)
  dg    = rism1t%rfft%ggrid(2) - rism1t%rfft%ggrid(1)
  pidg  = dg / pi
  !
  ! ... calculate list of gs
  DO igxy = 1, rismlt%lfft%nglxy
    ggxy = tpiba2 * rismlt%lfft%glxy(igxy)
!$omp parallel do default(shared) private(igz, gz, ggz, gs)
    DO igz = 1, rism1t%rfft%ngrid
      gz  = rism1t%rfft%ggrid(igz)
      ggz = gz * gz
      gs  = SQRT(ggxy + ggz)
      gs_t(igz, igxy) = gs
    END DO
!$omp end parallel do
  END DO
  !
  ! ... calculate cos(gz*rz)
  cosgz = 0.0_DP
  irank = mp_rank(rismlt%intra_comm)
  nproc = mp_size(rismlt%intra_comm)
  !
  DO irz = 1, rismlt%lfft%nrz
    IF (irank /= MOD(irz - 1, nproc)) THEN
      CYCLE
    END IF
    !
    rz = alat * DBLE(irz - 1) * rismlt%lfft%zstep
!$omp parallel do default(shared) private(igz, gz)
    DO igz = 1, rism1t%rfft%ngrid
      gz = rism1t%rfft%ggrid(igz)
      cosgz(igz, irz) = COS(gz * rz)
    END DO
!$omp end parallel do
    cosgz(1, irz) = 0.5_DP * cosgz(1, irz)
  END DO
  !
  CALL mp_sum(cosgz, rismlt%intra_comm)
  !
  ! ... calculate susceptibility
  DO iq1 = 1, nq
    ! ... properties of unique site1
    iv1 = iuniq_to_isite(1, iq1)
    !
    DO iq2 = 1, nq
      ! ... properties of unique site2
      nv2 = iuniq_to_nsite(iq2)
      IF (rismlt%mp_site%isite_start <= iq2 .AND. iq2 <= rismlt%mp_site%isite_end) THEN
        iiq2 = iq2 - rismlt%mp_site%isite_start + 1
      ELSE
        iiq2 = 0
      END IF
      !
      IF (iiq2 > 0) THEN
        IF (rismlt%nsite > 0) THEN
          xg_0(iiq2, iq1) = 0.0_DP
        END IF
        IF ((rismlt%nrzl * rismlt%ngs) > 0) THEN
          xgs21(:) = 0.0_DP
        END IF
      END IF
      !
      DO iiv2 = 1, nv2
        ! ... properties of a site2
        iv2    = iuniq_to_isite(iiv2, iq2)
        isolV2 = isite_to_isolV(iv2)
        IF (lhand) THEN
          rho2 = solVs(isolV2)%density
        ELSE
          rho2 = solVs(isolV2)%subdensity
        END IF
        !
        iw1 = MAX(iv1, iv2)
        iw2 = MIN(iv1, iv2)
        ivv = iw1 * (iw1 - 1) / 2 + iw2
        !
        ! ... create h21(g) or x21(g) of 1D-RISM
        IF (rism1t%is_intra) THEN
          IF (iv1 == iv2) THEN
            xg_1d(:) = rho2 * rism1t%hg(:, ivv)
          ELSE
            xg_1d(:) = rism1t%wg(:, ivv) + rho2 * rism1t%hg(:, ivv)
          END IF
          CALL mp_gather(xg_1d, xg_spl, rism1t%mp_task%ilen_vecs, rism1t%mp_task%idis_vecs, &
                       & root_spline(iq2), rism1t%mp_task%itask_comm)
        END IF
        !
        IF (rank_map(2, iq2) /= rank_map(3, iq2)) THEN
          CALL mp_get(xg_spl, xg_spl, rank_map(1, iq2), &
          & rank_map(2, iq2), rank_map(3, iq2), ivv, rismlt%intra_comm)
        END IF
        !
        IF (iiq2 > 0) THEN
          ! ... prepare spline correction
          IF (rismlt%mp_site%me_sitg == rismlt%mp_site%root_sitg) THEN
            CALL suscept_g0(rism1t%mp_task%nvec, rism1t%rfft%ggrid, xg_spl, dxg0, ddxg0)
            CALL spline(rism1t%rfft%ggrid(1:rism1t%mp_task%nvec), xg_spl, dxg0, ddxg0, xg_d2y)
          END IF
          !
          CALL mp_bcast(xg_spl, rismlt%mp_site%root_sitg, rismlt%mp_site%intra_sitg_comm)
          CALL mp_bcast(xg_d2y, rismlt%mp_site%root_sitg, rismlt%mp_site%intra_sitg_comm)
          !
          ! ... perform spline correction fitting h21(g) or x21(g) from 1D-RISM to Laue-RISM
          DO igxy = 1, rismlt%lfft%nglxy
!$omp parallel do default(shared) private(igz, gs, xgs)
            DO igz = 1, rism1t%rfft%ngrid
              gs = gs_t(igz, igxy)
              IF (gs <= (gsmax + eps12)) THEN
                xgs = splint(rism1t%rfft%ggrid(1:rism1t%mp_task%nvec), xg_spl, xg_d2y, gs)
                xg_t(igz, igxy) = xgs
              ELSE
                xg_t(igz, igxy) = 0.0_DP
              END IF
            END DO
!$omp end parallel do
          END DO
          !
          ! ... calculate x21(rz,gxy)
          IF (rismlt%nsite > 0) THEN
            IF (iv1 == iv2) THEN
              xg_0(iiq2, iq1) = xg_0(iiq2, iq1) + xg_spl(1) + 1.0_DP
            ELSE
              xg_0(iiq2, iq1) = xg_0(iiq2, iq1) + xg_spl(1)
            END IF
          END IF
          !
          IF ((rismlt%nrzl * rismlt%ngs) > 0) THEN
            CALL dgemm('T', 'N', rismlt%lfft%nrz, rismlt%lfft%nglxy, rism1t%rfft%ngrid, &
                     & pidg, cosgz, rism1t%rfft%ngrid, xg_t, rism1t%rfft%ngrid, &
                     & 1.0_DP, xgs21, rismlt%nrzl)
          END IF
          !
          IF (iv1 == iv2) THEN
            DO igxy = 1, rismlt%lfft%nglxy
              jgxy = (igxy - 1) * rismlt%nrzl
              ggxy = tpiba2 * rismlt%lfft%glxy(igxy)
!$omp parallel do default(shared) private(irz, rz)
              DO irz = 1, rismlt%lfft%nrz
                rz = alat * DBLE(irz - 1) * rismlt%lfft%zstep
                xgs21(irz + jgxy) = xgs21(irz + jgxy) + &
                & EXP(-rz * rz / alpha / alpha - 0.25_DP * alpha * alpha * ggxy) / alpha / sqrtpi
              END DO
!$omp end parallel do
            END DO
          END IF
        END IF
        !
        CALL mp_barrier(rismlt%intra_comm)
        !
      END DO
      !
      IF (iiq2 > 0) THEN
        IF ((rismlt%nrzl * rismlt%ngs) > 0) THEN
          IF (lhand) THEN
            rismlt%xgs(:, iiq2, iq1) = xgs21(:)
          ELSE
            rismlt%ygs(:, iiq2, iq1) = xgs21(:)
          END IF
        END IF
      END IF
      !
    END DO
  END DO
  !
  ! ... renormalize at G = 0
  CALL renormalize_g0()
  !
  ! ... correct at Gxy = 0
  CALL correct_gxy0()
  !
  ! ... print data
  CALL print_x21()
  !
  ! ... deallocate working memory
  DEALLOCATE(rank_map)
  DEALLOCATE(root_spline)
  DEALLOCATE(xg_1d)
  DEALLOCATE(xg_spl)
  DEALLOCATE(xg_d2y)
  IF ((rism1t%rfft%ngrid * rismlt%lfft%nglxy) > 0) THEN
    DEALLOCATE(gs_t)
  END IF
  IF ((rism1t%rfft%ngrid * rismlt%lfft%nglxy) > 0) THEN
    DEALLOCATE(xg_t)
  END IF
  IF (rismlt%nsite > 0) THEN
    DEALLOCATE(xg_0)
  END IF
  IF ((rismlt%nrzl * rismlt%ngs) > 0) THEN
    DEALLOCATE(xgs21)
  END IF
  IF ((rism1t%rfft%ngrid * rismlt%lfft%nrz) > 0) THEN
    DEALLOCATE(cosgz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE renormalize_g0()
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE :: msol(:)
    REAL(DP), ALLOCATABLE :: qsol(:)
    REAL(DP)              :: qtot
    REAL(DP)              :: qsqr
    !
    ALLOCATE(msol(nsolV))
    ALLOCATE(qsol(nsolV))
    !
    DO iq1 = 1, nq
      !
      ! ... sum numbers and charges of solvent atoms in a molecule
      msol = 0.0_DP
      qsol = 0.0_DP
      !
      DO iq2 = rismlt%mp_site%isite_start, rismlt%mp_site%isite_end
        iiq2   = iq2 - rismlt%mp_site%isite_start + 1
        iv2    = iuniq_to_isite(1, iq2)
        nv2    = iuniq_to_nsite(iq2)
        isolV2 = isite_to_isolV(iv2)
        iatom2 = isite_to_iatom(iv2)
        qv2    = solVs(isolV2)%charge(iatom2)
        !
        msol(isolV2) = msol(isolV2) + xg_0(iiq2, iq1)
        qsol(isolV2) = qsol(isolV2) + DBLE(nv2) * qv2
      END DO
      !
      CALL mp_sum(msol, rismlt%mp_site%inter_sitg_comm)
      CALL mp_sum(qsol, rismlt%mp_site%inter_sitg_comm)
      !
      DO isolV2 = 1, nsolV
        IF (solVs(isolV2)%natom > 0) THEN
          msol(isolV2) = msol(isolV2) / DBLE(solVs(isolV2)%natom)
          qsol(isolV2) = qsol(isolV2) / DBLE(solVs(isolV2)%natom)
        ELSE
          msol(isolV2) = 0.0_DP
          qsol(isolV2) = 0.0_DP
        END IF
      END DO
      !
      ! ... renormalize: to correct stoichiometry
      DO iq2 = rismlt%mp_site%isite_start, rismlt%mp_site%isite_end
        iiq2   = iq2 - rismlt%mp_site%isite_start + 1
        iv2    = iuniq_to_isite(1, iq2)
        nv2    = iuniq_to_nsite(iq2)
        isolV2 = isite_to_isolV(iv2)
        !
        xg_0(iiq2, iq1) = DBLE(nv2) * msol(isolV2)
      END DO
      !
      ! ... total charge and square sum of charge
      qtot = 0.0_DP
      qsqr = 0.0_DP
      !
      DO iq2 = rismlt%mp_site%isite_start, rismlt%mp_site%isite_end
        iiq2   = iq2 - rismlt%mp_site%isite_start + 1
        iv2    = iuniq_to_isite(1, iq2)
        nv2    = iuniq_to_nsite(iq2)
        isolV2 = isite_to_isolV(iv2)
        !
        qtot = qtot + qsol(isolV2) * xg_0(iiq2, iq1)
        qsqr = qsqr + DBLE(nv2) * qsol(isolV2) * qsol(isolV2)
      END DO
      !
      CALL mp_sum(qtot, rismlt%mp_site%inter_sitg_comm)
      CALL mp_sum(qsqr, rismlt%mp_site%inter_sitg_comm)
      !
      ! ... renormalize: to correct total charge
      IF (ABS(qtot) > eps12) THEN
        IF (ABS(qsqr) <= eps12) THEN  ! will not be occurred
          CALL errore('renormalize_g0', 'qsqr is zero', 1)
        END IF
        !
        DO iq2 = rismlt%mp_site%isite_start, rismlt%mp_site%isite_end
          iiq2   = iq2 - rismlt%mp_site%isite_start + 1
          iv2    = iuniq_to_isite(1, iq2)
          nv2    = iuniq_to_nsite(iq2)
          isolV2 = isite_to_isolV(iv2)
          !
          xg_0(iiq2, iq1) = xg_0(iiq2, iq1) - DBLE(nv2) * qsol(isolV2) * qtot / qsqr
        END DO
      END IF
      !
    END DO
    !
    DEALLOCATE(msol)
    DEALLOCATE(qsol)
    !
  END SUBROUTINE renormalize_g0
  !
  SUBROUTINE correct_gxy0()
    IMPLICIT NONE
    REAL(DP) :: dz
    REAL(DP) :: xg_int
    REAL(DP) :: xg_scale
    !
    dz = alat * rismlt%lfft%zstep
    !
    DO iq1 = 1, nq
      DO iq2 = rismlt%mp_site%isite_start, rismlt%mp_site%isite_end
        iiq2 = iq2 - rismlt%mp_site%isite_start + 1
        !
        IF (rismlt%lfft%gxystart > 1) THEN
          !
          ! ... integrate at Gxy = 0
          IF (lhand) THEN
            xg_int = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:xg_int)
            DO irz = 2, rismlt%lfft%nrz
              xg_int = xg_int + 2.0_DP * dz * rismlt%xgs(irz, iiq2, iq1)
            END DO
!$omp end parallel do
            xg_int = xg_int + dz * rismlt%xgs(1, iiq2, iq1)
          ELSE
            xg_int = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:xg_int)
            DO irz = 2, rismlt%lfft%nrz
              xg_int = xg_int + 2.0_DP * dz * rismlt%ygs(irz, iiq2, iq1)
            END DO
!$omp end parallel do
            xg_int = xg_int + dz * rismlt%ygs(1, iiq2, iq1)
          END IF
          !
          ! ... rescale x21 at Gxy = 0
          IF (ABS(xg_int) > eps12) THEN
            xg_scale = xg_0(iiq2, iq1) / xg_int
            IF (lhand) THEN
!$omp parallel do default(shared) private(irz)
              DO irz = 1, rismlt%lfft%nrz
                rismlt%xgs(irz, iiq2, iq1) = rismlt%xgs(irz, iiq2, iq1) * xg_scale
              END DO
!$omp end parallel do
            ELSE
!$omp parallel do default(shared) private(irz)
              DO irz = 1, rismlt%lfft%nrz
                rismlt%ygs(irz, iiq2, iq1) = rismlt%ygs(irz, iiq2, iq1) * xg_scale
              END DO
!$omp end parallel do
            END IF
          END IF
          !
        END IF
      END DO
    END DO
    !
  END SUBROUTINE correct_gxy0
  !
  SUBROUTINE print_x21()
    IMPLICIT NONE
#if defined (__DEBUG_RISM)
    INTEGER                  :: ista
    INTEGER                  :: my_group_id
    INTEGER                  :: io_group_id
    INTEGER                  :: owner_group_id
    COMPLEX(DP), ALLOCATABLE :: xtmp(:)
    !
    ! ... get process info.
    my_group_id = mp_rank(rismlt%mp_site%inter_sitg_comm)
    !
    ! ... find the index of the group which includes ionode
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, rismlt%mp_site%intra_sitg_comm)
    CALL mp_sum(io_group_id, rismlt%mp_site%inter_sitg_comm)
    !
    ! ... init solvavg
    IF (my_group_id == io_group_id) THEN
      CALL solvavg_init(rismlt%lfft, rismlt%mp_site%intra_sitg_comm, .TRUE.)
    END IF
    !
    ! ... put data to solvavg
    ALLOCATE(xtmp(rismlt%nrzl * rismlt%ngs))
    !
    DO iq1 = 1, nq
      iv1    = iuniq_to_isite(1, iq1)
      isolV1 = isite_to_isolV(iv1)
      iatom1 = isite_to_iatom(iv1)
      satom1 = ADJUSTL(solVs(isolV1)%aname(iatom1))
      !
      DO iq2 = 1, nq
        iv2    = iuniq_to_isite(1, iq2)
        isolV2 = isite_to_isolV(iv2)
        iatom2 = isite_to_iatom(iv2)
        satom2 = ADJUSTL(solVs(isolV2)%aname(iatom2))
        !
        IF (rismlt%mp_site%isite_start <= iq2 .AND. iq2 <= rismlt%mp_site%isite_end) THEN
          owner_group_id = my_group_id
          IF (lhand) THEN
            xtmp = rismlt%xgs(:, iq2 - rismlt%mp_site%isite_start + 1, iq1)
          ELSE
            xtmp = rismlt%ygs(:, iq2 - rismlt%mp_site%isite_start + 1, iq1)
          END IF
        ELSE
          owner_group_id = 0
          xtmp = CMPLX(0.0_DP, 0.0_DP, kind=DP)
        END IF
        !
        CALL mp_sum(owner_group_id, rismlt%mp_site%inter_sitg_comm)
        CALL mp_get(xtmp, xtmp, my_group_id, io_group_id, &
                  & owner_group_id, iq2, rismlt%mp_site%inter_sitg_comm)
        !
        IF (my_group_id == io_group_id) THEN
          CALL solvavg_put('x_' // TRIM(satom2) // ':' // TRIM(satom1), .FALSE., xtmp, rismlt%nrzl, .TRUE.)
        END IF
        !
        CALL mp_barrier(rismlt%mp_site%inter_sitg_comm)
      END DO
    END DO
    !
    DEALLOCATE(xtmp)
    !
    ! ... print solvavg
    IF (my_group_id == io_group_id) THEN
      IF (lhand) THEN
        CALL solvavg_print(TRIM(tmp_dir) // TRIM(prefix) // '.rism1_x21', 'solvent susceptibility', ista)
      ELSE
        CALL solvavg_print(TRIM(tmp_dir) // TRIM(prefix) // '.rism1_y21', 'solvent susceptibility', ista)
      END IF
      ista = ABS(ista)
    ELSE
      ista = 0
    END IF
    !
    IF (ista /= 0) THEN
      CALL errore('print_x21', 'cannot write file', ista)
    END IF
    !
    ! ... finalize solvavg
    IF (my_group_id == io_group_id) THEN
      CALL solvavg_clear()
    END IF
    !
#endif
  END SUBROUTINE print_x21
  !
END SUBROUTINE suscept_laue
