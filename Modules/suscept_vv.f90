!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE suscept_vv(rism1t, rism3t, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... create inter-site susceptibility for 3D-RISM (from 1D-RISM)
  ! ...   x21(g) = w21(g) + rho2 * h21(g)
  ! ...
  ! ... (A.Kovalenko, F.Hirata, Chem. Phys. Lett. 1998, 290, 237-244)
  !
  USE cell_base, ONLY : tpiba
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE, &
                      & IERR_RISM_1DRISM_IS_NOT_AVAIL
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_rank, mp_max, mp_get, mp_gather, mp_bcast, mp_barrier
  USE rism,      ONLY : rism_type, ITYPE_1DRISM, ITYPE_3DRISM
  USE solvmol,   ONLY : solVs, get_nsite_in_solVs, get_nuniq_in_solVs, &
                      & iuniq_to_nsite, iuniq_to_isite, isite_to_isolV
  USE splinelib, ONLY : spline, splint
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)    :: rism1t
  TYPE(rism_type), INTENT(INOUT) :: rism3t
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: nv
  INTEGER               :: nq
  INTEGER               :: iq1, iq2
  INTEGER               :: iv1, iv2
  INTEGER               :: iw1, iw2
  INTEGER               :: iiq2
  INTEGER               :: nv2, iiv2
  INTEGER               :: isolV2
  INTEGER               :: ivv
  INTEGER               :: igs
  INTEGER,  ALLOCATABLE :: rank_map(:,:)
  INTEGER,  ALLOCATABLE :: root_spline(:)
  REAL(DP)              :: rho2
  REAL(DP)              :: gs
  REAL(DP)              :: xgs
  REAL(DP)              :: dxg0
  REAL(DP)              :: ddxg0
  REAL(DP), ALLOCATABLE :: xg_1d(:)
  REAL(DP), ALLOCATABLE :: xg_spl(:)
  REAL(DP), ALLOCATABLE :: xg_d2y(:)
  REAL(DP), ALLOCATABLE :: gs_t(:)
  !
  ! ... number of sites in solvents
  nv = get_nsite_in_solVs()
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type of 1D-RISM
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
  ! ... check data type of 3D-RISM
  IF (rism3t%itype /= ITYPE_3DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rism3t%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rism3t%ngs < rism3t%gvec%ngl) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate working memory
  ALLOCATE(rank_map(3, nq))
  ALLOCATE(root_spline(nq))
  ALLOCATE(xg_1d(rism1t%ng))
  ALLOCATE(xg_spl(rism1t%mp_task%nvec))
  ALLOCATE(xg_d2y(rism1t%mp_task%nvec))
  IF (rism3t%gvec%ngl > 0) THEN
    ALLOCATE(gs_t(rism3t%gvec%ngl))
  END IF
  !
  ! ... setup roots to prepare spline
  DO iq1 = 1, nq
    rank_map(1, iq1) = mp_rank(rism3t%intra_comm)
    rank_map(2, iq1) = 0
    rank_map(3, iq1) = 0
    !
    root_spline(iq1) = 0
    IF (rism3t%mp_site%isite_start <= iq1 .AND. iq1 <= rism3t%mp_site%isite_end) THEN
      IF (rism3t%mp_site%me_sitg == rism3t%mp_site%root_sitg) THEN
        rank_map(2, iq1) = rank_map(1, iq1) + 1
        IF (rism1t%is_intra) THEN
          root_spline(iq1) = rism1t%mp_task%me_task + 1
        END IF
      END IF
    END IF
    !
    CALL mp_max(rank_map(2, iq1), rism3t%intra_comm)
    rank_map(2, iq1) = rank_map(2, iq1) - 1
    !
    CALL mp_max(root_spline(iq1), rism3t%intra_comm)
    root_spline(iq1) = root_spline(iq1) - 1
    IF (root_spline(iq1) < 0) THEN
      root_spline(iq1) = rism1t%mp_task%root_task
    END IF
    !
    IF (rism1t%is_intra .AND. rism1t%mp_task%me_task == root_spline(iq1)) THEN
      rank_map(3, iq1) = rank_map(1, iq1) + 1
    END IF
    !
    CALL mp_max(rank_map(3, iq1), rism3t%intra_comm)
    rank_map(3, iq1) = rank_map(3, iq1) - 1
  END DO
  !
  ! ... calculate list of gs
  DO igs = 1, rism3t%gvec%ngl
    gs = tpiba * SQRT(rism3t%gvec%gl(igs))
    gs_t(igs) = gs
  END DO
  !
  ! ... calculate susceptibility
  DO iq1 = 1, nq
    ! ... properties of unique site1
    iv1 = iuniq_to_isite(1, iq1)
    !
    DO iq2 = 1, nq
      ! ... properties of unique site2
      nv2 = iuniq_to_nsite(iq2)
      IF (rism3t%mp_site%isite_start <= iq2 .AND. iq2 <= rism3t%mp_site%isite_end) THEN
        iiq2 = iq2 - rism3t%mp_site%isite_start + 1
      ELSE
        iiq2 = 0
      END IF
      !
      IF (iiq2 > 0) THEN
        IF (rism3t%ngs > 0) THEN
          rism3t%xgs(:, iiq2, iq1) = 0.0_DP
        END IF
      END IF
      !
      DO iiv2 = 1, nv2
        ! ... properties of a site2
        iv2    = iuniq_to_isite(iiv2, iq2)
        isolV2 = isite_to_isolV(iv2)
        rho2   = solVs(isolV2)%density
        !
        iw1 = MAX(iv1, iv2)
        iw2 = MIN(iv1, iv2)
        ivv = iw1 * (iw1 - 1) / 2 + iw2
        !
        ! ... create x21(g) of 1D-RISM
        IF (rism1t%is_intra) THEN
          xg_1d(:) = rism1t%wg(:, ivv) + rho2 * rism1t%hg(:, ivv)
          CALL mp_gather(xg_1d, xg_spl, rism1t%mp_task%ilen_vecs, rism1t%mp_task%idis_vecs, &
                       & root_spline(iq2), rism1t%mp_task%itask_comm)
        END IF
        !
        IF (rank_map(2, iq2) /= rank_map(3, iq2)) THEN
          CALL mp_get(xg_spl, xg_spl, rank_map(1, iq2), &
          & rank_map(2, iq2), rank_map(3, iq2), ivv, rism3t%intra_comm)
        END IF
        !
        IF (iiq2 > 0) THEN
          ! ... prepare spline correction
          IF (rism3t%mp_site%me_sitg == rism3t%mp_site%root_sitg) THEN
            CALL suscept_g0(rism1t%mp_task%nvec, rism1t%rfft%ggrid, xg_spl, dxg0, ddxg0)
            CALL spline(rism1t%rfft%ggrid(1:rism1t%mp_task%nvec), xg_spl, dxg0, ddxg0, xg_d2y)
          END IF
          !
          CALL mp_bcast(xg_spl, rism3t%mp_site%root_sitg, rism3t%mp_site%intra_sitg_comm)
          CALL mp_bcast(xg_d2y, rism3t%mp_site%root_sitg, rism3t%mp_site%intra_sitg_comm)
          !
          ! ... perform spline correction fitting x21(g) from 1D-RISM to 3D-RISM
!$omp parallel do default(shared) private(igs, gs, xgs)
          DO igs = 1, rism3t%gvec%ngl
            gs  = gs_t(igs)
            xgs = splint(rism1t%rfft%ggrid(1:rism1t%mp_task%nvec), xg_spl, xg_d2y, gs)
            rism3t%xgs(igs, iiq2, iq1) = rism3t%xgs(igs, iiq2, iq1) + xgs
          END DO
!$omp end parallel do
        END IF
        !
        CALL mp_barrier(rism3t%intra_comm)
        !
      END DO
      !
    END DO
  END DO
  !
  ! ... deallocate working memory
  DEALLOCATE(rank_map)
  DEALLOCATE(root_spline)
  DEALLOCATE(xg_1d)
  DEALLOCATE(xg_spl)
  DEALLOCATE(xg_d2y)
  IF (rism3t%gvec%ngl > 0) THEN
    DEALLOCATE(gs_t)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE suscept_vv
