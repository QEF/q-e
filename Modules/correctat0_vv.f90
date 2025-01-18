!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE correctat0_vv(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... correct correlations of 1D-RISM at g = 0 or r = 0.
  ! ...
  ! ... for g = 0,
  ! ...                       / inf
  ! ...   A(g = 0) = 4 * pi * | dr r^2 * A(r)
  ! ...                       / 0
  ! ...
  ! ... for r = 0,
  ! ...                  1        / inf
  ! ...   A(r = 0) = ---------- * | dg g^2 * A(g)
  ! ...               2 * pi^2    / 0
  ! ...
  !
  USE constants, ONLY : pi, tpi, fpi
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_1DRISM
  USE solvmol,   ONLY : get_nsite_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: nv
  INTEGER  :: ivv
  INTEGER  :: ir, iir, jr
  INTEGER  :: ig, iig, jg
  REAL(DP) :: dr
  REAL(DP) :: rfac
  REAL(DP) :: r, rr
  REAL(DP) :: dg
  REAL(DP) :: gfac
  REAL(DP) :: g, gg
  REAL(DP) :: csr0
  REAL(DP) :: csg0
  REAL(DP) :: hr0
  REAL(DP) :: hg0
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
  IF (rismt%nr /= rismt%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nsite < (nv * (nv + 1) / 2)) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set variables
  dr   = rismt%rfft%rgrid(2) - rismt%rfft%rgrid(1)
  dg   = rismt%rfft%ggrid(2) - rismt%rfft%ggrid(1)
  rfac = fpi * dr
  gfac = dg / tpi / pi
  !
  ! ... calculate correlations at G = 0 and R = 0
  DO ivv = 1, rismt%nsite
    !
    ! ... for G = 0
    csg0 = 0.0_DP
    hg0  = 0.0_DP
    !
    IF (rismt%mp_task%ivec_start == 1) THEN
      jr = 2
    ELSE
      jr = 1
    END IF
    !
    DO ir = jr, rismt%nr
      iir  = rismt%mp_task%ivec_start + ir - 1
      r    = rismt%rfft%rgrid(iir)
      rr   = r * r
      csg0 = csg0 + rfac * rr * rismt%csr(ir, ivv)
      hg0  = hg0  + rfac * rr * rismt%hr( ir, ivv)
    END DO
    !
    CALL mp_sum(csg0, rismt%mp_task%itask_comm)
    CALL mp_sum(hg0,  rismt%mp_task%itask_comm)
    !
    IF (rismt%mp_task%ivec_start == 1) THEN
      rismt%csg(1, ivv) = csg0
      rismt%hg( 1, ivv) = hg0
    END IF
    !
    ! ... for R = 0
    csr0 = 0.0_DP
    hr0  = 0.0_DP
    !
    IF (rismt%mp_task%ivec_start == 1) THEN
      jg = 2
    ELSE
      jg = 1
    END IF
    !
    DO ig = jg, rismt%ng
      iig = rismt%mp_task%ivec_start + ig - 1
      g    = rismt%rfft%ggrid(iig)
      gg   = g * g
      csr0 = csr0 + gfac * gg * rismt%csg(ig, ivv)
      hr0  = hr0  + gfac * gg * rismt%hg( ig, ivv)
    END DO
    !
    CALL mp_sum(csr0, rismt%mp_task%itask_comm)
    CALL mp_sum(hr0,  rismt%mp_task%itask_comm)
    !
    IF (rismt%mp_task%ivec_start == 1) THEN
      rismt%csr(1, ivv) = csr0
      rismt%hr( 1, ivv) = hr0
      !rismt%hr( 1, ivv) = -1.0_DP
    END IF
    !
  END DO
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE correctat0_vv
