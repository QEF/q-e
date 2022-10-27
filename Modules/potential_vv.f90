!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE potential_vv(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... create potentials for 1D-RISM
  ! ...
  ! ... short-range R-space potential:
  ! ...   Lennard-Jones:
  ! ...               [ / sig \12     / sig \6 ]
  ! ...     4 * esp * [ |-----|   -   |-----|  ]
  ! ...               [ \  r  /       \  r  /  ]
  ! ...   Coulomb:
  ! ...                      erfc(r / tau)
  ! ...     e^2 * q1 * q2 * ---------------
  ! ...                            r
  ! ... long-range R-space potential:
  ! ...   Coulomb:
  ! ...                      erf(r / tau)
  ! ...     e^2 * q1 * q2 * --------------
  ! ...                            r
  ! ... long-range G-space potential:
  ! ...   Coulomb:
  ! ...                                exp(-g^2 * tau^2 / 4)
  ! ...     e^2 * 4 * pi * q1 * q2 * ------------------------
  ! ...                                         g^2
  !
  USE constants, ONLY : fpi, sqrtpi, ele2 => e2
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_1DRISM
  USE solvmol,   ONLY : solVs, get_nsite_in_solVs, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: nv
  INTEGER  :: iv1, iv2
  INTEGER  :: ivv
  INTEGER  :: isolV1, isolV2
  INTEGER  :: iatom1, iatom2
  INTEGER  :: ir, iir, jr
  INTEGER  :: ig, iig, jg
  REAL(DP) :: q1, q2, q12
  REAL(DP) :: e1, e2, e12
  REAL(DP) :: s1, s2, s12
  REAL(DP) :: r
  REAL(DP) :: tau
  REAL(DP) :: g
  REAL(DP) :: sr, sr2, sr6, sr12
  REAL(DP) :: ulj_r
  REAL(DP) :: utc_r
  REAL(DP) :: usc_r
  REAL(DP) :: ulc_r
  REAL(DP) :: ulc_g
  REAL(DP) :: erf0
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
  ! ... coulomb smearing radius
  tau = rismt%tau
  !
  ! ... calculate potentials
  DO iv1 = 1, nv
    ! ... properties of site1
    isolV1 = isite_to_isolV(iv1)
    iatom1 = isite_to_iatom(iv1)
    q1 = solVs(isolV1)%charge(iatom1)
    e1 = solVs(isolV1)%ljeps(iatom1)
    s1 = solVs(isolV1)%ljsig(iatom1)
    !
    DO iv2 = 1, iv1
      ! ... properties of site2
      isolV2 = isite_to_isolV(iv2)
      iatom2 = isite_to_iatom(iv2)
      q2 = solVs(isolV2)%charge(iatom2)
      e2 = solVs(isolV2)%ljeps(iatom2)
      s2 = solVs(isolV2)%ljsig(iatom2)
      !
      ! ... properties of pair of site1 and 2
      ivv = iv1 * (iv1 - 1) / 2 + iv2
      q12  = q1 * q2
      e12  = SQRT(e1 * e2)
      s12  = 0.5_DP * (s1 + s2)
      !
      ! ... short-range potential in R-space
      ! ...   in case R = 0
      IF (rismt%mp_task%ivec_start == 1) THEN
        jr = 2
        rismt%usr(1, ivv) = 0.0_DP
        rismt%ulr(1, ivv) = ele2 * q12 * 2.0_DP / sqrtpi / tau
      ELSE
        jr = 1
      END IF
      ! ...   in case R /= 0
!$omp parallel do default(shared) private(ir, iir, r, sr, sr2, sr6, sr12, &
!$omp          ulj_r, utc_r, usc_r, ulc_r, erf0)
      DO ir = jr, rismt%nr
        iir = rismt%mp_task%ivec_start + ir - 1
        r = rismt%rfft%rgrid(iir)
        sr = s12 / r
        sr2 = sr * sr
        sr6 = sr2 * sr2 * sr2
        sr12 = sr6 * sr6
        ulj_r = 4.0_DP * e12 * (sr12 - sr6)
        utc_r = ele2 * q12 / r
        erf0  = erf(r / tau)
        usc_r = utc_r * (1.0_DP - erf0)
        ulc_r = utc_r * erf0
        rismt%usr(ir, ivv) = ulj_r + usc_r
        rismt%ulr(ir, ivv) = ulc_r
      END DO
!$omp end parallel do
      !
      ! ... long-range potential in G-space
      ! ...   in case G = 0
      IF (rismt%mp_task%ivec_start == 1) THEN
        jg = 2
        rismt%ulg(1, ivv) = 0.0_DP
      ELSE
        jg = 1
      END IF
      ! ...   in case G /= 0
!$omp parallel do default(shared) private(ig, iig, g, ulc_g)
      DO ig = jg, rismt%ng
        iig = rismt%mp_task%ivec_start + ig - 1
        g = rismt%rfft%ggrid(iig)
        ulc_g = ele2 * fpi * q12 * EXP(-0.25_DP * g * g * tau * tau) / g / g
        rismt%ulg(ig, ivv) = ulc_g
      END DO
!$omp end parallel do
      !
    END DO
  END DO
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE potential_vv
