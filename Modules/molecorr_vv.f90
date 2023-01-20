!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE molecorr_vv(rismt, alpha, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... create intra-molecular correlation for 1D-RISM
  ! ...
  ! ... if site1 and 2 are in a molecule
  ! ...
  ! ...   1) in case that alpha is zero or negative.
  ! ...
  ! ...             sin(g * r12)
  ! ...     w(g) = --------------
  ! ...               g * r12
  ! ...
  ! ...   2) in case that alpha is positive.
  ! ...
  ! ...                        1                        [ r-r12 ]^2          [ r+r12 ]^2
  ! ...     w(r) = -------------------------- * [ exp(- [-------]  ) - exp(- [-------]  ) ]
  ! ...             4*(pi)^(3/2)*alpha*r*r12            [ alpha ]            [ alpha ]
  ! ...
  ! ...             sin(g * r12)
  ! ...     w(g) = -------------- * exp(- g^2 * alpha^2 / 4)
  ! ...               g * r12
  ! ...
  ! ... if site1 is in another molecule from site2
  ! ...
  ! ...   w(g) = 0
  ! ...
  ! ... (F.Hirata et al., Chem. Phys. Lett. 1981, 83, 329-334)
  ! ...
  !
  USE constants, ONLY : eps8
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_1DRISM
  USE solvmol,   ONLY : solVs, get_nsite_in_solVs, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: alpha
  INTEGER,         INTENT(OUT)   :: ierr
  !
  REAL(DP), PARAMETER :: RMIN = eps8
  !
  INTEGER  :: nv
  INTEGER  :: iv1, iv2
  INTEGER  :: ivv
  INTEGER  :: isolV1, isolV2
  INTEGER  :: iatom1, iatom2
  INTEGER  :: ig, iig, jg
  REAL(DP) :: x1, x2
  REAL(DP) :: y1, y2
  REAL(DP) :: z1, z2
  REAL(DP) :: rr, r
  REAL(DP) :: g
  REAL(DP) :: exp0
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
  ! ... calculate correlation functions
  DO iv1 = 1, nv
    ! ... properties of site1
    isolV1 = isite_to_isolV(iv1)
    iatom1 = isite_to_iatom(iv1)
    x1 = solVs(isolV1)%coord(1, iatom1)
    y1 = solVs(isolV1)%coord(2, iatom1)
    z1 = solVs(isolV1)%coord(3, iatom1)
    !
    DO iv2 = 1, iv1
      ! ... properties of site2
      isolV2 = isite_to_isolV(iv2)
      iatom2 = isite_to_iatom(iv2)
      x2 = solVs(isolV2)%coord(1, iatom2)
      y2 = solVs(isolV2)%coord(2, iatom2)
      z2 = solVs(isolV2)%coord(3, iatom2)
      !
      ivv = iv1 * (iv1 - 1) / 2 + iv2
      !
      ! ... site1 and 2 are in a solvent molecule
      IF (isolV1 == isolV2) THEN
        ! ... in case G = 0
        IF (rismt%mp_task%ivec_start == 1) THEN
          jg = 2
          rismt%wg(1, ivv) = 1.0_DP
        ELSE
          jg = 1
        END IF
        !
        rr = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)
        ! ... in case G /= 0, but R = 0
        IF (rr < (RMIN * RMIN)) THEN
          IF (alpha <= 0.0_DP) THEN
            DO ig = jg, rismt%ng
              rismt%wg(ig, ivv) = 1.0_DP
            END DO
            !
          ELSE !IF (alpha > 0.0_DP) THEN
            DO ig = jg, rismt%ng
              iig = rismt%mp_task%ivec_start + ig - 1
              g = rismt%rfft%ggrid(iig)
              exp0 = EXP(-0.25_DP * g * g * alpha * alpha)
              rismt%wg(ig, ivv) = exp0
            END DO
          END IF
          !
        ! ... in case G /= 0, and R /= 0
        ELSE
          r = SQRT(rr)
          !
          IF (alpha <= 0.0_DP) THEN
            DO ig = jg, rismt%ng
              iig = rismt%mp_task%ivec_start + ig - 1
              g = rismt%rfft%ggrid(iig)
              rismt%wg(ig, ivv) = SIN(g * r) / g / r
            END DO
            !
          ELSE !IF (alpha > 0.0_DP) THEN
            DO ig = jg, rismt%ng
              iig = rismt%mp_task%ivec_start + ig - 1
              g = rismt%rfft%ggrid(iig)
              exp0 = EXP(-0.25_DP * g * g * alpha * alpha)
              rismt%wg(ig, ivv) = SIN(g * r) / g / r * exp0
            END DO
          END IF
        END IF
        !
      ! ... site1 and 2 are in different solvent molecules
      ELSE
        rismt%wg(:, ivv) = 0.0_DP
      END IF
      !
    END DO
  END DO
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE molecorr_vv
