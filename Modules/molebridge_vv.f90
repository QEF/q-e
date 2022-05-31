!
! Copyright (C) 2018 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE molebridge_vv(rismt, epsr, tau, lhand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... create intra-molecular dielectric bridge for DRISM
  ! ...
  ! ... if site1 and 2 are in a molecule
  ! ...
  ! ...   z(g) = j0(g*x1) j0(g*y1) * j1(g*z1) * hc(g) * j0(g*x2) * j0(g*y2) * j1(g*z2)
  ! ...
  ! ... , where
  ! ...
  ! ...              sin(g * r)
  ! ...   j0(g*r) = ------------ ,
  ! ...                g * r
  ! ...
  ! ...              sin(g * r)       cos(g * r)
  ! ...   j1(g*r) = ------------  -  ------------  and
  ! ...               (g * r)^2         g * r
  ! ...
  ! ...
  ! ...            eps/y - 3
  ! ...   hc(g) = ----------- * exp(- g^2 * tau^2 / 4)
  ! ...               rho
  ! ...
  ! ... if site1 is in another molecule from site2
  ! ...
  ! ...   z(g) = 0
  ! ...
  ! ... (J.S.Perkyns and B.M.Pettitt, CPL 1992, 190, 626)
  ! ...
  !
  USE constants, ONLY : K_BOLTZMANN_RY, fpi, e2, eps8, eps32
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_1DRISM
  USE solvmol,   ONLY : solVs, nsolV, get_nsite_in_solVs, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: epsr   ! dielectric constant
  REAL(DP),        INTENT(IN)    :: tau    ! size of molecule
  LOGICAL,         INTENT(IN)    :: lhand  ! if true, right-hand. if false, left-hand.
  INTEGER,         INTENT(OUT)   :: ierr
  !
  REAL(DP), PARAMETER :: RMIN = eps8
  !
  INTEGER               :: nv
  INTEGER               :: iv1, iv2
  INTEGER               :: ivv
  INTEGER               :: isolV1, isolV2
  INTEGER               :: iatom1, iatom2
  INTEGER               :: ig, iig, jg
  REAL(DP)              :: beta
  REAL(DP)              :: x1, x2
  REAL(DP)              :: y1, y2
  REAL(DP)              :: z1, z2
  REAL(DP)              :: jx1, jx2
  REAL(DP)              :: jy1, jy2
  REAL(DP)              :: jz1, jz2
  REAL(DP)              :: g
  REAL(DP)              :: den, dip
  REAL(DP)              :: y0, rho0
  REAL(DP)              :: a0, exp0
  REAL(DP), ALLOCATABLE :: hc(:)
  !
  ! ... set zero, if not DRISM
  IF (epsr <= 0.0_DP .OR. tau <= 0.0_DP .OR. (.NOT. ANY(solVs(:)%is_polar))) THEN
    !
    rismt%zg = 0.0_DP
    !
    ierr = IERR_RISM_NULL
    RETURN
  END IF
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
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... allocate memory
  ALLOCATE(hc(rismt%ng))
  !
  ! ... calculate envelope function (hc)
  rho0 = 0.0_DP
  y0   = 0.0_DP
  !
  DO isolV1 = 1, nsolV
    IF (solVs(isolV1)%is_polar) THEN
      IF (lhand) THEN
        den = solVs(isolV1)%density
      ELSE
        den = solVs(isolV1)%subdensity
      END IF
      dip  = solVs(isolV1)%dipole
      rho0 = rho0 + den
      y0   = y0   + den * dip * dip
    END IF
  END DO
  !
  IF (ABS(rho0) < eps32 .OR. ABS(y0) < eps32) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  y0 = fpi * beta / 9.0_DP * y0
  !
  a0 = ((epsr / fpi / e2) / y0 - 3.0_DP) / rho0
  !
  DO ig = 1, rismt%ng
    iig = rismt%mp_task%ivec_start + ig - 1
    g = rismt%rfft%ggrid(iig)
    !
    exp0 = EXP(-0.25_DP * g * g * tau * tau)
    hc(ig) = a0 * exp0
  END DO
  !
  ! ... calculate dielectric bridge functions
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
      ! ... site1 and 2 are in a polar solvent molecule
      IF (isolV1 == isolV2 .AND. solVs(isolV1)%is_polar) THEN
        ! ... in case G = 0
        IF (rismt%mp_task%ivec_start == 1) THEN
          jg = 2
          rismt%zg(1, ivv) = 0.0_DP
        ELSE
          jg = 1
        END IF
        !
        ! ... in case G /= 0
        DO ig = jg, rismt%ng
          iig = rismt%mp_task%ivec_start + ig - 1
          g = rismt%rfft%ggrid(iig)
          !
          jx1 = bessel0(g, x1)
          jx2 = bessel0(g, x2)
          jy1 = bessel0(g, y1)
          jy2 = bessel0(g, y2)
          jz1 = bessel1(g, z1)
          jz2 = bessel1(g, z2)
          !
          rismt%zg(ig, ivv) = jx1 * jx2 * jy1 * jy2 * jz1 * jz2 * hc(ig)
        END DO
        !
      ! ... site1 and 2 are in different solvent molecules,
      ! ... or not polar solvent molecule
      ELSE
        rismt%zg(:, ivv) = 0.0_DP
      END IF
      !
    END DO
  END DO
  !
  ! ... deallocate memory
  DEALLOCATE(hc)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  FUNCTION bessel0(g, r) RESULT(j0)
    !
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: g
    REAL(DP), INTENT(IN) :: r
    !
    REAL(DP) :: j0
    REAL(DP) :: gr
    !
    IF (ABS(r) < RMIN) THEN
      j0 = 1.0_DP
    ELSE
      gr = g * r
      j0 = SIN(gr) / gr
    END IF
    !
  END FUNCTION bessel0
  !
  FUNCTION bessel1(g, r) RESULT(j1)
    !
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: g
    REAL(DP), INTENT(IN) :: r
    !
    REAL(DP) :: j1
    REAL(DP) :: gr
    !
    IF (ABS(r) < RMIN) THEN
      j1 = 0.0_DP
    ELSE
      gr = g * r
      j1 = SIN(gr) / gr / gr - COS(gr) / gr
    END IF
    !
  END FUNCTION bessel1
  !
END SUBROUTINE molebridge_vv
