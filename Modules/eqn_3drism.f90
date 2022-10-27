!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_3drism(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve 3D-RISM equation, which is defined as
  ! ...   h1(g) = c2(g) * x21(g)
  ! ...
  ! ... (A.Kovalenko, F.Hirata, Chem. Phys. Lett. 1998, 290, 237-244)
  !
  USE constants, ONLY : K_BOLTZMANN_RY
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_3DRISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iq1, iq2
  INTEGER                  :: iiq1, iiq2
  INTEGER                  :: ig
  INTEGER                  :: igs
  REAL(DP)                 :: beta
  REAL(DP)                 :: xg21
  COMPLEX(DP)              :: cgz2
  COMPLEX(DP), ALLOCATABLE :: hgz1(:)
#if defined (__DEBUG_RISM)
  !
  CALL start_clock('3DRISM_eqn')
#endif
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngs < rismt%gvec%ngl) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ng < rismt%gvec%ngm) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... allocate working memory
  IF (rismt%gvec%ngm > 0) THEN
    ALLOCATE(hgz1(rismt%gvec%ngm))
  END IF
  !
  ! ... 3D-RISM equation
  DO iq1 = 1, nq
    ! ... properties of site1
    IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
    ELSE
      iiq1 = 0
    END IF
    !
    IF (rismt%gvec%ngm > 0) THEN
      hgz1 = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    END IF
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      ! ... properties of site2
      iiq2 = iq2 - rismt%mp_site%isite_start + 1
      !
      ! ... solve 3D-RISM equation for each ig
!$omp parallel do default(shared) private(ig, igs, xg21, cgz2)
      DO ig = rismt%gvec%gstart, rismt%gvec%ngm
        igs      = rismt%gvec%igtongl(ig)
        xg21     = rismt%xgs(igs, iiq2, iq1)
        cgz2     = rismt%csgz(ig, iiq2) - beta * rismt%ulgz(ig, iiq2)
        hgz1(ig) = hgz1(ig) + cgz2 * xg21
      END DO
!$omp end parallel do
      !
    END DO
    !
    IF (rismt%gvec%ngm > 0) THEN
      CALL mp_sum(hgz1, rismt%mp_site%inter_sitg_comm)
    END IF
    !
    ! ... canonical condition
    IF (rismt%gvec%gstart > 1) THEN
      IF (rismt%gvec%ngm > 0) THEN
        hgz1(1) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
    END IF
    !
    IF (iiq1 > 0) THEN
      IF (rismt%gvec%ngm > 0) THEN
        rismt%hgz(1:rismt%gvec%ngm, iiq1) = hgz1
      END IF
    END IF
    !
  END DO
  !
  ! ... deallocate working memory
  IF (rismt%gvec%ngm > 0) THEN
    DEALLOCATE(hgz1)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
#if defined (__DEBUG_RISM)
  !
  CALL stop_clock('3DRISM_eqn')
#endif
  !
END SUBROUTINE eqn_3drism
