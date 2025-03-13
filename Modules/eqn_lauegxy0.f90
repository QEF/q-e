!
! Copyright (C) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_lauegxy0(rismt, lboth, expand, long, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve short-range part of Laue-RISM equation at Gxy = 0, which is defined as
  ! ...
  ! ...               /+inf
  ! ...   hs1(0,z1) = | dz2 cs2(0,z2) * x21(0,z2-z1)
  ! ...               /-inf
  ! ...
  ! ... total correlations are calculated around the unit cell or the expanded cell.
  ! ... optionally, long-range part can be added.
  ! ...
  !
  USE cell_base, ONLY : alat
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lboth  ! both-hands calculation, or not
  LOGICAL,         INTENT(IN)    :: expand ! expand-cell(.TRUE.) or unit-cell(.FALSE.)
  LOGICAL,         INTENT(IN)    :: long   ! add long-range part, or not
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: nq
  INTEGER               :: iq1, iq2
  INTEGER               :: iiq1, iiq2
  INTEGER               :: iz1, iz2
  INTEGER               :: iiz2
  INTEGER               :: nzint1
  INTEGER               :: izint1
  INTEGER               :: nzint2
  INTEGER               :: izint2
  INTEGER               :: izdelt
  INTEGER               :: nzright1
  INTEGER               :: izright1_sta
  INTEGER               :: izright1_end
  INTEGER               :: nzright2
  INTEGER               :: izright2_sta
  INTEGER               :: izright2_end
  INTEGER               :: nzleft1
  INTEGER               :: izleft1_sta
  INTEGER               :: izleft1_end
  INTEGER               :: nzleft2
  INTEGER               :: izleft2_sta
  INTEGER               :: izleft2_end
  REAL(DP), ALLOCATABLE :: xgt(:)
  REAL(DP), ALLOCATABLE :: ygt(:)
  REAL(DP)              :: zstep
  REAL(DP), ALLOCATABLE :: x21(:,:)
  REAL(DP), ALLOCATABLE :: cs2(:)
  REAL(DP), ALLOCATABLE :: hs1(:)
  !
  EXTERNAL :: dgemv
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
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set dz (in a.u.)
  zstep = alat * rismt%lfft%zstep
  !
  ! ... set integral regions as index of long Z-stick (i.e. expanded cell)
  IF (expand) THEN
    izright1_sta = rismt%lfft%izright_gedge
    izright1_end = rismt%lfft%nrz
    izleft1_sta  = 1
    izleft1_end  = rismt%lfft%izleft_gedge
  ELSE
    izright1_sta = rismt%lfft%izright_start0
    izright1_end = rismt%lfft%izright_end0
    izleft1_sta  = rismt%lfft%izleft_start0
    izleft1_end  = rismt%lfft%izleft_end0
  END IF
  !
  izright2_sta = rismt%lfft%izright_start0
  izright2_end = rismt%lfft%izright_end0
  izleft2_sta  = rismt%lfft%izleft_start0
  izleft2_end  = rismt%lfft%izleft_end0
  !
  ! ... count integral points along Z
  nzright1 = MAX(izright1_end - izright1_sta + 1, 0)
  nzleft1  = MAX(izleft1_end  - izleft1_sta  + 1, 0)
  nzint1   = nzright1 + nzleft1
  nzright2 = MAX(izright2_end - izright2_sta + 1, 0)
  nzleft2  = MAX(izleft2_end  - izleft2_sta  + 1, 0)
  nzint2   = nzright2 + nzleft2
  !
  ! ... allocate working memory
  IF (rismt%nrzl > 0) THEN
    ALLOCATE(xgt(rismt%nrzl))
    ALLOCATE(ygt(rismt%nrzl))
  END IF
  IF (nzint2 * nzint1 > 0) THEN
    ALLOCATE(x21(nzint2, nzint1))
  END IF
  IF (nzint2 > 0) THEN
    ALLOCATE(cs2(nzint2))
  END IF
  IF (nzint1 > 0) THEN
    ALLOCATE(hs1(nzint1))
  END IF
  !
  ! ... initialize hg0
  IF (.NOT. expand) THEN
    IF (rismt%nrzl * rismt%nsite > 0) THEN
      rismt%hg0 = 0.0_DP
    END IF
  END IF
  !
  ! ... Laue-RISM equation of short-range
  DO iq1 = 1, nq
    ! ... properties of site1
    IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
    ELSE
      iiq1 = 0
    END IF
    !
    IF (nzint1 > 0) THEN
      hs1 = 0.0_DP
    END IF
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      ! ... properties of site2
      iiq2 = iq2 - rismt%mp_site%isite_start + 1
      !
      ! ... solve Laue-RISM equation for Gxy = 0
      IF (rismt%lfft%gxystart > 1) THEN
        !
        ! ... x(z2-z1)
        xgt(1:rismt%nrzl) = rismt%xgs(1:rismt%nrzl, iiq2, iq1)
        IF (.NOT. lboth) THEN
          ygt(1:rismt%nrzl) = rismt%xgs(1:rismt%nrzl, iiq2, iq1)
        ELSE
          ygt(1:rismt%nrzl) = rismt%ygs(1:rismt%nrzl, iiq2, iq1)
        END IF
        !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
        DO iz1 = izleft1_sta, izleft1_end
          izint1 = iz1 - izleft1_sta + 1
          DO iz2 = izleft2_sta, izleft2_end
            izint2 = iz2 - izleft2_sta + 1
            izdelt = ABS(iz1 - iz2) + 1
            x21(izint2, izint1) = ygt(izdelt)
          END DO
          DO iz2 = izright2_sta, izright2_end
            izint2 = nzleft2 + iz2 - izright2_sta + 1
            izdelt = ABS(iz1 - iz2) + 1
            x21(izint2, izint1) = ygt(izdelt)
          END DO
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
        DO iz1 = izright1_sta, izright1_end
          izint1 = nzleft1 + iz1 - izright1_sta + 1
          DO iz2 = izleft2_sta, izleft2_end
            izint2 = iz2 - izleft2_sta + 1
            izdelt = ABS(iz1 - iz2) + 1
            x21(izint2, izint1) = xgt(izdelt)
          END DO
          DO iz2 = izright2_sta, izright2_end
            izint2 = nzleft2 + iz2 - izright2_sta + 1
            izdelt = ABS(iz1 - iz2) + 1
            x21(izint2, izint1) = xgt(izdelt)
          END DO
        END DO
!$omp end parallel do
        !
        ! ... cs(z2)
!$omp parallel do default(shared) private(iz2, izint2)
        DO iz2 = izleft2_sta, izleft2_end
          izint2 = iz2 - izleft2_sta + 1
          cs2(izint2) = rismt%csg0(iz2, iiq2)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(iz2, izint2)
        DO iz2 = izright2_sta, izright2_end
          izint2 = nzleft2 + iz2 - izright2_sta + 1
          cs2(izint2) = rismt%csg0(iz2, iiq2)
        END DO
!$omp end parallel do
        !
        ! ... hs(z1)
        IF (nzint2 * nzint1 > 0) THEN
          CALL dgemv('T', nzint2, nzint1, zstep, x21, nzint2, cs2, 1, 1.0_DP, hs1(1), 1)
        END IF
        !
      END IF
      !
    END DO
    !
    IF (nzint1 > 0) THEN
      CALL mp_sum(hs1, rismt%mp_site%inter_sitg_comm)
    END IF
    !
    IF (iiq1 > 0 .AND. rismt%lfft%gxystart > 1) THEN
      !
      IF (expand) THEN
        ! ... copy hs1 -> hsgz (expand-cell)
        rismt%hsgz(1:rismt%lfft%nrz, iiq1) = CMPLX(-1.0_DP, 0.0_DP, kind=DP)
        !
        IF (long) THEN
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izleft1_sta, izleft1_end
            izint1 = iz1 - izleft1_sta + 1
            rismt%hsgz(iz1, iiq1) = rismt%hlgz(iz1, iiq1) + CMPLX(hs1(izint1), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izright1_sta, izright1_end
            izint1 = nzleft1 + iz1 - izright1_sta + 1
            rismt%hsgz(iz1, iiq1) = rismt%hlgz(iz1, iiq1) + CMPLX(hs1(izint1), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
        ELSE
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izleft1_sta, izleft1_end
            izint1 = iz1 - izleft1_sta + 1
            rismt%hsgz(iz1, iiq1) = CMPLX(hs1(izint1), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izright1_sta, izright1_end
            izint1 = nzleft1 + iz1 - izright1_sta + 1
            rismt%hsgz(iz1, iiq1) = CMPLX(hs1(izint1), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
        END IF
        !
      ELSE
        ! ... copy hs1 -> hg0 (unit-cell)
        IF (rismt%nrzl > 0) THEN
          rismt%hg0(:, iiq1) = -1.0_DP
        END IF
        !
        DO iz1 = 1, (izleft1_sta - 1)
          rismt%hg0(iz1, iiq1) = 0.0_DP
        END DO
        !
        DO iz1 = (izright1_end + 1), rismt%lfft%nrz
          rismt%hg0(iz1, iiq1) = 0.0_DP
        END DO
        !
        IF (long) THEN
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izleft1_sta, izleft1_end
            izint1 = iz1 - izleft1_sta + 1
            rismt%hg0(iz1, iiq1) = DBLE(rismt%hlgz(iz1, iiq1)) + hs1(izint1)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izright1_sta, izright1_end
            izint1 = nzleft1 + iz1 - izright1_sta + 1
            rismt%hg0(iz1, iiq1) = DBLE(rismt%hlgz(iz1, iiq1)) + hs1(izint1)
          END DO
!$omp end parallel do
          !
        ELSE
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izleft1_sta, izleft1_end
            izint1 = iz1 - izleft1_sta + 1
            rismt%hg0(iz1, iiq1) = hs1(izint1)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izright1_sta, izright1_end
            izint1 = nzleft1 + iz1 - izright1_sta + 1
            rismt%hg0(iz1, iiq1) = hs1(izint1)
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
  ! ... share hg0, for R-space calculation
  IF (.NOT. expand) THEN
    IF (rismt%nrzl * rismt%nsite > 0) THEN
      CALL mp_sum(rismt%hg0, rismt%mp_site%intra_sitg_comm)
    END IF
  END IF
  !
  ! ... deallocate working memory
  IF (rismt%nrzl > 0) THEN
    DEALLOCATE(xgt)
    DEALLOCATE(ygt)
  END IF
  IF (nzint2 * nzint1 > 0) THEN
    DEALLOCATE(x21)
  END IF
  IF (nzint2 > 0) THEN
    DEALLOCATE(cs2)
  END IF
  IF (nzint1 > 0) THEN
    DEALLOCATE(hs1)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE eqn_lauegxy0
