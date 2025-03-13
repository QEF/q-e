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
SUBROUTINE eqn_laueshort(rismt, lboth, lgxy0, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve short-range part of Laue-RISM equation, which is defined as
  ! ...
  ! ...                 /+inf
  ! ...   hs1(gxy,z1) = | dz2 cs2(gxy,z2) * x21(gxy,z2-z1)
  ! ...                 /-inf
  ! ...
  ! ... total correlations are calculated around the expanded cell.
  ! ... also, dipole-part and contribution from void-region are included.
  ! ... this subroutine will be performed when RISM's equation will be converged.
  ! ...
  !
  USE err_rism, ONLY : IERR_RISM_NULL
  USE kinds,    ONLY : DP
  USE rism,     ONLY : rism_type
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lboth  ! both-hands calculation, or not
  LOGICAL,         INTENT(IN)    :: lgxy0  ! only for Gxy=0, or not
  INTEGER,         INTENT(OUT)   :: ierr
  !
  ! ... Laue-RISM equation of short-range (Gxy /= 0)
  IF (.NOT. lgxy0) THEN
#if defined (__DEBUG_RISM)
    CALL start_clock('3DRISM_eqnSx')
    !
#endif
    CALL eqn_laueshort_x(rismt, lboth, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('3DRISM_eqnSx')
#endif
  END IF
  !
  ! ... Laue-RISM equation of short-range (Gxy = 0)
#if defined (__DEBUG_RISM)
  CALL start_clock('3DRISM_eqnS0')
  !
#endif
  CALL eqn_lauegxy0(rismt, lboth, .TRUE., .FALSE., ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... add dipole part of Laue-RISM (Gxy = 0)
  CALL eqn_lauedipole(rismt, .TRUE., .FALSE., ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... add contribution from void-region (Gxy = 0)
  CALL eqn_lauevoid(rismt, lboth, .TRUE., ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
#if defined (__DEBUG_RISM)
  !
  CALL stop_clock('3DRISM_eqnS0')
#endif
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE eqn_laueshort
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_laueshort_x(rismt, lboth, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve short-range part of Laue-RISM equation, for Gxy /= 0.
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
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iq1, iq2
  INTEGER                  :: iiq1, iiq2
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  INTEGER                  :: iglxy
  INTEGER                  :: iglxy_old
  INTEGER                  :: jglxy
  INTEGER                  :: iz1, iz2
  INTEGER                  :: iiz2
  INTEGER                  :: nzint1
  INTEGER                  :: izint1
  INTEGER                  :: nzint2
  INTEGER                  :: izint2
  INTEGER                  :: izdelt
  INTEGER                  :: nzright1
  INTEGER                  :: izright1_sta
  INTEGER                  :: izright1_end
  INTEGER                  :: nzright2
  INTEGER                  :: izright2_sta
  INTEGER                  :: izright2_end
  INTEGER                  :: nzleft1
  INTEGER                  :: izleft1_sta
  INTEGER                  :: izleft1_end
  INTEGER                  :: nzleft2
  INTEGER                  :: izleft2_sta
  INTEGER                  :: izleft2_end
  REAL(DP),    ALLOCATABLE :: xgt(:)
  REAL(DP),    ALLOCATABLE :: ygt(:)
  COMPLEX(DP)              :: zstep
  COMPLEX(DP), ALLOCATABLE :: x21(:,:)
  COMPLEX(DP), ALLOCATABLE :: cs2(:)
  COMPLEX(DP), ALLOCATABLE :: hs1(:,:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX( 0.0_DP, 0.0_DP, kind=DP)
  COMPLEX(DP), PARAMETER   :: C_ONE  = CMPLX( 1.0_DP, 0.0_DP, kind=DP)
  !
  EXTERNAL :: zgemv
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
  IF (rismt%ngs < rismt%lfft%nglxy) THEN
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
  ! ... set dz (in a.u.)
  zstep = CMPLX(alat * rismt%lfft%zstep, 0.0_DP, kind=DP)
  !
  ! ... set integral regions as index of long Z-stick (i.e. expanded cell)
  izright1_sta = rismt%lfft%izright_start
  izright1_end = rismt%lfft%nrz
  izleft1_sta  = 1
  izleft1_end  = rismt%lfft%izleft_end
  izright2_sta = rismt%lfft%izright_start
  izright2_end = rismt%lfft%izright_end
  izleft2_sta  = rismt%lfft%izleft_start
  izleft2_end  = rismt%lfft%izleft_end
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
  IF (nzint1 * rismt%lfft%ngxy > 0) THEN
    ALLOCATE(hs1(nzint1, rismt%lfft%ngxy))
  END IF
  !
  DO iq1 = 1, nq
    ! ... properties of site1
    IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
    ELSE
      iiq1 = 0
    END IF
    !
    IF (nzint1 * rismt%lfft%ngxy > 0) THEN
      hs1 = C_ZERO
    END IF
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      ! ... properties of site2
      iiq2 = iq2 - rismt%mp_site%isite_start + 1
      !
      ! ... solve Laue-RISM equation for each igxy
      iglxy_old = -1
      !
      ! ... loop for Gxy /= 0
      DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
        jgxy  = rismt%nrzs * (igxy - 1)
        iglxy = rismt%lfft%igtonglxy(igxy)
        jglxy = rismt%nrzl * (iglxy - 1)
        !
        ! ... x(z2-z1)
        IF (iglxy /= iglxy_old) THEN
          iglxy_old = iglxy
          !
          xgt(1:rismt%nrzl) = rismt%xgs((1 + jglxy):(rismt%nrzl + jglxy), iiq2, iq1)
          IF (.NOT. lboth) THEN
            ygt(1:rismt%nrzl) = rismt%xgs((1 + jglxy):(rismt%nrzl + jglxy), iiq2, iq1)
          ELSE
            ygt(1:rismt%nrzl) = rismt%ygs((1 + jglxy):(rismt%nrzl + jglxy), iiq2, iq1)
          END IF
          !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
          DO iz1 = izleft1_sta, izleft1_end
            izint1 = iz1 - izleft1_sta + 1
            DO iz2 = izleft2_sta, izleft2_end
              izint2 = iz2 - izleft2_sta + 1
              izdelt = ABS(iz1 - iz2) + 1
              x21(izint2, izint1) = CMPLX(ygt(izdelt), 0.0_DP, kind=DP)
            END DO
            DO iz2 = izright2_sta, izright2_end
              izint2 = nzleft2 + iz2 - izright2_sta + 1
              izdelt = ABS(iz1 - iz2) + 1
              x21(izint2, izint1) = CMPLX(ygt(izdelt), 0.0_DP, kind=DP)
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
              x21(izint2, izint1) = CMPLX(xgt(izdelt), 0.0_DP, kind=DP)
            END DO
            DO iz2 = izright2_sta, izright2_end
              izint2 = nzleft2 + iz2 - izright2_sta + 1
              izdelt = ABS(iz1 - iz2) + 1
              x21(izint2, izint1) = CMPLX(xgt(izdelt), 0.0_DP, kind=DP)
            END DO
          END DO
!$omp end parallel do
          !
        END IF
        !
        ! ... cs(z2)
!$omp parallel do default(shared) private(iz2, izint2, iiz2)
        DO iz2 = izleft2_sta, izleft2_end
          izint2 = iz2 - izleft2_sta + 1
          iiz2 = iz2 - rismt%lfft%izcell_start + 1
          cs2(izint2) = rismt%csgz(iiz2 + jgxy, iiq2)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(iz2, izint2, iiz2)
        DO iz2 = izright2_sta, izright2_end
          izint2 = nzleft2 + iz2 - izright2_sta + 1
          iiz2 = iz2 - rismt%lfft%izcell_start + 1
          cs2(izint2) = rismt%csgz(iiz2 + jgxy, iiq2)
        END DO
!$omp end parallel do
        !
        ! ... hs(z1)
        IF (nzint2 * nzint1 > 0) THEN
          CALL zgemv('T', nzint2, nzint1, zstep, x21, nzint2, cs2, 1, C_ONE, hs1(1, igxy), 1)
        END IF
        !
      END DO
      !
    END DO
    !
    IF (nzint1 * rismt%lfft%ngxy > 0) THEN
      CALL mp_sum(hs1, rismt%mp_site%inter_sitg_comm)
    END IF
    !
    IF (iiq1 > 0) THEN
      ! ... copy hs1 -> hsgz
      IF (rismt%nrzl * rismt%ngxy > 0) THEN
        rismt%hsgz(:, iiq1) = C_ZERO
      END IF
      !
      ! ... loop for Gxy /= 0
      DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
        jgxy = rismt%nrzl * (igxy - 1)
        !
!$omp parallel do default(shared) private(iz1, izint1)
        DO iz1 = izleft1_sta, izleft1_end
          izint1 = iz1 - izleft1_sta + 1
          rismt%hsgz(iz1 + jgxy, iiq1) = hs1(izint1, igxy)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(iz1, izint1)
        DO iz1 = izright1_sta, izright1_end
          izint1 = nzleft1 + iz1 - izright1_sta + 1
          rismt%hsgz(iz1 + jgxy, iiq1) = hs1(izint1, igxy)
        END DO
!$omp end parallel do
      END DO
    END IF
    !
  END DO
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
  IF (nzint1 * rismt%lfft%ngxy > 0) THEN
    DEALLOCATE(hs1)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE eqn_laueshort_x
