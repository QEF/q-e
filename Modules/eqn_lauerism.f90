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
SUBROUTINE eqn_lauerism(rismt, lboth, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve Laue-RISM equation, which is defined as
  ! ...
  ! ...                /+inf
  ! ...   h1(gxy,z1) = | dz2 cs2(gxy,z2) * x21(gxy,z2-z1) + hl1(gxy,z1)
  ! ...                /-inf
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
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER :: iq
  INTEGER :: iiq
  INTEGER :: irz
  INTEGER :: iirz
  !
  ! ... Laue-RISM equation (Gxy /= 0)
#if defined (__DEBUG_RISM)
  CALL start_clock('3DRISM_eqnx')
  !
#endif
  CALL eqn_lauerism_x(rismt, lboth, ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
#if defined (__DEBUG_RISM)
  !
  CALL stop_clock('3DRISM_eqnx')
#endif
  !
  ! ... Laue-RISM equation of short-range and long-range (Gxy = 0)
#if defined (__DEBUG_RISM)
  CALL start_clock('3DRISM_eqn0')
  !
#endif
  CALL eqn_lauegxy0(rismt, lboth, .FALSE., .TRUE., ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... add dipole part of Laue-RISM (Gxy = 0)
  CALL eqn_lauedipole(rismt, .FALSE., .FALSE., ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... add contribution from void-region (Gxy = 0)
  CALL eqn_lauevoid(rismt, lboth, .FALSE., ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... set hgz at Gxy = 0
  IF (rismt%lfft%gxystart > 1) THEN
    !
    DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq = iq - rismt%mp_site%isite_start + 1
      !
!$omp parallel do default(shared) private(irz, iirz)
      DO irz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
        iirz = irz - rismt%lfft%izcell_start + 1
        rismt%hgz(iirz, iiq) = CMPLX(rismt%hg0(irz, iiq), 0.0_DP, kind=DP)
      END DO
!$omp end parallel do
      !
    END DO
    !
  END IF
#if defined (__DEBUG_RISM)
  !
  CALL stop_clock('3DRISM_eqn0')
#endif
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE eqn_lauerism
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_lauerism_x(rismt, lboth, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve Laue-RISM equation, for Gxy /= 0.
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
  INTEGER                  :: kgxy
  INTEGER                  :: iglxy
  INTEGER                  :: iglxy_old
  INTEGER                  :: jglxy
  INTEGER                  :: iz1, iz2
  INTEGER                  :: iiz1
  INTEGER                  :: nzint
  INTEGER                  :: izint1
  INTEGER                  :: izint2
  INTEGER                  :: izdelt
  INTEGER                  :: nzright
  INTEGER                  :: izright_sta
  INTEGER                  :: izright_end
  INTEGER                  :: nzleft
  INTEGER                  :: izleft_sta
  INTEGER                  :: izleft_end
  REAL(DP),    ALLOCATABLE :: xgt(:)
  REAL(DP),    ALLOCATABLE :: ygt(:)
  COMPLEX(DP)              :: zstep
  COMPLEX(DP), ALLOCATABLE :: x21(:,:)
  COMPLEX(DP), ALLOCATABLE :: cs2(:)
  COMPLEX(DP), ALLOCATABLE :: hs1(:,:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX( 0.0_DP, 0.0_DP, kind=DP)
  COMPLEX(DP), PARAMETER   :: C_ONE  = CMPLX( 1.0_DP, 0.0_DP, kind=DP)
  COMPLEX(DP), PARAMETER   :: C_MONE = CMPLX(-1.0_DP, 0.0_DP, kind=DP)
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
  ! ... set integral regions as index of short Z-stick (i.e. unit cell)
  izright_sta = rismt%lfft%izright_start - rismt%lfft%izcell_start + 1
  izright_end = rismt%lfft%izright_end   - rismt%lfft%izcell_start + 1
  izleft_sta  = rismt%lfft%izleft_start  - rismt%lfft%izcell_start + 1
  izleft_end  = rismt%lfft%izleft_end    - rismt%lfft%izcell_start + 1
  !
  ! ... count integral points along Z
  nzright = MAX(izright_end - izright_sta + 1, 0)
  nzleft  = MAX(izleft_end  - izleft_sta  + 1, 0)
  nzint   = nzright + nzleft
  !
  ! ... allocate working memory
  IF (rismt%nrzl > 0) THEN
    ALLOCATE(xgt(rismt%nrzl))
  END IF
  IF (lboth .AND. rismt%nrzl > 0) THEN
    ALLOCATE(ygt(rismt%nrzl))
  END IF
  IF (nzint > 0) THEN
    ALLOCATE(x21(nzint, nzint))
    ALLOCATE(cs2(nzint))
  END IF
  IF (nzint * rismt%lfft%ngxy > 0) THEN
    ALLOCATE(hs1(nzint, rismt%lfft%ngxy))
  END IF
  !
  ! ... Laue-RISM equation of short-range (Gxy /= 0)
  DO iq1 = 1, nq
    ! ... properties of site1
    IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
    ELSE
      iiq1 = 0
    END IF
    !
    IF (nzint * rismt%lfft%ngxy > 0) THEN
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
          IF (.NOT. lboth) THEN
            ! ... susceptibility matrix of one-hand calculation, which is symmetric
            xgt(1:rismt%nrzl) = rismt%xgs((1 + jglxy):(rismt%nrzl + jglxy), iiq2, iq1)
            !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
            DO iz1 = izleft_sta, izleft_end
              izint1 = iz1 - izleft_sta + 1
              DO iz2 = izleft_sta, iz1
                izint2 = iz2 - izleft_sta + 1
                izdelt = iz1 - iz2 + 1
                x21(izint2, izint1) = CMPLX(xgt(izdelt), 0.0_DP, kind=DP)
              END DO
            END DO
!$omp end parallel do
            !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
            DO iz1 = izright_sta, izright_end
              izint1 = nzleft + iz1 - izright_sta + 1
              DO iz2 = izright_sta, iz1
                izint2 = nzleft + iz2 - izright_sta + 1
                izdelt = iz1 - iz2 + 1
                x21(izint2, izint1) = CMPLX(xgt(izdelt), 0.0_DP, kind=DP)
              END DO
            END DO
!$omp end parallel do
            !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
            DO iz1 = izright_sta, izright_end
              izint1 = nzleft + iz1 - izright_sta + 1
              DO iz2 = izleft_sta, izleft_end
                izint2 = iz2 - izleft_sta + 1
                izdelt = iz1 - iz2 + 1
                x21(izint2, izint1) = CMPLX(xgt(izdelt), 0.0_DP, kind=DP)
              END DO
            END DO
!$omp end parallel do
            !
!$omp parallel do default(shared) private(izint1, izint2)
            DO izint1 = 1, nzint
              DO izint2 = 1, (izint1 - 1)
                x21(izint1, izint2) = x21(izint2, izint1)
              END DO
            END DO
!$omp end parallel do
            !
          ELSE !IF (lboth) THEN
            ! ... susceptibility matrix of both-hands calculation
            xgt(1:rismt%nrzl) = rismt%xgs((1 + jglxy):(rismt%nrzl + jglxy), iiq2, iq1)
            ygt(1:rismt%nrzl) = rismt%ygs((1 + jglxy):(rismt%nrzl + jglxy), iiq2, iq1)
            !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
            DO iz1 = izleft_sta, izleft_end
              izint1 = iz1 - izleft_sta + 1
              DO iz2 = izleft_sta, izleft_end
                izint2 = iz2 - izleft_sta + 1
                izdelt = ABS(iz1 - iz2) + 1
                x21(izint2, izint1) = CMPLX(ygt(izdelt), 0.0_DP, kind=DP)
              END DO
            END DO
!$omp end parallel do
            !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
            DO iz1 = izright_sta, izright_end
              izint1 = nzleft + iz1 - izright_sta + 1
              DO iz2 = izright_sta, izright_end
                izint2 = nzleft + iz2 - izright_sta + 1
                izdelt = ABS(iz1 - iz2) + 1
                x21(izint2, izint1) = CMPLX(xgt(izdelt), 0.0_DP, kind=DP)
              END DO
            END DO
!$omp end parallel do
            !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
            DO iz1 = izleft_sta, izleft_end
              izint1 = iz1 - izleft_sta + 1
              DO iz2 = izright_sta, izright_end
                izint2 = nzleft + iz2 - izright_sta + 1
                izdelt = iz2 - iz1 + 1
                x21(izint2, izint1) = CMPLX(ygt(izdelt), 0.0_DP, kind=DP)
              END DO
            END DO
!$omp end parallel do
            !
!$omp parallel do default(shared) private(iz1, izint1, iz2, izint2, izdelt)
            DO iz1 = izright_sta, izright_end
              izint1 = nzleft + iz1 - izright_sta + 1
              DO iz2 = izleft_sta, izleft_end
                izint2 = iz2 - izleft_sta + 1
                izdelt = iz1 - iz2 + 1
                x21(izint2, izint1) = CMPLX(xgt(izdelt), 0.0_DP, kind=DP)
              END DO
            END DO
!$omp end parallel do
            !
          END IF
          !
        END IF
        !
        ! ... cs(z2)
!$omp parallel do default(shared) private(iz2, izint2)
        DO iz2 = izleft_sta, izleft_end
          izint2 = iz2 - izleft_sta + 1
          cs2(izint2) = rismt%csgz(iz2 + jgxy, iiq2)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(iz2, izint2)
        DO iz2 = izright_sta, izright_end
          izint2 = nzleft + iz2 - izright_sta + 1
          cs2(izint2) = rismt%csgz(iz2 + jgxy, iiq2)
        END DO
!$omp end parallel do
        !
        ! ... hs(z1)
        IF (nzint > 0) THEN
          CALL zgemv('T', nzint, nzint, zstep, x21, nzint, cs2, 1, C_ONE, hs1(1, igxy), 1)
        END IF
        !
      END DO
      !
    END DO
    !
    IF (nzint * rismt%lfft%ngxy > 0) THEN
      CALL mp_sum(hs1, rismt%mp_site%inter_sitg_comm)
    END IF
    !
    IF (iiq1 > 0) THEN
      IF (rismt%nrzs * rismt%ngxy > 0) THEN
        rismt%hgz(:, iiq1) = C_ZERO
      END IF
      !
      IF (rismt%lfft%gxystart > 1) THEN
        rismt%hgz(1:rismt%dfft%nr3, iiq1) = C_MONE
      END IF
      !
      ! ... loop for Gxy /= 0
      DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
        jgxy  = rismt%nrzs * (igxy - 1)
        !
!$omp parallel do default(shared) private(iz1, izint1)
        DO iz1 = izleft_sta, izleft_end
          izint1 = iz1 - izleft_sta + 1
          rismt%hgz(iz1 + jgxy, iiq1) = hs1(izint1, igxy)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(iz1, izint1)
        DO iz1 = izright_sta, izright_end
          izint1 = nzleft + iz1 - izright_sta + 1
          rismt%hgz(iz1 + jgxy, iiq1) = hs1(izint1, igxy)
        END DO
!$omp end parallel do
      END DO
    END IF
    !
  END DO
  !
  ! ... add long-range correlation (Gxy /= 0)
  DO iq1 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq1 = iq1 - rismt%mp_site%isite_start + 1
    !
    ! ... loop for Gxy /= 0
    DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
      jgxy = rismt%nrzs * (igxy - 1)
      kgxy = rismt%nrzl * (igxy - 1)
      !
!$omp parallel do default(shared) private(iz1, iiz1)
      DO iz1 = izleft_sta, izleft_end
        iiz1 = iz1 + rismt%lfft%izcell_start - 1
        rismt%hgz(iz1 + jgxy, iiq1) = rismt%hgz(iz1 + jgxy, iiq1) + rismt%hlgz(iiz1 + kgxy, iiq1)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(iz1, iiz1)
      DO iz1 = izright_sta, izright_end
        iiz1 = iz1 + rismt%lfft%izcell_start - 1
        rismt%hgz(iz1 + jgxy, iiq1) = rismt%hgz(iz1 + jgxy, iiq1) + rismt%hlgz(iiz1 + kgxy, iiq1)
      END DO
!$omp end parallel do
    END DO
    !
  END DO
  !
  ! ... set zero to hgz at Gxy = 0
  IF (rismt%lfft%gxystart > 1) THEN
    !
    DO iq1 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq1 = iq1 - rismt%mp_site%isite_start + 1
      rismt%hgz(1:rismt%nrzs, iiq1) = C_ZERO
    END DO
    !
  END IF
  !
  ! ... deallocate working memory
  IF (rismt%nrzl > 0) THEN
    DEALLOCATE(xgt)
  END IF
  IF (lboth .AND. rismt%nrzl > 0) THEN
    DEALLOCATE(ygt)
  END IF
  IF (nzint > 0) THEN
    DEALLOCATE(x21)
    DEALLOCATE(cs2)
  END IF
  IF (nzint * rismt%lfft%ngxy > 0) THEN
    DEALLOCATE(hs1)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE eqn_lauerism_x
