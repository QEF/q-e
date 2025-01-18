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
SUBROUTINE corrgxy0_laue(rismt, lextract, ar, ag0, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... sum or extract Gxy=0 terms of correlations, in Laue-RISM calculation.
  ! ...
  ! ... Variables:
  ! ...   lextract: if .TRUE.  extract Gxy=0 terms of correlations
  ! ...             if .FALSE. sum A(r) + A(gxy=0,z)
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE fft_types, ONLY : fft_index_to_3d
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)    :: rismt
  LOGICAL,         INTENT(IN)    :: lextract
  REAL(DP),        INTENT(INOUT) :: ar(rismt%nr, 1:*)
  REAL(DP),        INTENT(INOUT) :: ag0(rismt%nrzl, 1:*)
  INTEGER,         INTENT(OUT)   :: ierr
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (lextract) THEN
    !
    ! ... extract Gxy=0 terms
    CALL extract_gxy0(ar, ag0)
    !
  ELSE
    !
    ! ... sum A(r) + A(gxy=0,z)
    CALL sum_gxy0(ar, ag0)
    !
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE extract_gxy0(ar, ag0)
    IMPLICIT NONE
    REAL(DP), INTENT(IN)    :: ar(rismt%nr, 1:*)
    REAL(DP), INTENT(INOUT) :: ag0(rismt%nrzl, 1:*)
    !
    INTEGER               :: ir
    INTEGER               :: i1, i2, i3
    INTEGER               :: iz
    INTEGER               :: iiz
    INTEGER               :: isite
    LOGICAL               :: offrange
    REAL(DP), ALLOCATABLE :: bg0(:,:)
#if defined(_OPENMP)
    REAL(DP), ALLOCATABLE :: bg1(:,:)
#endif
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    ALLOCATE(bg0(rismt%dfft%nr3, rismt%nsite))
    bg0 = 0.0_DP
    !
!$omp parallel default(shared) private(ir, i1, i2, i3, iz, isite, offrange, bg1)
#if defined(_OPENMP)
    ALLOCATE(bg1(rismt%dfft%nr3, rismt%nsite))
    bg1 = 0.0_DP
#endif
!$omp do
    DO ir = 1, rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
      !
      CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
      IF (offrange) THEN
        CYCLE
      END IF
      !
      IF (i3 < (rismt%dfft%nr3 - (rismt%dfft%nr3 / 2))) THEN
        iz = i3 + (rismt%dfft%nr3 / 2)
      ELSE
        iz = i3 - rismt%dfft%nr3 + (rismt%dfft%nr3 / 2)
      END IF
      iz = iz + 1
      !
      DO isite = 1, rismt%nsite
#if defined(_OPENMP)
        bg1(iz, isite) = bg1(iz, isite) + ar(ir, isite)
#else
        bg0(iz, isite) = bg0(iz, isite) + ar(ir, isite)
#endif
      END DO
      !
    END DO
!$omp end do
#if defined(_OPENMP)
!$omp critical
    bg0 = bg0 + bg1
!$omp end critical
    DEALLOCATE(bg1)
#endif
!$omp end parallel
    !
    CALL mp_sum(bg0, rismt%mp_site%intra_sitg_comm)
    !
    bg0 = bg0 / DBLE(rismt%dfft%nr1 * rismt%dfft%nr2)
    !
    DO isite = 1, rismt%nsite
      DO iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
        iiz = iz - rismt%lfft%izcell_start + 1
        ag0(iz, isite) = bg0(iiz, isite)
      END DO
    END DO
    !
    DEALLOCATE(bg0)
    !
  END SUBROUTINE extract_gxy0
  !
  SUBROUTINE sum_gxy0(ar, ag0)
    IMPLICIT NONE
    REAL(DP), INTENT(INOUT) :: ar(rismt%nr, 1:*)
    REAL(DP), INTENT(IN)    :: ag0(rismt%nrzl, 1:*)
    !
    INTEGER :: ir
    INTEGER :: i1, i2, i3
    INTEGER :: iz
    INTEGER :: isite
    LOGICAL :: offrange
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
!$omp parallel do default(shared) private(ir, i1, i2, i3, iz, isite, offrange)
    DO ir = 1, rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p
      !
      CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
      IF (offrange) THEN
        CYCLE
      END IF
      !
      IF (i3 < (rismt%dfft%nr3 - (rismt%dfft%nr3 / 2))) THEN
        iz = i3 + (rismt%dfft%nr3 / 2)
      ELSE
        iz = i3 - rismt%dfft%nr3 + (rismt%dfft%nr3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      DO isite = 1, rismt%nsite
        ar(ir, isite) = ar(ir, isite) + ag0(iz, isite)
      END DO
      !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE sum_gxy0
  !
END SUBROUTINE corrgxy0_laue
