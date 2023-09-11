!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE allocate_lauefft_rz(lauefft0, dzright, dzleft)
  !--------------------------------------------------------------------------
  !
  USE cell_base,   ONLY : at
  USE fft_support, ONLY : good_fft_order, good_fft_dimension
  USE kinds,       ONLY : DP
  USE lauefft,     ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  REAL(DP),           INTENT(IN)    :: dzright
  REAL(DP),           INTENT(IN)    :: dzleft
  !
  INTEGER  :: nzright
  INTEGER  :: nzleft
  INTEGER  :: nzrest
  INTEGER  :: nzres1
  INTEGER  :: nzres2
  REAL(DP) :: z0
  REAL(DP) :: hz
  !
  ! ... check lauefft0%dfft
  IF (lauefft0%dfft%nr3 < 1) THEN
    CALL errore(' allocate_lauefft_rz ', ' lauefft0%dfft%nr3 is not positive ', 1)
  END IF
  !
  ! ... size of unit cell
  z0 = 0.5_DP * at(3, 3)
  !
  ! ... step of R-space grid
  hz = 2.0_DP * z0 / DBLE(lauefft0%dfft%nr3)
  !
  ! ... count R-space grids for right
  nzright = 0
  IF (dzright > 0.0_DP) THEN
    nzright = INT(dzright / hz) + 1
  END IF
  !
  ! ... count R-space grids for left
  nzleft = 0
  IF (dzleft > 0.0_DP) THEN
    nzleft = INT(dzleft / hz) + 1
  END IF
  !
  ! ... define dimension of grids
  lauefft0%nrz     = lauefft0%dfft%nr3 + nzright + nzleft
  lauefft0%nrz     = good_fft_order(    lauefft0%nrz)
  lauefft0%nrzx    = good_fft_dimension(lauefft0%nrz)
  lauefft0%zstep   = hz
#if defined (__ESM_NOT_SYMMETRIC)
  lauefft0%zoffset = 0.5_DP * DBLE(MOD(lauefft0%dfft%nr3, 2)) * hz ! grid on z=0 must be
#else
  lauefft0%zoffset = 0.5_DP * hz ! grids shall be symmetric
#endif
  !
  ! ... check and correct number of grids
  IF (nzright < 1 .AND. nzleft < 1) THEN
    ! ... case that right or left are not.
    nzright = 0
    nzleft  = 0
    !
  ELSE IF (nzright < 1) THEN
    ! ... case that right is not and left is.
    nzright = 0
    nzleft  = lauefft0%nrz - lauefft0%dfft%nr3
    IF (nzleft < 1) THEN
      CALL errore(' allocate_lauefft_rz ', ' nzleft is not positive ', 1)
    END IF
    !
  ELSE IF (nzleft < 1) THEN
    ! ... case that right is and left is not.
    nzright = lauefft0%nrz - lauefft0%dfft%nr3
    nzleft  = 0
    IF (nzright < 1) THEN
      CALL errore(' allocate_lauefft_rz ', ' nzright is not positive ', 1)
    END IF
    !
  ELSE
    ! ... case that right and left are.
    nzrest  = lauefft0%nrz - lauefft0%dfft%nr3 - nzright - nzleft
    nzres1  = nzrest / 2
    nzres2  = nzrest  - nzres1
    nzright = nzright + nzres1
    nzleft  = nzleft  + nzres2
    IF (nzright < 1) THEN
      CALL errore(' allocate_lauefft_rz ', ' nzright is not positive ', 1)
    END IF
    IF (nzleft < 1) THEN
      CALL errore(' allocate_lauefft_rz ', ' nzleft is not positive ', 1)
    END IF
  END IF
  !
  IF (lauefft0%nrz /= (lauefft0%dfft%nr3 + nzright + nzleft)) THEN
    CALL errore(' allocate_lauefft_rz ', ' lauefft0%nrz is not consistent ', 1)
  END IF
  !
  ! ... set properties of unit cell
  IF (nzleft > 0) THEN
    lauefft0%izcell_start = nzleft + 1
    lauefft0%izcell_end   = nzleft + lauefft0%dfft%nr3
  ELSE
    lauefft0%izcell_start = 1
    lauefft0%izcell_end   = lauefft0%dfft%nr3
  END IF
  !
  ! ... set properties of right
  IF (nzright > 0) THEN
    lauefft0%xright = .TRUE.
    lauefft0%zright = z0 + DBLE(nzright) * hz
    lauefft0%izright_start = lauefft0%izcell_start + lauefft0%dfft%nr3 / 2
    lauefft0%izright_end   = lauefft0%izcell_end
    IF (lauefft0%izright_start > lauefft0%izright_end) THEN
      CALL errore(' allocate_lauefft_rz ', ' izright_start > izright_end ', 1)
    END IF
    !
  ELSE
    lauefft0%xright        = .FALSE.
    lauefft0%zright        = z0
    lauefft0%izright_start = lauefft0%izcell_end + 1
    lauefft0%izright_end   = lauefft0%izcell_end
  END IF
  !
  ! ... set properties of left
  IF (nzleft > 0) THEN
    lauefft0%xleft = .TRUE.
    lauefft0%zleft = -z0 - DBLE(nzleft) * hz
    lauefft0%izleft_start = lauefft0%izcell_start
    IF (lauefft0%xright) THEN
      lauefft0%izleft_end = lauefft0%izcell_start + lauefft0%dfft%nr3 / 2 - 1
    ELSE
      lauefft0%izleft_end = lauefft0%izcell_end - lauefft0%dfft%nr3 / 2
    END IF
    IF (lauefft0%izleft_start > lauefft0%izleft_end) THEN
      CALL errore(' allocate_lauefft_rz ', ' izleft_start > izleft_end ', 1)
    END IF
    !
  ELSE
    lauefft0%xleft        = .FALSE.
    lauefft0%zleft        = -z0
    lauefft0%izleft_start = lauefft0%izcell_start
    lauefft0%izleft_end   = lauefft0%izcell_start - 1
  END IF
  !
  ! ... check expanded cell
  IF ((.NOT. lauefft0%xright) .AND. (.NOT. lauefft0%xleft)) THEN
    CALL errore(' allocate_lauefft_rz ', ' expanded cell is not defined ', 1)
  END IF
  !
  ! ... set domain for Gxy = 0
  lauefft0%izright_start0 = lauefft0%izright_start
  lauefft0%izright_end0   = lauefft0%izright_end
  lauefft0%izleft_start0  = lauefft0%izleft_start
  lauefft0%izleft_end0    = lauefft0%izleft_end
  !
  ! ... set position of barriers
  lauefft0%izright_gedge = lauefft0%izright_start
  lauefft0%izleft_gedge  = lauefft0%izleft_end
  !
END SUBROUTINE allocate_lauefft_rz
!
!--------------------------------------------------------------------------
SUBROUTINE set_lauefft_offset_x(lauefft0, wright, wleft)
  !--------------------------------------------------------------------------
  !
  ! ... set offsets,
  ! ... i.e. set izright_start and izleft_end.
  !
  USE cell_base, ONLY : alat
  USE constants, ONLY : eps6
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  REAL(DP),           INTENT(IN)    :: wright
  REAL(DP),           INTENT(IN)    :: wleft
  !
  REAL(DP) :: wright_
  REAL(DP) :: wleft_
  REAL(DP) :: dwright
  REAL(DP) :: dwleft
  !
  ! ... check zstep
  IF (lauefft0%zstep <= 0.0_DP) THEN
    RETURN
  END IF
  !
  ! ... modify values of offsets
  wright_ = wright
  wleft_  = wleft
  !
  IF (lauefft0%xright .AND. lauefft0%xleft) THEN
    IF (wright < wleft) THEN
      wright_ = 0.5_DP * (wright + wleft)
      wleft_  = wright_
    END IF
  END IF
  !
  ! ... set offset of right
  IF (lauefft0%xright) THEN
    dwright = wright_ - lauefft0%zleft - lauefft0%zoffset + eps6 / alat
    lauefft0%izright_start = NINT(dwright / lauefft0%zstep) + 1
    lauefft0%izright_start = MAX(lauefft0%izright_start, lauefft0%izcell_start)
    !
    IF (lauefft0%izright_start > lauefft0%izright_end) THEN
      CALL errore(' set_lauefft_offset_x ', ' izright_start > izright_end ', 1)
    END IF
    !
    ! ... update offset for Gxy = 0
    lauefft0%izright_start0 = MIN(lauefft0%izright_start0, lauefft0%izright_start)
    !
    ! ... update offset barrier
    lauefft0%izright_gedge = MAX(lauefft0%izright_gedge, lauefft0%izright_start)
  END IF
  !
  ! ... set offset of left
  IF (lauefft0%xleft) THEN
    dwleft = wleft_ - lauefft0%zleft - lauefft0%zoffset - eps6 / alat
    lauefft0%izleft_end = NINT(dwleft / lauefft0%zstep) + 1
    lauefft0%izleft_end = MIN(lauefft0%izleft_end, lauefft0%izcell_end)
    !
    IF (lauefft0%izleft_end == lauefft0%izright_start) THEN
      lauefft0%izleft_end = lauefft0%izright_start - 1
    END IF
    !
    IF (lauefft0%izleft_start > lauefft0%izleft_end) THEN
      CALL errore(' set_lauefft_offset_x ', ' izleft_start > izleft_end ', 1)
    END IF
    !
    ! ... update offset for Gxy = 0
    lauefft0%izleft_end0 = MAX(lauefft0%izleft_end0, lauefft0%izleft_end)
    !
    ! ... update offset barrier
    lauefft0%izleft_gedge = MIN(lauefft0%izleft_gedge, lauefft0%izleft_end)
  END IF
  !
  IF (lauefft0%izleft_end >= lauefft0%izright_start) THEN
    CALL errore(' set_lauefft_offset_x ', ' izleft_end >= izright_start ', 1)
  END IF
  !
END SUBROUTINE set_lauefft_offset_x
!
!--------------------------------------------------------------------------
SUBROUTINE set_lauefft_offset0_x(lauefft0, wright1, wright2, wleft1, wleft2)
  !--------------------------------------------------------------------------
  !
  ! ... set offsets for Gxy = 0,
  ! ... i.e. set izright_start0 izright_end izleft_start0 and izleft_end0.
  !
  USE cell_base, ONLY : alat, at
  USE constants, ONLY : eps6
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  REAL(DP),           INTENT(IN)    :: wright1  ! solute-side
  REAL(DP),           INTENT(IN)    :: wright2  ! solvent-side
  REAL(DP),           INTENT(IN)    :: wleft1   ! solute-side
  REAL(DP),           INTENT(IN)    :: wleft2   ! solvent-side
  !
  REAL(DP) :: wright1_
  REAL(DP) :: wleft1_
  REAL(DP) :: dwright
  REAL(DP) :: dwleft
  !
  ! ... check zstep
  IF (lauefft0%zstep <= 0.0_DP) THEN
    RETURN
  END IF
  !
  ! ... modify values of solute-side offsets
  wright1_ = wright1
  wleft1_  = wleft1
  !
  IF (lauefft0%xright .AND. lauefft0%xleft) THEN
    IF (wright1 < wleft1) THEN
      wright1_ = 0.5_DP * (wright1 + wleft1)
      wleft1_  = wright1_
    END IF
  END IF
  !
  ! ... set offset of right
  IF (lauefft0%xright) THEN
    dwright = wright1_ - lauefft0%zleft - lauefft0%zoffset + eps6 / alat
    lauefft0%izright_start0 = NINT(dwright / lauefft0%zstep) + 1
    lauefft0%izright_start0 = MAX(lauefft0%izright_start0, lauefft0%izcell_start)
    lauefft0%izright_start0 = MIN(lauefft0%izright_start0, lauefft0%izright_start)
    !
    dwright = wright2 - lauefft0%zleft - lauefft0%zoffset
    lauefft0%izright_end0 = NINT(dwright / lauefft0%zstep) + 1
    lauefft0%izright_end0 = MIN(lauefft0%izright_end0, lauefft0%nrz)
    !
    IF (lauefft0%izright_start0 > lauefft0%izright_start) THEN
      CALL errore(' set_lauefft_offset0_x ', ' izright_start0 > izright_start ', 1)
    END IF
    IF (lauefft0%izright_end0 < lauefft0%izright_end) THEN
      CALL errore(' set_lauefft_offset0_x ', ' izright_end0 < izright_end ', 1)
    END IF
  END IF
  !
  ! ... set offset of left
  IF (lauefft0%xleft) THEN
    dwleft = wleft1_ - lauefft0%zleft - lauefft0%zoffset - eps6 / alat
    lauefft0%izleft_end0 = NINT(dwleft / lauefft0%zstep) + 1
    lauefft0%izleft_end0 = MIN(lauefft0%izleft_end0, lauefft0%izcell_end)
    lauefft0%izleft_end0 = MAX(lauefft0%izleft_end0, lauefft0%izleft_end)
    !
    IF (lauefft0%izleft_end0 == lauefft0%izright_start0) THEN
      lauefft0%izleft_end0 = lauefft0%izright_start0 - 1
    END IF
    !
    dwleft = wleft2 - lauefft0%zleft - lauefft0%zoffset
    lauefft0%izleft_start0 = NINT(dwleft / lauefft0%zstep) + 1
    lauefft0%izleft_start0 = MAX(lauefft0%izleft_start0, 1)
    !
    IF (lauefft0%izleft_start0 > lauefft0%izleft_start) THEN
      CALL errore(' set_lauefft_offset0_x ', ' izleft_start0 > izleft_start ', 1)
    END IF
    IF (lauefft0%izleft_end0 < lauefft0%izleft_end) THEN
      CALL errore(' set_lauefft_offset0_x ', ' izleft_end0 < izleft_end ', 1)
    END IF
  END IF
  !
  IF (lauefft0%izleft_end0 >= lauefft0%izright_start0) THEN
    CALL errore(' set_lauefft_offset0_x ', ' izleft_end0 >= izright_start0 ', 1)
  END IF
  !
END SUBROUTINE set_lauefft_offset0_x
!
!--------------------------------------------------------------------------
SUBROUTINE set_lauefft_barrier_x(lauefft0, wright, wleft)
  !--------------------------------------------------------------------------
  !
  ! ... set offset of barriers where g(z) = 0,
  ! ... i.e. set izright_gedge and izleft_gedge.
  !
  USE cell_base, ONLY : alat
  USE constants, ONLY : eps6
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  REAL(DP),           INTENT(IN)    :: wright
  REAL(DP),           INTENT(IN)    :: wleft
  !
  REAL(DP) :: dwright
  REAL(DP) :: dwleft
  !
  ! ... check zstep
  IF (lauefft0%zstep <= 0.0_DP) THEN
    RETURN
  END IF
  !
  ! ... set offset of right
  IF (lauefft0%xright) THEN
    dwright = wright - lauefft0%zleft - lauefft0%zoffset + eps6 / alat
    lauefft0%izright_gedge = NINT(dwright / lauefft0%zstep) + 1
    !
    IF (lauefft0%izright_gedge > lauefft0%izright_end) THEN
      CALL errore(' set_lauefft_barrier_x ', ' izright_gedge > izright_end ', 1)
    END IF
    IF (lauefft0%izright_gedge < lauefft0%izright_start) THEN
      CALL errore(' set_lauefft_barrier_x ', ' izright_gedge < izright_start ', 1)
    END IF
  END IF
  !
  ! ... set offset of left
  IF (lauefft0%xleft) THEN
    dwleft = wleft - lauefft0%zleft - lauefft0%zoffset - eps6 / alat
    lauefft0%izleft_gedge = NINT(dwleft / lauefft0%zstep) + 1
    !
    IF (lauefft0%izleft_gedge == lauefft0%izright_gedge) THEN
      lauefft0%izleft_gedge = lauefft0%izright_gedge - 1
    END IF
    !
    IF (lauefft0%izleft_start > lauefft0%izleft_gedge) THEN
      CALL errore(' set_lauefft_barrier_x ', ' izleft_start > izleft_gedge ', 1)
    END IF
    IF (lauefft0%izleft_end < lauefft0%izleft_gedge) THEN
      CALL errore(' set_lauefft_barrier_x ', ' izleft_end < izleft_gedge ', 1)
    END IF
  END IF
  !
END SUBROUTINE set_lauefft_barrier_x
!
!--------------------------------------------------------------------------
SUBROUTINE allocate_lauefft_gz(lauefft0, ngmt, ig1t, ig2t, ig3t, gt)
  !--------------------------------------------------------------------------
  !
  USE constants,     ONLY : eps4, tpi
  USE control_flags, ONLY : gamma_only
  USE kinds,         ONLY : DP
  USE lauefft,       ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  INTEGER,            INTENT(IN)    :: ngmt
  INTEGER,            INTENT(IN)    :: ig1t(1:*)
  INTEGER,            INTENT(IN)    :: ig2t(1:*)
  INTEGER,            INTENT(IN)    :: ig3t(1:*)
  REAL(DP),           INTENT(IN)    :: gt(3,1:*)
  !
  INTEGER               :: ig
  INTEGER               :: igz
  INTEGER               :: iigz
  INTEGER               :: isign
  INTEGER               :: mx, my, mz
  REAL(DP)              :: tz
  REAL(DP)              :: tt
  REAL(DP)              :: phi
  !
  REAL(DP), ALLOCATABLE :: gz(:)
  INTEGER,  ALLOCATABLE :: millz(:)
  LOGICAL,  ALLOCATABLE :: has_gz(:)
  INTEGER,  ALLOCATABLE :: igzmap(:)
  !
  INTEGER,  ALLOCATABLE :: igzsrt(:)
  REAL(DP), ALLOCATABLE :: gz2sort(:)
  REAL(DP), ALLOCATABLE :: gz_unsorted(:)
  INTEGER,  ALLOCATABLE :: millz_unsorted(:)
  !
  ! ... check lauefft0%dfft
  IF (lauefft0%dfft%nr3 < 1) THEN
    CALL errore(' allocate_lauefft_gz ', ' lauefft0%dfft%nr3 is not positive ', 1)
  END IF
  !
  ! ... alloc memory
  ALLOCATE(gz(    lauefft0%dfft%nr3))
  ALLOCATE(millz( lauefft0%dfft%nr3))
  ALLOCATE(has_gz(lauefft0%dfft%nr3))
  ALLOCATE(igzmap(lauefft0%dfft%nr3))
  gz     = 0.0_DP
  millz  = 0
  has_gz = .FALSE.
  igzmap = 0
  !
  ! ... count Gz-vectors
  DO ig = 1, ngmt
    mx = ig1t(ig)
    my = ig2t(ig)
    !
    isign = 1
 10 CONTINUE
    !
    mz  = isign * ig3t(ig)
    igz = mz + 1
    IF (igz < 1) THEN
      igz = igz + lauefft0%dfft%nr3
    END IF
    !
    IF (igz < 1 .OR. igz > lauefft0%dfft%nr3) THEN
      CALL errore(' allocate_lauefft_gz ', ' incorrect igz ', ig)
    END IF
    !
    IF (.NOT. has_gz(igz)) THEN
      tz = DBLE(isign) * gt(3, ig)
      gz(    igz) = tz
      millz( igz) = mz
      has_gz(igz) = .TRUE.
    END IF
    !
    IF (gamma_only .AND. isign > 0 .AND. mx == 0 .AND. my == 0) THEN
      isign = -1
      GOTO 10
    END IF
  END DO
  !
  lauefft0%ngz = COUNT(has_gz)
  IF (lauefft0%ngz < 1) THEN
    CALL errore(' allocate_lauefft_gz ', ' ngz is not positive ', 1)
  END IF
  !
  ! ... sort Gz-vectors
  ALLOCATE(igzsrt(        lauefft0%ngz))
  ALLOCATE(gz2sort(       lauefft0%ngz))
  ALLOCATE(gz_unsorted(   lauefft0%ngz))
  ALLOCATE(millz_unsorted(lauefft0%ngz))
  !
  iigz = 0
  DO igz = 1, lauefft0%dfft%nr3
    IF (has_gz(igz)) THEN
      iigz = iigz + 1
      tz = gz(   igz)
      mz = millz(igz)
      !
      IF (ABS(tz) > eps4) THEN
        gz2sort(iigz) = tz
      ELSE
        gz2sort(iigz) = 0.0_DP
      END IF
      !
      gz_unsorted(   iigz) = tz
      millz_unsorted(iigz) = mz
    END IF
  END DO
  !
  igzsrt(1) = 0
  CALL hpsort_eps(lauefft0%ngz, gz2sort, igzsrt, eps4)
  !
  ! ... define Gz-vectors
  ALLOCATE(lauefft0%nlz(  lauefft0%ngz))
  ALLOCATE(lauefft0%gz(   lauefft0%ngz))
  ALLOCATE(lauefft0%millz(lauefft0%ngz))
  !
  lauefft0%gz(   :) = gz_unsorted(   igzsrt(:))
  lauefft0%millz(:) = millz_unsorted(igzsrt(:))
  !
  lauefft0%gzzero = -1
  DO iigz = 1, lauefft0%ngz
    mz  = lauefft0%millz(iigz)
    igz = mz + 1
    IF (igz < 1) THEN
      igz = igz + lauefft0%dfft%nr3
    END IF
    !
    ! ... detect Gz = 0
    IF (ABS(lauefft0%gz(iigz)) <= eps4) THEN
      lauefft0%gzzero = iigz
    END IF
    !
    ! ... mapping iigz -> igz
    lauefft0%nlz(iigz) = igz
    !
    ! ... mapping igz -> iigz
    igzmap(igz) = iigz
    !
  END DO
  !
  IF (lauefft0%gzzero < 1) THEN
    CALL errore(' allocate_lauefft_gz ', ' gzzero was not detected ', 1)
  END IF
  !
  ! ... define index G -> Gz
  ALLOCATE(lauefft0%igtoigz(2, ngmt))
  !
  DO ig = 1, ngmt
    mx = ig1t(ig)
    my = ig2t(ig)
    !
    isign = 1
 11 CONTINUE
    !
    mz  = isign * ig3t(ig)
    igz = mz + 1
    IF (igz < 1) THEN
      igz = igz + lauefft0%dfft%nr3
    END IF
    !
    IF (isign > 0) THEN
      lauefft0%igtoigz(1, ig) = igzmap(igz)
      lauefft0%igtoigz(2, ig) = -1
    ELSE
      lauefft0%igtoigz(2, ig) = igzmap(igz)
    END IF
    !
    IF (gamma_only .AND. isign > 0 .AND. mx == 0 .AND. my == 0) THEN
      isign = -1
      GOTO 11
    END IF
  END DO
  !
  ! ... evaluate phase factor of Z <--> Gz
  ALLOCATE(lauefft0%zphase(lauefft0%ngz))
#if defined (__ESM_NOT_SYMMETRIC)
  !
  lauefft0%zphase = CMPLX(1.0_DP, 0.0_DP, KIND=DP)
  !
#else
  !
  IF (MOD(lauefft0%dfft%nr3, 2) == 1) THEN
    ! ... if odd number of meshs, do not shift
    lauefft0%zphase = CMPLX(1.0_DP, 0.0_DP, KIND=DP)
    !
  ELSE
    ! ... if even number of meshs, shift a half mesh
    lauefft0%zphase = CMPLX(0.0_DP, 0.0_DP, KIND=DP)
    !
    DO iigz = 1, lauefft0%ngz
      tz  = lauefft0%gz(iigz)
      phi = tpi * tz * 0.5_DP * lauefft0%zstep
      lauefft0%zphase(iigz) = CMPLX(COS(phi), -SIN(phi), KIND=DP)
    END DO
    !
  END IF
  !
#endif
  ! ... dealloc memory
  DEALLOCATE(gz)
  DEALLOCATE(millz)
  DEALLOCATE(has_gz)
  DEALLOCATE(igzmap)
  DEALLOCATE(igzsrt)
  DEALLOCATE(gz2sort)
  DEALLOCATE(gz_unsorted)
  DEALLOCATE(millz_unsorted)
  !
END SUBROUTINE allocate_lauefft_gz
!
!--------------------------------------------------------------------------
SUBROUTINE allocate_lauefft_gz_exp(lauefft0, gcutm)
  !--------------------------------------------------------------------------
  !
  USE constants, ONLY : tpi
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  REAL(DP),           INTENT(IN)    :: gcutm
  !
  INTEGER               :: i
  INTEGER               :: ni
  INTEGER               :: igz
  INTEGER               :: iigz
  INTEGER               :: mz
  REAL(DP)              :: bg
  REAL(DP)              :: tz
  REAL(DP)              :: tt
  REAL(DP)              :: phi
  !
  REAL(DP), ALLOCATABLE :: gz(:)
  INTEGER,  ALLOCATABLE :: millz(:)
  !
  ! ... check lauefft0%dfft
  IF (lauefft0%nrz < 1) THEN
    CALL errore(' allocate_lauefft_gz_exp ', ' lauefft0%nrz is not positive ', 1)
  END IF
  !
  IF (lauefft0%nrzx < 1) THEN
    CALL errore(' allocate_lauefft_gz_exp ', ' lauefft0%nrzx is not positive ', 1)
  END IF
  !
  ! ... set variable
  ni = (lauefft0%nrz - 1) / 2
  bg = 1.0_DP / (lauefft0%zright - lauefft0%zleft)
  !
  ! ... alloc memory
  ALLOCATE(gz(   2 * ni + 1))
  ALLOCATE(millz(2 * ni + 1))
  !
  ! ... count Gz-vectors
  lauefft0%ngz_x = 0
  !
  DO i = -ni, ni
    tz = DBLE(i) * bg
    tt = tz * tz
    IF (tt <= gcutm) THEN
      lauefft0%ngz_x = lauefft0%ngz_x + 1
      gz(   lauefft0%ngz_x) = tz
      millz(lauefft0%ngz_x) = i
    END IF
  END DO
  !
  ! ... define Gz-vectors
  lauefft0%gzzero_x = -1
  ALLOCATE(lauefft0%nlz_x(  lauefft0%ngz_x))
  ALLOCATE(lauefft0%gz_x(   lauefft0%ngz_x))
  ALLOCATE(lauefft0%millz_x(lauefft0%ngz_x))
  !
  DO iigz = 1, lauefft0%ngz_x
    mz  = millz(iigz)
    igz = mz + 1
    IF (igz < 1) THEN
      igz = igz + lauefft0%nrz
    END IF
    !
    IF (mz == 0) THEN
      lauefft0%gzzero_x = iigz
    END IF
    !
    lauefft0%nlz_x(  iigz) = igz
    lauefft0%gz_x(   iigz) = gz(iigz)
    lauefft0%millz_x(iigz) = mz
  END DO
  !
  IF (lauefft0%gzzero_x < 1) THEN
    CALL errore(' allocate_lauefft_gz_exp ', ' gzzero_x was not detected ', 1)
  END IF
  !
  ! ... evaluate phase factor of Z <--> Gz
  ALLOCATE(lauefft0%zphase_x(lauefft0%ngz_x))
#if defined (__ESM_NOT_SYMMETRIC)
  !
  lauefft0%zphase_x = CMPLX(1.0_DP, 0.0_DP, KIND=DP)
  !
#else
  !
  IF (MOD(lauefft0%dfft%nr3, 2) == 1) THEN
    ! ... if odd number of meshs, do not shift
    lauefft0%zphase_x = CMPLX(1.0_DP, 0.0_DP, KIND=DP)
    !
  ELSE
    ! ... if even number of meshs, shift a half mesh
    lauefft0%zphase_x = CMPLX(0.0_DP, 0.0_DP, KIND=DP)
    !
    DO iigz = 1, lauefft0%ngz_x
      tz  = lauefft0%gz_x(iigz)
      phi = tpi * tz * 0.5_DP * lauefft0%zstep
      lauefft0%zphase_x(iigz) = CMPLX(COS(phi), -SIN(phi), KIND=DP)
    END DO
    !
  END IF
  !
#endif
  ! ... dealloc memory
  DEALLOCATE(gz)
  DEALLOCATE(millz)
  !
END SUBROUTINE allocate_lauefft_gz_exp
!
!--------------------------------------------------------------------------
SUBROUTINE allocate_lauefft_gxy(lauefft0, ngmt, ig1t, ig2t, gt, comm)
  !--------------------------------------------------------------------------
  !
  USE constants,     ONLY : eps8
  USE control_flags, ONLY : gamma_only
  USE kinds,         ONLY : DP
  USE lauefft,       ONLY : lauefft_type
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  INTEGER,            INTENT(IN)    :: ngmt
  INTEGER,            INTENT(IN)    :: ig1t(1:*)
  INTEGER,            INTENT(IN)    :: ig2t(1:*)
  REAL(DP),           INTENT(IN)    :: gt(3,1:*)
  INTEGER,            INTENT(IN)    :: comm
  !
  INTEGER               :: ig
  INTEGER               :: igx, igy
  INTEGER               :: igxy
  INTEGER               :: mx, my
  REAL(DP)              :: tx, ty
  REAL(DP)              :: tt
  !
  REAL(DP), ALLOCATABLE :: gxy(:,:,:)
  INTEGER,  ALLOCATABLE :: millxy(:,:,:)
  LOGICAL,  ALLOCATABLE :: has_gxy(:,:)
  INTEGER,  ALLOCATABLE :: igxymap(:,:)
  !
  INTEGER,  ALLOCATABLE :: igxysrt(:)
  REAL(DP), ALLOCATABLE :: gxy2sort(:)
  REAL(DP), ALLOCATABLE :: gxy_unsorted(:,:)
  REAL(DP), ALLOCATABLE :: gnxy_unsorted(:)
  REAL(DP), ALLOCATABLE :: ggxy_unsorted(:)
  INTEGER,  ALLOCATABLE :: millxy_unsorted(:,:)
  !
  ! ... check lauefft0%dfft
  IF (lauefft0%dfft%nr1 < 1) THEN
    CALL errore(' allocate_lauefft_gxy ', ' lauefft0%dfft%nr1 is not positive ', 1)
  END IF
  !
  IF (lauefft0%dfft%nr2 < 1) THEN
    CALL errore(' allocate_lauefft_gxy ', ' lauefft0%dfft%nr2 is not positive ', 1)
  END IF
  !
  ! ... alloc memory
  ALLOCATE(gxy(   2, lauefft0%dfft%nr1, lauefft0%dfft%nr2))
  ALLOCATE(millxy(2, lauefft0%dfft%nr1, lauefft0%dfft%nr2))
  ALLOCATE(has_gxy(  lauefft0%dfft%nr1, lauefft0%dfft%nr2))
  ALLOCATE(igxymap(  lauefft0%dfft%nr1, lauefft0%dfft%nr2))
  gxy     = 0.0_DP
  millxy  = 0
  has_gxy = .FALSE.
  igxymap = 0
  !
  ! ... count Gxy-vectors
  DO ig = 1, ngmt
    mx  = ig1t(ig)
    igx = mx + 1
    IF (igx < 1) THEN
      igx = igx + lauefft0%dfft%nr1
    END IF
    !
    IF (igx < 1 .OR. igx > lauefft0%dfft%nr1) THEN
      CALL errore(' allocate_lauefft_gxy ', ' incorrect igx ', ig)
    END IF
    !
    my  = ig2t(ig)
    igy = my + 1
    IF (igy < 1) THEN
      igy = igy + lauefft0%dfft%nr2
    END IF
    !
    IF (igy < 1 .OR. igy > lauefft0%dfft%nr2) THEN
      CALL errore(' allocate_lauefft_gxy ', ' incorrect igy ', ig)
    END IF
    !
    IF (.NOT. has_gxy(igx, igy)) THEN
      tx = gt(1, ig)
      ty = gt(2, ig)
      gxy(   1, igx, igy) = tx
      gxy(   2, igx, igy) = ty
      millxy(1, igx, igy) = mx
      millxy(2, igx, igy) = my
      has_gxy(  igx, igy) = .TRUE.
    END IF
  END DO
  !
  lauefft0%ngxy = COUNT(has_gxy)
  IF (lauefft0%ngxy < 1) THEN
    CALL errore(' allocate_lauefft_gxy ', ' ngxy is not positive ', 1)
  END IF
  !
  lauefft0%ngxy_g = lauefft0%ngxy
  CALL mp_sum(lauefft0%ngxy_g, comm)
  !
  ! ... sort Gxy-vectors
  ALLOCATE(igxysrt(           lauefft0%ngxy))
  ALLOCATE(gxy2sort(          lauefft0%ngxy))
  ALLOCATE(gxy_unsorted(   2, lauefft0%ngxy))
  ALLOCATE(gnxy_unsorted(     lauefft0%ngxy))
  ALLOCATE(ggxy_unsorted(     lauefft0%ngxy))
  ALLOCATE(millxy_unsorted(2, lauefft0%ngxy))
  !
  igxy = 0
  DO igx = 1, lauefft0%dfft%nr1
    DO igy = 1, lauefft0%dfft%nr2
      IF (has_gxy(igx, igy)) THEN
        igxy = igxy + 1
        tx = gxy(1, igx, igy)
        ty = gxy(2, igx, igy)
        tt = tx * tx + ty * ty
        mx = millxy(1, igx, igy)
        my = millxy(2, igx, igy)
        !
        IF (tt > eps8) THEN
          gxy2sort(igxy) = tt
        ELSE
          gxy2sort(igxy) = 0.0_DP
        END IF
        !
        gxy_unsorted(   1, igxy) = tx
        gxy_unsorted(   2, igxy) = ty
        gnxy_unsorted(     igxy) = SQRT(tt)
        ggxy_unsorted(     igxy) = tt
        millxy_unsorted(1, igxy) = mx
        millxy_unsorted(2, igxy) = my
      END IF
    END DO
  END DO
  !
  igxysrt(1) = 0
  CALL hpsort_eps(lauefft0%ngxy, gxy2sort, igxysrt, eps8)
  !
  ! ... define Gxy-vectors
  ALLOCATE(lauefft0%nlxy(     lauefft0%ngxy))
  ALLOCATE(lauefft0%nlmxy(    lauefft0%ngxy))
  ALLOCATE(lauefft0%gxy(   2, lauefft0%ngxy))
  ALLOCATE(lauefft0%gnxy(     lauefft0%ngxy))
  ALLOCATE(lauefft0%ggxy(     lauefft0%ngxy))
  ALLOCATE(lauefft0%millxy(2, lauefft0%ngxy))
  !
  lauefft0%gxy(   1, :) = gxy_unsorted(   1, igxysrt(:))
  lauefft0%gxy(   2, :) = gxy_unsorted(   2, igxysrt(:))
  lauefft0%gnxy(     :) = gnxy_unsorted(     igxysrt(:))
  lauefft0%ggxy(     :) = ggxy_unsorted(     igxysrt(:))
  lauefft0%millxy(1, :) = millxy_unsorted(1, igxysrt(:))
  lauefft0%millxy(2, :) = millxy_unsorted(2, igxysrt(:))
  !
  IF (lauefft0%ggxy(1) <= eps8) THEN
    lauefft0%gxystart = 2
  ELSE
    lauefft0%gxystart = 1
  END IF
  !
  DO igxy = 1, lauefft0%ngxy
    mx  = lauefft0%millxy(1, igxy)
    igx = mx + 1
    IF (igx < 1) THEN
      igx = igx + lauefft0%dfft%nr1
    END IF
    !
    my  = lauefft0%millxy(2, igxy)
    igy = my + 1
    IF (igy < 1) THEN
      igy = igy + lauefft0%dfft%nr2
    END IF
    !
    ! ... 2D-FFT index
    IF (lauefft0%dfft%lpara) THEN
      lauefft0%nlxy(igxy) = &
      & (lauefft0%dfft%isind(igx + (igy - 1) * lauefft0%dfft%nr1x) - 1) * lauefft0%dfft%nr3x
    ELSE
      lauefft0%nlxy(igxy) = igx + (igy - 1) * lauefft0%dfft%nr1x
    END IF
    !
    ! ... mapping (igx,igy) -> igxy
    igxymap(igx, igy) = igxy
    !
  END DO
  !
  IF (gamma_only) THEN
    DO igxy = 1, lauefft0%ngxy
      mx  = -lauefft0%millxy(1, igxy)
      igx = mx + 1
      IF (igx < 1) THEN
        igx = igx + lauefft0%dfft%nr1
      END IF
      !
      my  = -lauefft0%millxy(2, igxy)
      igy = my + 1
      IF (igy < 1) THEN
        igy = igy + lauefft0%dfft%nr2
      END IF
      !
      ! ... 2D-FFT index (0 > Gxy)
      IF (lauefft0%dfft%lpara) THEN
        lauefft0%nlmxy(igxy) = &
        & (lauefft0%dfft%isind(igx + (igy - 1) * lauefft0%dfft%nr1x) - 1) * lauefft0%dfft%nr3x
      ELSE
        lauefft0%nlmxy(igxy) = igx + (igy - 1) * lauefft0%dfft%nr1x
      END IF
      !
    END DO
  END IF
  !
  ! ... define index G -> Gxy
  ALLOCATE(lauefft0%igtoigxy(ngmt))
  !
  DO ig = 1, ngmt
    mx  = ig1t(ig)
    igx = mx + 1
    IF (igx < 1) THEN
      igx = igx + lauefft0%dfft%nr1
    END IF
    !
    my  = ig2t(ig)
    igy = my + 1
    IF (igy < 1) THEN
      igy = igy + lauefft0%dfft%nr2
    END IF
    !
    lauefft0%igtoigxy(ig) = igxymap(igx, igy)
  END DO
  !
  ! ... dealloc memory
  DEALLOCATE(gxy)
  DEALLOCATE(millxy)
  DEALLOCATE(has_gxy)
  DEALLOCATE(igxymap)
  DEALLOCATE(igxysrt)
  DEALLOCATE(gxy2sort)
  DEALLOCATE(gxy_unsorted)
  DEALLOCATE(gnxy_unsorted)
  DEALLOCATE(ggxy_unsorted)
  DEALLOCATE(millxy_unsorted)
  !
END SUBROUTINE allocate_lauefft_gxy
!
!--------------------------------------------------------------------------
SUBROUTINE gxyshells(lauefft0, vc)
  !--------------------------------------------------------------------------
  !
  USE constants, ONLY : eps8
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : lauefft_type
  !
  IMPLICIT NONE
  !
  TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
  LOGICAL,            INTENT(IN)    :: vc
  !
  INTEGER :: ng, igl
  !
  ! ... deallocate memory, if needed
  IF (ASSOCIATED(lauefft0%glxy))      DEALLOCATE(lauefft0%glxy)
  IF (ASSOCIATED(lauefft0%igtonglxy)) DEALLOCATE(lauefft0%igtonglxy)
  !
  IF (vc) THEN
    !
    ! ... in case of a variable cell run each Gxy vector has its shell
    ALLOCATE(lauefft0%glxy(     lauefft0%ngxy))
    ALLOCATE(lauefft0%igtonglxy(lauefft0%ngxy))
    !
    lauefft0%nglxy = lauefft0%ngxy
    lauefft0%glxy  = lauefft0%ggxy
    DO ng = 1, lauefft0%ngxy
       lauefft0%igtonglxy(ng) = ng
    END DO
    !
  ELSE
    !
    ! ... Gxy vectors are grouped in shells with the same norm
    ALLOCATE(lauefft0%igtonglxy(lauefft0%ngxy))
    !
    lauefft0%nglxy        = 1
    lauefft0%igtonglxy(1) = 1
    DO ng = 2, lauefft0%ngxy
      IF (lauefft0%ggxy(ng) > lauefft0%ggxy(ng - 1) + eps8) THEN
        lauefft0%nglxy = lauefft0%nglxy + 1
      END IF
      lauefft0%igtonglxy(ng) = lauefft0%nglxy
    END DO
    !
    ALLOCATE(lauefft0%glxy(lauefft0%nglxy))
    !
    lauefft0%glxy(1) = lauefft0%ggxy(1)
    igl = 1
    DO ng = 2, lauefft0%ngxy
      IF (lauefft0%ggxy(ng) > lauefft0%ggxy(ng - 1) + eps8) THEN
        igl = igl + 1
        lauefft0%glxy(igl) = lauefft0%ggxy(ng)
      END IF
    END DO
    !
    IF (igl /= lauefft0%nglxy) THEN
      CALL errore(' gxyshells ', ' igl <> ngl ', lauefft0%nglxy)
    END IF
  END IF
  !
END SUBROUTINE gxyshells
