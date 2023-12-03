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
MODULE solvavg
  !--------------------------------------------------------------------------
  !
  ! ... planar average of solvents.
  !
  USE constants, ONLY : BOHR_RADIUS_ANGS
  USE cell_base, ONLY : at, alat
  USE fft_types, ONLY : fft_type_descriptor, fft_index_to_3d
  USE io_global, ONLY : ionode
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : lauefft_type
  USE mp,        ONLY : mp_sum
  USE parallel_include
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
  INTEGER, PARAMETER :: LEN_LABEL    = 20
  !
#if defined (__DEBUG_RISM)
  INTEGER, PARAMETER :: MAX_NUM_DATA = 256
#else
  INTEGER, PARAMETER :: MAX_NUM_DATA = 64
#endif
  !
  ! ... related 3D-FFT
  TYPE(fft_type_descriptor), POINTER     :: dfft => NULL()
  !
  ! ... related Laue-FFT
  TYPE(lauefft_type),        POINTER     :: lfft => NULL()
  !
  ! ... MPI-communicator
  INTEGER                                :: intra_group_comm = MPI_COMM_NULL
  !
  ! ... radial data, or not
  LOGICAL                                :: radial = .FALSE.
  !
  ! ... number of data
  INTEGER                                :: ndata = 0
  !
  ! ... labels of data
  CHARACTER(LEN=LEN_LABEL),  ALLOCATABLE :: label(:)
  !
  ! ... data of planar average
  REAL(DP),                  ALLOCATABLE :: rdata(:,:)
  !
  ! ... public components
  PUBLIC :: solvavg_size
  PUBLIC :: solvavg_init
  PUBLIC :: solvavg_clear
  PUBLIC :: solvavg_put
  PUBLIC :: solvavg_add
  PUBLIC :: solvavg_print
  !
  ! ... overload subroutines
  INTERFACE solvavg_init
    MODULE PROCEDURE solvavg_init_3d
    MODULE PROCEDURE solvavg_init_laue
  END INTERFACE
  !
  INTERFACE solvavg_put
    MODULE PROCEDURE solvavg_put_real
    MODULE PROCEDURE solvavg_put_laue
  END INTERFACE
  !
  INTERFACE solvavg_add
    MODULE PROCEDURE solvavg_add_real
    MODULE PROCEDURE solvavg_add_laue
  END INTERFACE
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  INTEGER FUNCTION solvavg_size()
    !--------------------------------------------------------------------------
    !
    ! ... get ndata
    !
    IMPLICIT NONE
    !
    solvavg_size = ndata
    !
  END FUNCTION solvavg_size
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_init_3d(dfft_, comm, is_radial)
    !--------------------------------------------------------------------------
    !
    ! ... initialize this module with 3D-FFT.
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), TARGET, INTENT(IN) :: dfft_
    INTEGER,                           INTENT(IN) :: comm
    LOGICAL,                           INTENT(IN) :: is_radial
    !
    dfft => dfft_
    !
    intra_group_comm = comm
    !
    radial = is_radial
    !
    ndata = 0
    ALLOCATE(label(          MAX_NUM_DATA))
    ALLOCATE(rdata(dfft%nr3, MAX_NUM_DATA))
    !
  END SUBROUTINE solvavg_init_3d
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_init_laue(lfft_, comm, is_radial)
    !--------------------------------------------------------------------------
    !
    ! ... initialize this module with Laue-FFT.
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), TARGET, INTENT(IN) :: lfft_
    INTEGER,                    INTENT(IN) :: comm
    LOGICAL,                    INTENT(IN) :: is_radial
    !
    lfft => lfft_
    !
    intra_group_comm = comm
    !
    radial = is_radial
    !
    ndata = 0
    ALLOCATE(label(          MAX_NUM_DATA))
    ALLOCATE(rdata(lfft%nrz, MAX_NUM_DATA))
    !
  END SUBROUTINE solvavg_init_laue
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_clear()
    !--------------------------------------------------------------------------
    !
    ! ... finalize this module.
    !
    IMPLICIT NONE
    !
    IF (ASSOCIATED(dfft)) THEN
      NULLIFY(dfft)
    END IF
    !
    IF (ASSOCIATED(lfft)) THEN
      NULLIFY(lfft)
    END IF
    !
    intra_group_comm = MPI_COMM_NULL
    !
    ndata = 0
    IF (ALLOCATED(label)) THEN
      DEALLOCATE(label)
    END IF
    IF (ALLOCATED(rdata)) THEN
      DEALLOCATE(rdata)
    END IF
    !
  END SUBROUTINE solvavg_clear
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_put_real(title, integrate, zuv)
    !--------------------------------------------------------------------------
    !
    ! ... put solvent data as R-space.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: title
    LOGICAL,          INTENT(IN) :: integrate
    REAL(DP),         INTENT(IN) :: zuv(:)  !dimension(nnr)
    !
    IF (ndata < MAX_NUM_DATA) THEN
      !
      ndata = ndata + 1
      label(   ndata) = TRIM(title)
      rdata(:, ndata) = 0.0_DP
      !
      CALL solvavg_add_real(ndata, integrate, zuv)
      !
    END IF
    !
  END SUBROUTINE solvavg_put_real
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_put_laue(title, integrate, zuv, nrz, expanded, igxy)
    !--------------------------------------------------------------------------
    !
    ! ... put solvent data as Laue-rep.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*),  INTENT(IN) :: title
    LOGICAL,           INTENT(IN) :: integrate
    COMPLEX(DP),       INTENT(IN) :: zuv(:)  !dimension(nrz*ngxy)
    INTEGER,           INTENT(IN) :: nrz
    LOGICAL,           INTENT(IN) :: expanded
    INTEGER, OPTIONAL, INTENT(IN) :: igxy
    !
    IF (ndata < MAX_NUM_DATA) THEN
      !
      ndata = ndata + 1
      label(   ndata) = TRIM(title)
      rdata(:, ndata) = 0.0_DP
      !
      IF (PRESENT(igxy)) THEN
        CALL solvavg_add_laue(ndata, integrate, zuv, nrz, expanded, igxy)
      ELSE
        CALL solvavg_add_laue(ndata, integrate, zuv, nrz, expanded)
      END IF
      !
    END IF
    !
  END SUBROUTINE solvavg_put_laue
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_add_real(idata, integrate, zuv)
    !--------------------------------------------------------------------------
    !
    ! ... add solvent data as R-space, at idata.
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: idata
    LOGICAL,  INTENT(IN) :: integrate
    REAL(DP), INTENT(IN) :: zuv(:)  !dimension(nnr)
    !
    LOGICAL               :: laue
    INTEGER               :: ir, nr
    INTEGER               :: i1, i2, i3
    INTEGER               :: n1, n2, n3
    INTEGER               :: nrz
    INTEGER               :: irz
    INTEGER               :: irz0
    LOGICAL               :: offrange
    REAL(DP)              :: area_xy
    REAL(DP), ALLOCATABLE :: ztmp(:)
    !
    ! ... check dfft, lfft
    IF (.NOT. (ASSOCIATED(dfft) .OR. ASSOCIATED(lfft))) THEN
      RETURN
    END IF
    !
    laue = .FALSE.
    IF (ASSOCIATED(lfft)) THEN
      laue = .TRUE.
    END IF
    !
    ! ... FFT box
    IF (laue) THEN
      n1   = lfft%dfft%nr1
      n2   = lfft%dfft%nr2
      n3   = lfft%dfft%nr3
      nr   = lfft%dfft%nr1x * lfft%dfft%my_nr2p * lfft%dfft%my_nr3p
      nrz  = lfft%nrz
      irz0 = lfft%izcell_start
    ELSE
      n1   = dfft%nr1
      n2   = dfft%nr2
      n3   = dfft%nr3
      nr   = dfft%nr1x * dfft%my_nr2p * dfft%my_nr3p
      nrz  = dfft%nr3
      irz0 = 1
    END IF
    !
    ! ... allocate memory
    ALLOCATE(ztmp(nrz))
    !
    ztmp = 0.0_DP
    !
    ! ... calculate planar average
    DO ir = 1, nr
      !
      IF (laue) THEN
        CALL fft_index_to_3d(ir, lfft%dfft, i1, i2, i3, offrange)
      ELSE
        CALL fft_index_to_3d(ir, dfft, i1, i2, i3, offrange)
      END IF
      !
      IF (offrange) THEN
        CYCLE
      END IF
      !
      IF (i3 >= (n3 - (n3 / 2))) THEN
        irz = i3 - n3
      ELSE
        irz = i3
      END IF
      irz = irz + (n3 / 2)
      irz = irz + irz0
      !
      ztmp(irz) = ztmp(irz) + zuv(ir)
      !
    END DO
    !
    CALL mp_sum(ztmp, intra_group_comm)
    !
    IF (integrate) THEN
      area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
      ztmp = ztmp * (area_xy / DBLE(n1 * n2))
    ELSE
      ztmp = ztmp / DBLE(n1 * n2)
    END IF
    !
    ! ... add to data
    IF (1 <= idata .AND. idata <= ndata) THEN
      DO irz = 1, nrz
        rdata(irz, idata) = rdata(irz, idata) + ztmp(irz)
      END DO
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(ztmp)
    !
  END SUBROUTINE solvavg_add_real
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_add_laue(idata, integrate, zuv, nrz, expanded, igxy)
    !--------------------------------------------------------------------------
    !
    ! ... add solvent data as Laue-rep., at idata.
    !
    IMPLICIT NONE
    !
    INTEGER,           INTENT(IN) :: idata
    LOGICAL,           INTENT(IN) :: integrate
    COMPLEX(DP),       INTENT(IN) :: zuv(:)  !dimension(nrz*ngxy)
    INTEGER,           INTENT(IN) :: nrz
    LOGICAL,           INTENT(IN) :: expanded
    INTEGER, OPTIONAL, INTENT(IN) :: igxy
    !
    INTEGER                  :: irz
    INTEGER                  :: jgxy
    REAL(DP)                 :: area_xy
    REAL(DP)                 :: zr
    REAL(DP)                 :: zi
    COMPLEX(DP), ALLOCATABLE :: ztmp(:)
    !
    ! ... check lfft
    IF (.NOT. ASSOCIATED(lfft)) THEN
      RETURN
    END IF
    !
    ! ... check nrz
    IF (expanded) THEN
      IF (nrz < lfft%nrz) THEN
        RETURN
      END IF
    ELSE
      IF (nrz < lfft%dfft%nr3) THEN
        RETURN
      END IF
    END IF
    !
    ! ... check igxy
    IF (PRESENT(igxy)) THEN
      jgxy = igxy
    ELSE
      jgxy = -1
    END IF
    !
    ! ... allocate memory
    ALLOCATE(ztmp(lfft%nrz))
    !
    ztmp = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    !
    IF (jgxy < 1) THEN
      !
      ! ... in case Gxy = 0
      IF (lfft%gxystart > 1) THEN
        IF (expanded) THEN
          DO irz = 1, lfft%nrz
            ztmp(irz) = zuv(irz)
          END DO
        ELSE
          DO irz = lfft%izcell_start, lfft%izcell_end
            ztmp(irz) = zuv(irz - lfft%izcell_start + 1)
          END DO
        END IF
      END IF
      !
    ELSE
      !
      ! ... in case Gxy /= 0
      IF (ionode .AND. (jgxy <= lfft%ngxy)) THEN
        IF (expanded) THEN
          DO irz = 1, lfft%nrz
            ztmp(irz) = zuv(irz + nrz * (jgxy - 1))
          END DO
        ELSE
          DO irz = lfft%izcell_start, lfft%izcell_end
            ztmp(irz) = zuv(irz - lfft%izcell_start + 1 + nrz * (jgxy - 1))
          END DO
        END IF
      END IF
      !
    END IF
    !
    CALL mp_sum(ztmp, intra_group_comm)
    !
    IF (integrate) THEN
      area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
      ztmp = ztmp * area_xy
    END IF
    !
    ! ... add to data
    IF (1 <= idata .AND. idata <= ndata) THEN
      DO irz = 1, lfft%nrz
        zr = DBLE( ztmp(irz))
        zi = AIMAG(ztmp(irz))
        !rdata(irz, idata) = rdata(irz, idata) + SQRT(zr * zr + zi * zi)
        rdata(irz, idata) = rdata(irz, idata) + zr
      END DO
    END IF
    !
    ! ... deallocate memory
    DEALLOCATE(ztmp)
    !
  END SUBROUTINE solvavg_add_laue
  !
  !--------------------------------------------------------------------------
  SUBROUTINE solvavg_print(filename, comment, ierr)
    !--------------------------------------------------------------------------
    !
    ! ... print planar average of solvents.
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    CHARACTER(LEN=*), INTENT(IN)  :: comment
    INTEGER,          INTENT(OUT) :: ierr
    !
    INTEGER :: iun
    !
    ierr = 0
    !
    CALL open_file(iun, filename, ierr)
    IF (ierr /= 0) THEN
      RETURN
    END IF
    !
    IF (ionode .AND. LEN_TRIM(comment) > 0) THEN
      WRITE(iun, '("#",A)') TRIM(comment)
    END IF
    !
    CALL write_data(iun)
    !
    CALL close_file(iun)
    !
  END SUBROUTINE solvavg_print
  !
  SUBROUTINE open_file(iun, filename, ierr)
    IMPLICIT NONE
    INTEGER,          INTENT(OUT) :: iun
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER,          INTENT(OUT) :: ierr
    !
    INTEGER :: ios
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    iun = find_free_unit()
    IF (ionode) THEN
      OPEN(unit=iun, file=filename, status='unknown', form='formatted', action='write', iostat=ios)
      ios = ABS(ios)
    ELSE
      ios = 0
    END IF
    !
    CALL mp_sum(ios, intra_group_comm)
    ierr = ios
  END SUBROUTINE open_file
  !
  SUBROUTINE close_file(iun)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iun
    !
    IF (ionode) THEN
      CLOSE(unit=iun)
    END IF
  END SUBROUTINE close_file
  !
  SUBROUTINE write_data(iun)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: iun
    !
    LOGICAL                  :: laue
    INTEGER                  :: idata
    INTEGER                  :: irz
    INTEGER                  :: nrz
    REAL(DP)                 :: rz
    REAL(DP)                 :: dz
    REAL(DP)                 :: z0
    REAL(DP)                 :: c
#if defined (__DEBUG_RISM)
    INTEGER                  :: ilabel
    INTEGER                  :: nlabel
#endif
    CHARACTER(LEN=LEN_LABEL) :: label1
    !
    IF (.NOT. ionode) THEN
      RETURN
    END IF
    !
    ! ... check dfft, lfft
    IF (.NOT. (ASSOCIATED(dfft) .OR. ASSOCIATED(lfft))) THEN
      RETURN
    END IF
    !
    laue = .FALSE.
    IF (ASSOCIATED(lfft)) THEN
      laue = .TRUE.
    END IF
    !
    ! ... FFT box
    IF (laue) THEN
      nrz = lfft%nrz
      dz  = lfft%zstep
      IF (.NOT. radial) THEN
        z0 = lfft%zleft + lfft%zoffset
      ELSE
        z0 = 0.0_DP
      END IF
    ELSE
      c   = SQRT(at(1, 3) * at(1, 3) + at(2, 3) * at(2, 3) + at(3, 3) * at(3, 3))
      nrz = dfft%nr3
      dz  = c / DBLE(dfft%nr3)
      IF (.NOT. radial) THEN
        z0 = -0.5_DP * c + dz * MOD(nrz, 2)
      ELSE
        z0 = 0.0_DP
      END IF
    END IF
    !
    ! ... write labels
#if defined (__DEBUG_RISM)
    WRITE(iun, '("#__z_(A) ")', advance='no')
#else
    WRITE(iun, '("#  z (A) ")', advance='no')
#endif
    !
    DO idata = 1, ndata
      label1 = label(idata)
#if defined (__DEBUG_RISM)
      nlabel = LEN_TRIM(label1)
      DO ilabel = 1, nlabel
        IF (label1(ilabel:ilabel) == ' ') THEN
          label1(ilabel:ilabel) = '_'
        END IF
      END DO
#endif
      WRITE(iun, '(2X,A20)', advance='no') label1
    END DO
    WRITE(iun, '()', advance='yes')
    !
    ! ... write data
    DO irz = 1, nrz
      rz = z0 + dz * DBLE(irz - 1)
      rz = rz * alat
      rz = rz * BOHR_RADIUS_ANGS
      WRITE(iun, '(F8.3)', advance='no') rz
      !
      DO idata = 1, ndata
        WRITE(iun, '(E22.12e3)', advance='no') rdata(irz, idata)
      END DO
      WRITE(iun, '()', advance='yes')
    END DO
    !
  END SUBROUTINE write_data
  !
END MODULE solvavg
