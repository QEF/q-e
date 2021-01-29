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
MODULE lauefft
  !--------------------------------------------------------------------------
  !
  ! ... this module perform 1D-FFT to convert between G-space and Laue-rep.,
  ! ...                 and 2D-FFT to convert between R-space and Laue-rep.
  ! ... in Laue-rep., X and Y are in G-space, and Z is in R-space.
  ! ... the unit cell(-z0<z<z0) is expanded along Z, as
  !
  !                                   izleft_end0
  !                                   |
  !                  izleft_start0    |    izright_start0    izright_end0
  !                  |                |    |                 |
  !           zleft  v  -z0           v    v             z0  v      zright
  !     ------|----------|-------------------------------|----------|--------> z
  !                      ^  ^   ^  ^           ^  ^   ^  ^
  !                      |  |   |  |           |  |   |  |
  !                      |  |   |  izleft_end  izright_start
  !                      |  |   |                 |   |  |
  !                      |  |   izleft_gedge      izright_gedge
  !                      |  |                         |  |
  !                      |  izleft_start              izright_end
  !                      |                               |
  !                      izcell_start                    izcell_end
  !
  ! ... , where it has to be izright_end <= izcell_end and izleft_start >= izcell_start.
  ! ... in [izleft_gedge,izleft_end] or [izright_start,izright_gedge], g(z) has to be 0.
  !
  USE control_flags,  ONLY : gamma_only
  USE fft_scalar,     ONLY : cft_1z, cft_2xy
  USE fft_scatter,    ONLY : fft_scatter_xy, fft_scatter_yz
  USE fft_scatter_2d, ONLY : fft_scatter2x1 => fft_scatter
  USE fft_types,      ONLY : fft_type_descriptor
  USE kinds,          ONLY : DP
  USE mp,             ONLY : mp_sum, mp_rank, mp_size
  USE parallel_include
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define variables
  TYPE lauefft_type
    !
    ! ... referenced 3D-FFT
    TYPE(fft_type_descriptor), POINTER :: dfft
    !
    ! ... properties about Z-stick (R-space, corresponding to expanded cell)
    INTEGER  :: nrz             ! 1D-FFT dimension
    INTEGER  :: nrzx            ! 1D-FFT grid leading dimension
    LOGICAL  :: xright          ! expand cell for right(z>z0), or not.
    LOGICAL  :: xleft           ! expand cell for left(z<-z0), or not.
    REAL(DP) :: zstep           ! R-space step of grid (in alat units)
    REAL(DP) :: zoffset         ! R-space offset of grid (in alat units)
    REAL(DP) :: zright          ! cell is expanded to zright(>z0) (in alat units)
    REAL(DP) :: zleft           ! cell is expanded to zleft(<-z0) (in alat units)
    INTEGER  :: izcell_start    ! starting index of unit cell.
    INTEGER  :: izcell_end      ! ending index of unit cell.
    INTEGER  :: izright_start   ! starting index of integral region for right.
    INTEGER  :: izright_end     ! ending index of integral region for right.
    INTEGER  :: izright_start0  ! starting index of integral region for right. (for Gxy = 0)
    INTEGER  :: izright_end0    ! ending index of integral region for right. (for Gxy = 0)
    INTEGER  :: izright_gedge   ! in [izright_start, izright_gedge], g(z) = 0.
    INTEGER  :: izleft_start    ! starting index of integral region for left.
    INTEGER  :: izleft_end      ! ending index of integral region for left.
    INTEGER  :: izleft_start0   ! starting index of integral region for left. (for Gxy = 0)
    INTEGER  :: izleft_end0     ! ending index of integral region for left. (for Gxy = 0)
    INTEGER  :: izleft_gedge    ! in [izleft_gedge, izleft_end], g(z) = 0.
    !
    ! ... properties about Z-stick (G-space, corresponding to unit cell)
    INTEGER              :: ngz           ! number of Gz-vectors
    INTEGER              :: gzzero        ! index, where Gz = 0
    INTEGER,     POINTER :: nlz(:)        ! 1D-FFT index for Gz-vectors
    REAL(DP),    POINTER :: gz(:)         ! Gz-vectors (in units of tpiba=(2pi/a))
    INTEGER,     POINTER :: millz(:)      ! miller index
    INTEGER,     POINTER :: igtoigz(:,:)  ! index of G to index of Gz
    COMPLEX(DP), POINTER :: zphase(:)     ! phase factor of 1D-FFT
    !
    ! ... properties about Z-stick (G-space, corresponding to expanded cell)
    INTEGER              :: ngz_x        ! number of Gz-vectors
    INTEGER              :: gzzero_x     ! index, where Gz = 0
    INTEGER,     POINTER :: nlz_x(:)     ! 1D-FFT index for Gz-vectors
    REAL(DP),    POINTER :: gz_x(:)      ! Gz-vectors (in units of tpiba=(2pi/a))
    INTEGER,     POINTER :: millz_x(:)   ! miller index
    COMPLEX(DP), POINTER :: zphase_x(:)  ! phase factor of 1D-FFT
    !
    ! ... properties about XY-plane (G-space)
    INTEGER           :: ngxy          ! number of Gxy-vectors
    INTEGER           :: ngxy_g        ! number of Gxy-vectors (in global)
    INTEGER           :: nglxy         ! number of Gxy-vector shells
    INTEGER           :: gxystart      ! index of the first Gxy-vectors
    INTEGER,  POINTER :: nlxy(:)       ! 2D-FFT index for Gxy-vectors (0 < Gxy)
    INTEGER,  POINTER :: nlmxy(:)      ! 2D-FFT index for Gxy-vectors (0 > Gxy)
    REAL(DP), POINTER :: gxy(:,:)      ! Gxy-vectors (in units of tpiba=(2pi/a))
    REAL(DP), POINTER :: gnxy(:)       ! |Gxy| (in units of tpiba=(2pi/a))
    REAL(DP), POINTER :: ggxy(:)       ! Gxy^2 (in units of tpiba2=(2pi/a)^2)
    INTEGER,  POINTER :: millxy(:,:)   ! miller index
    REAL(DP), POINTER :: glxy(:)       ! Gxy^2 for each shell (in units of tpiba2=(2pi/a)^2)
    INTEGER,  POINTER :: igtonglxy(:)  ! shell index for Gxy-vectors
    INTEGER,  POINTER :: igtoigxy(:)   ! index of G to index of Gxy
    !
  END TYPE lauefft_type
  !
  ! ... public components
  PUBLIC :: lauefft_type
  PUBLIC :: allocate_lauefft
  PUBLIC :: deallocate_lauefft
  PUBLIC :: set_lauefft_offset
  PUBLIC :: set_lauefft_offset0
  PUBLIC :: set_lauefft_barrier
  PUBLIC :: fw_lauefft_1z
  PUBLIC :: inv_lauefft_1z
  PUBLIC :: fw_lauefft_1z_exp
  PUBLIC :: inv_lauefft_1z_exp
  PUBLIC :: fw_lauefft_2xy
  PUBLIC :: inv_lauefft_2xy
  PUBLIC :: gather_lauefft
  PUBLIC :: gather_lauefft_to_real
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_lauefft(lauefft0, dfft_, dzright, dzleft, ngmt, ig1t, ig2t, ig3t, gt, gcutm, comm)
    !--------------------------------------------------------------------------
    !
    ! ... initialize lauefft_type
    ! ...
    ! ... Variables:
    ! ...   dfft_:   target of referenced 3D-FFT, which must be initialized.
    ! ...   dzright: |zright - z0| (in alat units)
    ! ...   dzleft:  |zleft  + z0| (in alat units)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type),                INTENT(INOUT) :: lauefft0
    TYPE(fft_type_descriptor), TARGET, INTENT(IN)    :: dfft_
    REAL(DP),                          INTENT(IN)    :: dzright
    REAL(DP),                          INTENT(IN)    :: dzleft
    INTEGER,                           INTENT(IN)    :: ngmt
    INTEGER,                           INTENT(IN)    :: ig1t(:)
    INTEGER,                           INTENT(IN)    :: ig2t(:)
    INTEGER,                           INTENT(IN)    :: ig3t(:)
    REAL(DP),                          INTENT(IN)    :: gt(:,:)
    REAL(DP),                          INTENT(IN)    :: gcutm
    INTEGER,                           INTENT(IN)    :: comm
    !
    ! ... set referenced 3D-FFT
    lauefft0%dfft => dfft_
    !
    ! ... set Z-stick (R-space of expanded cell)
    CALL allocate_lauefft_rz(lauefft0, dzright, dzleft)
    !
    ! ... set Z-stick (G-space of unit cell)
    CALL allocate_lauefft_gz(lauefft0, ngmt, ig1t, ig2t, ig3t, gt)
    !
    ! ... set Z-stick (G-space of expanded cell)
    CALL allocate_lauefft_gz_exp(lauefft0, gcutm)
    !
    ! ... set XY-plane (G-space)
    CALL allocate_lauefft_gxy(lauefft0, ngmt, ig1t, ig2t, gt, comm)
    CALL gxyshells(lauefft0, .FALSE.) ! lmovecell must be .FALSE.
    !
  END SUBROUTINE allocate_lauefft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_lauefft(lauefft0)
    !--------------------------------------------------------------------------
    !
    ! ... finalize lauefft_type
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
    !
    ! ... clear referenced 3D-FFT
    IF (ASSOCIATED(lauefft0%dfft)) THEN
      NULLIFY(lauefft0%dfft)
    END IF
    !
    ! ... clear properties about Z-stick (R-space)
    lauefft0%nrz            = 0
    lauefft0%nrzx           = 0
    lauefft0%xright         = .FALSE.
    lauefft0%xleft          = .FALSE.
    lauefft0%zstep          = 0.0_DP
    lauefft0%zoffset        = 0.0_DP
    lauefft0%zright         = 0.0_DP
    lauefft0%zleft          = 0.0_DP
    lauefft0%izcell_start   = 0
    lauefft0%izcell_end     = 0
    lauefft0%izright_start  = 0
    lauefft0%izright_end    = 0
    lauefft0%izright_start0 = 0
    lauefft0%izright_end0   = 0
    lauefft0%izright_gedge  = 0
    lauefft0%izleft_start   = 0
    lauefft0%izleft_end     = 0
    lauefft0%izleft_start0  = 0
    lauefft0%izleft_end0    = 0
    lauefft0%izleft_gedge   = 0
    !
    ! ... clear properties about Z-stick (G-space of unit cell)
    lauefft0%ngz    = 0
    lauefft0%gzzero = 0
    IF (ASSOCIATED(lauefft0%nlz))     DEALLOCATE(lauefft0%nlz)
    IF (ASSOCIATED(lauefft0%gz))      DEALLOCATE(lauefft0%gz)
    IF (ASSOCIATED(lauefft0%millz))   DEALLOCATE(lauefft0%millz)
    IF (ASSOCIATED(lauefft0%igtoigz)) DEALLOCATE(lauefft0%igtoigz)
    IF (ASSOCIATED(lauefft0%zphase))  DEALLOCATE(lauefft0%zphase)
    !
    ! ... clear properties about Z-stick (G-space of expanded cell)
    lauefft0%ngz_x    = 0
    lauefft0%gzzero_x = 0
    IF (ASSOCIATED(lauefft0%nlz_x))    DEALLOCATE(lauefft0%nlz_x)
    IF (ASSOCIATED(lauefft0%gz_x))     DEALLOCATE(lauefft0%gz_x)
    IF (ASSOCIATED(lauefft0%millz_x))  DEALLOCATE(lauefft0%millz_x)
    IF (ASSOCIATED(lauefft0%zphase_x)) DEALLOCATE(lauefft0%zphase_x)
    !
    ! ... clear properties about XY-plane (G-space)
    lauefft0%ngxy     = 0
    lauefft0%ngxy_g   = 0
    lauefft0%nglxy    = 0
    lauefft0%gxystart = 0
    IF (ASSOCIATED(lauefft0%nlxy))      DEALLOCATE(lauefft0%nlxy)
    IF (ASSOCIATED(lauefft0%nlmxy))     DEALLOCATE(lauefft0%nlmxy)
    IF (ASSOCIATED(lauefft0%gxy))       DEALLOCATE(lauefft0%gxy)
    IF (ASSOCIATED(lauefft0%gnxy))      DEALLOCATE(lauefft0%gnxy)
    IF (ASSOCIATED(lauefft0%ggxy))      DEALLOCATE(lauefft0%ggxy)
    IF (ASSOCIATED(lauefft0%millxy))    DEALLOCATE(lauefft0%millxy)
    IF (ASSOCIATED(lauefft0%glxy))      DEALLOCATE(lauefft0%glxy)
    IF (ASSOCIATED(lauefft0%igtonglxy)) DEALLOCATE(lauefft0%igtonglxy)
    IF (ASSOCIATED(lauefft0%igtoigxy))  DEALLOCATE(lauefft0%igtoigxy)
    !
  END SUBROUTINE deallocate_lauefft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_lauefft_offset(lauefft0, wright, wleft)
    !--------------------------------------------------------------------------
    !
    ! ... set offsets from z=0.
    ! ...
    ! ... Variables:
    ! ...   wright: offset on right-hand side (in alat units)
    ! ...   wleft:  offset on left-hand side (in alat units)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
    REAL(DP),           INTENT(IN)    :: wright
    REAL(DP),           INTENT(IN)    :: wleft
    !
    CALL set_lauefft_offset_x(lauefft0, wright, wleft)
    !
  END SUBROUTINE set_lauefft_offset
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_lauefft_offset0(lauefft0, wright1, wright2, wleft1, wleft2)
    !--------------------------------------------------------------------------
    !
    ! ... set offsets from z=0, for Gxy = 0.
    ! ...
    ! ... Variables:
    ! ...   wright1: offset of solute-side on right-hand side (in alat units)
    ! ...   wright2: offset of solvent-side on right-hand side (in alat units)
    ! ...   wleft1:  offset of solute-side on left-hand side (in alat units)
    ! ...   wleft2:  offset of solvent-side on left-hand side (in alat units)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
    REAL(DP),           INTENT(IN)    :: wright1
    REAL(DP),           INTENT(IN)    :: wright2
    REAL(DP),           INTENT(IN)    :: wleft1
    REAL(DP),           INTENT(IN)    :: wleft2
    !
    CALL set_lauefft_offset0_x(lauefft0, wright1, wright2, wleft1, wleft2)
    !
  END SUBROUTINE set_lauefft_offset0
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_lauefft_barrier(lauefft0, wright, wleft)
    !--------------------------------------------------------------------------
    !
    ! ... set offset of barrier, where g(z) = 0.
    ! ...
    ! ... Variables:
    ! ...   wright: offset of barrier on right-hand side (in alat units)
    ! ...   wleft:  offset of barrier on left-hand side (in alat units)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(INOUT) :: lauefft0
    REAL(DP),           INTENT(IN)    :: wright
    REAL(DP),           INTENT(IN)    :: wleft
    !
    CALL set_lauefft_barrier_x(lauefft0, wright, wleft)
    !
  END SUBROUTINE set_lauefft_barrier
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fw_lauefft_1z(lauefft0, cl, nrz, irz_start, cg)
    !--------------------------------------------------------------------------
    !
    ! ... FFT Rz,Gy,Gx -> Gz,Gy,Gx (in unit cell)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    COMPLEX(DP),        INTENT(IN)  :: cl(1:*)    ! Laue-rep., dimension(nrz*ngxy)
    INTEGER,            INTENT(IN)  :: nrz        ! leading dimension of Z
    INTEGER,            INTENT(IN)  :: irz_start  ! starting index of Z
    COMPLEX(DP),        INTENT(OUT) :: cg(1:*)    ! G-space,   dimension(nnr)
    !
    INTEGER                  :: irz
    INTEGER                  :: jrz1
    INTEGER                  :: jrz2
    INTEGER                  :: igz
    INTEGER                  :: jgz1
    INTEGER                  :: jgz2
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: nx1
    INTEGER                  :: nx2
    INTEGER                  :: nx3, n3
    COMPLEX(DP), ALLOCATABLE :: cinp(:)
    COMPLEX(DP), ALLOCATABLE :: cout(:)
    !
    n3  = lauefft0%dfft%nr3
    nx1 = lauefft0%dfft%nr1x
    nx2 = lauefft0%dfft%nr2x
    nx3 = lauefft0%dfft%nr3x
    !
    ALLOCATE(cinp(nx3 * lauefft0%ngxy))
    ALLOCATE(cout(nx3 * lauefft0%ngxy))
    !
    cinp = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * nrz
      jgxy2 = (igxy - 1) * nx3
!$omp parallel do default(shared) private(irz, jrz1, jrz2)
      DO irz = 1, n3
        jrz1 = irz + irz_start - 1
        IF (irz <= (n3 / 2)) THEN
          jrz2 = irz + (n3 - (n3 / 2))
        ELSE
          jrz2 = irz - (n3 / 2)
        END IF
        cinp(jrz2 + jgxy2) = cl(jrz1 + jgxy1)
      END DO
!$omp end parallel do
    END DO
    !
    CALL cft_1z(cinp, lauefft0%ngxy, n3, nx3, -1, cout)
    !
    cg(1:lauefft0%dfft%nnr) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = lauefft0%nlxy(igxy)
      jgxy2 = (igxy - 1) * nx3
      IF (lauefft0%dfft%lpara) THEN
!$omp parallel do default(shared) private(igz)
        DO igz = 1, lauefft0%ngz
          cg(lauefft0%nlz(igz) + jgxy1) = &
          & cout(lauefft0%nlz(igz) + jgxy2) * lauefft0%zphase(igz)
        END DO
!$omp end parallel do
      ELSE
!$omp parallel do default(shared) private(igz)
        DO igz = 1, lauefft0%ngz
          cg(nx1 * nx2 * (lauefft0%nlz(igz) - 1) + jgxy1) = &
          & cout(lauefft0%nlz(igz) + jgxy2) * lauefft0%zphase(igz)
        END DO
!$omp end parallel do
      END IF
    END DO
    !
    IF (gamma_only) THEN
      DO igxy = lauefft0%gxystart, lauefft0%ngxy
        jgxy1 = lauefft0%nlxy( igxy)
        jgxy2 = lauefft0%nlmxy(igxy)
        IF (lauefft0%dfft%lpara) THEN
!$omp parallel do default(shared) private(igz, jgz1, jgz2)
          DO igz = 1, lauefft0%ngz
            jgz1 = lauefft0%nlz(igz)
            jgz2 = lauefft0%nlz(lauefft0%ngz - igz + 1)
            cg(jgz2 + jgxy2) = CONJG(cg(jgz1 + jgxy1))
          END DO
!$omp end parallel do
        ELSE
!$omp parallel do default(shared) private(igz, jgz1, jgz2)
          DO igz = 1, lauefft0%ngz
            jgz1 = nx1 * nx2 * (lauefft0%nlz(igz) - 1)
            jgz2 = nx1 * nx2 * (lauefft0%nlz(lauefft0%ngz - igz + 1) - 1)
            cg(jgz2 + jgxy2) = CONJG(cg(jgz1 + jgxy1))
          END DO
!$omp end parallel do
        END IF
      END DO
    END IF
    !
    DEALLOCATE(cinp)
    DEALLOCATE(cout)
    !
  END SUBROUTINE fw_lauefft_1z
  !
  !--------------------------------------------------------------------------
  SUBROUTINE inv_lauefft_1z(lauefft0, cg, cl, nrz, irz_start)
    !--------------------------------------------------------------------------
    !
    ! ... FFT Gz,Gy,Gx -> Rz,Gy,Gx (in unit cell)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    COMPLEX(DP),        INTENT(IN)  :: cg(1:*)    ! G-space,   dimension(nnr)
    COMPLEX(DP),        INTENT(OUT) :: cl(1:*)    ! Laue-rep., dimension(nrz*ngxy)
    INTEGER,            INTENT(IN)  :: nrz        ! leading dimension of Z
    INTEGER,            INTENT(IN)  :: irz_start  ! starting index of Z
    !
    INTEGER                  :: irz
    INTEGER                  :: jrz1
    INTEGER                  :: jrz2
    INTEGER                  :: igz
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: nx1
    INTEGER                  :: nx2
    INTEGER                  :: nx3, n3
    COMPLEX(DP), ALLOCATABLE :: cinp(:)
    COMPLEX(DP), ALLOCATABLE :: cout(:)
    !
    n3  = lauefft0%dfft%nr3
    nx1 = lauefft0%dfft%nr1x
    nx2 = lauefft0%dfft%nr2x
    nx3 = lauefft0%dfft%nr3x
    !
    ALLOCATE(cinp(nx3 * lauefft0%ngxy))
    ALLOCATE(cout(nx3 * lauefft0%ngxy))
    !
    cinp = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = lauefft0%nlxy(igxy)
      jgxy2 = (igxy - 1) * nx3
      IF (lauefft0%dfft%lpara) THEN
!$omp parallel do default(shared) private(igz)
        DO igz = 1, lauefft0%ngz
          cinp(lauefft0%nlz(igz) + jgxy2) = &
          & cg(lauefft0%nlz(igz) + jgxy1) * CONJG(lauefft0%zphase(igz))
        END DO
!$omp end parallel do
      ELSE
!$omp parallel do default(shared) private(igz)
        DO igz = 1, lauefft0%ngz
          cinp(lauefft0%nlz(igz) + jgxy2) = &
          & cg(nx1 * nx2 * (lauefft0%nlz(igz) - 1) + jgxy1) * CONJG(lauefft0%zphase(igz))
        END DO
!$omp end parallel do
      END IF
    END DO
    !
    CALL cft_1z(cinp, lauefft0%ngxy, n3, nx3, +1, cout)
    !
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * nrz
      jgxy2 = (igxy - 1) * nx3
!$omp parallel do default(shared) private(irz, jrz1, jrz2)
      DO irz = 1, n3
        jrz1 = irz + irz_start - 1
        IF (irz <= (n3 / 2)) THEN
          jrz2 = irz + (n3 - (n3 / 2))
        ELSE
          jrz2 = irz - (n3 / 2)
        END IF
        cl(jrz1 + jgxy1) = cout(jrz2 + jgxy2)
      END DO
!$omp end parallel do
    END DO
    !
    DEALLOCATE(cinp)
    DEALLOCATE(cout)
    !
  END SUBROUTINE inv_lauefft_1z
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fw_lauefft_1z_exp(lauefft0, cl, nrz, cg, ngz)
    !--------------------------------------------------------------------------
    !
    ! ... FFT Rz,Gy,Gx -> Gz,Gy,Gx (in expanded cell)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    COMPLEX(DP),        INTENT(IN)  :: cl(1:*)    ! Laue-rep., dimension(nrz*ngxy)
    INTEGER,            INTENT(IN)  :: nrz        ! leading dimension of Z
    COMPLEX(DP),        INTENT(OUT) :: cg(1:*)    ! G-space,   dimension(ngz*ngxy)
    INTEGER,            INTENT(IN)  :: ngz        ! leading dimension of Gz
    !
    INTEGER                  :: irz
    INTEGER                  :: irz0
    INTEGER                  :: jrz1
    INTEGER                  :: jrz2
    INTEGER                  :: igz
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: n3, nx3
    COMPLEX(DP), ALLOCATABLE :: cinp(:)
    COMPLEX(DP), ALLOCATABLE :: cout(:)
    !
    n3   = lauefft0%nrz
    nx3  = lauefft0%nrzx
    irz0 = (lauefft0%dfft%nr3 / 2) + lauefft0%izcell_start - 1
    !
    ALLOCATE(cinp(nx3 * lauefft0%ngxy))
    ALLOCATE(cout(nx3 * lauefft0%ngxy))
    !
    cinp = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * nrz
      jgxy2 = (igxy - 1) * nx3
!$omp parallel do default(shared) private(irz, jrz1, jrz2)
      DO irz = 1, n3
        jrz1 = irz
        IF (irz <= irz0) THEN
          jrz2 = irz + (n3 - irz0)
        ELSE
          jrz2 = irz - irz0
        END IF
        cinp(jrz2 + jgxy2) = cl(jrz1 + jgxy1)
      END DO
!$omp end parallel do
    END DO
    !
    CALL cft_1z(cinp, lauefft0%ngxy, n3, nx3, -1, cout)
    !
    cg(1:(ngz * lauefft0%ngxy)) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * ngz
      jgxy2 = (igxy - 1) * nx3
!$omp parallel do default(shared) private(igz)
      DO igz = 1, lauefft0%ngz_x
        cg(igz + jgxy1) = &
        & cout(lauefft0%nlz_x(igz) + jgxy2) * lauefft0%zphase_x(igz)
      END DO
!$omp end parallel do
    END DO
    !
    DEALLOCATE(cinp)
    DEALLOCATE(cout)
    !
  END SUBROUTINE fw_lauefft_1z_exp
  !
  !--------------------------------------------------------------------------
  SUBROUTINE inv_lauefft_1z_exp(lauefft0, cg, ngz, cl, nrz)
    !--------------------------------------------------------------------------
    !
    ! ... FFT Gz,Gy,Gx -> Rz,Gy,Gx (in expanded cell)
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    COMPLEX(DP),        INTENT(IN)  :: cg(1:*)    ! G-space,   dimension(ngz*ngxy)
    INTEGER,            INTENT(IN)  :: ngz        ! leading dimension of Gz
    COMPLEX(DP),        INTENT(OUT) :: cl(1:*)    ! Laue-rep., dimension(nrz*ngxy)
    INTEGER,            INTENT(IN)  :: nrz        ! leading dimension of Z
    !
    INTEGER                  :: irz
    INTEGER                  :: irz0
    INTEGER                  :: jrz1
    INTEGER                  :: jrz2
    INTEGER                  :: igz
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: n3, nx3
    COMPLEX(DP), ALLOCATABLE :: cinp(:)
    COMPLEX(DP), ALLOCATABLE :: cout(:)
    !
    n3   = lauefft0%nrz
    nx3  = lauefft0%nrzx
    irz0 = (lauefft0%dfft%nr3 / 2) + lauefft0%izcell_start - 1
    !
    ALLOCATE(cinp(nx3 * lauefft0%ngxy))
    ALLOCATE(cout(nx3 * lauefft0%ngxy))
    !
    cinp = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * ngz
      jgxy2 = (igxy - 1) * nx3
!$omp parallel do default(shared) private(igz)
      DO igz = 1, lauefft0%ngz_x
        cinp(lauefft0%nlz_x(igz) + jgxy2) = &
        & cg(igz + jgxy1) * CONJG(lauefft0%zphase_x(igz))
      END DO
!$omp end parallel do
    END DO
    !
    CALL cft_1z(cinp, lauefft0%ngxy, n3, nx3, +1, cout)
    !
    cl(1:(nrz * lauefft0%ngxy)) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * nrz
      jgxy2 = (igxy - 1) * nx3
!$omp parallel do default(shared) private(irz, jrz1, jrz2)
      DO irz = 1, n3
        jrz1 = irz
        IF (irz <= irz0) THEN
          jrz2 = irz + (n3 - irz0)
        ELSE
          jrz2 = irz - irz0
        END IF
        cl(jrz1 + jgxy1) = cout(jrz2 + jgxy2)
      END DO
!$omp end parallel do
    END DO
    !
    DEALLOCATE(cinp)
    DEALLOCATE(cout)
    !
  END SUBROUTINE inv_lauefft_1z_exp
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fw_lauefft_2xy(lauefft0, cr, cl, nrz, irz_start, i3mask)
    !--------------------------------------------------------------------------
    !
    ! ... FFT Rx,Ry,Rz -> Rz,Gy,Gx
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    REAL(DP),           INTENT(IN)  :: cr(1:*)      ! R-space,   dimension(nnr)
    COMPLEX(DP),        INTENT(OUT) :: cl(1:*)      ! Laue-rep., dimension(nrz*ngxy)
    INTEGER,            INTENT(IN)  :: nrz          ! leading dimension of Z
    INTEGER,            INTENT(IN)  :: irz_start    ! starting index of Z
    LOGICAL, OPTIONAL,  INTENT(IN)  :: i3mask(1:*)  ! mask of i3-axis, dimension(n3)
    !
    INTEGER                  :: ir
    INTEGER                  :: irz
    INTEGER                  :: jrz1
    INTEGER                  :: jrz2
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: n1, nx1, np1
    INTEGER                  :: n2, nx2, np2
    INTEGER                  :: n3, nx3, np3
    INTEGER                  :: i3
    INTEGER                  :: i3min, i3max
    INTEGER                  :: i3sta, i3end
    INTEGER                  :: j3sta, j3end
    LOGICAL                  :: do_fft
    COMPLEX(DP), ALLOCATABLE :: cinp(:)
    COMPLEX(DP), ALLOCATABLE :: cout(:)
    !
    n1  = lauefft0%dfft%nr1
    n2  = lauefft0%dfft%nr2
    n3  = lauefft0%dfft%nr3
    nx1 = lauefft0%dfft%nr1x
    nx2 = lauefft0%dfft%nr2x
    nx3 = lauefft0%dfft%nr3x
    np1 = lauefft0%dfft%nr1p(lauefft0%dfft%mype2 + 1)
    np2 = lauefft0%dfft%my_nr2p
    np3 = lauefft0%dfft%my_nr3p
    !
    ALLOCATE(cinp(lauefft0%dfft%nnr))
    ALLOCATE(cout(lauefft0%dfft%nnr))
    !
!$omp parallel do default(shared) private(ir)
    DO ir = 1, lauefft0%dfft%nnr
      cinp(ir) = CMPLX(cr(ir), 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
    IF (np2 == nx2) THEN
      !
      ! ... Ry-axis is NOT distributed
      !
      IF (PRESENT(i3mask)) THEN
        !
        ! ... 2D-FFT for specified planes
        !
        i3min = lauefft0%dfft%my_i0r3p
        i3max = MIN(np3 + i3min, n3)
        i3sta = i3min
        i3end = i3min
        !
        DO i3 = (i3min + 1), i3max
          IF (i3mask(i3)) THEN
            i3sta = i3
          ELSE
            i3end = i3
          END IF
          !
          do_fft = (.NOT. i3mask(i3))
          IF (i3 < i3max) THEN
            do_fft = do_fft .AND. i3mask(i3 + 1)
          END IF
          !
          IF (do_fft .AND. i3sta < i3end) THEN
            !
            j3sta = nx1 * nx2 * (i3sta - i3min) + 1
            j3end = nx1 * nx2 * (i3end - i3min)
            !
            ! Rx, Ry, Rz
            CALL cft_2xy(cinp(j3sta:j3end), i3end - i3sta, &
                       & n1, n2, nx1, nx2, -1, lauefft0%dfft%iplp)
            ! Gx, Gy, Rz
            !
          END IF
        END DO
        !
      ELSE
        !
        ! ... 2D-FFT for all planes
        !
        ! Rx, Ry, Rz
        CALL cft_2xy(cinp, np3, n1, n2, nx1, nx2, -1, lauefft0%dfft%iplp)
        ! Gx, Gy, Rz
        !
      END IF
      !
      IF (lauefft0%dfft%lpara .AND. lauefft0%dfft%use_pencil_decomposition) THEN
        !
        ! Gx, Gy, Rz
        CALL fft_scatter_xy(lauefft0%dfft, cout, cinp, lauefft0%dfft%nnr, -1)
        ! Gy, Gx, Rz
        CALL fft_scatter_yz(lauefft0%dfft, cinp, cout, lauefft0%dfft%nnr, -1)
        ! Rz, Gy, Gx
        !
      ELSE IF (lauefft0%dfft%lpara) THEN
        !
        ! Gx, Gy, Rz
        CALL fft_scatter2x1(lauefft0%dfft, &
        & cout, nx3, lauefft0%dfft%nnr, cinp, lauefft0%dfft%nsp, lauefft0%dfft%nr3p, -1)
        ! Rz, Gx, Gy
        !
      END IF
      !
    ELSE
      !
      ! ... Ry-axis is distributed
      !
      IF (.NOT. lauefft0%dfft%lpara) THEN
        CALL errore('fw_lauefft_2xy', 'my_nr2p != nr2x, but not parallel', 1)
      END IF
      !
      IF (.NOT. lauefft0%dfft%use_pencil_decomposition) THEN
        CALL errore('fw_lauefft_2xy', 'my_nr2p != nr2x, but not pencil-decomposed', 1)
      END IF
      !
      ! Rx, Ry, Rz
      CALL cft_1z(cinp, np2 * np3, n1, nx1, -1, cout)
      ! Gx, Ry, Rz
      CALL fft_scatter_xy(lauefft0%dfft, cinp, cout, lauefft0%dfft%nnr, -1)
      ! Ry, Gx, Rz
      CALL cft_1z(cinp, np1 * np3, n2, nx2, -1, cout)
      ! Gy, Gx, Rz
      CALL fft_scatter_yz(lauefft0%dfft, cinp, cout, lauefft0%dfft%nnr, -1)
      ! Rz, Gy, Gx
      !
    END IF
    !
    cout = cinp
    !
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * nrz
      jgxy2 = lauefft0%nlxy(igxy)
!$omp parallel do default(shared) private(irz, jrz1, jrz2)
      DO irz = 1, n3
        jrz1 = irz + irz_start - 1
        IF (irz <= (n3 / 2)) THEN
          jrz2 = irz + (n3 - (n3 / 2))
        ELSE
          jrz2 = irz - (n3 / 2)
        END IF
        IF (lauefft0%dfft%lpara) THEN
          cl(jrz1 + jgxy1) = cout(jrz2 + jgxy2)
        ELSE
          cl(jrz1 + jgxy1) = cout(nx1 * nx2 * (jrz2 - 1) + jgxy2)
        END IF
      END DO
!$omp end parallel do
    END DO
    !
    DEALLOCATE(cinp)
    DEALLOCATE(cout)
    !
  END SUBROUTINE fw_lauefft_2xy
  !
  !--------------------------------------------------------------------------
  SUBROUTINE inv_lauefft_2xy(lauefft0, cl, nrz, irz_start, cr, i3mask)
    !--------------------------------------------------------------------------
    !
    ! ... FFT Rz,Gy,Gx -> Rx,Ry,Rz
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    COMPLEX(DP),        INTENT(IN)  :: cl(1:*)      ! Laue-rep., dimension(nrz*ngxy)
    INTEGER,            INTENT(IN)  :: nrz          ! leading dimension of Z
    INTEGER,            INTENT(IN)  :: irz_start    ! starting index of Z
    REAL(DP),           INTENT(OUT) :: cr(1:*)      ! R-space,   dimension(nnr)
    LOGICAL, OPTIONAL,  INTENT(IN)  :: i3mask(1:*)  ! mask of i3-axis, dimension(n3)
    !
    INTEGER                  :: ir
    INTEGER                  :: irz
    INTEGER                  :: jrz1
    INTEGER                  :: jrz2
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: n1, nx1, np1
    INTEGER                  :: n2, nx2, np2
    INTEGER                  :: n3, nx3, np3
    INTEGER                  :: i3
    INTEGER                  :: i3min, i3max
    INTEGER                  :: i3sta, i3end
    INTEGER                  :: j3sta, j3end
    LOGICAL                  :: do_fft
    COMPLEX(DP), ALLOCATABLE :: cinp(:)
    COMPLEX(DP), ALLOCATABLE :: cout(:)
    !
    n1  = lauefft0%dfft%nr1
    n2  = lauefft0%dfft%nr2
    n3  = lauefft0%dfft%nr3
    nx1 = lauefft0%dfft%nr1x
    nx2 = lauefft0%dfft%nr2x
    nx3 = lauefft0%dfft%nr3x
    np1 = lauefft0%dfft%nr1p(lauefft0%dfft%mype2 + 1)
    np2 = lauefft0%dfft%my_nr2p
    np3 = lauefft0%dfft%my_nr3p
    !
    ALLOCATE(cinp(lauefft0%dfft%nnr))
    ALLOCATE(cout(lauefft0%dfft%nnr))
    !
    cinp = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    DO igxy = 1, lauefft0%ngxy
      jgxy1 = (igxy - 1) * nrz
      jgxy2 = lauefft0%nlxy(igxy)
!$omp parallel do default(shared) private(irz, jrz1, jrz2)
      DO irz = 1, n3
        jrz1 = irz + irz_start - 1
        IF (irz <= (n3 / 2)) THEN
          jrz2 = irz + (n3 - (n3 / 2))
        ELSE
          jrz2 = irz - (n3 / 2)
        END IF
        IF (lauefft0%dfft%lpara) THEN
          cinp(jrz2 + jgxy2) = cl(jrz1 + jgxy1)
        ELSE
          cinp(nx1 * nx2 * (jrz2 - 1) + jgxy2) = cl(jrz1 + jgxy1)
        END IF
      END DO
!$omp end parallel do
    END DO
    !
    IF (gamma_only) THEN
      DO igxy = lauefft0%gxystart, lauefft0%ngxy
        jgxy1 = lauefft0%nlxy( igxy)
        jgxy2 = lauefft0%nlmxy(igxy)
        IF (lauefft0%dfft%lpara) THEN
!$omp parallel do default(shared) private(irz)
          DO irz = 1, n3
            cinp(irz + jgxy2) = CONJG(cinp(irz + jgxy1))
          END DO
!$omp end parallel do
        ELSE
!$omp parallel do default(shared) private(irz)
          DO irz = 1, n3
            cinp(nx1 * nx2 * (irz - 1) + jgxy2) = CONJG(cinp(nx1 * nx2 * (irz - 1) + jgxy1))
          END DO
!$omp end parallel do
        END IF
      END DO
    END IF
    !
    cout = cinp
    !
    IF (np2 == nx2) THEN
      !
      ! ... Ry-axis is NOT distributed
      !
      IF (lauefft0%dfft%lpara .AND. lauefft0%dfft%use_pencil_decomposition) THEN
        !
        ! Rz, Gy, Gx
        CALL fft_scatter_yz(lauefft0%dfft, cout, cinp, lauefft0%dfft%nnr, +1)
        ! Gy, Gx, Rz
        CALL fft_scatter_xy(lauefft0%dfft, cinp, cout, lauefft0%dfft%nnr, +1)
        ! Gx, Gy, Rz
        !
      ELSE IF (lauefft0%dfft%lpara) THEN
        !
        ! Rz, Gy, Gx
        CALL fft_scatter2x1(lauefft0%dfft, &
        & cinp, nx3, lauefft0%dfft%nnr, cout, lauefft0%dfft%nsp, lauefft0%dfft%nr3p, +1)
        ! Gy, Gx, Rz
        !
      END IF
      !
      IF (PRESENT(i3mask)) THEN
        !
        ! ... 2D-FFT for specified planes
        !
        i3min = lauefft0%dfft%my_i0r3p
        i3max = MIN(np3 + i3min, n3)
        i3sta = i3min
        i3end = i3min
        !
        DO i3 = (i3min + 1), i3max
          IF (i3mask(i3)) THEN
            i3sta = i3
          ELSE
            i3end = i3
          END IF
          !
          do_fft = (.NOT. i3mask(i3))
          IF (i3 < i3max) THEN
            do_fft = do_fft .AND. i3mask(i3 + 1)
          END IF
          !
          IF (do_fft .AND. i3sta < i3end) THEN
            !
            j3sta = nx1 * nx2 * (i3sta - i3min) + 1
            j3end = nx1 * nx2 * (i3end - i3min)
            !
            ! Gx, Gy, Rz
            CALL cft_2xy(cout(j3sta:j3end), i3end - i3sta, &
                       & n1, n2, nx1, nx2, +1, lauefft0%dfft%iplp)
            ! Rx, Ry, Rz
            !
          END IF
        END DO
        !
      ELSE
        !
        ! ... 2D-FFT for all planes
        !
        ! Gx, Gy, Rz
        CALL cft_2xy(cout, np3, n1, n2, nx1, nx2, +1, lauefft0%dfft%iplp)
        ! Rx, Ry, Rz
        !
      END IF
      !
    ELSE
      !
      ! ... Ry-axis is distributed
      !
      IF (.NOT. lauefft0%dfft%lpara) THEN
        CALL errore('inv_lauefft_2xy', 'my_nr2p != nr2x, but not parallel', 1)
      END IF
      !
      IF (.NOT. lauefft0%dfft%use_pencil_decomposition) THEN
        CALL errore('inv_lauefft_2xy', 'my_nr2p != nr2x, but not pencil-decomposed', 1)
      END IF
      !
      ! Rz, Gy, Gx
      CALL fft_scatter_yz(lauefft0%dfft, cout, cinp, lauefft0%dfft%nnr, +1)
      ! Gy, Gx, Rz
      CALL cft_1z(cinp, np1 * np3, n2, nx2, +1, cout)
      ! Ry, Gx, Rz
      CALL fft_scatter_xy(lauefft0%dfft, cout, cinp, lauefft0%dfft%nnr, +1)
      ! Gx, Ry, Rz
      CALL cft_1z(cinp, np2 * np3, n1, nx1, +1, cout)
      ! Rx, Ry, Rz
      !
    END IF
    !
!$omp parallel do default(shared) private(ir)
    DO ir = 1, lauefft0%dfft%nnr
      cr(ir) = DBLE(cout(ir))
    END DO
!$omp end parallel do
    !
    DEALLOCATE(cinp)
    DEALLOCATE(cout)
    !
  END SUBROUTINE inv_lauefft_2xy
  !
  !--------------------------------------------------------------------------
  SUBROUTINE gather_lauefft(lauefft0, cl, nrz, cltot, lall)
    !--------------------------------------------------------------------------
    !
    ! ... gathers a distributed Laue-rep. FFT grid to lauefft0%dfft%root, that is,
    ! ... the first processor of input descriptor lauefft0%dfft
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    COMPLEX(DP),        INTENT(IN)  :: cl(1:*)     ! distributed variable (nrz*ngxy)  ; stick major
    INTEGER,            INTENT(IN)  :: nrz         ! leading dimension of Z
    COMPLEX(DP),        INTENT(OUT) :: cltot(1:*)  ! gathered variable (nr1x*nr2x*nrz); plane major
    LOGICAL, OPTIONAL,  INTENT(IN)  :: lall        ! bcast to all ranks, or not
    !
    LOGICAL                  :: lall_
    INTEGER                  :: ierr
    INTEGER                  :: isign
    INTEGER                  :: irz
    INTEGER                  :: igxy
    INTEGER                  :: jgxy1
    INTEGER                  :: jgxy2
    INTEGER                  :: igx, igy
    INTEGER                  :: mx, my
    INTEGER                  :: n1, nx1
    INTEGER                  :: n2, nx2
    INTEGER                  :: n3
    INTEGER                  :: ntot
    REAL(DP)                 :: creal
    REAL(DP)                 :: cimag
    COMPLEX(DP), ALLOCATABLE :: cltmp(:)
    !
    lall_ = .FALSE.
    IF (PRESENT(lall)) THEN
      lall_ = lall
    END IF
    !
    n1   = lauefft0%dfft%nr1
    n2   = lauefft0%dfft%nr2
    n3   = lauefft0%nrz
    nx1  = lauefft0%dfft%nr1x
    nx2  = lauefft0%dfft%nr2x
    ntot = nx1 * nx2 * n3
    !
    ALLOCATE(cltmp(ntot))
    cltmp(1:ntot) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    !
    DO igxy = 1, lauefft0%ngxy
      isign = 1
   10 CONTINUE
      !
      mx  = isign * lauefft0%millxy(1, igxy)
      igx = mx + 1
      IF (igx < 1) THEN
        igx = igx + n1
      END IF
      !
      my  = isign * lauefft0%millxy(2, igxy)
      igy = my + 1
      IF (igy < 1) THEN
        igy = igy + n2
      END IF
      !
      jgxy1 = nrz * (igxy - 1)
      jgxy2 = igx + (igy - 1) * nx1
!$omp parallel do default(shared) private(irz, creal, cimag)
      DO irz = 1, n3
        creal = DBLE( cl(irz + jgxy1))
        cimag = AIMAG(cl(irz + jgxy1))
        cltmp(nx1 * nx2 * (irz - 1) + jgxy2) = CMPLX(creal, DBLE(isign) * cimag, kind=DP)
      END DO
!$omp end parallel do
      !
      IF (gamma_only .AND. isign > 0 .AND. igxy >= lauefft0%gxystart) THEN
        isign = -1
        GOTO 10
      END IF
    END DO
    !
#if defined (__MPI)
    IF (lall_) THEN
      CALL mp_sum(cltmp, lauefft0%dfft%comm)
      cltot(1:ntot) = cltmp(1:ntot)
      !
    ELSE
      CALL MPI_REDUCE(cltmp(1), cltot(1), ntot, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                    & lauefft0%dfft%root, lauefft0%dfft%comm, ierr)
      IF (ierr /= MPI_SUCCESS) THEN
        CALL errore('gather_lauefft', 'error at MPI_REDUCE', 1)
      END IF
    END IF
    !
#else
    cltot(1:ntot) = cltmp(1:ntot)
    !
#endif
    DEALLOCATE(cltmp)
    !
  END SUBROUTINE gather_lauefft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE gather_lauefft_to_real(lauefft0, cl, nrz, crtot)
    !--------------------------------------------------------------------------
    !
    ! ... gathers a distributed Laue-rep. FFT grid to lauefft0%dfft%root, that is,
    ! ... the first processor of input descriptor lauefft0%dfft, as R-space
    !                                                               ^^^^^^^
    !
    IMPLICIT NONE
    !
    TYPE(lauefft_type), INTENT(IN)  :: lauefft0
    COMPLEX(DP),        INTENT(IN)  :: cl(1:*)     ! distributed variable (nrz*ngxy)  ; stick major
    INTEGER,            INTENT(IN)  :: nrz         ! leading dimension of Z
    REAL(DP),           INTENT(OUT) :: crtot(1:*)  ! gathered variable (nr1x*nr2x*nrz); plane major
    !
    INTEGER                  :: ierr
    INTEGER                  :: irz
    INTEGER                  :: nsta
    INTEGER                  :: nend
    INTEGER                  :: n1, nx1
    INTEGER                  :: n2, nx2
    INTEGER                  :: n3
    INTEGER                  :: ntot
    INTEGER                  :: irank
    INTEGER                  :: nproc
    REAL(DP),    ALLOCATABLE :: crtmp(:)
    COMPLEX(DP), ALLOCATABLE :: cltot(:)
    !
    n1   = lauefft0%dfft%nr1
    n2   = lauefft0%dfft%nr2
    n3   = lauefft0%nrz
    nx1  = lauefft0%dfft%nr1x
    nx2  = lauefft0%dfft%nr2x
    ntot = nx1 * nx2 * n3
    !
    ALLOCATE(crtmp(ntot))
    ALLOCATE(cltot(ntot))
    !
    crtmp = 0.0_DP
    CALL gather_lauefft(lauefft0, cl, nrz, cltot, lall=.TRUE.)
    !
    irank = mp_rank(lauefft0%dfft%comm)
    nproc = mp_size(lauefft0%dfft%comm)
    !
    DO irz = 1, n3
      IF (irank /= MOD(irz - 1, nproc)) THEN
        CYCLE
      END IF
      !
      nsta = nx1 * nx2 * (irz - 1) + 1
      nend = nx1 * nx2 * irz
      CALL cft_2xy(cltot(nsta:nend), 1, n1, n2, nx1, nx2, +1)
      crtmp(nsta:nend) = DBLE(cltot(nsta:nend))
    END DO
    !
#if defined (__MPI)
    CALL MPI_REDUCE(crtmp(1), crtot(1), ntot, MPI_DOUBLE_PRECISION, MPI_SUM, &
                  & lauefft0%dfft%root, lauefft0%dfft%comm, ierr)
    IF (ierr /= MPI_SUCCESS) THEN
      CALL errore('gather_lauefft_to_real', 'error at MPI_REDUCE', 1)
    END IF
    !
#else
    crtot(1:ntot) = crtmp(1:ntot)
    !
#endif
    DEALLOCATE(crtmp)
    DEALLOCATE(cltot)
    !
  END SUBROUTINE gather_lauefft_to_real
  !
END MODULE lauefft
