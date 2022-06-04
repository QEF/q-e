!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE radfft
  !--------------------------------------------------------------------------
  !
  ! ... this module perform radial FFT, which is defined as
  ! ...                       / inf
  ! ...   g * A(g) = 4 * pi * | dr sin(g * r) * r * A(r) ,
  ! ...                       / 0
  ! ... and
  ! ...                  1        / inf
  ! ...   r * A(r) = ---------- * | dg sin(g * r) * g * A(g) .
  ! ...               2 * pi^2    / 0
  !
  USE constants,   ONLY : tpi
  USE fft_scalar,  ONLY : cft_1z
  USE fft_support, ONLY : good_fft_order
  USE kinds,       ONLY : DP
  USE mp,          ONLY : mp_sum
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define variables
  TYPE radfft_type
    INTEGER           :: ngrid        ! number of grids
    INTEGER           :: mgrid        ! number of grids for FFT-box
    INTEGER           :: lgrid        ! modified mgrid
    INTEGER           :: igrid_start  ! starting index of grids (for MPI)
    INTEGER           :: igrid_end    ! ending index of grids (for MPI)
    INTEGER           :: igrid_len    ! length of grids (for MPI)
    INTEGER           :: comm         ! MPI-communicator
    LOGICAL           :: lmpi         ! use MPI, or not ?
    REAL(DP), POINTER :: rgrid(:)     ! grids in R-space
    REAL(DP), POINTER :: ggrid(:)     ! grids in G-space
    REAL(DP), POINTER :: singr(:,:)   ! sin(g*r) for Fourier Transform, with BLAS level-3
  END TYPE radfft_type
  !
  ! ... public components
  PUBLIC :: radfft_type
  PUBLIC :: allocate_radfft
  PUBLIC :: init_mpi_radfft
  PUBLIC :: deallocate_radfft
  PUBLIC :: fw_radfft
  PUBLIC :: inv_radfft
  PUBLIC :: fw_mpi_radfft
  PUBLIC :: inv_mpi_radfft
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_radfft(radfft0, nr, rmax)
    !--------------------------------------------------------------------------
    !
    ! ... initialize radfft_type
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(INOUT) :: radfft0
    INTEGER,           INTENT(IN)    :: nr
    REAL(DP),          INTENT(IN)    :: rmax
    !
    INTEGER  :: igrid
    INTEGER  :: ngrid
    REAL(DP) :: hr
    REAL(DP) :: hg
    !
    ! check nr
    IF (nr < 2) THEN
      CALL errore(' allocate_radfft ', ' too small number of grids ', 1)
    END IF
    !
    ! number of grids
    ngrid         = nr
    radfft0%ngrid = ngrid
    radfft0%mgrid = 2 * ngrid - 1
    radfft0%lgrid = -1 ! do it after
    !
    ! not MPI (as default)
    radfft0%lmpi = .FALSE.
    !
    ! R-space
    ALLOCATE(radfft0%rgrid(ngrid))
    hr = rmax / DBLE(ngrid)
    DO igrid = 1, ngrid
      radfft0%rgrid(igrid) = hr * DBLE(igrid - 1)
    END DO
    !
    ! G-space
    ALLOCATE(radfft0%ggrid(ngrid))
    hg = (tpi / rmax) * (DBLE(ngrid) / DBLE(2 * ngrid - 1))
    DO igrid = 1, ngrid
      radfft0%ggrid(igrid) = hg * DBLE(igrid - 1)
    END DO
    !
  END SUBROUTINE allocate_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE init_mpi_radfft(radfft0, comm_, igrid1, igrid2)
    !--------------------------------------------------------------------------
    !
    ! ... initialize radfft_type about MPI
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(INOUT) :: radfft0
    INTEGER,           INTENT(IN)    :: comm_
    INTEGER,           INTENT(IN)    :: igrid1
    INTEGER,           INTENT(IN)    :: igrid2
    !
    INTEGER  :: ir, ig
    INTEGER  :: iir
    REAL(DP) :: r, g
    !
    ! set MPI
    radfft0%lmpi = .TRUE.
    radfft0%comm = comm_
    !
    ! size of grids for MPI
    radfft0%igrid_start = MAX(igrid1, 1)
    radfft0%igrid_end   = MIN(igrid2, radfft0%ngrid)
    radfft0%igrid_len   = radfft0%igrid_end - radfft0%igrid_start + 1
    !
    IF (radfft0%igrid_len < 1) THEN
      RETURN
    END IF
    !
    ! calculate sin(g*r)
    ALLOCATE(radfft0%singr(radfft0%ngrid, radfft0%igrid_len))
    !
!$omp parallel do default(shared) private(ir, iir, r, ig, g)
    DO ir = radfft0%igrid_start, radfft0%igrid_end
      iir = ir - radfft0%igrid_start + 1
      r = radfft0%rgrid(ir)
      DO ig = 1, radfft0%ngrid
        g = radfft0%ggrid(ig)
        radfft0%singr(ig, iir) = SIN(g * r)
      END DO
    END DO
!$omp end parallel do
    !
  END SUBROUTINE init_mpi_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_radfft(radfft0)
    !--------------------------------------------------------------------------
    !
    ! ... finalize radfft_type
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(INOUT) :: radfft0
    !
    radfft0%ngrid       = 0
    radfft0%mgrid       = 0
    radfft0%lgrid       = 0
    radfft0%igrid_start = 0
    radfft0%igrid_end   = 0
    radfft0%igrid_len   = 0
    radfft0%comm        = 0
    radfft0%lmpi        = .FALSE.
    IF (ASSOCIATED(radfft0%rgrid)) DEALLOCATE(radfft0%rgrid)
    IF (ASSOCIATED(radfft0%ggrid)) DEALLOCATE(radfft0%ggrid)
    IF (ASSOCIATED(radfft0%singr)) DEALLOCATE(radfft0%singr)
    !
  END SUBROUTINE deallocate_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fw_radfft(radfft0, cr, cg)
    !--------------------------------------------------------------------------
    !
    ! ... FFT R -> G
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(INOUT) :: radfft0
    REAL(DP),          INTENT(IN)    :: cr(:)
    REAL(DP),          INTENT(OUT)   :: cg(:)
    !
    INTEGER                  :: igrid
    REAL(DP)                 :: dr
    REAL(DP)                 :: fac
    COMPLEX(DP), ALLOCATABLE :: crr(:)
    COMPLEX(DP), ALLOCATABLE :: cgg(:)
    !
    ! init lgrid
    IF (radfft0%lgrid < 1) THEN
      radfft0%lgrid = good_fft_order(radfft0%mgrid)
    END IF
    !
    ! allocate memory
    ALLOCATE(crr(radfft0%lgrid))
    ALLOCATE(cgg(radfft0%lgrid))
    !
    ! cr -> crr (ungerade)
    dr  = radfft0%rgrid(2) - radfft0%rgrid(1)
    fac = dr * tpi
    DO igrid = 1, radfft0%ngrid
      crr(igrid) = CMPLX(0.0_DP, fac * radfft0%rgrid(igrid) * cr(igrid), kind=DP)
    END DO
    DO igrid = (radfft0%ngrid + 1), radfft0%mgrid
      crr(igrid) = -crr(2 * radfft0%ngrid - igrid + 1)
    END DO
    !
    ! 1D-FFT (crr -> cgg)
    CALL cft_1z(crr, 1, radfft0%mgrid, radfft0%lgrid, -1, cgg)
    !
    ! cgg -> cg (only sin)
    cg(1) = 0.0_DP
    DO igrid = 2, radfft0%ngrid
      cg(igrid) = DBLE(cgg(igrid)) / radfft0%ggrid(igrid) * DBLE(radfft0%mgrid)
    END DO
    !
    ! deallocate memory
    DEALLOCATE(crr)
    DEALLOCATE(cgg)
    !
  END SUBROUTINE fw_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE inv_radfft(radfft0, cg, cr)
    !--------------------------------------------------------------------------
    !
    ! ... FFT G -> R
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(INOUT) :: radfft0
    REAL(DP),          INTENT(IN)    :: cg(:)
    REAL(DP),          INTENT(OUT)   :: cr(:)
    !
    INTEGER                  :: igrid
    REAL(DP)                 :: dg
    REAL(DP)                 :: fac
    COMPLEX(DP), ALLOCATABLE :: cgg(:)
    COMPLEX(DP), ALLOCATABLE :: crr(:)
    !
    ! init lgrid
    IF (radfft0%lgrid < 1) THEN
      radfft0%lgrid = good_fft_order(radfft0%mgrid)
    END IF
    !
    ! allocate memory
    ALLOCATE(cgg(radfft0%lgrid))
    ALLOCATE(crr(radfft0%lgrid))
    !
    ! cg -> cgg (only sin)
    dg  = radfft0%ggrid(2) - radfft0%ggrid(1)
    fac = -dg / tpi / tpi
    DO igrid = 1, radfft0%ngrid
      cgg(igrid) = CMPLX(0.0_DP, fac * radfft0%ggrid(igrid) * cg(igrid), kind=DP)
    END DO
    DO igrid = (radfft0%ngrid + 1), radfft0%mgrid
      cgg(igrid) = -cgg(2 * radfft0%ngrid - igrid + 1)
    END DO
    !
    ! 1D-FFT (cgg -> crr)
    CALL cft_1z(cgg, 1, radfft0%mgrid, radfft0%lgrid, +1, crr)
    !
    ! crr -> cr (ungerade)
    cr(1) = 0.0_DP
    DO igrid = 2, radfft0%ngrid
      cr(igrid) = DBLE(crr(igrid)) / radfft0%rgrid(igrid)
    END DO
    !
    ! deallocate memory
    DEALLOCATE(cgg)
    DEALLOCATE(crr)
    !
  END SUBROUTINE inv_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fw_mpi_radfft(radfft0, cr, cg, mult)
    !--------------------------------------------------------------------------
    !
    ! ... FFT R -> G
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(IN)  :: radfft0
    REAL(DP),          INTENT(IN)  :: cr(1:*)
    REAL(DP),          INTENT(OUT) :: cg(1:*)
    INTEGER,           INTENT(IN)  :: mult
    !
    INTEGER               :: i
    INTEGER               :: igrid
    INTEGER               :: jgrid
    INTEGER               :: mgrid
    INTEGER               :: iigrid_start
    REAL(DP)              :: dr
    REAL(DP)              :: fac
    REAL(DP), ALLOCATABLE :: crr(:,:)
    REAL(DP), ALLOCATABLE :: cgg(:,:)
    !
    IF (mult < 1) RETURN
    !
    ! allocate memory
    ALLOCATE(crr(radfft0%igrid_len, mult))
    ALLOCATE(cgg(radfft0%ngrid,     mult))
    !
    cgg = 0.0_DP
    !
    IF (radfft0%igrid_len > 0) THEN
      ! cr -> crr
      DO i = 1, mult
        mgrid = (i - 1) * radfft0%igrid_len
!$omp parallel do default(shared) private(igrid, jgrid)
        DO igrid = radfft0%igrid_start, radfft0%igrid_end
          jgrid = igrid - radfft0%igrid_start + 1
          crr(jgrid, i) = cr(jgrid + mgrid) * radfft0%rgrid(igrid)
        END DO
!$omp end parallel do
      END DO
      !
      ! perform integration
      dr  = radfft0%rgrid(2) - radfft0%rgrid(1)
      fac = 2.0_DP * dr * tpi
      CALL dgemm('N', 'N', radfft0%ngrid, mult, radfft0%igrid_len, &
               & fac, radfft0%singr, radfft0%ngrid, crr, radfft0%igrid_len, &
               & 0.0_DP, cgg, radfft0%ngrid)
    END IF
    !
    CALL mp_sum(cgg, radfft0%comm)
    !
    IF (radfft0%igrid_len > 0) THEN
      ! cgg -> cg
      DO i = 1, mult
        mgrid = (i - 1) * radfft0%igrid_len
        iigrid_start = radfft0%igrid_start
        IF (iigrid_start == 1) THEN
          iigrid_start = 2
          cg(1 + mgrid) = 0.0_DP
        END IF
        !
!$omp parallel do default(shared) private(igrid, jgrid)
        DO igrid = iigrid_start, radfft0%igrid_end
          jgrid = igrid - radfft0%igrid_start + 1
          cg(jgrid + mgrid) = cgg(igrid, i) / radfft0%ggrid(igrid)
        END DO
!$omp end parallel do
      END DO
    END IF
    !
    ! deallocate memory
    DEALLOCATE(crr)
    DEALLOCATE(cgg)
    !
  END SUBROUTINE fw_mpi_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE inv_mpi_radfft(radfft0, cg, cr, mult)
    !--------------------------------------------------------------------------
    !
    ! ... FFT G -> R
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(IN)  :: radfft0
    REAL(DP),          INTENT(IN)  :: cg(1:*)
    REAL(DP),          INTENT(OUT) :: cr(1:*)
    INTEGER,           INTENT(IN)  :: mult
    !
    INTEGER               :: i
    INTEGER               :: igrid
    INTEGER               :: jgrid
    INTEGER               :: mgrid
    INTEGER               :: iigrid_start
    REAL(DP)              :: dg
    REAL(DP)              :: fac
    REAL(DP), ALLOCATABLE :: cgg(:,:)
    REAL(DP), ALLOCATABLE :: crr(:,:)
    !
    IF (mult < 1) RETURN
    !
    ! allocate memory
    ALLOCATE(cgg(radfft0%ngrid,     mult))
    ALLOCATE(crr(radfft0%igrid_len, mult))
    !
    cgg = 0.0_DP
    !
    IF (radfft0%igrid_len > 0) THEN
      ! cg -> cgg
      DO i = 1, mult
        mgrid = (i - 1) * radfft0%igrid_len
!$omp parallel do default(shared) private(igrid, jgrid)
        DO igrid = radfft0%igrid_start, radfft0%igrid_end
          jgrid = igrid - radfft0%igrid_start + 1
          cgg(igrid, i) = cg(jgrid + mgrid) * radfft0%ggrid(igrid)
        END DO
!$omp end parallel do
      END DO
    END IF
    !
    CALL mp_sum(cgg, radfft0%comm)
    !
    IF (radfft0%igrid_len > 0) THEN
      ! perform integration
      dg  = radfft0%ggrid(2) - radfft0%ggrid(1)
      fac = 2.0_DP * dg / tpi / tpi
      CALL dgemm('T', 'N', radfft0%igrid_len, mult, radfft0%ngrid, &
               & fac, radfft0%singr, radfft0%ngrid, cgg, radfft0%ngrid, &
               & 0.0_DP, crr, radfft0%igrid_len)
      !
      ! crr -> cr
      DO i = 1, mult
        mgrid = (i - 1) * radfft0%igrid_len
        iigrid_start = radfft0%igrid_start
        IF (iigrid_start == 1) THEN
          iigrid_start = 2
          cr(1 + mgrid) = 0.0_DP
        END IF
        !
!$omp parallel do default(shared) private(igrid, jgrid)
        DO igrid = iigrid_start, radfft0%igrid_end
          jgrid = igrid - radfft0%igrid_start + 1
          cr(jgrid + mgrid) = crr(jgrid, i) / radfft0%rgrid(igrid)
        END DO
!$omp end parallel do
      END DO
    END IF
    !
    ! deallocate memory
    DEALLOCATE(cgg)
    DEALLOCATE(crr)
    !
  END SUBROUTINE inv_mpi_radfft
  !
END MODULE radfft
