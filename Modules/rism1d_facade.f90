!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE rism1d_facade
  !--------------------------------------------------------------------------
  !
  ! ... Facade (or Interface) of 1D-RISM's library.
  ! ... External codes, which utilize 1D-RISM, can access data and subroutines
  ! ... throught this module.
  !
  USE err_rism,    ONLY : stop_by_err_rism, IERR_RISM_NULL, IERR_RISM_NOT_CONVERGED
  USE io_global,   ONLY : stdout, ionode
  USE io_rism_xml, ONLY : write_1drism, read_1drism
  USE kinds,       ONLY : DP
  USE mp,          ONLY : mp_rank, mp_size, mp_sum, mp_barrier, mp_comm_split
  USE mp_images,   ONLY : intra_image_comm
  USE rism,        ONLY : rism_type, clean_rism_data, allocate_1drism, deallocate_rism
  USE solvmol,     ONLY : get_nsite_in_solVs
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
  INTEGER, PARAMETER     :: LEN_STR = 30
  !
  ! ... define variables
  LOGICAL                :: lrism1d       = .FALSE.  ! to calculate 1D-RISM, or not
  INTEGER                :: nproc_sub     = 0        ! total size of MPI's subspace
  INTEGER                :: nproc_switch  = 0        ! size of MPI's subspace to switch algorithm
  LOGICAL                :: has_any_corr  = .FALSE.  ! has nonzero correlations
  CHARACTER(LEN=LEN_STR) :: starting_corr = ''       ! initial correlations: 'zero', 'file', 'fix'
  INTEGER                :: niter         = 0        ! maximum number of iteration
  REAL(DP)               :: epsv          = 0.0_DP   ! convergence threshold
  REAL(DP)               :: bond_width    = 0.0_DP   ! gaussian width of bonds (in bohr)
  INTEGER                :: mdiis_size    = 0        ! size of MDIIS
  REAL(DP)               :: mdiis_step    = 0.0_DP   ! step of MDIIS
  REAL(DP)               :: dielectric    = 0.0_DP   ! dielectric constant (for DRISM)
  REAL(DP)               :: molesize      = 0.0_DP   ! size of molecule (in bohr, for DRISM)
  !
  ! ... define 1D-RISM's main data
  TYPE(rism_type), POINTER :: rism1t => NULL()
  TYPE(rism_type), TARGET  :: rism1t_right
  TYPE(rism_type), TARGET  :: rism1t_left
  LOGICAL                  :: init_rism1t_right = .FALSE.
  LOGICAL                  :: init_rism1t_left  = .FALSE.
  !
  ! ... public components
  PUBLIC :: lrism1d
  PUBLIC :: nproc_sub
  PUBLIC :: nproc_switch
  PUBLIC :: starting_corr
  PUBLIC :: niter
  PUBLIC :: epsv
  PUBLIC :: bond_width
  PUBLIC :: mdiis_size
  PUBLIC :: mdiis_step
  PUBLIC :: dielectric
  PUBLIC :: molesize
  !
  PUBLIC :: rism1t
  PUBLIC :: rism1d_iosys
  PUBLIC :: rism1d_summary
  PUBLIC :: rism1d_initialize
  PUBLIC :: rism1d_finalize
  PUBLIC :: rism1d_is_avail
  PUBLIC :: rism1d_activate_right
  PUBLIC :: rism1d_activate_left
  PUBLIC :: rism1d_prepare
  PUBLIC :: rism1d_run
  PUBLIC :: rism1d_write_to_show
  PUBLIC :: rism1d_write_to_restart
  PUBLIC :: rism1d_read_to_restart
  PUBLIC :: rism1d_print_clock
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_iosys(trism, laue)
    !----------------------------------------------------------------------------
    !
    ! ... set variables and initialize this module
    !
    IMPLICIT NONE
    !
    LOGICAL,           INTENT(IN) :: trism
    LOGICAL, OPTIONAL, INTENT(IN) :: laue
    !
    LOGICAL :: laue_
    !
    lrism1d = trism
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    IF (PRESENT(laue)) THEN
      laue_ = laue
    ELSE
      laue_ = .FALSE.
    END IF
    !
    CALL iosys_1drism(laue_)
    !
  END SUBROUTINE rism1d_iosys
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_summary()
    !----------------------------------------------------------------------------
    !
    ! ... print conditions
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    CALL summary_1drism()
    !
  END SUBROUTINE rism1d_summary
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_initialize(ngrid, rmax, lboth)
    !----------------------------------------------------------------------------
    !
    ! ... initialize this module
    ! ...
    ! ... Variables:
    ! ...   ngrid: number of grids
    ! ...   rmax:  maximum radius of R-space (in bohr)
    ! ...   lboth: initialize the both rism1t_right and rism1t_left, or not
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)  :: ngrid
    REAL(DP), INTENT(IN)  :: rmax
    LOGICAL,  INTENT(IN)  :: lboth
    !
    INTEGER :: nv
    INTEGER :: rism1d_comm
    INTEGER :: rism1d_root
    LOGICAL :: rism1d_primal
    LOGICAL :: mpi_radfft
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    nv = get_nsite_in_solVs()
    !
    ! ... initialize MPI's information
    CALL rism1d_mpi_init(rism1d_comm, rism1d_root, rism1d_primal, mpi_radfft)
    !
    ! ... initialize rism1t_right
    init_rism1t_right = .TRUE.
    CALL allocate_1drism(rism1t_right, nv, ngrid, rmax, mpi_radfft, &
    & intra_image_comm, rism1d_root, rism1d_primal, rism1d_comm)
    !
    IF (lboth) THEN
      ! ... initialize rism1t_left
      init_rism1t_left = .TRUE.
      CALL allocate_1drism(rism1t_left, nv, ngrid, rmax, mpi_radfft, &
      & intra_image_comm, rism1d_root, rism1d_primal, rism1d_comm)
    END IF
    !
    ! ... activate rism1t_right
    CALL rism1d_activate_right()
    !
  END SUBROUTINE rism1d_initialize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_mpi_init(rism1d_comm, rism1d_root, rism1d_primal, mpi_radfft)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(OUT) :: rism1d_comm
    INTEGER, INTENT(OUT) :: rism1d_root
    LOGICAL, INTENT(OUT) :: rism1d_primal
    LOGICAL, INTENT(OUT) :: mpi_radfft
    !
    INTEGER :: color
    INTEGER :: irank
    INTEGER :: nproc
    !
    irank = mp_rank(intra_image_comm)
    nproc = mp_size(intra_image_comm)
    !
    IF (nproc <= nproc_sub) THEN
      rism1d_comm   = intra_image_comm
      rism1d_root   = 0
      rism1d_primal = .TRUE.
      !
    ELSE
      color = irank / nproc_sub
      rism1d_root = 0
      IF (ionode) THEN
        color = 0
        rism1d_root = irank
      END IF
      CALL mp_barrier(intra_image_comm)
      CALL mp_comm_split(intra_image_comm, color, irank, rism1d_comm)
      CALL mp_sum(rism1d_root, intra_image_comm)
      rism1d_primal = (color == 0)
    END IF
    !
    IF (MIN(nproc, nproc_sub) > nproc_switch) THEN
      mpi_radfft = .TRUE.
    ELSE
      mpi_radfft = .FALSE.
    END IF
    !
  END SUBROUTINE rism1d_mpi_init
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_finalize()
    !----------------------------------------------------------------------------
    !
    ! ... finalize this module
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    ! ... finalize rism1t
    IF (ASSOCIATED(rism1t)) THEN
      NULLIFY(rism1t)
    END IF
    !
    ! ... finalize rism1t_right
    IF (init_rism1t_right) THEN
      init_rism1t_right = .FALSE.
      CALL deallocate_rism(rism1t_right, .TRUE.)
    END IF
    !
    ! ... finalize rism1t_left
    IF (init_rism1t_left) THEN
      init_rism1t_left = .FALSE.
      CALL deallocate_rism(rism1t_left, .TRUE.)
    END IF
    !
  END SUBROUTINE rism1d_finalize
  !
  !----------------------------------------------------------------------------
  FUNCTION rism1d_is_avail() RESULT(avail)
    !----------------------------------------------------------------------------
    !
    ! ... result of 1D-RISM is available, or not
    !
    IMPLICIT NONE
    LOGICAL :: avail
    !
    IF (.NOT. lrism1d) THEN
      avail = .FALSE.
      RETURN
    END IF
    !
    IF (.NOT. (init_rism1t_right .OR. init_rism1t_left)) THEN
      avail = .FALSE.
      RETURN
    END IF
    !
    avail = .TRUE.
    !
    IF (init_rism1t_right) THEN
      avail = avail .AND. rism1t_right%avail
    END IF
    !
    IF (init_rism1t_left) THEN
      avail = avail .AND. rism1t_left%avail
    END IF
    !
  END FUNCTION rism1d_is_avail
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_activate_right()
    !----------------------------------------------------------------------------
    !
    ! ... activate rism1t_right
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    !IF (.NOT. init_rism1t_right) THEN
    !  CALL stop_by_err_rism('rism1d_activate_right', IERR_RISM_1DRISM_IS_NOT_AVAIL)
    !END IF
    !
    rism1t => rism1t_right
    !
  END SUBROUTINE rism1d_activate_right
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_activate_left()
    !----------------------------------------------------------------------------
    !
    ! ... activate rism1t_left
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    !IF (.NOT. init_rism1t_left) THEN
    !  CALL stop_by_err_rism('rism1d_activate_left', IERR_RISM_1DRISM_IS_NOT_AVAIL)
    !END IF
    !
    rism1t => rism1t_left
    !
  END SUBROUTINE rism1d_activate_left
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_prepare()
    !----------------------------------------------------------------------------
    !
    ! ... prepare 1D-RISM's iterative calculation
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    CALL start_clock('1DRISM_pre')
    !
    ! ... initial calculation (potential, intra-molecular correlation, dielectric bridge)
    IF (init_rism1t_right) THEN
      CALL init_1drism(rism1t_right, bond_width, dielectric, molesize, .TRUE., ierr)
      !
      IF (ierr /= IERR_RISM_NULL) THEN
        CALL stop_by_err_rism('rism1d_prepare', ierr)
      END IF
    END IF
    !
    IF (init_rism1t_left) THEN
      CALL init_1drism(rism1t_left, bond_width, dielectric, molesize, .FALSE., ierr)
      !
      IF (ierr /= IERR_RISM_NULL) THEN
        CALL stop_by_err_rism('rism1d_prepare', ierr)
      END IF
    END IF
    !
    ! ... create initial correlation
    IF (TRIM(starting_corr) == 'file' .OR. TRIM(starting_corr) == 'fix') THEN
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Correlation function is read from file")')
      WRITE(stdout, '()')
      !
      IF (init_rism1t_right) THEN
        CALL clean_rism_data(rism1t_right)
      END IF
      IF (init_rism1t_left) THEN
        CALL clean_rism_data(rism1t_left)
      END IF
      !
      CALL rism1d_read_to_restart()
      has_any_corr = .TRUE.
      !
      IF (TRIM(starting_corr) == 'fix') THEN
        rism1t%avail = .TRUE.
      END IF
      !
    ELSE !IF (TRIM(starting_corr) == 'zero') THEN
      IF (init_rism1t_right) THEN
        CALL clean_rism_data(rism1t_right)
      END IF
      IF (init_rism1t_left) THEN
        CALL clean_rism_data(rism1t_left)
      END IF
      !
      has_any_corr = .FALSE.
    END IF
    !
    CALL stop_clock('1DRISM_pre')
    !
  END SUBROUTINE rism1d_prepare
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_run(lconv)
    !----------------------------------------------------------------------------
    !
    ! ... perform 1D-RISM's iterative calculation
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(OUT) :: lconv
    !
    INTEGER           :: ierr
    INTEGER           :: ihand
    CHARACTER(LEN=64) :: title
    !
    IF (.NOT. lrism1d) THEN
      lconv = .FALSE.
      RETURN
    END IF
    !
    ! ... does not calculate, if 'fix'
    IF (TRIM(starting_corr) == 'fix') THEN
      lconv = .FALSE.
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Correlation function is fixed")')
      WRITE(stdout, '()')
      RETURN
    END IF
    !
    lconv = .TRUE.
    !
    DO ihand = 1, 2
      IF (ihand == 1) THEN
        IF (.NOT. init_rism1t_right) CYCLE
      ELSE
        IF (.NOT. init_rism1t_left)  CYCLE
      END IF
      !
      CALL start_clock('1DRISM_run')
      !
      ! ... calculate 1D-RISM
      IF (ihand == 1) THEN
        title = ''
        IF (init_rism1t_right .AND. init_rism1t_left) THEN
          title = 'the right-hand side'
        END IF
        CALL do_1drism(rism1t_right, niter, epsv, mdiis_size, mdiis_step, bond_width, &
                     & .TRUE., .NOT. has_any_corr, TRIM(ADJUSTL(title)), ierr)
        !
      ELSE
        title = ''
        IF (init_rism1t_right .AND. init_rism1t_left) THEN
          title = 'the left-hand side'
        END IF
        CALL do_1drism(rism1t_left,  niter, epsv, mdiis_size, mdiis_step, bond_width, &
                     & .FALSE., .NOT. has_any_corr, TRIM(ADJUSTL(title)), ierr)
      END IF
      !
      ! ... 1D-RISM has been converged ?
      IF (ierr == IERR_RISM_NOT_CONVERGED) THEN
        lconv = lconv .AND. .FALSE.
        !
      ELSE IF (ierr == IERR_RISM_NULL) THEN
        lconv = lconv .AND. .TRUE.
        !
      ELSE ! an error has been occurred
        lconv = lconv .AND. .FALSE.
        CALL stop_by_err_rism('rism1d_run', ierr)
      END IF
      !
      CALL stop_clock('1DRISM_run')
      !
    END DO
    !
    ! ... here, correlation is nonzero.
    has_any_corr = .TRUE.
    !
  END SUBROUTINE rism1d_run
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_write_to_show()
    !----------------------------------------------------------------------------
    !
    ! ... write 1D-RISM's data (human readable)
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    IF (init_rism1t_right .AND. init_rism1t_left) THEN
      CALL print_corr_vv(rism1t_right, '#right', ierr)
      CALL print_corr_vv(rism1t_left,  '#left', ierr)
      !
    ELSE IF (init_rism1t_right) THEN
      CALL print_corr_vv(rism1t_right, '', ierr)
      !
    ELSE IF (init_rism1t_right) THEN
      CALL print_corr_vv(rism1t_left,  '', ierr)
      !
    END IF
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism1d_write_to_show', ierr)
    END IF
    !
    CALL mp_barrier(rism1t%super_comm)
    !
  END SUBROUTINE rism1d_write_to_show
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_write_to_restart(ext)
    !----------------------------------------------------------------------------
    !
    ! ... write 1D-RISM's data (for restart calculation)
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ext
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    IF (init_rism1t_right) THEN
      IF (PRESENT(ext)) THEN
        CALL write_1drism(rism1t_right, '1.' // TRIM(ADJUSTL(ext)))
      ELSE
        CALL write_1drism(rism1t_right, '1')
      END IF
    END IF
    !
    IF (init_rism1t_left) THEN
      IF (PRESENT(ext)) THEN
        CALL write_1drism(rism1t_left, '2.' // TRIM(ADJUSTL(ext)))
      ELSE
        CALL write_1drism(rism1t_left, '2')
      END IF
    END IF
    !
    CALL mp_barrier(rism1t%super_comm)
    !
  END SUBROUTINE rism1d_write_to_restart
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_read_to_restart(ext)
    !----------------------------------------------------------------------------
    !
    ! ... read 1D-RISM's data (for restart calculation)
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ext
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    IF (init_rism1t_right) THEN
      IF (PRESENT(ext)) THEN
        CALL read_1drism(rism1t_right, '1.' // TRIM(ADJUSTL(ext)))
      ELSE
        CALL read_1drism(rism1t_right, '1')
      END IF
    END IF
    !
    IF (init_rism1t_left) THEN
      IF (PRESENT(ext)) THEN
        CALL read_1drism(rism1t_left, '2.' // TRIM(ADJUSTL(ext)))
      ELSE
        CALL read_1drism(rism1t_left, '2')
      END IF
    END IF
    !
  END SUBROUTINE rism1d_read_to_restart
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism1d_print_clock()
    !----------------------------------------------------------------------------
    !
    ! ... print clock for 1D-RISM
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism1d) THEN
      RETURN
    END IF
    !
    CALL print_clock('1DRISM_pre')
    CALL print_clock('1DRISM_run')
#if defined (__DEBUG_RISM)
    CALL print_clock('1DRISM_eqn')
    CALL print_clock('1DRISM_clos')
    CALL print_clock('1DRISM_fft')
    CALL print_clock('1DRISM_mdiis')
#endif
    !
  END SUBROUTINE rism1d_print_clock
  !
END MODULE rism1d_facade

