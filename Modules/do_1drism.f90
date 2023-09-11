!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE do_1drism(rismt, maxiter, rmsconv, nbox, eta, gbond, lhand, cool, title, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... perform 1D-RISM method
  ! ...
  ! ... Variables:
  ! ...   rismt:   data structure
  ! ...   maxiter: maximum number of iterations
  ! ...   rmsconv: RMS of residual vectors to check convergence
  ! ...   nbox:    box size of MDIIS
  ! ...   eta:     step radius of MDIIS
  ! ...   gbond:   gaussian width of bonds
  ! ...   lhand:   if true, right-hand. if false, left-hand
  ! ...   cool:    apply cooling, or not
  ! ...   title:   subtitle of calculation
  ! ...   ierr:    status of calculation
  !
  USE check_stop,    ONLY : check_stop_now, stopped_by_user
  USE control_flags, ONLY : iverbosity
  USE err_rism,      ONLY : merge_ierr_rism, IERR_RISM_NULL, &
                          & IERR_RISM_INCORRECT_DATA_TYPE, IERR_RISM_NOT_CONVERGED
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE mdiis,         ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis, reset_mdiis
  USE mp,            ONLY : mp_bcast
  USE radfft,        ONLY : fw_radfft, inv_radfft, fw_mpi_radfft, inv_mpi_radfft
  USE rism,          ONLY : rism_type, ITYPE_1DRISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type),  INTENT(INOUT) :: rismt
  INTEGER,          INTENT(IN)    :: maxiter
  REAL(DP),         INTENT(IN)    :: rmsconv
  INTEGER,          INTENT(IN)    :: nbox
  REAL(DP),         INTENT(IN)    :: eta
  REAL(DP),         INTENT(IN)    :: gbond
  CHARACTER(LEN=*), INTENT(IN)    :: title
  LOGICAL,          INTENT(IN)    :: lhand
  LOGICAL,          INTENT(IN)    :: cool
  INTEGER,          INTENT(OUT)   :: ierr
  !
  INTEGER               :: iter
  INTEGER               :: msite
  LOGICAL               :: lconv
  REAL(DP)              :: rmscurr
  REAL(DP)              :: rmssave
  REAL(DP)              :: gmax
  REAL(DP)              :: temporg
  REAL(DP), ALLOCATABLE :: work(:,:)
  REAL(DP), ALLOCATABLE :: dcsrr(:,:)
  REAL(DP), ALLOCATABLE :: csrr(:,:)
  REAL(DP), ALLOCATABLE :: csr_(:,:)
  TYPE(mdiis_type)      :: mdiist
  ! if mdiist is an automatic variable,
  ! pointers in mdiis_type may not work well.
  SAVE                  :: mdiist
  !
  INTEGER,  PARAMETER   :: NPRINT     = 100
  INTEGER,  PARAMETER   :: MDIIS_EXT  = 5
  REAL(DP), PARAMETER   :: RMS_LARGE  = 1.2_DP
  REAL(DP), PARAMETER   :: GMAX_SCALE = 2.0_DP
  !
  REAL(DP), PARAMETER   :: TEMP_MAX   = 3000.0_DP
  REAL(DP), PARAMETER   :: TEMP_SCALE = 8.0_DP
  REAL(DP), PARAMETER   :: TEMP_DENOM = 1.5_DP
  REAL(DP), PARAMETER   :: TEMP_EXPON = 0.75_DP
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_1DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr /= rismt%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  msite = rismt%mp_site%isite_end - rismt%mp_site%isite_start + 1
  IF (msite > 0) THEN
    ALLOCATE(work(rismt%mp_task%nvec, msite))
  ELSE
    ALLOCATE(work(1, 1))
  END IF
  ALLOCATE(dcsrr(rismt%nr, rismt%nsite))
  ALLOCATE(csrr(rismt%nr, rismt%nsite))
  !
  CALL allocate_mdiis(mdiist, nbox, rismt%nr * rismt%nsite, eta, MDIIS_EXT)
  !
  ! ... reset conditions
  lconv       = .FALSE.
  rismt%avail = .FALSE.
  rmssave     = 1.0E+99_DP
  rmscurr     = 1.0E+99_DP
  !
  IF (gbond > 0.0_DP) THEN
    gmax = GMAX_SCALE / gbond
  ELSE
    gmax = 0.0_DP
  END IF
  !
  temporg = rismt%temp
  IF (cool) THEN
    rismt%temp = MIN(TEMP_MAX, TEMP_SCALE * rismt%temp)
  END IF
  !
  ! ... start 1D-RISM iteration
  WRITE(stdout, '()')
  IF (LEN_TRIM(title) > 0) THEN
    WRITE(stdout, '(5X,"1D-RISM Calculation (",A,")")') TRIM(title)
  ELSE
    WRITE(stdout, '(5X,"1D-RISM Calculation")')
  END IF
  WRITE(stdout, '()')
#if defined (__DEBUG_RISM)
  !
  IF (iverbosity > 0) THEN
    CALL write_rism_type(rismt)
  END IF
#endif
  !
  DO iter = 1, maxiter
    !
    ! ... stop by user
    IF (check_stop_now()) THEN
      EXIT
    END IF
    !
    ! ... FFT: Cs(r) -> Cs(g)
    IF (rismt%is_intra) THEN
      CALL fft_csr_to_csg()
    END IF
    !
    ! ... 1D-RISM eq.: Cs(g) -> H(g)
    IF (rismt%is_intra) THEN
      !CALL eqn_1drism(rismt, gmax, lhand, ierr)
      CALL eqn_1drism(rismt, -1.0_DP, lhand, ierr)
    ELSE
      ierr = IERR_RISM_NULL
    END IF
    CALL merge_ierr_rism(ierr, rismt%super_comm)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... FFT: H(g) -> H(r)
    IF (rismt%is_intra) THEN
      CALL fft_hg_to_hr()
    END IF
    !
    ! ... Closure: H(r) -> G(r)
#if defined (__DEBUG_RISM)
    CALL start_clock('1DRISM_clos')
    !
#endif
    IF (rismt%is_intra) THEN
      CALL closure(rismt, ierr)
    ELSE
      ierr = IERR_RISM_NULL
    END IF
    CALL merge_ierr_rism(ierr, rismt%super_comm)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('1DRISM_clos')
#endif
    !
    ! ... Residual: G(r) -> dCs(r) * r
    IF (rismt%is_intra) THEN
      CALL make_residual()
    END IF
    !
    ! ... calculate RMS
    IF (rismt%is_intra) THEN
      CALL rms_residual(rismt%mp_task%nvec * rismt%nsite, &
      & rismt%nr * rismt%nsite, dcsrr, rmscurr, rismt%intra_comm)
      rmscurr = rmscurr / rismt%rfft%rgrid(rismt%rfft%ngrid) * SQRT(3.0_DP)
      !
      IF (rismt%temp <= temporg .AND. rmscurr < rmsconv) THEN
        lconv = .TRUE.
      END IF
      IF (cool .AND. rmscurr < (rmsconv ** TEMP_EXPON)) THEN
        rismt%temp = MAX(rismt%temp / TEMP_DENOM, temporg)
      END IF
    END IF
    !
    CALL mp_bcast(lconv,      rismt%super_root, rismt%super_comm)
    CALL mp_bcast(rismt%temp, rismt%super_root, rismt%super_comm)
    !
    ! ... write data
    IF (iverbosity > 0 .OR. MOD(iter - 1, NPRINT) == 0 .OR. lconv .OR. iter == maxiter) THEN
      WRITE(stdout, '(5X,"iter. #",I6,"  RMS(g-h-1)=",1PE10.3,"  nbox=",I3,"  T=",0PF8.2,"K")') &
      & iter, rmscurr, mdiist%nbox, rismt%temp
      FLUSH(stdout)
    END IF
#if defined (__DEBUG_RISM)
    !
    IF (iverbosity > 0) THEN
      CALL write_rism_type(rismt)
    END IF
#endif
    !
    ! ... converged ?
    IF (lconv) THEN
      EXIT
    END IF
    !
    ! ... check RMS to reset MDIIS
    IF (rismt%is_intra) THEN
      IF (rmscurr > (RMS_LARGE * rmssave)) THEN
        CALL reset_mdiis(mdiist)
        rmssave = rmscurr
      ELSE
        rmssave = MIN(rmscurr, rmssave)
      END IF
    END IF
    !
    ! ... MDIIS: dCs(r) * r -> Cs(r)
    IF (rismt%is_intra) THEN
      CALL perform_mdiis()
    END IF
    !
  ! ... end 1D-RISM iteration
  END DO
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"End of 1D-RISM calculation")')
  WRITE(stdout, '()')
  !
  ! ... iteration has been converged ?
  IF (lconv) THEN
    ! ... write convergence message
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"convergence has been achieved in ",I5," iterations")') iter
    WRITE(stdout, '()')
    FLUSH(stdout)
    !
    ! ... remove large `g' components from Cs
    IF (rismt%is_intra) THEN
      IF (gmax > 0.0_DP) THEN
        ALLOCATE(csr_(rismt%nr, rismt%nsite))
        csr_ = rismt%csr
        CALL remove_glarge_csr()
      END IF
    END IF
    !
    ! ... correct correlations at G = 0 or R = 0
    IF (rismt%is_intra) THEN
      CALL correctat0_vv(rismt, ierr)
    ELSE
      ierr = IERR_RISM_NULL
    END IF
    CALL merge_ierr_rism(ierr, rismt%super_comm)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... calculate chemical potential
    IF (rismt%is_intra) THEN
      CALL chempot(rismt, ierr)
    ELSE
      ierr = IERR_RISM_NULL
    END IF
    CALL merge_ierr_rism(ierr, rismt%super_comm)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... print chemical potential
    IF (rismt%is_intra) THEN
      CALL print_chempot_vv(rismt, lhand, ierr)
    ELSE
      ierr = IERR_RISM_NULL
    END IF
    CALL merge_ierr_rism(ierr, rismt%super_comm)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... restore Cs
    IF (rismt%is_intra) THEN
      IF (gmax > 0.0_DP) THEN
#if !defined (__RISM_SMOOTH_CSVV)
        IF (rismt%mp_task%ivec_start == 1) THEN
          rismt%csr(2:rismt%nr, :) = csr_(2:rismt%nr, :)
        ELSE
          rismt%csr = csr_
        END IF
#endif
        DEALLOCATE(csr_)
      END IF
    END IF
    !
    ! ... set conditions
    ierr = IERR_RISM_NULL
    rismt%avail = .TRUE.
    rismt%temp  = temporg
    !
  ELSE
    ! ... write NOT convergence message
    WRITE(stdout, '()')
    IF (stopped_by_user) THEN
      WRITE(stdout, '(5X,"convergence NOT achieved: stopped by user")')
    ELSE
      WRITE(stdout, '(5X,"convergence NOT achieved")')
    END IF
    WRITE(stdout, '()')
    FLUSH(stdout)
    !
    ! ... set conditions
    ierr = IERR_RISM_NOT_CONVERGED
    rismt%avail = .FALSE.
    rismt%temp  = temporg
  END IF
#if defined (__DEBUG_RISM)
  !
  IF (iverbosity > 0) THEN
    CALL write_rism_type(rismt)
  END IF
#endif
  !
  ! ... deallocate memory
100 CONTINUE
  !
  DEALLOCATE(work)
  DEALLOCATE(dcsrr)
  DEALLOCATE(csrr)
  !
  CALL deallocate_mdiis(mdiist)
  !
CONTAINS
  !
  SUBROUTINE make_residual()
    IMPLICIT NONE
    INTEGER  :: isite
    INTEGER  :: ir
    INTEGER  :: iir
    REAL(DP) :: r
    !
    DO isite = 1, rismt%nsite
      DO ir = 1, rismt%nr
        iir = rismt%mp_task%ivec_start + ir - 1
        r = rismt%rfft%rgrid(iir)
        dcsrr(ir, isite) = (rismt%gr(ir, isite) - rismt%hr(ir, isite) - 1.0_DP) * r
      END DO
    END DO
  END SUBROUTINE make_residual
  !
  SUBROUTINE perform_mdiis()
    IMPLICIT NONE
    INTEGER  :: isite
    INTEGER  :: ir
    INTEGER  :: iir
    INTEGER  :: jr
    REAL(DP) :: r
    !
#if defined (__DEBUG_RISM)
    CALL start_clock('1DRISM_mdiis')
    !
#endif
    ! ... Cs(r) -> Cs(r) * r
    DO isite = 1, rismt%nsite
      DO ir = 1, rismt%nr
        iir = rismt%mp_task%ivec_start + ir - 1
        r = rismt%rfft%rgrid(iir)
        csrr(ir, isite) = rismt%csr(ir, isite) * r
      END DO
    END DO
    !
    ! ... update Cs(r) * r
    CALL update_by_mdiis(mdiist, csrr, dcsrr, rismt%intra_comm)
    !
    ! ... Cs(r) * r -> Cs(r)
    DO isite = 1, rismt%nsite
      IF (rismt%mp_task%ivec_start == 1) THEN
        jr = 2
        rismt%csr(1, isite) = 0.0_DP
      ELSE
        jr = 1
      END IF
      DO ir = jr, rismt%nr
        iir = rismt%mp_task%ivec_start + ir - 1
        r = rismt%rfft%rgrid(iir)
        rismt%csr(ir, isite) = csrr(ir, isite) / r
      END DO
    END DO
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('1DRISM_mdiis')
#endif
    !
  END SUBROUTINE perform_mdiis
  !
  SUBROUTINE fft_csr_to_csg()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    !
#if defined (__DEBUG_RISM)
    CALL start_clock('1DRISM_fft')
    !
#endif
    IF (rismt%rfft%lmpi) THEN
      ! ... Fourier Transform, with BLAS level-3
      CALL fw_mpi_radfft(rismt%rfft, rismt%csr, rismt%csg, rismt%nsite)
      !
    ELSE
      ! ... csr -> work
      CALL mp_swap_ax_rism(rismt%mp_site, rismt%mp_task, &
      & rismt%mp_task%nvec, work, rismt%nr, rismt%csr, +1)
      !
      ! ... FFT work
      DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        jsite = isite - rismt%mp_site%isite_start + 1
        CALL fw_radfft(rismt%rfft, work(:, jsite), work(:, jsite))
      END DO
      !
      ! ... work -> csg
      CALL mp_swap_ax_rism(rismt%mp_site, rismt%mp_task, &
      & rismt%mp_task%nvec, work, rismt%ng, rismt%csg, -1)
      !
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('1DRISM_fft')
#endif
    !
  END SUBROUTINE fft_csr_to_csg
  !
  SUBROUTINE fft_hg_to_hr()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    !
#if defined (__DEBUG_RISM)
    CALL start_clock('1DRISM_fft')
    !
#endif
    IF (rismt%rfft%lmpi) THEN
      ! ... Fourier Transform, with BLAS level-3
      CALL inv_mpi_radfft(rismt%rfft, rismt%hg, rismt%hr, rismt%nsite)
      !
    ELSE
      ! ... hg -> work
      CALL mp_swap_ax_rism(rismt%mp_site, rismt%mp_task, &
      & rismt%mp_task%nvec, work(1, 1), rismt%ng, rismt%hg(1, 1), +1)
      !
      ! ... FFT work
      DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        jsite = isite - rismt%mp_site%isite_start + 1
        CALL inv_radfft(rismt%rfft, work(:, jsite), work(:, jsite))
      END DO
      !
      ! ... work -> hr
      CALL mp_swap_ax_rism(rismt%mp_site, rismt%mp_task, &
      & rismt%mp_task%nvec, work(1, 1), rismt%nr, rismt%hr(1, 1), -1)
      !
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('1DRISM_fft')
#endif
    !
  END SUBROUTINE fft_hg_to_hr
  !
  SUBROUTINE remove_glarge_csr()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    INTEGER :: ig
    INTEGER :: iig
    !
    ! ... remove large `g' components
    DO isite = 1, rismt%nsite
      DO ig = 1, rismt%ng
        iig = rismt%mp_task%ivec_start + ig - 1
        IF (rismt%rfft%ggrid(iig) > gmax) THEN
          rismt%csg(ig, isite) = 0.0_DP
        END IF
      END DO
    END DO
    !
    IF (rismt%rfft%lmpi) THEN
      ! ... Fourier Transform, with BLAS level-3
      CALL inv_mpi_radfft(rismt%rfft, rismt%csg, rismt%csr, rismt%nsite)
      !
    ELSE
      ! ... csg -> work
      CALL mp_swap_ax_rism(rismt%mp_site, rismt%mp_task, &
      & rismt%mp_task%nvec, work(1, 1), rismt%ng, rismt%csg(1, 1), +1)
      !
      ! ... FFT work
      DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        jsite = isite - rismt%mp_site%isite_start + 1
        CALL inv_radfft(rismt%rfft, work(:, jsite), work(:, jsite))
      END DO
      !
      ! ... work -> csr
      CALL mp_swap_ax_rism(rismt%mp_site, rismt%mp_task, &
      & rismt%mp_task%nvec, work(1, 1), rismt%nr, rismt%csr(1, 1), -1)
      !
    END IF
    !
  END SUBROUTINE remove_glarge_csr
  !
END SUBROUTINE do_1drism
