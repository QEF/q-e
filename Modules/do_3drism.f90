!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE do_3drism(rismt, maxiter, rmsconv, nbox, eta, title, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... perform 3D-RISM method
  ! ...
  ! ... Variables:
  ! ...   rismt:   data structure
  ! ...   maxiter: maximum number of iterations
  ! ...   rmsconv: RMS of residual vectors to check convergence
  ! ...   nbox:    box size of MDIIS
  ! ...   eta:     step radius of MDIIS
  ! ...   title:   subtitle of calculation
  ! ...   ierr:    status of calculation
  !
  USE check_stop,     ONLY : check_stop_now, stopped_by_user
  USE control_flags,  ONLY : iverbosity, gamma_only
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE, IERR_RISM_NOT_CONVERGED
  USE fft_interfaces, ONLY : fwfft, invfft
  USE fft_types,      ONLY : fft_index_to_3d
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE mdiis,          ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis, reset_mdiis
  USE rism,           ONLY : rism_type, ITYPE_3DRISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, get_nsite_in_solVs, nsolV, solVs, &
                           & iuniq_to_nsite, iuniq_to_isite, isite_to_isolV
  !
  IMPLICIT NONE
  !
  TYPE(rism_type),  INTENT(INOUT) :: rismt
  INTEGER,          INTENT(IN)    :: maxiter
  REAL(DP),         INTENT(IN)    :: rmsconv
  INTEGER,          INTENT(IN)    :: nbox
  CHARACTER(LEN=*), INTENT(IN)    :: title
  REAL(DP),         INTENT(IN)    :: eta
  INTEGER,          INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iter
  INTEGER                  :: ngrid
  INTEGER                  :: nsite
  LOGICAL                  :: lconv
  REAL(DP)                 :: rmscurr
  REAL(DP)                 :: rmssave
  INTEGER                  :: rmswarn
  REAL(DP),    ALLOCATABLE :: csr (:,:)
  REAL(DP),    ALLOCATABLE :: dcsr(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  TYPE(mdiis_type)         :: mdiist
  ! if mdiist is an automatic variable,
  ! pointers in mdiis_type may not work well.
  SAVE                     :: mdiist
  REAL(DP)                 :: csr_ (1, 1)
  REAL(DP)                 :: dcsr_(1, 1)
  !
  INTEGER,     PARAMETER   :: NPRINT      = 10
  INTEGER,     PARAMETER   :: MDIIS_EXT   = 3
  INTEGER,     PARAMETER   :: RMSWARN_MAX = 16
  REAL(DP),    PARAMETER   :: RMS_SMALL   = 0.95_DP
  REAL(DP),    PARAMETER   :: RMS_LARGE   = 2.00_DP
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ng < rismt%gvec%ngm) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  nq = get_nuniq_in_solVs()
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
    ALLOCATE(csr( rismt%dfft%nnr, rismt%nsite))
    ALLOCATE(dcsr(rismt%dfft%nnr, rismt%nsite))
  END IF
  IF (rismt%dfft%nnr > 0) THEN
    ALLOCATE(aux(rismt%dfft%nnr))
  END IF
  !
  CALL allocate_mdiis(mdiist, nbox, rismt%dfft%nnr * rismt%nsite, eta, MDIIS_EXT)
  !
  ! ... reset conditions
  lconv       = .FALSE.
  rismt%avail = .FALSE.
  rmssave     = 1.0E+99_DP
  rmswarn     = 0
  !
  ! ... start 3D-RISM iteration
  WRITE(stdout, '()')
  IF (LEN_TRIM(title) > 0) THEN
    WRITE(stdout, '(5X,"3D-RISM Calculation (",A,")")') TRIM(title)
  ELSE
    WRITE(stdout, '(5X,"3D-RISM Calculation")')
  END IF
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"convergence threshold    =",1PE10.3)') rmsconv
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
    CALL fft_csr_to_csg()
    !
    ! ... 3D-RISM eq.: Cs(g) -> H(g)
    CALL eqn_3drism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... FFT: H(g) -> H(r)
    CALL fft_hg_to_hr()
    !
    ! ... Closure: H(r) -> G(r)
#if defined (__DEBUG_RISM)
    CALL start_clock('3DRISM_clos')
    !
#endif
    CALL closure(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('3DRISM_clos')
#endif
    !
    ! ... Residual: G(r) -> dCs(r)
    CALL prepare_csr_dcsr()
    !
    ! ... clean date out of physical range
    CALL clean_out_of_range()
    !
    ! ... calculate RMS
    ngrid = rismt%dfft%nr1 * rismt%dfft%nr2 * rismt%dfft%nr3
    nsite = get_nsite_in_solVs()
    !
    IF (rismt%nr * rismt%nsite > 0) THEN
      CALL rms_residual(ngrid * nsite, rismt%nr * rismt%nsite, &
                      & dcsr,  rmscurr, rismt%intra_comm)
    ELSE
      CALL rms_residual(ngrid * nsite, rismt%nr * rismt%nsite, &
                      & dcsr_, rmscurr, rismt%intra_comm)
    END IF
    !
    IF (rmscurr < rmsconv) THEN
      lconv = .TRUE.
    END IF
    !
    ! ... write data
    IF (iverbosity > 0 .OR. MOD(iter - 1, NPRINT) == 0 .OR. lconv .OR. iter == maxiter) THEN
      WRITE(stdout, '(5X,"iter. #",I6,"  RMS(g-h-1)=",1PE10.3,"  nbox=",I3)') &
      & iter, rmscurr, mdiist%nbox
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
    IF (rmscurr > (RMS_SMALL * rmssave)) THEN
      rmswarn = rmswarn + 1
    ELSE
      rmswarn = 0
    END IF
    IF (rmswarn >= RMSWARN_MAX) THEN
      CALL reset_mdiis(mdiist, .TRUE.)
      rmssave = rmscurr
      rmswarn = 0
    END IF
    !
    IF (rmscurr > (RMS_LARGE * rmssave)) THEN
      CALL reset_mdiis(mdiist, .TRUE.)
      rmssave = rmscurr
      rmswarn = 0
    ELSE
      rmssave = MIN(rmscurr, rmssave)
    END IF
    !
    ! ... MDIIS: dCs(r) -> Cs(r)
#if defined (__DEBUG_RISM)
    CALL start_clock('3DRISM_mdiis')
    !
#endif
    IF (rismt%nr * rismt%nsite > 0) THEN
      CALL update_by_mdiis(mdiist, csr,  dcsr,  rismt%intra_comm)
    ELSE
      CALL update_by_mdiis(mdiist, csr_, dcsr_, rismt%intra_comm)
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('3DRISM_mdiis')
#endif
    !
    CALL restore_csr()
    !
  ! ... end 3D-RISM iteration
  END DO
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"End of 3D-RISM calculation")')
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
    ! ... calculate chemical potential
#if defined (__DEBUG_RISM)
    CALL start_clock('3DRISM_chem')
    !
#endif
    CALL chempot(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('3DRISM_chem')
#endif
    !
    ! ... calculate solvation's charge density, potential and energy
#if defined (__DEBUG_RISM)
    CALL start_clock('3DRISM_solva')
    !
#endif
    CALL solvation_3drism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('3DRISM_solva')
#endif
    !
    ! ... print chemical potential
    CALL print_chempot_3drism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... set conditions
    ierr = IERR_RISM_NULL
    rismt%avail = .TRUE.
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
  IF (rismt%dfft%nnr * rismt%nsite > 0) THEN
    DEALLOCATE(csr)
    DEALLOCATE(dcsr)
  END IF
  IF (rismt%dfft%nnr > 0) THEN
    DEALLOCATE(aux)
  END IF
  !
  CALL deallocate_mdiis(mdiist)
  !
CONTAINS
  !
  SUBROUTINE fft_csr_to_csg()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    INTEGER :: ir
    INTEGER :: ig
    !
#if defined (__DEBUG_RISM)
    CALL start_clock('3DRISM_fft')
    !
#endif
    DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      jsite = isite - rismt%mp_site%isite_start + 1
      IF (rismt%gvec%ngm > 0) THEN
        rismt%csgz(:, jsite) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      !
!$omp parallel do default(shared) private(ir)
      DO ir = 1, rismt%dfft%nnr
        aux(ir) = CMPLX(rismt%csr(ir, jsite), 0.0_DP, kind=DP)
      END DO
!$omp end parallel do
      IF (rismt%dfft%nnr > 0) THEN
        CALL fwfft('Rho', aux, rismt%dfft)
      END IF
      !
!$omp parallel do default(shared) private(ig)
      DO ig = 1, rismt%gvec%ngm
        rismt%csgz(ig, jsite) = aux(rismt%dfft%nl(ig))
      END DO
!$omp end parallel do
    END DO
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('3DRISM_fft')
#endif
    !
  END SUBROUTINE fft_csr_to_csg
  !
  SUBROUTINE fft_hg_to_hr()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    INTEGER :: ir
    INTEGER :: ig
    !
#if defined (__DEBUG_RISM)
    CALL start_clock('3DRISM_fft')
    !
#endif
    DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      jsite = isite - rismt%mp_site%isite_start + 1
      IF (rismt%dfft%nnr > 0) THEN
        rismt%hr(:, jsite) = 0.0_DP
      END IF
      !
      IF (rismt%dfft%nnr > 0) THEN
        aux = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
!$omp parallel do default(shared) private(ig)
      DO ig = 1, rismt%gvec%ngm
        aux(rismt%dfft%nl(ig)) = rismt%hgz(ig, jsite)
      END DO
!$omp end parallel do
      IF (gamma_only) THEN
!$omp parallel do default(shared) private(ig)
        DO ig = 1, rismt%gvec%ngm
          aux(rismt%dfft%nlm(ig)) = CONJG(rismt%hgz(ig, jsite))
        END DO
!$omp end parallel do
      END IF
      IF (rismt%dfft%nnr > 0) THEN
        CALL invfft('Rho', aux, rismt%dfft)
      END IF
      !
!$omp parallel do default(shared) private(ir)
      DO ir = 1, rismt%dfft%nnr
        rismt%hr(ir, jsite) = DBLE(aux(ir))
      END DO
!$omp end parallel do
    END DO
#if defined (__DEBUG_RISM)
    !
    CALL stop_clock('3DRISM_fft')
#endif
    !
  END SUBROUTINE fft_hg_to_hr
  !
  SUBROUTINE clean_out_of_range()
    IMPLICIT NONE
    INTEGER :: ir
    INTEGER :: i1, i2, i3
    LOGICAL :: offrange
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
!$omp parallel do default(shared) private(ir, i1, i2, i3, offrange)
    DO ir = 1, rismt%dfft%nnr
      !
      IF (ir <= rismt%dfft%nr1x * rismt%dfft%my_nr2p * rismt%dfft%my_nr3p) THEN
        CALL fft_index_to_3d(ir, rismt%dfft, i1, i2, i3, offrange)
      ELSE
        offrange = .TRUE.
      END IF
      !
      IF (offrange) THEN
        rismt%csr(ir, :) = 0.0_DP
        rismt%hr (ir, :) = 0.0_DP
        rismt%gr (ir, :) = 0.0_DP
        csr      (ir, :) = 0.0_DP
        dcsr     (ir, :) = 0.0_DP
      END IF
      !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE clean_out_of_range
  !
  SUBROUTINE prepare_csr_dcsr()
    IMPLICIT NONE
    INTEGER  :: iq
    INTEGER  :: iiq
    INTEGER  :: iv
    INTEGER  :: nv
    INTEGER  :: isolV
    INTEGER  :: natom
    INTEGER  :: nr
    REAL(DP) :: rhov
    REAL(DP) :: rhovt
    REAL(DP) :: vscale
    !
    rhovt = 0.0_DP
    DO isolV = 1, nsolV
      natom = solVs(isolV)%natom
      rhov  = solVs(isolV)%density
      rhovt = rhovt + rhov * DBLE(natom)
    END DO
    !
    IF (rhovt <= 0.0_DP) THEN ! will not be occurred
      CALL errore('do_3drism', 'rhovt is not positive', 1)
    END IF
    !
    nr = rismt%dfft%nnr
    !
    DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq    = iq - rismt%mp_site%isite_start + 1
      iv     = iuniq_to_isite(1, iq)
      nv     = iuniq_to_nsite(iq)
      isolV  = isite_to_isolV(iv)
      rhov   = solVs(isolV)%density
      vscale = SQRT(DBLE(nv))
      !
      IF (nr > 0) THEN
        csr( 1:nr, iiq) = vscale * rismt%csr(1:nr, iiq)
        dcsr(1:nr, iiq) = vscale * (rhov / rhovt) &
                      & * (rismt%gr(1:nr, iiq) - rismt%hr(1:nr, iiq) - 1.0_DP)
      END IF
    END DO
    !
  END SUBROUTINE prepare_csr_dcsr
  !
  SUBROUTINE restore_csr()
    IMPLICIT NONE
    INTEGER  :: iq
    INTEGER  :: iiq
    INTEGER  :: nv
    INTEGER  :: nr
    REAL(DP) :: vscale
    !
    nr = rismt%dfft%nnr
    !
    DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq    = iq - rismt%mp_site%isite_start + 1
      nv     = iuniq_to_nsite(iq)
      vscale = SQRT(DBLE(nv))
      !
      IF (nr > 0) THEN
        rismt%csr(1:nr, iiq) = csr(1:nr, iiq) / vscale
      END IF
    END DO
    !
  END SUBROUTINE restore_csr
  !
END SUBROUTINE do_3drism
