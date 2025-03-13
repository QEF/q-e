!
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!----------------------------------------------------------------------
MODULE wannier2bloch_opt
!----------------------------------------------------------------------
!!
!! Modules that contains all the routines that transforms quantities from Wannier
!! space to Bloch space.
!!
!! Optimizd for q-points on a grid.
!!
!! Reference: J. Kaye et al., SciPost Phys. 15, 062 (2023) (See Eqs. (8, 9))
!!
!
USE kinds, ONLY : DP
!
IMPLICIT NONE
!
LOGICAL, SAVE :: xq_is_initialized = .FALSE.
!! Set to .TRUE. once xq1_save and xq2_save are initialized
!
REAL(KIND = DP), SAVE :: xq1_save
!! First component of the q-point saved in epmatwp_23
INTEGER, SAVE :: nrr_g_23
!! Number of el-ph WS points projected to the yz plane
INTEGER, SAVE, ALLOCATABLE :: irvec_g_23(:, :)
!! (2, nrr_g_23)
!! Coordinates of real space vector for el-ph projected to the yz plane
COMPLEX(KIND = DP), SAVE, ALLOCATABLE :: epmatwp_23(:, :, :, :, :)
!! (nbndsub, nbndsub, nrr_k, nmodes, nrr_g_23)
!! Partially Fourier-transformed el-ph matrix elements.
!! Allocated only if etf_mem == 0.
!
REAL(KIND = DP), SAVE :: xq2_save
!! Second component of the q-point saved in epmatwp_3
INTEGER, SAVE :: nrr_g_3
!! Number of el-ph WS points projected to the z axis
INTEGER, SAVE, ALLOCATABLE :: irvec_g_3(:)
!! (nrr_g_23)
!! Coordinates of real space vector for el-ph projected to the z axis
COMPLEX(KIND = DP), SAVE, ALLOCATABLE :: epmatwp_3(:, :, :, :, :)
!! (nbndsub, nbndsub, nrr_k, nmodes, nrr_g_3)
!! Partially Fourier-transformed el-ph matrix elements.
!! Allocated only if etf_mem == 0.
!
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ephwan2blochp_opt(nmodes, xq, irvec_g, nrr_g, epmatf, nbnd, nrr_k)
  !--------------------------------------------------------------------------
  !!
  !! Optimized Fourier transformation of el-ph matrix from phonon Wannier to
  !! phonon Bloch basis.
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE mp_pools,         ONLY : me_pool
  USE io_global,        ONLY : stdout
  USE ep_constants,     ONLY : eps8
  USE input,            ONLY : etf_mem
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nmodes
  !! Total number of modes
  INTEGER, INTENT(in) :: nrr_g
  !! Number of phononic WS points
  INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
  !! Coordinates of WS points
  INTEGER, INTENT(in) :: nbnd
  !! Number of bands
  INTEGER, INTENT(in) ::  nrr_k
  !! Number of electronic WS points
  REAL(KIND = DP), INTENT(in) :: xq(3)
  !! q-point vector for the Fourier transformation
  COMPLEX(KIND = DP), INTENT(out) :: epmatf(nbnd, nbnd, nrr_k, nmodes)
  !! el-ph matrix in Bloch representation, fine grid
  !
  INTEGER :: ierr
  !! Error status
  !
  CALL start_clock('ephW2Bp_opt')
  !
  IF (.NOT. xq_is_initialized) THEN
    !
    ! At the first call, the first and second components are updated.
    !
    CALL ephwan2blochp_reset23(xq(1), nmodes, irvec_g, nrr_g, nbnd, nrr_k)
    CALL ephwan2blochp_reset3(xq(2), nmodes, nbnd, nrr_k)
    !
    WRITE(stdout, '(5x, a, I6)') 'nrr_g_3 = ', nrr_g_3
    xq_is_initialized = .TRUE.
    !
  ELSEIF (ABS(xq(1) - xq1_save) > eps8) THEN
    !
    ! If xq(1) /= xq1_save, the first and second components are updated.
    !
    DEALLOCATE(epmatwp_23, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_opt', 'Error deallocating epmatwp_23', 1)
    DEALLOCATE(epmatwp_3, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_opt', 'Error deallocating epmatwp_3', 1)
    DEALLOCATE(irvec_g_23, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_opt', 'Error deallocating irvec_g_23', 1)
    DEALLOCATE(irvec_g_3, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_opt', 'Error deallocating irvec_g_3', 1)
    !
    CALL ephwan2blochp_reset23(xq(1), nmodes, irvec_g, nrr_g, nbnd, nrr_k)
    CALL ephwan2blochp_reset3(xq(2), nmodes, nbnd, nrr_k)
    !
  ELSEIF (ABS(xq(2) - xq2_save) > eps8 ) THEN
    !
    ! The first component is the same, but the second component is updated.
    !
    DEALLOCATE(epmatwp_3, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_opt', 'Error deallocating epmatwp_3', 1)
    DEALLOCATE(irvec_g_3, STAT = ierr)
    IF (ierr /= 0) CALL errore('ephwan2blochp_opt', 'Error deallocating irvec_g_3', 1)
    !
    CALL ephwan2blochp_reset3(xq(2), nmodes, nbnd, nrr_k)
    !
  ENDIF
  !
  ! Fourier transform epmatwp_23 to epmatf.
  !
  CALL ephwan2blochp_get3(xq(3), epmatf, nmodes, nbnd, nrr_k)
  !
  CALL stop_clock('ephW2Bp_opt')
  !
  END SUBROUTINE ephwan2blochp_opt
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ephwan2blochp_get3(xq3, epmat, nmodes, nbnd, nrr_k)
  !--------------------------------------------------------------------------
  !!
  !! Fourier transform along xq3 direction to compute epmat(xq1, xq2, xq3).
  !!
  !! If etf_mem == 1, read epmatwp_3(R_p3) from file.
  !!
  !! epmat = sum_R_p3 e^{i * xq3 * R_p3} epmatwp_3(R_p3)
  !!
  !--------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : inter_pool_comm
  USE ep_constants,     ONLY : twopi, ci, czero
  USE io,               ONLY : rwepmatw
  USE parallelism,      ONLY : para_bounds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nmodes
  !! Total number of modes
  INTEGER, INTENT(in) :: nbnd
  !! Number of bands
  INTEGER, INTENT(in) ::  nrr_k
  !! Number of electronic WS points
  REAL(KIND = DP), INTENT(in) :: xq3
  !! y and z component of the q-point vector for the Fourier transformation
  COMPLEX(KIND = DP), INTENT(out) :: epmat(nbnd, nbnd, nrr_k, nmodes)
  !! el-ph matrix at given phonon q-point, in Cartesian basis
  !
  ! Local variables
  INTEGER :: ir
  !! Real space WS index
  INTEGER :: irn
  !! Combined WS and atom index
  INTEGER :: imode
  !! Counter for modes
  INTEGER :: irn3_start
  !! Starting irn3 for this core
  INTEGER :: irn3_stop
  !! Stoping irn3 for this core
  REAL(KIND = DP) :: rdotk
  !! Exponential for the Fourier transformation
  COMPLEX(KIND = DP) :: cfac
  !! Factor for the Fourier transformation
  !
  CALL start_clock('ephW2Bp_g3')
  !
  epmat = czero
  !
  ! Because nrr_g_3 can be quite small, we do a combined parallelization
  ! on WS vectors and modes.
  !
  CALL para_bounds(irn3_start, irn3_stop, nrr_g_3 * nmodes)
  !
  DO irn = irn3_start, irn3_stop
    !
    ir = (irn - 1) / nmodes + 1
    imode = MOD(irn - 1, nmodes) + 1
    !
    ! Set Fourier transform factor cfac = e^{i * xq23 * (R_p23)}.
    !
    rdotk = twopi * xq3 * REAL(irvec_g_3(ir), KIND=DP)
    cfac = EXP(ci * rdotk)
    !
    ! Fourier transform epmatwp along the third direction to get eptmp_23
    !
    CALL ZAXPY(nbnd * nbnd * nrr_k, cfac, epmatwp_3(:, :, :, imode, ir), 1, epmat(:, :, :, imode), 1)
    !
  ENDDO ! irn
  !
  CALL mp_sum(epmat, inter_pool_comm)
  !
  CALL stop_clock('ephW2Bp_g3')
  !
  !---------------------------------------------------------------------------
  END SUBROUTINE ephwan2blochp_get3
  !---------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ephwan2blochp_reset23(xq1, nmodes, irvec_g, nrr_g, nbnd, nrr_k)
  !--------------------------------------------------------------------------
  !!
  !! Fourier transform along xq3 direction to compute epmat_23(xq1, R2, R3).
  !!
  !! If etf_mem == 1, write epmat_23(xq1, R2, R3) to file. The file is read
  !! in subroutine ephwan2blochp_get23.
  !!
  !! g_23(R_p2, R_p3) = sum_R_p1 e^{i * xq1 * R_p1} g(R_p)
  !!
  !! The degeneracy factor is absorbed into the definition of g_23.
  !!
  !--------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : inter_pool_comm
  USE io_global,        ONLY : stdout
  USE global_var,       ONLY : epmatwp, epmatwp_dist
  USE ep_constants,     ONLY : twopi, ci, czero
  USE parallelism,      ONLY : para_bounds
  USE io,               ONLY : rwepmatw
  USE input,            ONLY : epw_memdist
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nmodes
  !! Total number of modes
  INTEGER, INTENT(in) :: nrr_g
  !! Number of phononic WS points
  INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
  !! Coordinates of WS points
  INTEGER, INTENT(in) :: nbnd
  !! Number of bands
  INTEGER, INTENT(in) ::  nrr_k
  !! Number of electronic WS points
  REAL(KIND = DP), INTENT(in) :: xq1
  !! First component of the q-point vector for the Fourier transformation
  !
  ! Local variables
  LOGICAL :: found
  !! True if irvec(2:3, ir) is already in the list
  INTEGER :: ir
  !! Real space WS index
  INTEGER :: ir_23
  !! Real space WS index
  INTEGER :: ir_23_tmp
  !! Real space WS index
  INTEGER :: jr
  !! Real space WS index
  INTEGER :: irn
  !! Combined WS and atom index
  INTEGER :: ir_start
  !! Starting ir for this cores
  INTEGER :: ir_stop
  !! Ending ir for this pool
  INTEGER :: na
  !! Atom index
  INTEGER :: imode
  !! Number of modes
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: rdotk
  !! Exponential for the FT
  COMPLEX(KIND = DP) :: phase
  !! phase factor for FT
  INTEGER, ALLOCATABLE :: irvec_g_23_temp(:, :)
  !! Temporary storage for irvec_g_23
  !
  CALL start_clock('ephW2Bp_s23')
  !
  xq1_save = xq1
  !
  ! Project the real-space WS points to the yz plane
  ! nrr_g_23: Number of inequivalent WS points in the yz plane
  ! irvec_g_23: List of inequivalent WS points in the yz plane
  !
  nrr_g_23 = 0
  ALLOCATE(irvec_g_23_temp(2, nrr_g), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset23', 'Error allocating irvec_g_23_temp', 1)
  !
  DO ir = 1, nrr_g
    !
    IF (nrr_g_23 == 0) THEN
      nrr_g_23 = nrr_g_23 + 1
      irvec_g_23_temp(:, nrr_g_23) = irvec_g(2:3, ir)
    ELSE ! nrr_g_23 > 0
      !
      ! Find whether irvec_g(2:3, ir) already exists in irvec_g_23_temp
      found = .FALSE.
      DO jr = 1, nrr_g_23
        IF (ALL(irvec_g_23_temp(:, jr) == irvec_g(2:3, ir))) found = .TRUE.
      ENDDO
      !
      ! If not found, append irvec_g to irvec_g_23_temp
      IF (.NOT. found) THEN
        nrr_g_23 = nrr_g_23 + 1
        irvec_g_23_temp(:, nrr_g_23) = irvec_g(2:3, ir)
      ENDIF
      !
    ENDIF ! nrr_g_23
    !
  ENDDO ! ir
  !
  ! Copy irvec_g_23_temp to irvec_g_23 with the correct size
  !
  ALLOCATE(irvec_g_23(2, nrr_g_23), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset23', 'Error allocating irvec_g_23', 1)
  !
  irvec_g_23 = irvec_g_23_temp(:, 1:nrr_g_23)
  !
  DEALLOCATE(irvec_g_23_temp, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset23', 'Error deallocating irvec_g_23_temp', 1)
  !
  ALLOCATE(epmatwp_23(nbnd, nbnd, nrr_k, nmodes, nrr_g_23), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset23', 'Error allocating epmatwp_23', 1)
  epmatwp_23 = czero
  !
  ! Because nrr_g_23 can be quite small, we do a combined parallelization
  ! on WS vectors and modes.
  !
  CALL para_bounds(ir_start, ir_stop, nrr_g * nmodes)
  !
  DO irn = ir_start, ir_stop
    !
    ir = (irn - 1) / nmodes + 1
    imode = MOD(irn - 1, nmodes) + 1
    na = (imode - 1) / 3 + 1
    !
    ! Find ir_23 such that irvec_g_23(:, ir_23) == irvec_g(2:3, ir)
    !
    ir_23 = -1
    DO ir_23_tmp = 1, nrr_g_23
      IF (ALL(irvec_g_23(:, ir_23_tmp) == irvec_g(2:3, ir))) THEN
        ir_23 = ir_23_tmp
        EXIT
      ENDIF
    ENDDO
    !
    IF (ir_23 == -1) THEN
      WRITE(stdout, '(5x, a, I8, a, I8)') "ir_23 not found for irn = ", irn, ", ir = ", ir
      CALL errore("ephwan2blochp_reset23", "ir_23 not found", 1)
    ENDIF
    !
    ! Set Fourier transform factor phase = e^{i * xq1 * R_p1}
    !
    rdotk = twopi * xq1 * REAL(irvec_g(1, ir), KIND=DP)
    phase = EXP(ci * rdotk)
    !
    ! Fourier transform epmatwp along the third direction to get epmatwp_23
    !
    IF (epw_memdist) THEN
      CALL ZAXPY(nbnd * nbnd * nrr_k, phase, epmatwp_dist(:, :, :, irn - ir_start + 1), 1, epmatwp_23(:, :, :, imode, ir_23), 1)
    ELSE
      CALL ZAXPY(nbnd * nbnd * nrr_k, phase, epmatwp(:, :, :, imode, ir), 1, epmatwp_23(:, :, :, imode, ir_23), 1)
    ENDIF
    !
  ENDDO ! irn
  !
  CALL mp_sum(epmatwp_23, inter_pool_comm)
  !
  CALL stop_clock('ephW2Bp_s23')
  !
  !---------------------------------------------------------------------------
  END SUBROUTINE ephwan2blochp_reset23
  !---------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ephwan2blochp_reset3(xq2, nmodes, nbnd, nrr_k)
  !--------------------------------------------------------------------------
  !!
  !! Fourier transform along xq3 direction to compute epmat_23(xq1, R2, R3).
  !!
  !! If etf_mem == 1, write epmat_23(xq1, R2, R3) to file. The file is read
  !! in subroutine ephwan2blochp_get23.
  !!
  !! g_23(R_p2, R_p3) = sum_R_p1 e^{i * xq1 * R_p1} g(R_p)
  !!
  !! The degeneracy factor is absorbed into the definition of g_23.
  !!
  !--------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : inter_pool_comm
  USE ep_constants,     ONLY : twopi, ci, czero
  USE parallelism,      ONLY : para_bounds
  USE io,               ONLY : rwepmatw
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nmodes
  !! Total number of modes
  INTEGER, INTENT(in) :: nbnd
  !! Number of bands
  INTEGER, INTENT(in) ::  nrr_k
  !! Number of electronic WS points
  REAL(KIND = DP), INTENT(in) :: xq2
  !! y component of the q-point vector for the Fourier transformation
  !
  ! Local variables
  LOGICAL :: found
  !! True if irvec(3, ir) is already in the list
  INTEGER :: ir
  !! Real space WS index
  INTEGER :: ir_3
  !! Real space WS index
  INTEGER :: jr
  !! Real space WS index
  INTEGER :: irn
  !! Combined WS and atom index
  INTEGER :: imode
  !! Number of modes
  INTEGER :: irn3_start
  !! Starting irn3 for this core
  INTEGER :: irn3_stop
  !! Stoping irn3 for this core
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: rdotk
  !! Exponential for the FT
  COMPLEX(KIND = DP) :: cfac
  !! Factor for the FT
  INTEGER, ALLOCATABLE :: irvec_g_3_temp(:)
  !! Temporary storage for irvec_g_3
  COMPLEX(KIND = DP), ALLOCATABLE :: eptmp_3(:, :, :)
  !! Temporary matrix to store el-ph
  !
  CALL start_clock('ephW2Bp_s3')
  !
  ALLOCATE(eptmp_3(nbnd, nbnd, nrr_k), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset3', 'Error allocating eptmp_3', 1)
  !
  xq2_save = xq2
  !
  ! Project the real-space WS points to the z axis
  ! nrr_g_3: Number of inequivalent WS points in the z axis
  ! irvec_g_3: List of inequivalent WS points in the z axis
  !
  nrr_g_3 = 0
  ALLOCATE(irvec_g_3_temp(nrr_g_23), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset3', 'Error allocating irvec_g_3_temp', 1)
  !
  DO ir = 1, nrr_g_23
    !
    IF (nrr_g_3 == 0) THEN
      nrr_g_3 = nrr_g_3 + 1
      irvec_g_3_temp(nrr_g_3) = irvec_g_23(2, ir)
    ELSE ! nrr_g_3 > 0
      !
      ! Find whether irvec_g_23(2, ir) already exists in irvec_g_3_temp
      found = .FALSE.
      DO jr = 1, nrr_g_3
        IF (irvec_g_3_temp(jr) == irvec_g_23(2, ir)) found = .TRUE.
      ENDDO
      !
      ! If not found, append irvec_g to irvec_g_3_temp
      IF (.NOT. found) THEN
        nrr_g_3 = nrr_g_3 + 1
        irvec_g_3_temp(nrr_g_3) = irvec_g_23(2, ir)
      ENDIF
      !
    ENDIF ! nrr_g_3
    !
  ENDDO ! ir
  !
  ! Copy irvec_g_3_temp to irvec_g_3 with the correct size
  !
  ALLOCATE(irvec_g_3(nrr_g_3), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset3', 'Error allocating irvec_g_3', 1)
  !
  irvec_g_3 = irvec_g_3_temp(1:nrr_g_3)
  !
  DEALLOCATE(irvec_g_3_temp, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset3', 'Error deallocating irvec_g_3_temp', 1)
  !
  ALLOCATE(epmatwp_3(nbnd, nbnd, nrr_k, nmodes, nrr_g_3), STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset3', 'Error allocating epmatwp_3', ierr)
  epmatwp_3 = czero
  !
  ! Because nrr_g_3 can be quite small, we do a combined parallelization
  ! on WS vectors and modes.
  !
  CALL para_bounds(irn3_start, irn3_stop, nrr_g_3 * nmodes)
  !
  DO irn = irn3_start, irn3_stop
    !
    eptmp_3(:, :, :) = czero
    !
    ir_3 = (irn - 1) / nmodes + 1
    imode = MOD(irn - 1, nmodes) + 1
    !
    DO ir = 1, nrr_g_23
      !
      ! Loop over ir such that irvec_g_3(ir_23) == irvec_g_23(2, ir)
      IF (irvec_g_3(ir_3) /= irvec_g_23(2, ir)) CYCLE
      !
      ! Set Fourier transform factor cfac = e^{i * xq2 * R_p2}
      ! R_p2 = irvec_g_23(1, ir). The x component is already Fourier
      ! transformed in ephwan2blochp_reset23.
      !
      rdotk = twopi * xq2 * REAL(irvec_g_23(1, ir), KIND=DP)
      !
      cfac = EXP(ci * rdotk)
      !
      ! Fourier transform epmatwp along the third direction to get eptmp_23
      !
      CALL ZAXPY(nbnd * nbnd * nrr_k, cfac, epmatwp_23(:, :, :, imode, ir), 1, eptmp_3, 1)
      !
    ENDDO ! ir
    !
    ! Since eptmp_3 is not distributed, we should not call mp_sum here.
    !
    epmatwp_3(:, :, :, imode, ir_3) = eptmp_3
    !
  ENDDO ! irn
  !
  CALL mp_sum(epmatwp_3, inter_pool_comm)
  !
  DEALLOCATE(eptmp_3, STAT = ierr)
  IF (ierr /= 0) CALL errore('ephwan2blochp_reset3', 'Error deallocating eptmp_3', 1)
  !
  CALL stop_clock('ephW2Bp_s3')
  !
  !---------------------------------------------------------------------------
  END SUBROUTINE ephwan2blochp_reset3
  !---------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------
  SUBROUTINE wan2bloch_opt_finalize()
  !---------------------------------------------------------------------------
  !! Deallocate variables used for optimized Fourier transformation
  !---------------------------------------------------------------------------
  IMPLICIT NONE
  !
  INTEGER :: ierr
  !! Error status
  !
  IF (xq_is_initialized) THEN
    DEALLOCATE(irvec_g_23, STAT = ierr)
    IF (ierr /= 0) CALL errore('wan2bloch_opt_finalize', 'Error allocating irvec_g_23', 1)
    DEALLOCATE(epmatwp_23, STAT = ierr)
    IF (ierr /= 0) CALL errore('wan2bloch_opt_finalize', 'Error allocating epmatwp_23', 1)
    DEALLOCATE(irvec_g_3, STAT = ierr)
    IF (ierr /= 0) CALL errore('wan2bloch_opt_finalize', 'Error allocating irvec_g_3', 1)
    DEALLOCATE(epmatwp_3, STAT = ierr)
    IF (ierr /= 0) CALL errore('wan2bloch_opt_finalize', 'Error allocating epmatwp_3', 1)
  ENDIF
  !
  !---------------------------------------------------------------------------
  END SUBROUTINE wan2bloch_opt_finalize
  !---------------------------------------------------------------------------
  !
!--------------------------------------------------------------------------
END MODULE wannier2bloch_opt
!--------------------------------------------------------------------------
