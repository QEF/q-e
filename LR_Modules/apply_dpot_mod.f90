!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
MODULE apply_dpot_mod
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL :: is_allocated = .FALSE.
  !! Check if temporary storages are allocated
  COMPLEX(DP), ALLOCATABLE :: psi_r(:, :)
  !! Temporary storage for a real-space wavefunction
  COMPLEX(DP), ALLOCATABLE :: tg_dv(:, :)
  !! Task groups: temporary storage for potential * wfct
  COMPLEX(DP), ALLOCATABLE :: tg_psic(:, :)
  !! Task groups: temporary storage for wavefunctions
  !
  CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE apply_dpot_allocate()
    !! Allocate temporary storages
    USE kinds,             ONLY : DP
    USE fft_base,          ONLY : dffts
    USE noncollin_module,  ONLY : npol, nspin_mag
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (is_allocated) RETURN
    is_allocated = .TRUE.
    !
    ALLOCATE(psi_r(dffts%nnr, npol), STAT=ierr)
    IF (ierr /= 0) CALL errore('apply_dpot_allocate', 'Error allocating psi_r', 1)
    !
    !$acc enter data create(psi_r(1:dffts%nnr, 1:npol))
    !
    IF (dffts%has_task_groups) THEN
      ALLOCATE(tg_dv(dffts%nnr_tg, nspin_mag), STAT=ierr)
      IF (ierr /= 0) CALL errore('apply_dpot_allocate', 'Error allocating tg_dv', 1)
      ALLOCATE(tg_psic(dffts%nnr_tg, npol), STAT=ierr)
      IF (ierr /= 0) CALL errore('apply_dpot_allocate', 'Error allocating tg_psic', 1)
    ENDIF
    !
  END SUBROUTINE apply_dpot_allocate
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE apply_dpot_deallocate()
    !! Deallocate temporary storages
    USE kinds,             ONLY : DP
    USE fft_base,          ONLY : dffts
    USE noncollin_module,  ONLY : npol, nspin_mag
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (.NOT. is_allocated) RETURN
    is_allocated = .FALSE.
    !
    !$acc exit data delete(psi_r)
    !
    DEALLOCATE(psi_r, STAT=ierr)
    IF (ierr /= 0) CALL errore('apply_dpot_deallocate', 'Error deallocating psi_r', 1)
    !
    IF (dffts%has_task_groups) THEN
      DEALLOCATE(tg_dv, STAT=ierr)
      IF (ierr /= 0) CALL errore('apply_dpot_deallocate', 'Error deallocating tg_dv', 1)
      DEALLOCATE(tg_psic, STAT=ierr)
      IF (ierr /= 0) CALL errore('apply_dpot_deallocate', 'Error deallocating tg_psic', 1)
    ENDIF
    !
  END SUBROUTINE apply_dpot_deallocate
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE apply_dpot_bands(ik, nbnd, dv, psi, dvpsi)
    !--------------------------------------------------------------------------
    !! Calculate dvpsi = dv * psi in G space, for nbnd bands.
    !! 1. inverse FFT to real space
    !! 2. Multiply the potential with the wavefunctions
    !! 3. FFT back to G-space
    !--------------------------------------------------------------------------
    !
    USE kinds,             ONLY : DP
    USE fft_base,          ONLY : dffts
    USE wvfct,             ONLY : npwx
    USE noncollin_module,  ONLY : noncolin, domag, npol, nspin_mag
    USE lsda_mod,          ONLY : current_spin
    USE fft_helper_subroutines, ONLY : fftx_ntgrp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik
    !! k-point index.
    INTEGER, INTENT(IN) :: nbnd
    !! Number of bands to compute dvpsi.
    COMPLEX(KIND=DP), INTENT(IN) :: dv(dffts%nnr, nspin_mag)
    !! Potential in real space.
    COMPLEX(KIND=DP), INTENT(IN) :: psi(npwx*npol, nbnd)
    !! Wavefunction to be multiplied.
    COMPLEX(KIND=DP), INTENT(INOUT) :: dvpsi(npwx*npol, nbnd)
    !! Output. Wavefunction multiplied by potential.
    !
    INTEGER :: ibnd
    !! Counter for bands
    INTEGER :: ipol
    !! Counter for polarization
    INTEGER :: incr
    !! Step size for loop over bands.
    INTEGER :: tg_v_siz
    !! Task groups: size of the potential.
    !
    CALL start_clock("apply_dpot_b")
    !
    !$acc enter data copyin(psi)
    !$acc update device(dv(1:dffts%nnr, 1:nspin_mag))
    !
    IF (.NOT. is_allocated) CALL apply_dpot_allocate()
    !
    incr = 1
    !
    ! Setup for task groups
    !
    IF (dffts%has_task_groups) THEN
      tg_v_siz = dffts%nnr_tg
      incr = fftx_ntgrp(dffts)
      !
      IF (noncolin) THEN
        CALL tg_cgather(dffts, dv(:, 1), tg_dv(:, 1))
        IF (domag) THEN
          DO ipol = 2, 4
            CALL tg_cgather(dffts, dv(:, ipol), tg_dv(:, ipol))
          ENDDO
        ENDIF
      ELSE
        CALL tg_cgather(dffts, dv(:, current_spin), tg_dv(:,1))
      ENDIF ! noncolin
    ENDIF ! has_task_groups
    !
    !$acc kernels present(dvpsi)
    dvpsi = (0.0_DP, 0.0_DP)
    !$acc end kernels
    !
    DO ibnd = 1, nbnd, incr
      IF (dffts%has_task_groups) THEN
        CALL cft_wave_tg(ik, psi, tg_psic, 1, tg_v_siz, ibnd, nbnd)
        CALL apply_dpot(tg_v_siz, tg_psic, tg_dv, 1)
        CALL cft_wave_tg(ik, dvpsi, tg_psic, -1, tg_v_siz, ibnd, nbnd)
      ELSE
        CALL cft_wave(ik, psi(:, ibnd), psi_r, +1)
        CALL apply_dpot(dffts%nnr, psi_r, dv, current_spin)
        CALL cft_wave(ik, dvpsi(:, ibnd), psi_r, -1)
      ENDIF ! has_task_groups
    ENDDO ! ibnd
    !
    !$acc update self(dvpsi(1:npwx*npol, 1:nbnd))
    !$acc exit data delete(psi)
    !
    CALL stop_clock("apply_dpot_b")
    !
  END SUBROUTINE apply_dpot_bands
  !
END MODULE apply_dpot_mod
