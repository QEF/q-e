!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! TODO
! TODO When stress tensor is implemented, following subroutines must be modified:
! TODO 1). rism_module#rism_check -> to check condition
! TODO 2). solvation_stress_ion   -> to calculate short-range local potential
! TODO 3). solvation_esm_stress   -> to calculate long-range local potential
! TODO 4). lj_get_stress_x        -> to calculate Lennard-Jones potential
! TODO
!--------------------------------------------------------------------------
MODULE rism_module
  !--------------------------------------------------------------------------
  !
  ! ... this module is a wrapper of 1D- and 3D-RISM,
  ! ... to perform 3D-RISM-SCF (H.Sato et al., J. Chem. Phys. 2000, 112, 9463).
  ! ...
  ! ... if ESM(BC1) is applied, ESM-RISM is performed.
  !
  USE cell_base,        ONLY : at, omega, iforceh
  USE cellmd,           ONLY : at_old, lmovecell
  USE check_stop,       ONLY : stopped_by_user
  USE constants,        ONLY : eps14
  USE control_flags,    ONLY : gamma_only, tr2, tstress
  USE esm,              ONLY : do_comp_esm, esm_bc
  USE exx_base,         ONLY : x_gamma_extrapolation
  USE fft_base,         ONLY : dfftp
  USE fft_interfaces,   ONLY : invfft
  USE xc_lib,           ONLY : exx_is_active
  USE gvect,            ONLY : ngm, gstart
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE ions_base,        ONLY : nat, ityp, zv, tau
  USE kinds,            ONLY : DP
  USE klist,            ONLY : nkstot, xk, nelec
  USE lsda_mod,         ONLY : lsda, nspin
  USE mp,               ONLY : mp_sum, mp_bcast
  USE mp_images,        ONLY : intra_image_comm
  USE noncollin_module, ONLY : nspin_lsda
  USE relax,            ONLY : starting_scf_threshold
  USE rism1d_facade,    ONLY : lrism1d, rism1d_finalize, rism1d_is_avail, &
                             & rism1d_iosys, rism1d_summary, rism1d_prepare, rism1d_run, &
                             & rism1d_write_to_restart, rism1d_write_to_show, rism1d_print_clock, &
                             & starting_1d => starting_corr
  USE rism3d_facade,    ONLY : lrism3d, epsv, starting_epsv, &
                             & rism3t, rism3d_initialize, rism3d_finalize, rism3d_is_laue, &
                             & rism3d_iosys, rism3d_summary, rism3d_prepare, rism3d_reprepare, &
                             & rism3d_run, rism3d_update_solute, rism3d_potential, &
                             & rism3d_force, rism3d_stress, rism3d_printpot, &
                             & rism3d_print_clock, starting_3d => starting_corr
  USE solute,           ONLY : deallocate_solU
  USE solvmol,          ONLY : deallocate_solVs
  USE vlocal,           ONLY : vloc
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define variables
  LOGICAL               :: lrism    = .FALSE.  ! to calculate RISM, or not
  LOGICAL               :: llaue    = .FALSE.  ! to apply Laue-RISM, or not
  REAL(DP), ALLOCATABLE :: vltot(:)            ! local potential in realspace
  !
  ! ... public components
  PUBLIC :: lrism
  PUBLIC :: rism_check
  PUBLIC :: deallocate_rism
  PUBLIC :: rism_iosys
  PUBLIC :: rism_calc1d
  PUBLIC :: rism_alloc3d
  PUBLIC :: rism_init3d
  PUBLIC :: rism_reinit3d
  PUBLIC :: rism_update_pos
  PUBLIC :: rism_calc3d
  PUBLIC :: rism_pot3d
  PUBLIC :: force_rism
  PUBLIC :: stres_rism
  PUBLIC :: rism_tobe_alive
  PUBLIC :: rism_set_restart
  PUBLIC :: rism_setlocal
  PUBLIC :: rism_new_conv_thr
  PUBLIC :: rism_printpot
  PUBLIC :: rism_print_clock
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_check()
    !----------------------------------------------------------------------------
    !
    ! ... check conditions,
    ! ... and stop program if incorrect conditions.
    !
    IMPLICIT NONE
    !
    INTEGER :: ia
    INTEGER :: ik, nk_
    !
    ! ...
    ! ... check ESM for Laue-RISM
    ! ...
    IF (do_comp_esm) THEN
      !
      ! ... BC1 only
      IF (TRIM(esm_bc) /= 'bc1' .AND. TRIM(esm_bc) /= 'pbc') THEN
        CALL errore('rism_check', 'Laue-RISM only supports ESM-BC1', 1)
      END IF
      !
      ! ... correct cell shape ?
      IF (ABS(at(1, 3)) > eps14 .OR. ABS(at(3, 1)) > eps14 .OR. &
        & ABS(at(2, 3)) > eps14 .OR. ABS(at(3, 2)) > eps14) THEN
        CALL errore('rism_check', 'incorrect unit cell for Laue-RISM', 1)
      END IF
      !
      ! ... correct atomic positions ?
      DO ia = 1, nat
        IF (tau(3, ia) <= (-0.5_DP * at(3, 3)) .OR. (0.5_DP * at(3, 3)) <= tau(3, ia)) THEN
          CALL errore('rism_check', 'incorrect atomic position for Laue-RISM', ia)
        END IF
      END DO
      !
      ! ... correct k-points ?
      IF (lsda) THEN
        nk_ = nkstot / 2
      ELSE
        nk_ = nkstot
      END IF
      !
      DO ik = 1, nk_
        IF (ABS(xk(3, ik)) > eps14) THEN
          CALL errore('rism_check', 'incorrect k-point for Laue-RISM', ik)
        END IF
      END DO
      !
      ! ... correct Vexx(G=0) ?
      IF (exx_is_active() .AND. (.NOT. x_gamma_extrapolation)) THEN
        CALL errore('rism_check', 'Laue-RISM requires Vexx(G=0)', 1)
      END IF
      !
    END IF
    !
    ! ...
    ! ... check Variable Cell
    ! ...
    IF (llaue) THEN
      !
#if defined (__RISM_STRESS)
      ! ... Laue-RISM only supports 2Dxy
      IF (lmovecell) THEN
        IF (iforceh(3, 1) /= 0 .OR. iforceh(3, 2) /= 0 .OR. iforceh(3, 3) /= 0 .OR. &
          & iforceh(1, 3) /= 0 .OR. iforceh(2, 3) /= 0) THEN
          CALL errore('rism_check', 'Laue-RISM only supports cell_dofree = "2Dxy"', 1)
        END IF
      END IF
#else
      ! ... Laue-RISM does not support storess tensor
      IF (tstress) THEN
        CALL errore('rism_check', 'Laue-RISM does not support stress tensor', 1)
      END IF
      !
      ! ... Laue-RISM does not support Variable Cell
      IF (lmovecell) THEN
        CALL errore('rism_check', 'Laue-RISM does not support variable cell', 1)
      END IF
#endif
      !
    ELSE
      !
      ! ... 3D-RISM does not support storess tensor
      IF (tstress) THEN
        CALL errore('rism_check', '3D-RISM does not support stress tensor', 1)
      END IF
      !
      ! ... 3D-RISM does not support Variable Cell
      IF (lmovecell) THEN
        CALL errore('rism_check', '3D-RISM does not support variable cell', 1)
      END IF
      !
    END IF
    !
  END SUBROUTINE rism_check
  !
  !----------------------------------------------------------------------------
  SUBROUTINE deallocate_rism(lall)
    !----------------------------------------------------------------------------
    !
    ! ... if lall=.TRUE., deallocate all data.
    ! ... if lall=.FALSE., deallocate data, which depend on FFT box.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lall
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (ALLOCATED(vltot)) THEN
      DEALLOCATE(vltot)
    END IF
    !
    IF (lall) THEN
      CALL rism1d_finalize()
      CALL deallocate_solVs()
    END IF
    !
    CALL rism3d_finalize(lall)
    CALL deallocate_solU(lall)
    !
  END SUBROUTINE deallocate_rism
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_iosys(trism)
    !----------------------------------------------------------------------------
    !
    ! ... set variables of 1D- and 3D-RISM,
    ! ... and allocate data, which are independent of FFT box.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: trism
    !
    lrism = trism
    llaue = (do_comp_esm .AND. (TRIM(esm_bc) == 'bc1'))
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    CALL rism1d_iosys(lrism, llaue)
    CALL rism3d_iosys(lrism, llaue, init=.FALSE.)
    !
  END SUBROUTINE rism_iosys
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_calc1d(lmust)
    !----------------------------------------------------------------------------
    !
    ! ... calculate 1D-RISM.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN), OPTIONAL :: lmust
    !
    LOGICAL :: lmust_
    LOGICAL :: conv_rism1d
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (.NOT. lrism1d) THEN
      CALL errore('rism_calc1d', '1D-RISM is not ready', 1)
    END IF
    !
    CALL rism_check()
    !
    lmust_ = .FALSE.
    IF (PRESENT(lmust)) THEN
      lmust_ = lmust
    END IF
    !
    IF ((.NOT. lmust_) .AND. rism1d_is_avail()) THEN
      CALL rism1d_write_to_restart()
      RETURN
    END IF
    !
    CALL rism1d_summary()
    CALL rism1d_prepare()
    CALL rism1d_run(conv_rism1d)
    CALL rism1d_write_to_restart()
    IF (conv_rism1d) THEN
      CALL rism1d_write_to_show()
    END IF
    !
    IF (.NOT. rism1d_is_avail()) THEN
      CALL errore('rism_calc1d', 'result of 1D-RISM calculation is not avairable', 1)
    END IF
    !
  END SUBROUTINE rism_calc1d
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_alloc3d()
    !----------------------------------------------------------------------------
    !
    ! ... allocate memory for 3D-RISM.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (.NOT. lrism3d) THEN
      CALL errore('rism_init3d', '3D-RISM is not ready', 1)
    END IF
    !
    CALL rism_check()
    !
    CALL rism3d_initialize(llaue)
    CALL rism3d_summary()
    !
  END SUBROUTINE rism_alloc3d
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_init3d()
    !----------------------------------------------------------------------------
    !
    ! ... initialize 3D-RISM.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (.NOT. lrism3d) THEN
      CALL errore('rism_init3d', '3D-RISM is not ready', 1)
    END IF
    !
    CALL rism_check()
    !
    CALL rism3d_prepare()
    !
  END SUBROUTINE rism_init3d
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_reinit3d()
    !----------------------------------------------------------------------------
    !
    ! ... re-initialize 3D-RISM for Variable Cell.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (.NOT. lrism3d) THEN
      CALL errore('rism_reinit3d', '3D-RISM is not ready', 1)
    END IF
    !
    CALL rism_check()
    !
    CALL rism3d_reprepare(at_old)
    !
  END SUBROUTINE rism_reinit3d
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_update_pos()
    !----------------------------------------------------------------------------
    !
    ! ... update positions of solute system.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    CALL rism_check()
    !
    CALL rism3d_update_solute()
    !
  END SUBROUTINE rism_update_pos
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_calc3d(rhog, esol, vsol, vr, dr2)
    !----------------------------------------------------------------------------
    !
    ! ... calculate 3D-RISM.
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN)    :: rhog(ngm)
    REAL(DP),    INTENT(OUT)   :: esol
    REAL(DP),    INTENT(OUT)   :: vsol
    REAL(DP),    INTENT(INOUT) :: vr(dfftp%nnr, nspin)
    REAL(DP),    INTENT(IN)    :: dr2  ! electronic SCF's error
    !
    INTEGER               :: is
    REAL(DP)              :: epsv_new
    REAL(DP)              :: tr2_
    REAL(DP), ALLOCATABLE :: vpot(:)
    LOGICAL               :: conv_rism3d
    !
    REAL(DP), PARAMETER   :: TR2_SCALE  = 10.0_DP
    REAL(DP), PARAMETER   :: TR2_EXPON  = 0.55_DP
    REAL(DP), PARAMETER   :: EPSV_EXPON = 0.50_DP
    REAL(DP), PARAMETER   :: EPSV_CEILING = 1.0E-2_DP
    !
    IF (.NOT. lrism) THEN
      esol = 0.0_DP
      vsol = 0.0_DP
      RETURN
    END IF
    !
    IF (.NOT. lrism3d) THEN
      CALL errore('rism_calc3d', '3D-RISM is not ready', 1)
    END IF
    !
    CALL rism_check()
    !
    ALLOCATE(vpot(dfftp%nnr))
    !
    ! ... make electronic potential and charge
    CALL solute_pot(rhog, vpot)
    !
    ! ... set convergence threshold
    IF (ionode) THEN
      epsv_new = epsv
      tr2_ = tr2 * DBLE(nelec) / TR2_SCALE
      tr2_ = tr2_ ** TR2_EXPON
      !
      IF (0.0_DP < epsv .AND. 0.0_DP < tr2_) THEN
        IF (0.0_DP < dr2) THEN
          ! ... dr2 is positive
          IF (dr2 >= tr2_) THEN
            epsv_new = 10.0_DP ** (LOG10(epsv) * LOG10(dr2) / LOG10(tr2_))
            epsv_new = MAX(epsv_new, epsv)
            epsv_new = MIN(epsv_new, MAX(EPSV_CEILING, epsv ** EPSV_EXPON))
          END IF
          !
        ELSE
          ! ... dr2 is negative or zero
          epsv_new = MAX(EPSV_CEILING, epsv ** EPSV_EXPON)
        END IF
      END IF
    END IF
    !
    CALL mp_bcast(epsv_new, ionode_id, intra_image_comm)
    !
    ! ... run 3D-RISM calculation
    CALL rism3d_run(vpot, rhog, conv_rism3d, epsv_new)
    !
    IF (.NOT. conv_rism3d) THEN
      stopped_by_user = .TRUE.
      esol = 0.0_DP
      vsol = 0.0_DP
      GOTO 10
    END IF
    !
    ! ... make solvation potential and energy
    CALL solvation_pot(vpot)
    CALL solvation_erg(esol, vsol, rhog)
    !
    DO is = 1, nspin_lsda
      vr(:, is) = vr(:, is) + vpot(:)
    END DO
    !
 10 CONTINUE
    DEALLOCATE(vpot)
    !
  END SUBROUTINE rism_calc3d
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_pot3d(rhog, vr)
    !----------------------------------------------------------------------------
    !
    ! ... create solvation potential of 3D-RISM, if possible.
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN)    :: rhog(ngm)
    REAL(DP),    INTENT(INOUT) :: vr(dfftp%nnr, nspin)
    !
    INTEGER               :: is
    REAL(DP), ALLOCATABLE :: vpot(:)
    !
    ! ... if 3D-RISM's data are kept,
    ! ... one can calculate solvation potential.
    IF (.NOT. lrism3d) THEN
      CALL errore('rism_pot3d', '3D-RISM is not ready', 1)
    END IF
    !
    ALLOCATE(vpot(dfftp%nnr))
    !
    CALL solute_pot(rhog, vpot)
    !
    CALL rism3d_potential(vpot, rhog)
    !
    CALL solvation_pot(vpot)
    !
    DO is = 1, nspin_lsda
      vr(:, is) = vr(:, is) + vpot(:)
    END DO
    !
    DEALLOCATE(vpot)
    !
  END SUBROUTINE rism_pot3d
  !
  !----------------------------------------------------------------------------
  SUBROUTINE solute_pot(rhog, vpot)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN)  :: rhog(ngm)
    REAL(DP),    INTENT(OUT) :: vpot(dfftp%nnr)
    !
    INTEGER               :: is
    REAL(DP)              :: ehart
    REAL(DP)              :: charge
    REAL(DP), ALLOCATABLE :: vh(:,:)
    !
    IF (.NOT. ALLOCATED(vltot)) THEN
      CALL errore('solute_pot', 'vltot is null', 1)
    END IF
    !
    ALLOCATE(vh(dfftp%nnr, nspin))
    !
    vh = 0.0_DP
    CALL v_h_without_esm(rhog, ehart, charge, vh)
    !
    vpot(1:dfftp%nnr) = vltot(1:dfftp%nnr)
    !
    DO is = 1, nspin_lsda
      vpot(:) = vpot(:) + vh(:, is) / DBLE(nspin_lsda)
    END DO
    !
    DEALLOCATE(vh)
    !
  END SUBROUTINE solute_pot
  !
  !----------------------------------------------------------------------------
  SUBROUTINE solvation_pot(vpot)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: vpot(dfftp%nnr)
    !
    INTEGER                  :: ig
    INTEGER                  :: ir
    COMPLEX(DP), ALLOCATABLE :: aux(:)
    !
    CALL start_clock('3DRISM_vsol')
    !
    ALLOCATE(aux(dfftp%nnr))
    !
    aux = CMPLX(0.0_DP, 0.0_DP, kind=DP)
    !
    IF (llaue) THEN
!$omp parallel do default(shared) private(ig)
      DO ig = 1, rism3t%gvec%ngm
        aux(dfftp%nl(ig)) = rism3t%vpot_pbc(ig)
      END DO
!$omp end parallel do
    ELSE
!$omp parallel do default(shared) private(ig)
      DO ig = 1, rism3t%gvec%ngm
        aux(dfftp%nl(ig)) = rism3t%vpot(ig)
      END DO
!$omp end parallel do
    END IF
    !
    IF (gamma_only) THEN
!$omp parallel do default(shared) private(ig)
      DO ig = 1, rism3t%gvec%ngm
        aux(dfftp%nlm(ig)) = CONJG(aux(dfftp%nl(ig)))
      END DO
!$omp end parallel do
    END IF
    !
    CALL invfft('Rho', aux, dfftp)
    !
!$omp parallel do default(shared) private(ir)
    DO ir = 1, dfftp%nnr
      vpot(ir) = -DBLE(aux(ir))  ! electron has negative charge
    END DO
!$omp end parallel do
    !
    DEALLOCATE(aux)
    !
    CALL stop_clock('3DRISM_vsol')
    !
  END SUBROUTINE solvation_pot
  !
  !----------------------------------------------------------------------------
  SUBROUTINE solvation_erg(esol, vsol, rhog)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP),    INTENT(OUT) :: esol
    REAL(DP),    INTENT(OUT) :: vsol
    COMPLEX(DP), INTENT(IN)  :: rhog(ngm)
    !
    INTEGER  :: is
    REAL(DP) :: rho0
    REAL(DP) :: ionic_charge
    !
    ! ... solvation energy
    esol = rism3t%esol
    !
    IF (llaue) THEN
       ! ... potential shifting energy (for Laue-RISM)
       rho0 = 0.0_DP
       !
       IF (gstart > 1) THEN
         rho0 = rhog(1)
       END IF
       !
       CALL mp_sum(rho0, rism3t%mp_site%intra_sitg_comm)
       !
       ionic_charge = SUM(zv(ityp(1:nat)))
       !
       vsol = -0.5_DP * rism3t%vsol * (ionic_charge - rho0 * omega)
       !
    ELSE
       ! ... zero (for 3D-RISM)
       vsol = 0.0_DP
    END IF
    !
  END SUBROUTINE solvation_erg
  !
  !----------------------------------------------------------------------------
  SUBROUTINE force_rism(forcesol)
    !----------------------------------------------------------------------------
    !
    ! ... calculate force from solvent.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: forcesol(3, nat)
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (.NOT. lrism3d) THEN
      CALL errore('force_rism', '3D-RISM is not ready', 1)
    END IF
    !
    IF (.NOT. rism3t%avail) THEN
      CALL errore('force_rism', 'result of 3D-RISM calculation is not avairable', 1)
    END IF
    !
    CALL rism_check()
    !
    forcesol = 0.0_DP
    CALL rism3d_force(forcesol, vloc)
    !
  END SUBROUTINE force_rism
  !
  !----------------------------------------------------------------------------
  SUBROUTINE stres_rism(sigmasol)
    !----------------------------------------------------------------------------
    !
    ! ... calculate stress tensor from solvent.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigmasol(3, 3)
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (.NOT. lrism3d) THEN
      CALL errore('stres_rism', '3D-RISM is not ready', 1)
    END IF
    !
    IF (.NOT. rism3t%avail) THEN
      CALL errore('stres_rism', 'result of 3D-RISM calculation is not avairable', 1)
    END IF
    !
    IF (.NOT. llaue) THEN
      CALL errore('stres_rism', 'you cannot calculate stress tensor of 3D-RISM', 1)
    END IF
    !
    CALL rism_check()
    !
    sigmasol = 0.0_DP
    CALL rism3d_stress(sigmasol)
    !
  END SUBROUTINE stres_rism
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_tobe_alive()
    !----------------------------------------------------------------------------
    !
    ! ... restore status from rism3d_facade.
    !
    IMPLICIT NONE
    !
    lrism = lrism3d
    llaue = rism3d_is_laue()
    !
  END SUBROUTINE rism_tobe_alive
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_set_restart()
    !----------------------------------------------------------------------------
    !
    ! ... set restart-mode.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    starting_1d = 'fix'
    starting_3d = 'file'
    !
  END SUBROUTINE rism_set_restart
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_setlocal(vltot_)
    !----------------------------------------------------------------------------
    !
    ! ... set local potential of realspace (vltot).
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: vltot_(dfftp%nnr)
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (ALLOCATED(vltot)) THEN
      DEALLOCATE(vltot)
    END IF
    !
    ALLOCATE(vltot(dfftp%nnr))
    !
    vltot = vltot_
    !
  END SUBROUTINE rism_setlocal
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_new_conv_thr()
    !----------------------------------------------------------------------------
    !
    ! ... update convergence threshold.
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: TR2_EXPON = 0.50_DP * LOG10(2.0_DP)
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      IF (epsv > 0.0_DP .AND. starting_epsv          > 0.0_DP .AND. &
          tr2  > 0.0_DP .AND. starting_scf_threshold > 0.0_DP) THEN
        epsv = starting_epsv * ((tr2 / starting_scf_threshold) ** TR2_EXPON)
      ELSE
        epsv = starting_epsv
      END IF
    END IF
    !
    CALL mp_bcast(epsv, ionode_id, intra_image_comm)
    !
  END SUBROUTINE rism_new_conv_thr
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_printpot()
    !----------------------------------------------------------------------------
    !
    ! ... print out planar averaged potentials and distributions
    ! ... of 3D-RISM or Laue-RISM.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism) THEN
      RETURN
    END IF
    !
    IF (.NOT. lrism3d) THEN
      CALL errore('rism_printpot', '3D-RISM is not ready', 1)
    END IF
    !
    IF (.NOT. rism3t%avail) THEN
      CALL errore('rism_printpot', 'result of 3D-RISM calculation is not avairable', 1)
    END IF
    !
    CALL rism3d_printpot()
    !
  END SUBROUTINE rism_printpot
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism_print_clock()
    !----------------------------------------------------------------------------
    !
    ! ... print clock for 1D- and 3D-RISM.
    !
    IMPLICIT NONE
    !
    IF (lrism1d .OR. lrism3d .OR. lrism) THEN
      WRITE(stdout, '(/,5X,"RISM routines")')
    END IF
    !
    IF (lrism1d) THEN
      CALL rism1d_print_clock()
    END IF
    !
    IF (lrism3d) THEN
      CALL rism3d_print_clock()
    END IF
    !
    IF (lrism) THEN
      CALL print_clock('3DRISM_vsol')
    END IF
    !
  END SUBROUTINE rism_print_clock
  !
END MODULE rism_module

