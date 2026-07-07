!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
MODULE response_kernels
CONTAINS
SUBROUTINE sternheimer_kernel(first_iter, time_reversed, npert, lrdvpsi, iudvpsi, &
         thresh, dvscfins, all_conv, avg_iter, drhoout, dbecsum, dbecsum_nc, &
         exclude_hubbard)
   !----------------------------------------------------------------------------
   !! Compute the density response to the perturbation dV = dV_bare + dV_ind by the
   !! non-interacting susceptibility. Solve Sternheimer equation
   !! (H - e) * dpsi = dvpsi = -P_c+ * (dV_bare + dV_ind) * psi.
   !!
   !! Perfoms the following tasks:
   !!  a) reads the bare potential term Delta V | psi > from buffer iudvpsi.
   !!  b) adds to it the screening term Delta V_{SCF} | psi >.
   !!     If lda_plus_u=.true. compute also the SCF part
   !!     of the response Hubbard potential.
   !!  c) applies P_c^+ (orthogonalization to valence states).
   !!  d) calls cgsolve_all to solve the linear system.
   !!  e) returns the Delta rho, and if lda_plus_u=.true. also return dbecsum
   !!
   !! dV_bare * psi is read from buffer iudvpsi, so they must be already calculated.
   !! dV_ind is given by input dvscfins, and dV_ind * psi is calculated in apply_dpot_bands.
   !!
   !! For USPPs, adddvscf is called, so relevant arrays must be already calculated.
   !! For DFT+U, adddvhubscf is called, so relevant arrays must be already calculated.
   !!
   !! Results are added to drhoout, dbecsum, dbecsum_nc, so these output arrays should
   !! be initialized before calling this subroutine.
   !!
   !! Input:
   !!    - first_iter: true if the first iteration, where dvscfins = 0
   !!    - time_reversed: false for ordinary cases, true for time-reversed wave functions
   !!                     which is need in the noncolliner magnetic case
   !!    - npert: number of perturbations
   !!    - lrdvpsi: record length for the buffer storing dV_bare * psi
   !!    - iudvpsi: unit for the buffer storing dV_bare * psi
   !!    - thresh: threshold for solving Sternheimer equation
   !!    - dvscfins: dV_ind calculated in the previous iteration
   !!    - exclude_hubbard: If TRUE, do not add the response of the Hubbard potential.
   !!                       Used in hp.x (Default: FALSE)
   !!
   !! Output:
   !!    - avg_iter: average number of iterations for the linear equation solver
   !!    - drhoout: induced charge density (dffts, without augmentation term)
   !!    - dbecsum: becsum with dpsi
   !!    - dbecsum_nc: becsum with dpsi. Optional, used if noncolin is true.
   !!
   !----------------------------------------------------------------------------
   USE kinds,                 ONLY : DP
   USE io_global,             ONLY : stdout
   USE mp,                    ONLY : mp_sum
   USE mp_pools,              ONLY : inter_pool_comm
   USE buffers,               ONLY : get_buffer, save_buffer
   USE fft_base,              ONLY : dffts
   USE ions_base,             ONLY : nat
   USE klist,                 ONLY : xk, wk, ngk, igk_k
   USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                 ONLY : nbnd, npwx, et
   USE wavefunctions,         ONLY : evc
   USE noncollin_module,      ONLY : noncolin, domag, npol, nspin_mag
   USE uspp,                  ONLY : vkb
   USE uspp_param,            ONLY : nhm
   USE uspp_init,             ONLY : init_us_2
   USE ldaU,                  ONLY : lda_plus_u
   USE units_lr,              ONLY : iuwfc, lrwfc, lrdwf, iudwf
   USE control_lr,            ONLY : nbnd_occ, lgamma
   USE qpoint,                ONLY : nksq, ikks, ikqs
   USE qpoint_aux,            ONLY : ikmks, ikmkmqs, becpt
   USE eqv,                   ONLY : dpsi, dvpsi, evq
   USE apply_dpot_mod,        ONLY : apply_dpot_bands
   USE lr_nc_mag,             ONLY : lr_apply_time_reversal
   USE incdrhoscf_mod,        ONLY : incdrhoscf, incdrhoscf_nc
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: first_iter
   !! true if the first iteration.
   LOGICAL, INTENT(IN) :: time_reversed
   !! true if solving for time reversed wave functions
   LOGICAL, INTENT(IN), OPTIONAL :: exclude_hubbard
   !! true if ignoring the Hubbard response term
   INTEGER, INTENT(IN) :: npert
   !! number of perturbations
   INTEGER, INTENT(IN) :: lrdvpsi
   !! record length for the buffer storing dV_bare * psi
   INTEGER, INTENT(IN) :: iudvpsi
   !! unit for the buffer storing dV_bare * psi
   REAL(DP), INTENT(IN) :: thresh
   !! threshold for solving the linear equation
   LOGICAL, INTENT(OUT) :: all_conv
   !! True if converged at all k points and perturbations
   REAL(DP), INTENT(OUT) :: avg_iter
   !! average number of iterations for the linear equation solver
   COMPLEX(DP), INTENT(INOUT) :: dvscfins(dffts%nnr, nspin_mag, npert)
   !! dV_ind calculated in the previous iteration
   COMPLEX(DP), INTENT(INOUT) :: drhoout(dffts%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   !
   LOGICAL :: conv_root
   !! true if linear system is converged
   LOGICAL :: exclude_hubbard_
   !! Local variable to set the default of exclude_hubbard to false
   INTEGER :: ikk, ikq, npw, npwq, ipert, num_iter, ik, nrec, ikmk, ikmkmq
   !! counters
   INTEGER :: tot_num_iter
   !! total number of iterations in cgsolve_all
   INTEGER :: tot_cg_calls
   !! total number of cgsolve_all calls
   REAL(DP) :: anorm
   !! the norm of the error of the conjugate gradient solution
   REAL(DP) :: rsign
   !! sign of the term in the magnetization
   REAL(DP), ALLOCATABLE :: h_diag(:, :)
   !! diagonal part of the Hamiltonian, used for preconditioning
   COMPLEX(DP) , ALLOCATABLE :: aux2(:, :)
   !! temporary storage used in apply_dpot_bands
   !
   EXTERNAL ch_psi_all, cg_psi
   !! functions passed to cgsolve_all
   !
   ! Initialization
   !
   CALL start_clock("sth_kernel")
   !
   exclude_hubbard_ = .FALSE.
   IF (PRESENT(exclude_hubbard)) exclude_hubbard_ = exclude_hubbard
   !
   !  change the sign of the magnetic field if required
   !
   IF (time_reversed) CALL lr_apply_time_reversal(first_iter, 2, dvscfins)
   !
   ALLOCATE(h_diag(npwx*npol, nbnd))
   ALLOCATE(aux2(npwx*npol, nbnd))
   h_diag = 0.d0
   aux2 = (0.d0, 0.d0)
   !
   !$acc enter data create(aux2(1:npwx*npol, 1:nbnd))
   !
   all_conv = .TRUE.
   tot_num_iter = 0
   tot_cg_calls = 0
   !
   DO ik = 1, nksq
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)
      !
      ! Set time-reversed k and k+q points
      !
      IF (time_reversed) THEN
         ikmk = ikmks(ik)
         ikmkmq = ikmkmqs(ik)
         rsign = -1.0_DP
      ELSE
         ikmk = ikk
         ikmkmq = ikq
         rsign = 1.0_DP
      ENDIF
      !
      IF (lsda) current_spin = isk(ikk)
      !
      ! reads unperturbed wavefunctions psi_k in G_space, for all bands
      ! if q=0, evq is a pointer to evc
      !
      IF (nksq > 1 .OR. (noncolin .AND. domag)) THEN
         IF (lgamma) THEN
            !civn: in this case evq is a pointer to evc
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
            !$acc update device(evc)
         ELSE
            !civn: in this case evq is allocated separately and needs to be updated on device
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
            !$acc update device(evc)
            CALL get_buffer(evq, lrwfc, iuwfc, ikmkmq)
            !$acc update device(evq)
         ENDIF
      ENDIF
      !
      ! compute beta functions and kinetic energy for k-point ik
      ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
      !
      CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
      !$acc update host(vkb)
      !
      ! compute the kinetic energy g2kin: (k+q+G)^2
      !
      CALL g2_kin(ikq)
      !
      ! compute preconditioning matrix h_diag used by cgsolve_all
      !
      CALL h_prec(ik, evq, h_diag)
      !
      DO ipert = 1, npert
         !
         ! read P_c^+ x psi_kpoint into dvpsi.
         !
         nrec = (ipert - 1) * nksq + ik
         IF (time_reversed) nrec = nrec + npert * nksq
         !
         CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrec)
         !
         IF (.NOT. first_iter) THEN
            !
            ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
            ! dvscf_q from previous iteration (mix_potential)
            !
            CALL apply_dpot_bands(ik, nbnd_occ(ikk), dvscfins(:, :, ipert), evc, aux2)
            dvpsi = dvpsi + aux2
            !
            !  In the case of US pseudopotentials there is an additional
            !  self-consistent term which comes from the dependence of D on
            !  V_{eff} on the bare change of the potential
            !
            CALL adddvscf(ipert, ik, time_reversed)
            !
            ! DFPT+U: add to dvpsi the scf part of the response
            ! Hubbard potential dV_hub
            !
            IF (lda_plus_u .AND. (.NOT. exclude_hubbard_)) CALL adddvhubscf(ipert, ik)
            !
         ENDIF ! .NOT. first_iter
         !
         ! Orthogonalize dvpsi to valence states
         !
         CALL orthogonalize(dvpsi, evq, ikmk, ikmkmq, dpsi, npwq, .FALSE.)
         !
         ! Initial guess for dpsi
         !
         IF (first_iter) THEN
            !
            !  At the first iteration dpsi is set to zero
            !
            dpsi(:, :) = (0.d0,0.d0)
         ELSE
            !
            ! starting value for delta_psi is read from iudwf
            !
            CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
         ENDIF
         !
         ! iterative solution of the linear system (H-e)*dpsi=dvpsi
         ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
         !
         conv_root = .TRUE.
         !
         ! TODO: should nbnd_occ(ikk) be nbnd_occ(ikmk)?
         CALL cgsolve_all(ch_psi_all, cg_psi, et(1, ikmk), dvpsi, dpsi, h_diag, &
            npwx, npwq, thresh, ik, num_iter, conv_root, anorm, nbnd_occ(ikk), npol)
         !
         tot_num_iter = tot_num_iter + num_iter
         tot_cg_calls = tot_cg_calls + 1
         !
         IF (.NOT. conv_root) THEN
            all_conv = .FALSE.
            WRITE( stdout, "(5x, 'kpoint', i4, ' sternheimer_kernel: &
               &root not converged, thresh < ', es10.3)") ik, anorm
         ENDIF
         !
         ! writes delta_psi on iunit iudwf, k=kpoint,
         !
         CALL save_buffer(dpsi, lrdwf, iudwf, nrec)
         !
         ! calculates dvscf, sum over k => dvscf_q_ipert
         !
         IF (noncolin) THEN
            CALL incdrhoscf_nc(drhoout(1,1,ipert), wk(ikk), ik, &
                               dbecsum_nc(1,1,1,1,ipert), dpsi, rsign)
         ELSE
            CALL incdrhoscf(drhoout(1,current_spin,ipert), wk(ikk), &
                            ik, dbecsum(1,1,current_spin,ipert), dpsi)
         ENDIF
      ENDDO ! ipert
   ENDDO ! ik
   !
   CALL mp_sum(tot_num_iter, inter_pool_comm)
   CALL mp_sum(tot_cg_calls, inter_pool_comm)
   avg_iter = REAL(tot_num_iter, DP) / REAL(tot_cg_calls, DP)
   !
   !$acc exit data delete(aux2)
   !
   !  reset the original magnetic field if it was changed
   !
   IF (time_reversed) CALL lr_apply_time_reversal(first_iter, 1, dvscfins)
   !
   DEALLOCATE(aux2)
   DEALLOCATE(h_diag)
   !
   CALL stop_clock("sth_kernel")
   !
!----------------------------------------------------------------------------
END SUBROUTINE sternheimer_kernel
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE sternheimer_kernel_freq(first_iter, time_reversed, npert, lrdvpsi, iudvpsi, &
         thresh, omega, dvscfins, all_conv, avg_iter, drhoout, dbecsum, dbecsum_nc, &
         exclude_hubbard)
   !----------------------------------------------------------------------------
   !! Compute the density response to the perturbation dV = dV_bare + dV_ind by the
   !! non-interacting susceptibility. Solve Sternheimer equation
   !! (H - e +- omega) * dpsi = dvpsi = -P_c+ * (dV_bare + dV_ind) * psi.
   !!
   !! See the comments in sternheimer_kernel for details.
   !!
   !! The differences are the following:
   !! - Have the isolv loop (Sternheimer for +omega and -omega) inside sternheimer_kernel_freq.
   !!   For the zero-frequency case, this loop is done in dfpt_kernel. This is to be unified.
   !! - Use h_prec_freq instead of h_prec. h_diag is complex, not real.
   !! - Use orthogonalize_omega instead of orthogonalize.
   !! - Use ccgsolve_all instead of cgsolve_all.
   !!
   !----------------------------------------------------------------------------
   USE kinds,                 ONLY : DP
   USE io_global,             ONLY : stdout
   USE mp,                    ONLY : mp_sum
   USE mp_pools,              ONLY : inter_pool_comm
   USE buffers,               ONLY : get_buffer, save_buffer
   USE fft_base,              ONLY : dffts
   USE ions_base,             ONLY : nat
   USE klist,                 ONLY : xk, wk, ngk, igk_k
   USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                 ONLY : nbnd, npwx, et
   USE wavefunctions,         ONLY : evc
   USE noncollin_module,      ONLY : noncolin, domag, npol, nspin_mag
   USE uspp,                  ONLY : vkb
   USE uspp_param,            ONLY : nhm
   USE uspp_init,             ONLY : init_us_2
   USE ldaU,                  ONLY : lda_plus_u
   USE units_lr,              ONLY : iuwfc, lrwfc, lrdwf, iudwf
   USE control_lr,            ONLY : nbnd_occ, lgamma
   USE qpoint,                ONLY : nksq, ikks, ikqs
   USE qpoint_aux,            ONLY : ikmks, ikmkmqs, becpt
   USE eqv,                   ONLY : dpsi, dvpsi, evq
   USE apply_dpot_mod,        ONLY : apply_dpot_bands
   USE lr_nc_mag,             ONLY : lr_apply_time_reversal
   USE incdrhoscf_mod,        ONLY : incdrhoscf, incdrhoscf_nc
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: first_iter
   !! true if the first iteration.
   LOGICAL, INTENT(IN) :: time_reversed
   !! true if solving for time reversed wave functions
   LOGICAL, INTENT(IN), OPTIONAL :: exclude_hubbard
   !! true if ignoring the Hubbard response term (Default: .FALSE.)
   INTEGER, INTENT(IN) :: npert
   !! number of perturbations
   INTEGER, INTENT(IN) :: lrdvpsi
   !! record length for the buffer storing dV_bare * psi
   INTEGER, INTENT(IN) :: iudvpsi
   !! unit for the buffer storing dV_bare * psi
   REAL(DP), INTENT(IN) :: thresh
   !! threshold for solving the linear equation
   COMPLEX(DP), INTENT(IN) :: omega
   !! frequency of the perturbation
   LOGICAL, INTENT(OUT) :: all_conv
   !! True if converged at all k points and perturbations
   REAL(DP), INTENT(OUT) :: avg_iter
   !! average number of iterations for the linear equation solver
   COMPLEX(DP), INTENT(INOUT) :: dvscfins(dffts%nnr, nspin_mag, npert)
   !! dV_ind calculated in the previous iteration
   COMPLEX(DP), INTENT(INOUT) :: drhoout(dffts%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   !
   LOGICAL :: conv_root
   !! true if linear system is converged
   LOGICAL :: exclude_hubbard_
   !! Local variable to set the default of exclude_hubbard to false
   INTEGER :: ikk, ikq, npw, npwq, ipert, num_iter, ik, nrec, ikmk, ikmkmq
   !! counters
   INTEGER :: tot_num_iter
   !! total number of iterations in cgsolve_all
   INTEGER :: tot_cg_calls
   !! total number of cgsolve_all calls
   REAL(DP) :: anorm
   !! the norm of the error of the conjugate gradient solution
   REAL(DP) :: rsign
   !! sign of the term in the magnetization
   COMPLEX(DP), ALLOCATABLE :: h_diag1(:, :)
   !! diagonal part of the Hamiltonian, used for preconditioning for +omega
   COMPLEX(DP), ALLOCATABLE :: h_diag2(:, :)
   !! diagonal part of the Hamiltonian, used for preconditioning for -omega
   COMPLEX(DP), ALLOCATABLE :: dvpsi1(:, :)
   !! Bare perturbation times wavefunction, orthogonalized for +omega
   COMPLEX(DP), ALLOCATABLE :: dvpsi2(:, :)
   !! Bare perturbation times wavefunction, orthogonalized for -omega
   COMPLEX(DP), ALLOCATABLE :: dpsi1(:, :)
   !! Solution of the Sternheimer equation for +omega
   COMPLEX(DP), ALLOCATABLE :: dpsi2(:, :)
   !! Solution of the Sternheimer equation for -omega
   COMPLEX(DP) , ALLOCATABLE :: aux2(:, :)
   !! temporary storage used in apply_dpot_bands
   !
   EXTERNAL ch_psi_all_complex, ccg_psi
   !! functions passed to ccgsolve_all
   !
   ! Initialization
   !
   ! JML: For time_reversed, one needs to filp both frequency and magnetic fields.
   !      The code below only accounts for frequency. The magnetic field part is
   !      copied from sternheimer_kernel, but is not adapted.
   !      The main TODO is to move the isolv loop from dfpt_kernel to sternheimer_kernel.
   IF (time_reversed) CALL errore("sternheimer_kernel_freq", &
      "finite-frequency Sternheimer with magnetism is not implemented", 1)
   !
   !
   CALL start_clock("sth_kernel_freq")
   !
   exclude_hubbard_ = .FALSE.
   IF (PRESENT(exclude_hubbard)) exclude_hubbard_ = exclude_hubbard
   !
   !  change the sign of the magnetic field if required
   !
   IF (time_reversed) CALL lr_apply_time_reversal(first_iter, 2, dvscfins)
   !
   ALLOCATE(h_diag1(npwx*npol, nbnd))
   ALLOCATE(h_diag2(npwx*npol, nbnd))
   ALLOCATE(dvpsi1(npwx*npol, nbnd))
   ALLOCATE(dvpsi2(npwx*npol, nbnd))
   ALLOCATE(dpsi1(npwx*npol, nbnd))
   ALLOCATE(dpsi2(npwx*npol, nbnd))
   ALLOCATE(aux2(npwx*npol, nbnd))
   h_diag1 = (0.d0, 0.d0)
   h_diag2 = (0.d0, 0.d0)
   aux2 = (0.d0, 0.d0)
   !
   !$acc enter data create(aux2(1:npwx*npol, 1:nbnd))
   !
   all_conv = .TRUE.
   tot_num_iter = 0
   tot_cg_calls = 0
   !
   DO ik = 1, nksq
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)
      !
      ! Set time-reversed k and k+q points
      !
      IF (time_reversed) THEN
         ikmk = ikmks(ik)
         ikmkmq = ikmkmqs(ik)
         rsign = -1.0_DP
      ELSE
         ikmk = ikk
         ikmkmq = ikq
         rsign = 1.0_DP
      ENDIF
      !
      IF (lsda) current_spin = isk(ikk)
      !
      ! reads unperturbed wavefunctions psi_k in G_space, for all bands
      ! if q=0, evq is a pointer to evc
      !
      IF (nksq > 1 .OR. (noncolin .AND. domag)) THEN
         IF (lgamma) THEN
            !civn: in this case evq is a pointer to evc
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
            !$acc update device(evc)
         ELSE
            !civn: in this case evq is allocated separately and needs to be updated on device
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
            !$acc update device(evc)
            CALL get_buffer(evq, lrwfc, iuwfc, ikmkmq)
            !$acc update device(evq)
         ENDIF
      ENDIF
      !
      ! compute beta functions and kinetic energy for k-point ik
      ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
      !
      CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
      !$acc update host(vkb)
      !
      ! compute the kinetic energy g2kin: (k+q+G)^2
      !
      CALL g2_kin(ikq)
      !
      ! compute preconditioning matrix h_diag used by cgsolve_all
      !
      CALL h_prec_freq(ik, +omega, h_diag1)
      CALL h_prec_freq(ik, -omega, h_diag2)
      !
      DO ipert = 1, npert
         !
         ! read P_c^+ x psi_kpoint into dvpsi.
         ! For nonmagnetic finite-frequency case, dvpsi is the same for +omega and -omega.
         ! For the magnetic case, one would need to read two different dvpsi.
         !
         nrec = (ipert - 1) * nksq + ik
         ! IF (time_reversed) nrec = nrec + npert * nksq
         !
         CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrec)
         !
         IF (.NOT. first_iter) THEN
            !
            ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
            ! dvscf_q from previous iteration (mix_potential)
            !
            CALL apply_dpot_bands(ik, nbnd_occ(ikk), dvscfins(:, :, ipert), evc, aux2)
            dvpsi = dvpsi + aux2
            !
            !  In the case of US pseudopotentials there is an additional
            !  self-consistent term which comes from the dependence of D on
            !  V_{eff} on the bare change of the potential
            !
            CALL adddvscf(ipert, ik, time_reversed)
            !
            ! DFPT+U: add to dvpsi the scf part of the response
            ! Hubbard potential dV_hub
            !
            IF (lda_plus_u .AND. (.NOT. exclude_hubbard_)) CALL adddvhubscf(ipert, ik)
            !
         ENDIF ! .NOT. first_iter
         !
         ! Orthogonalize dvpsi to valence states
         !
         dvpsi1 = dvpsi
         dvpsi2 = dvpsi
         CALL orthogonalize_omega(dvpsi1, evq, ikmk, ikmkmq, dpsi, npwq, +omega)
         CALL orthogonalize_omega(dvpsi2, evq, ikmk, ikmkmq, dpsi, npwq, -omega)
         !
         ! Initial guess for dpsi
         !
         IF (first_iter) THEN
            !
            !  At the first iteration dpsi is set to zero
            !
            dpsi1(:, :) = (0.d0, 0.d0)
            dpsi2(:, :) = (0.d0, 0.d0)
            !
         ELSE
            !
            ! starting value for delta_psi is read from iudwf
            ! Read dpsi1 for +omega and dpsi2 for -omega
            !
            CALL get_buffer(dpsi1, lrdwf, iudwf, nrec)
            CALL get_buffer(dpsi2, lrdwf, iudwf, nrec + npert * nksq)
            !
         ENDIF
         !
         ! iterative solution of the linear system (H-e)*dpsi=dvpsi
         ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
         !
         ! First Sternheimer with +omega, input dvpsi1 and h_diag1, output dpsi1
         !
         conv_root = .TRUE.
         !
         CALL ccgsolve_all(ch_psi_all_complex, ccg_psi, et(1, ikmk), dvpsi1, dpsi1, &
                           h_diag1, npwx, npwq, thresh, ik, num_iter, conv_root, &
                           anorm, nbnd_occ(ikk), npol, +omega)
         !
         tot_num_iter = tot_num_iter + num_iter
         tot_cg_calls = tot_cg_calls + 1
         !
         IF (.NOT. conv_root) THEN
            all_conv = .FALSE.
            WRITE( stdout, "(5x, 'kpoint', i4, ' sternheimer_kernel: &
               &root not converged, thresh < ', es10.3)") ik, anorm
         ENDIF
         !
         ! Second Sternheimer with -omega, input dvpsi2 and h_diag2, output dpsi2
         !
         conv_root = .TRUE.
         !
         CALL ccgsolve_all(ch_psi_all_complex, ccg_psi, et(1, ikmk), dvpsi2, dpsi2, &
                           h_diag2, npwx, npwq, thresh, ik, num_iter, conv_root, &
                           anorm, nbnd_occ(ikk), npol, -omega)
         !
         tot_num_iter = tot_num_iter + num_iter
         tot_cg_calls = tot_cg_calls + 1
         !
         IF (.NOT. conv_root) THEN
            all_conv = .FALSE.
            WRITE( stdout, "(5x, 'kpoint', i4, ' sternheimer_kernel: &
               &root not converged, thresh < ', es10.3)") ik, anorm
         ENDIF
         !
         ! writes delta_psi on iunit iudwf, k=kpoint,
         !
         CALL save_buffer(dpsi1, lrdwf, iudwf, nrec)
         CALL save_buffer(dpsi2, lrdwf, iudwf, nrec + npert * nksq)
         !
         ! calculates dvscf, sum over k => dvscf_q_ipert
         !
         dpsi = dpsi1 + dpsi2
         !
         IF (noncolin) THEN
            CALL incdrhoscf_nc(drhoout(1,1,ipert), wk(ikk), ik, &
                               dbecsum_nc(1,1,1,1,ipert), dpsi, rsign)
         ELSE
            CALL incdrhoscf(drhoout(1,current_spin,ipert), wk(ikk), &
                            ik, dbecsum(1,1,current_spin,ipert), dpsi)
         ENDIF
      ENDDO ! ipert
   ENDDO ! ik
   !
   CALL mp_sum(tot_num_iter, inter_pool_comm)
   CALL mp_sum(tot_cg_calls, inter_pool_comm)
   avg_iter = REAL(tot_num_iter, DP) / REAL(tot_cg_calls, DP)
   !
   !$acc exit data delete(aux2)
   !
   !  reset the original magnetic field if it was changed
   !
   IF (time_reversed) CALL lr_apply_time_reversal(first_iter, 1, dvscfins)
   !
   DEALLOCATE(aux2)
   DEALLOCATE(h_diag1)
   DEALLOCATE(h_diag2)
   DEALLOCATE(dvpsi1)
   DEALLOCATE(dvpsi2)
   DEALLOCATE(dpsi1)
   DEALLOCATE(dpsi2)
   !
   CALL stop_clock("sth_kernel_freq")
   !
!----------------------------------------------------------------------------
END SUBROUTINE sternheimer_kernel_freq
!------------------------------------------------------------------------------
!
SUBROUTINE sternheimer_postprocess(nsolv, dfpt_data, dbecsum_nc, dbecsum_nc_trev)
   !----------------------------------------------------------------------------
   !! Postprocess the results after a Sternheimer calculation.
   !! Modifies the drhos, drhop, and dbecsum fields of dfpt_data.
   !----------------------------------------------------------------------------
   USE kinds,                 ONLY : DP
   USE mp,                    ONLY : mp_sum
   USE mp_pools,              ONLY : inter_pool_comm
   USE mp_bands,              ONLY : intra_bgrp_comm
   USE fft_base,              ONLY : dffts, dfftp
   USE fft_interfaces,        ONLY : fft_interpolate
   USE gvecs,                 ONLY : doublegrid
   USE ions_base,             ONLY : nat
   USE lsda_mod,              ONLY : nspin
   USE noncollin_module,      ONLY : noncolin, nspin_mag, domag
   USE uspp,                  ONLY : okvan
   USE uspp_param,            ONLY : nhm
   USE paw_variables,         ONLY : okpaw
   USE control_lr,            ONLY : lgamma_gamma
   USE dfpt_type,             ONLY : dfpt_data_type
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: nsolv
   !! Number of Sternheimer equations solved (1 for nonmagnetic/LSDA zero frequency,
   !! 2 for noncollinear magnetism or finite frequency)
   TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
   !! Data that describes linear response quantities
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc(nhm, nhm, nat, nspin, dfpt_data%npert)
   !! becsum with dpsi. Used if noncolin is true.
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc_trev(nhm, nhm, nat, nspin, dfpt_data%npert)
   !! becsum with dpsi with time-reversed Sternheimer. Used if noncolin .AND. domag is true.
   !
   INTEGER :: ipert
   !! counter for perturbations
   INTEGER :: is
   !! counter for spins
   INTEGER :: npert
   !! number of perturbations, shorthand for dfpt_data%npert
   COMPLEX(DP), ALLOCATABLE :: dbecsum_aux(:, :, :, :)
   !
   CALL start_clock("sth_postproc")
   !
   IF (noncolin .AND. domag) THEN
      IF (.NOT. PRESENT(dbecsum_nc_trev)) CALL errore("sternheimer_postprocess", &
         "dbecsum_nc_trev must be present if noncolin and domag are true", 1)
   ENDIF
   !
   npert = dfpt_data%npert
   !
   IF (nsolv == 2) THEN
      dfpt_data%drhos = dfpt_data%drhos / 2.0_dp
      IF (noncolin) THEN
         dbecsum_nc = dbecsum_nc / 2.0_dp
         IF (domag) dbecsum_nc_trev = dbecsum_nc_trev / 2.0_dp
      ELSE
         dfpt_data%dbecsum = dfpt_data%dbecsum / 2.0_dp
      ENDIF
   ENDIF
   !
   !  The calculation of dbecsum is distributed across processors
   !  (see addusdbec) - we sum over processors the contributions
   !  coming from each slice of bands
   !
   IF (noncolin) THEN
      CALL mp_sum(dbecsum_nc, intra_bgrp_comm)
      IF (domag) CALL mp_sum(dbecsum_nc_trev, intra_bgrp_comm)
   ELSE
      CALL mp_sum(dfpt_data%dbecsum, intra_bgrp_comm)
   ENDIF
   !
   IF (doublegrid) THEN
      DO is = 1, nspin_mag
         DO ipert = 1, npert
            CALL fft_interpolate(dffts, dfpt_data%drhos(:, is, ipert), &
                                 dfftp, dfpt_data%drhop(:, is, ipert))
         ENDDO
      ENDDO
   ELSE
      CALL zcopy(dffts%nnr * nspin_mag * npert, dfpt_data%drhos, 1, dfpt_data%drhop, 1)
   ENDIF
   !
   !  In the noncolinear, spin-orbit case rotate dbecsum
   !
   IF (noncolin .AND. okvan) THEN
      ! dbecsum_nc has 2 if domag, 1 if nonmagnetic
      CALL set_dbecsum_nc(dbecsum_nc, dfpt_data%dbecsum, npert)
      !
      IF (domag) THEN
         ALLOCATE(dbecsum_aux((nhm * (nhm + 1))/2, nat, nspin_mag, npert))
         dbecsum_aux = (0.0_DP, 0.0_DP)
         CALL set_dbecsum_nc(dbecsum_nc_trev, dbecsum_aux, npert)
         dfpt_data%dbecsum(:,:,1,:) = dfpt_data%dbecsum(:,:,1,:) + dbecsum_aux(:,:,1,:)
         dfpt_data%dbecsum(:,:,2:4,:) = dfpt_data%dbecsum(:,:,2:4,:) - dbecsum_aux(:,:,2:4,:)
         DEALLOCATE(dbecsum_aux)
      ENDIF
   ENDIF
   !
   ! Add augmentation charge contribution to drhop (for USPP/PAW)
   !
   IF (okvan) CALL lr_addusddens(npert, dfpt_data%dbecsum, dfpt_data%drhop)
   !
   ! Add Pulay correction to drhop if present
   !
   IF (ALLOCATED(dfpt_data%drhop_pulay)) THEN
      CALL zaxpy(dfftp%nnr * nspin_mag * npert, (1.d0, 0.d0), dfpt_data%drhop_pulay, 1, &
                 dfpt_data%drhop, 1)
   ENDIF
   !
   !   Reduce the delta rho across pools
   !   Postprocess dbecsum only if okpaw, since it is not used anymore otherwise.
   !
   CALL mp_sum(dfpt_data%drhos, inter_pool_comm)
   CALL mp_sum(dfpt_data%drhop, inter_pool_comm)
   IF (okpaw) CALL mp_sum(dfpt_data%dbecsum, inter_pool_comm)
   !
   CALL stop_clock("sth_postproc")
   !
END SUBROUTINE sternheimer_postprocess
!------------------------------------------------------------------------------
END MODULE response_kernels
!------------------------------------------------------------------------------
