!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
MODULE rho_response
CONTAINS
SUBROUTINE nonint_rho_response(first_iter, npert, lrdvpsi, iudvpsi, thresh, dvscfins, &
         avg_iter, drhoout, dbecsum, dbecsum_nc)
   !----------------------------------------------------------------------------
   !! Compute the density response to the perturbation dV = dV_bare + dV_ind by the
   !! non-interacting susceptibility. Solve Sternheimer equation
   !! (H - e) * dpsi = dvpsi = -P_c+ * (dV_bare + dV_ind) * psi.
   !!
   !! dV_bare * psi is read from buffer iudvpsi, so they must be already calculated.
   !! dV_ind is given by input dvscfins, and dV_ind * psi is calculated in apply_dpot_bands.
   !!
   !! For USPPs, adddvscf is called, so relevant arrays must be already calculated.
   !! For DFT+U, adddvhubscf is called, so relevant arrays must be already calculated.
   !!
   !! Input:
   !!    - first_iter: true if the first iteration, where dvscfins = 0
   !!    - npert: number of perturbations
   !!    - lrdvpsi: record length for the buffer storing dV_bare * psi
   !!    - iudvpsi: unit for the buffer storing dV_bare * psi
   !!    - thresh: threshold for solving Sternheimer equation
   !!    - dvscfins: dV_ind calculated in the previous iteration
   !!
   !! Output:
   !!    - avg_iter: average number of iterations for the linear equation solver
   !!    - drhoout: induced charge density
   !!    - dbecsum: becsum with dpsi
   !!    - dbecsum_nc: becsum with dpsi. Optional, used if noncolin is true.
   !----------------------------------------------------------------------------
   USE kinds,                 ONLY : DP
   USE io_global,             ONLY : stdout
   USE buffers,               ONLY : get_buffer, save_buffer
   USE fft_base,              ONLY : dfftp
   USE ions_base,             ONLY : nat
   USE klist,                 ONLY : xk, wk, ngk, igk_k
   USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                 ONLY : nbnd, npwx, et
   USE wavefunctions,         ONLY : evc
   USE noncollin_module,      ONLY : noncolin, npol, nspin_mag
   USE uspp,                  ONLY : vkb
   USE uspp_param,            ONLY : nhm
   USE ldaU,                  ONLY : lda_plus_u
   USE units_ph,              ONLY : lrdwf, iudwf
   USE units_lr,              ONLY : iuwfc, lrwfc
   USE control_lr,            ONLY : nbnd_occ, lgamma
   USE qpoint,                ONLY : nksq, ikks, ikqs
   USE eqv,                   ONLY : dpsi, dvpsi, evq
   USE apply_dpot_mod,        ONLY : apply_dpot_bands
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: first_iter
   !! true if the first iteration.
   INTEGER, INTENT(IN) :: npert
   !! number of perturbations
   INTEGER, INTENT(IN) :: lrdvpsi
   !! record length for the buffer storing dV_bare * psi
   INTEGER, INTENT(IN) :: iudvpsi
   !! unit for the buffer storing dV_bare * psi
   REAL(DP), INTENT(IN) :: thresh
   !! threshold for solving the linear equation
   REAL(DP), INTENT(OUT) :: avg_iter
   !! average number of iterations for the linear equation solver
   COMPLEX(DP), POINTER, INTENT(IN) :: dvscfins(:, :, :)
   !! dV_ind calculated in the previous iteration
   COMPLEX(DP), INTENT(INOUT) :: drhoout(dfftp%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   !
   LOGICAL :: conv_root
   !! true if linear system is converged
   INTEGER :: ikk, ikq, npw, npwq, ipert, num_iter, ik, nrec
   !! counters
   INTEGER :: tot_num_iter
   !! total number of iterations in cgsolve_all
   INTEGER :: tot_cg_calls
   !! total number of cgsolve_all calls
   REAL(DP) :: anorm
   !! the norm of the error of the conjugate gradient solution
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
   ALLOCATE(h_diag(npwx*npol, nbnd))
   ALLOCATE(aux2(npwx*npol, nbnd))
   !
   tot_num_iter = 0
   tot_cg_calls = 0
   drhoout = (0.d0, 0.d0)
   dbecsum = (0.d0, 0.d0)
   IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
   !
   DO ik = 1, nksq
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)
      !
      IF (lsda) current_spin = isk(ikk)
      !
      ! reads unperturbed wavefunctions psi_k in G_space, for all bands
      ! if q=0, evq is a pointer to evc
      !
      IF (nksq > 1) THEN
         IF (lgamma) THEN
            CALL get_buffer(evc, lrwfc, iuwfc, ikk)
         ELSE
            CALL get_buffer(evc, lrwfc, iuwfc, ikk)
            CALL get_buffer(evq, lrwfc, iuwfc, ikq)
         ENDIF
      ENDIF
      !
      ! compute beta functions and kinetic energy for k-point ik
      ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
      !
      CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb)
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
            !  selfconsist term which comes from the dependence of D on
            !  V_{eff} on the bare change of the potential
            !
            CALL adddvscf(ipert, ik)
            !
            ! DFPT+U: add to dvpsi the scf part of the response
            ! Hubbard potential dV_hub
            !
            IF (lda_plus_u) CALL adddvhubscf(ipert, ik)
            !
         ENDIF
         !
         ! Orthogonalize dvpsi to valence states
         !
         CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)
         !
         ! Initial guess for dpsi
         !
         IF (first_iter) THEN
            !
            !  At the first iteration dpsi is set to zero
            !
            dpsi(:, :) = (0.d0,0.d0)
         ELSE
            ! starting value for delta_psi is read from iudwf
            !
            nrec = (ipert - 1) * nksq + ik
            CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
         ENDIF
         !
         ! iterative solution of the linear system (H-e)*dpsi=dvpsi
         ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
         !
         conv_root = .TRUE.
         !
         CALL cgsolve_all(ch_psi_all, cg_psi, et(1, ikk), dvpsi, dpsi, h_diag, &
            npwx, npwq, thresh, ik, num_iter, conv_root, anorm, nbnd_occ(ikk), npol)
         !
         tot_num_iter = tot_num_iter + num_iter
         tot_cg_calls = tot_cg_calls + 1
         !
         IF (.NOT. conv_root) WRITE( stdout, "(5x, 'kpoint', i4, &
                & ' solve_e: root not converged ', es10.3)") ik, anorm
         !
         ! writes delta_psi on iunit iudwf, k=kpoint,
         !
         nrec = (ipert - 1) * nksq + ik
         CALL save_buffer(dpsi, lrdwf, iudwf, nrec)
         !
         ! calculates dvscf, sum over k => dvscf_q_ipert
         !
         IF (noncolin) THEN
            CALL incdrhoscf_nc(drhoout(1,1,ipert), wk(ikk), ik, &
                               dbecsum_nc(1,1,1,1,ipert), dpsi, 1.0d0)
         ELSE
            CALL incdrhoscf(drhoout(1,current_spin,ipert), wk(ikk), &
                            ik, dbecsum(1,1,current_spin,ipert), dpsi)
         ENDIF
      ENDDO ! ipert
   ENDDO ! ik
   !
   avg_iter = REAL(tot_num_iter, DP) / REAL(tot_cg_calls, DP)
   !
!----------------------------------------------------------------------------
END SUBROUTINE nonint_rho_response
!------------------------------------------------------------------------------
END MODULE rho_response
!------------------------------------------------------------------------------
