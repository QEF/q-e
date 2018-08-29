!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_solve_linear_system (na, iq)
  !-----------------------------------------------------------------------
  !
  ! This is a driver routine for the solution of the linear-response Kohn-Sham 
  ! equations (45) in Ref. [1]. The solution defines the change of Kohn-Sham 
  ! wavefunctions due change of occupations.
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !    
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat
  USE io_global,            ONLY : stdout
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions,        ONLY : evc
  USE cell_base,            ONLY : tpiba2
  USE klist,                ONLY : lgauss, ltetra, xk, wk, nelec, ngk, igk_k
  USE gvect,                ONLY : g
  USE gvecs,                ONLY : doublegrid
  USE scf,                  ONLY : rho
  USE fft_base,             ONLY : dfftp, dffts
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE wvfct,                ONLY : nbnd, npwx, g2kin, et
  USE uspp,                 ONLY : okvan, vkb, nkb
  USE uspp_param,           ONLY : nhm
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, becp
  USE buffers,              ONLY : save_buffer, get_buffer
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE paw_variables,        ONLY : okpaw
  USE paw_onecenter,        ONLY : paw_dpotential
  USE paw_symmetry,         ONLY : paw_dusymmetrize, paw_dumqsymmetrize
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE hp_efermi_shift,      ONLY : hp_ef_shift, def
  USE eqv,                  ONLY : dvpsi, dpsi, evq
  USE qpoint,               ONLY : nksq, ikks, ikqs, xq
  USE control_lr,           ONLY : lgamma, nbnd_occ
  USE units_lr,             ONLY : iuwfc, lrwfc
  USE lrus,                 ONLY : int3, int3_paw
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE fft_helper_subroutines
  USE fft_interfaces,       ONLY : fft_interpolate
  USE lr_symm_base,         ONLY : irotmq, minus_q, nsymq, rtau
  USE ldaU_hp,              ONLY : thresh_init, dnsscf, dns0, trace_dns_tot_old,  & 
                                   conv_thr_chi_best, iter_best, niter_max, nmix, &
                                   alpha_mix, iudwfc, lrdwfc, code
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na ! number of the perturbed atom
  INTEGER, INTENT(IN) :: iq ! number of the q point
  !
  REAL(DP),  ALLOCATABLE :: h_diag (:,:) ! diagonal part of the Hamiltonian
  !
  REAL(DP) :: thresh, & ! convergence threshold
              anorm,  & ! the norm of the error
              averlt, & ! average number of iterations
              dr2       ! self-consistency error
  !
  REAL(DP) :: dos_ef, &  ! density of states at the Fermi level
              weight, &  ! Misc variables for metals
              aux_avg(2) ! Misc variables for metals
  !
  REAL(DP), ALLOCATABLE :: becsum1(:,:,:)
  ! 
  COMPLEX(DP), ALLOCATABLE, TARGET :: dvscfin(:,:)
  ! change of the scf potential (input)
  !
  COMPLEX(DP), POINTER :: dvscfins(:,:)
  ! change of the scf potential (smooth part only)
  !
  COMPLEX(DP), ALLOCATABLE :: drhoscf  (:,:), &
                              drhoscfh (:,:), &
                              dvscfout (:,:)
  ! change of rho / scf potential (output)
  !
  COMPLEX(DP), ALLOCATABLE :: &
     ldos (:,:),             & ! local density of states at Ef
     ldoss (:,:),            & ! as above, without augmentation charges
     dbecsum (:,:,:,:),      & ! the derivative of becsum
     aux1 (:,:), aux2 (:,:), & ! auxiliary arrays
     mixin(:), mixout(:),    & ! auxiliary arrays for mixing of the response potential
     tg_dv(:,:),             & ! Task groups: auxiliary array for potential * wfct
     tg_psic(:,:)              ! Task groups: auxiliary array for wavefunctions
 
  COMPLEX(DP), ALLOCATABLE :: t(:,:,:,:), tmq(:,:,:)
  ! PAW: auxiliary arrays
 
  LOGICAL :: conv_root,  & ! true if linear system is converged
             exst,       & ! used to open the recover file
             lmetq0,     & ! true if xq=(0,0,0) in a metal
             convt,      & ! not needed for HP 
             convt_chi     ! used instead of convt to control the convergence

  REAL(DP), PARAMETER :: tr2 = 1.D-30 ! threshold parameter

  INTEGER :: ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             ltaver,     & ! average counter
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ig,         & ! counter on G vectors
             ndim,       &
             is,         & ! counter on spin polarizations
             nt,         & ! counter on types
             ios,        & ! integer variable for I/O control
             incr,       & ! used for task groups
             v_siz,      & ! size of the potential
             npw,        & ! number of plane waves at k 
             npwq          ! number of plane waves at k+q 

  REAL(DP) :: tcpu, get_clock ! timing variables
  CHARACTER(LEN=256) :: filename, &
                        flmixdpot = 'mixd'
  EXTERNAL ch_psi_all, cg_psi
  !
  CALL start_clock ('hp_solve_linear_system')
  !
  WRITE( stdout,*) "     =--------------------------------------------="
  WRITE( stdout, '(13x,"START SOLVING THE LINEAR SYSTEM")')
  WRITE( stdout,*) "     =--------------------------------------------="
  !
  ! Allocate arrays for the SCF density/potential
  !
  ALLOCATE (drhoscf (dffts%nnr, nspin_mag)) 
  ALLOCATE (drhoscfh(dfftp%nnr, nspin_mag))
  ALLOCATE (dvscfin (dfftp%nnr, nspin_mag))
  ALLOCATE (dvscfout(dfftp%nnr, nspin_mag))
  !
  IF (doublegrid) THEN
     ALLOCATE (dvscfins(dffts%nnr, nspin_mag))
  ELSE
     dvscfins => dvscfin
  ENDIF
  !
  ! USPP-specific allocations
  !
  IF (okvan) ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, 1))
  IF (okpaw) ALLOCATE (int3_paw ( nhm, nhm, nat, nspin_mag, 1))
  CALL allocate_bec_type (nkb, nbnd, becp)
  !
  ALLOCATE (dbecsum((nhm*(nhm+1))/2, nat, nspin_mag, 1))
  !
  IF (okpaw) THEN
     !
     ALLOCATE (mixin(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
     ALLOCATE (mixout(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
     mixin = (0.0_DP,0.0_DP)
     !
     ! Auxiliary unitary arrays
     ALLOCATE ( tmq(1,1,3*nat) )
     ALLOCATE ( t(1,1,48,3*nat) )
     t(:,:,:,:) = (1.0_DP, 0.0_DP)
     tmq(:,:,:) = (1.0_DP, 0.0_DP)
     !
  ENDIF
  !
  ALLOCATE (aux1(dffts%nnr, npol))
  ALLOCATE (aux2(npwx*npol, nbnd))
  ALLOCATE (h_diag(npwx*npol, nbnd))
  ALLOCATE (trace_dns_tot_old(nat))
  trace_dns_tot_old(:) = (0.d0, 0.d0)
  !
  convt     = .FALSE.
  convt_chi = .FALSE.
  !
  incr = 1
  IF ( dffts%has_task_groups ) THEN
     !
     v_siz =  dffts%nnr_tg
     ALLOCATE( tg_dv  ( v_siz, nspin_mag ) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = fftx_ntgrp(dffts)
     !
  ENDIF
  !
  ! If q=0 for a metal: allocate and compute local DOS and DOS at Ef
  !
  lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
  !
  IF (lmetq0) THEN
     ALLOCATE (ldos (dfftp%nnr, nspin_mag))
     ALLOCATE (ldoss(dffts%nnr, nspin_mag))
     ALLOCATE (becsum1 ( (nhm * (nhm + 1))/2, nat, nspin_mag))
     CALL localdos (ldos, ldoss, becsum1, dos_ef)
     IF (.NOT.okpaw) DEALLOCATE (becsum1)
  ENDIF
  !
  ! The loop of the linear-response calculation
  !
  DO iter = 1, niter_max
     !
     WRITE(stdout,'(/6x,"atom #",i3,3x,"q point #",i4,3x,"iter # ",i3)') na, iq, iter
     !
     ltaver = 0
     lintercall = 0
     !
     drhoscf(:,:)     = (0.d0, 0.d0)
     dvscfout(:,:)    = (0.d0, 0.d0)
     dbecsum(:,:,:,:) = (0.d0, 0.d0)
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!   START OF THE K LOOP   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
     DO ik = 1, nksq
        !
        ikk  = ikks(ik)
        ikq  = ikqs(ik)
        npw  = ngk(ikk)
        npwq = ngk(ikq)
        !
        IF (lsda) current_spin = isk(ikk)
        !
        ! Read unperturbed KS wavefuctions psi(k) and psi(k+q)
        !
        IF (nksq.gt.1) THEN
           IF (lgamma) THEN
              CALL get_buffer (evc, lrwfc, iuwfc, ikk)
           ELSE
              CALL get_buffer (evc, lrwfc, iuwfc, ikk)
              CALL get_buffer (evq, lrwfc, iuwfc, ikq)
           ENDIF
        ENDIF
        !
        ! USPP: Compute the projectors vkb at k+q
        !
        CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkb)
        !
        ! Compute the kinetic energy at k+q
        !
        CALL g2_kin (ikq)
        !
        ! Compute preconditioning matrix h_diag used by cgsolve_all
        !
        CALL h_prec (ik, evq, h_diag)
        !
        ! Computes (iter=1) or reads (iter>1) the action of the perturbing
        ! potential on the unperturbed KS wavefunctions: |dvpsi> = dV_pert * |evc>
        ! See Eq. (46) in Ref. [1] 
        !
        CALL hp_dvpsi_pert(ik)
        !
        IF ( iter > 1 ) THEN
           !
           ! Add the contribution of the self consistent term.
           ! Calculates dvscf_q*psi(k) in G-space, for all bands, k=ik
           ! dvscf_q from previous iteration (mix_potential)
           !
           CALL start_clock ('hp_vpsifft')
           IF ( dffts%has_task_groups ) &
                 & CALL tg_cgather( dffts, dvscfins(:,current_spin), tg_dv(:,1))
           aux2 = (0.0_DP, 0.0_DP)
           DO ibnd = 1, nbnd_occ(ikk), incr
              IF ( dffts%has_task_groups ) THEN
                 CALL cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
                 CALL apply_dpot(v_siz, tg_psic, tg_dv, 1)
                 CALL cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, nbnd_occ(ikk))
              ELSE
                 CALL cft_wave (ik, evc(:,ibnd), aux1, +1)
                 CALL apply_dpot(dffts%nnr, aux1, dvscfins, current_spin)
                 CALL cft_wave (ik, aux2(:,ibnd), aux1, -1)
              ENDIF
           ENDDO
           dvpsi = dvpsi + aux2
           CALL stop_clock ('hp_vpsifft')
           !
           ! USPP: there is an additional self-consistent term proportional to int3
           ! |dvpsi> = |dvpsi> + dV_HXC*|evc> + int3 * |beta><beta|evc>
           !
           IF (okvan) CALL adddvscf(1, ik)
           !
        ENDIF
        !
        ! Ortogonalize dvpsi to valence states: ps = <evq|dvpsi>
        ! Apply -P_c^+. See Eq. (A21) in Ref. [1]
        !
        CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .FALSE.)
        !
        IF ( iter == 1 ) THEN
           !
           ! At the first iteration dpsi and dvscfin are set to zero
           !
           dpsi(:,:)    = (0.d0, 0.d0)
           dvscfin(:,:) = (0.d0, 0.d0)
           !
           ! Starting threshold for iterative solution of the linear system.
           ! A strickt threshold for the first iteration is needed, 
           ! because we need dns0 to very high precision.
           !
           thresh = thresh_init * nelec
           !
        ELSE
           !
           ! Starting value for dpsi is read from file
           !
           CALL get_buffer( dpsi, lrdwfc, iudwfc, ik)
           !
           ! Threshold for iterative solution of the linear system.
           ! We start with not a strict threshold for iter=2, and then
           ! it decreases with iterations.
           !
           thresh = MIN (1.D-1 * SQRT(dr2), 1.D-2)
           !
        ENDIF
        !
        ! Iterative solution of the linear system:
        ! (H + Q - eS) * |dpsi> = |dvpsi>,
        ! where |dvpsi> = - P_c^+ (dV_HXC + dV_pert) * |evc>
        ! See Eq. (43) in Ref. [1]
        !
        CALL cgsolve_all (ch_psi_all, cg_psi, et(1,ikk), dvpsi, dpsi, h_diag, &
              & npwx, npwq, thresh, ik, lter, conv_root, anorm, nbnd_occ(ikk), npol )
        !
        ltaver = ltaver + lter
        !
        lintercall = lintercall + 1
        !
        IF (.NOT.conv_root) THEN
           WRITE( stdout, '(6x,"kpoint",i4,  &
             & " hp_solve_linear_system: root not converged, thresh < ",e10.3)') ik , anorm
           IF (iter == 1) WRITE( stdout, '(6x,"Try to increase thresh_init...")')
        ENDIF
        !
        ! Writes dpsi on file for a given k
        !
        CALL save_buffer (dpsi, lrdwfc, iudwfc, ik)
        !
        ! Setup the weight at point k (normalized by the number of k points)
        !
        weight = wk(ikk)
        !
        ! Calculates the response charge density (sum over k)
        ! See Eq. (48) in Ref. [1]
        ! 
        CALL incdrhoscf (drhoscf(:,current_spin), weight, ik, &
                           & dbecsum(:,:,current_spin,1), dpsi)
        !
     ENDDO ! k points 
     !
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!   END OF THE K LOOP   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !
#if defined (__MPI)
     !
     ! USPP: The calculation of dbecsum is distributed across processors (see addusdbec)
     ! Sum over processors the contributions coming from each slice of bands
     !
     CALL mp_sum ( dbecsum, intra_pool_comm )
     !
#endif
     !
     ! Copy/interpolate the response density drhoscf -> drhoscfh
     !
     IF (doublegrid) THEN
        do is = 1, nspin_mag
           CALL fft_interpolate (dffts, drhoscf(:,is), dfftp, drhoscfh(:,is))
        enddo
     ELSE
        CALL zcopy (nspin_mag*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     ENDIF
     !
     ! USPP: Compute the total response charge density (standard term + US term)
     !
     IF (okvan) CALL lr_addusddens (drhoscfh, dbecsum)
     !
#if defined (__MPI)
     CALL mp_sum ( drhoscfh, inter_pool_comm ) 
     IF (okpaw) CALL mp_sum ( dbecsum, inter_pool_comm )
#endif
     !
     ! PAW: the factor of 2 is due to the presence of the CC term
     ! (see first two terms in Eq.(9) in PRB 81, 075123 (2010))
     !
     IF (okpaw) dbecsum = 2.0_DP * dbecsum
     !
     ! Metallic case and q=0: add a correction to the response charge density
     ! due to the shift of the Fermi energy (see Eq.(75) in Rev. Mod. Phys. 73, 515 (2001)).
     ! This term is added to the response charge density (in order to obtain correct 
     ! response HXC potential) and to the response occupation matrices.
     ! 
     IF (lmetq0) THEN
        ! 
        IF (okpaw) THEN
           CALL hp_ef_shift (drhoscfh, ldos, ldoss, dos_ef, dbecsum, becsum1)
        ELSE
           CALL hp_ef_shift (drhoscfh, ldos, ldoss, dos_ef)
        ENDIF
        !
        ! Check that def is not too large (it is in Ry). 
        !
        IF ( ABS(DBLE(def)) > 5.0d0 ) THEN
           !
           WRITE( stdout, '(/6x,"WARNING: The Fermi energy shift is too big!")')
           WRITE( stdout, '(6x, "This may happen in two cases:")')
           WRITE( stdout, '(6x, "1. The DOS at the Fermi level is too small:")')
           WRITE( stdout, '(6x, "   DOS(E_Fermi) = ",1x,2e12.4)') dos_ef
           WRITE( stdout, '(6x, "   This means that most likely the system has a gap,")')
           WRITE( stdout, '(6x, "   and hence it should NOT be treated as a metal")')
           WRITE( stdout, '(6x, "   (otherwise numerical instabilities will appear).")')
           WRITE( stdout, '(6x, "2. Numerical instabilities due to too low cutoff")')
           WRITE( stdout, '(6x, "   for hard pseudopotentials.")')
           WRITE( stdout, '(/6x,"Stopping...")')
           !
           CALL hp_stop_smoothly (.FALSE.)
           !
        ENDIF    
        !
     ENDIF
     !
     ! Symmetrization of the response charge density.
     !
     CALL hp_psymdvscf (drhoscfh)
     !
     ! Symmetrize dbecsum
     !
     IF (okpaw) THEN
        IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,1,1,1,irotmq,rtau,xq,tmq)
        CALL PAW_dusymmetrize(dbecsum,1,1,1,nsymq,rtau,xq,t)
     ENDIF
     !
     ! Copy drhoscfh to dvscfout 
     !
     CALL zcopy (nspin_mag*dfftp%nnr, drhoscfh, 1, dvscfout, 1)
     !
     ! Compute the response potential dV_HXC from the response charge density.
     ! See Eq. (47) in Ref. [1]
     !
     CALL dv_of_drho (dvscfout, .FALSE.)
     !
     ! Mix the new HXC response potential (dvscfout) with the old one (dvscfin).
     ! Note: dvscfin = 0 for iter = 1 (so there is no mixing).
     ! Output: dvscfin becomes a mixed potential
     !         dvscfout contains the difference vout-vin (it is not needed)
     ! PAW: mix also dbecsum
     !
     IF (okpaw) THEN
        CALL setmixout(dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag)/2, &
                  & mixout, dvscfout, dbecsum, ndim, -1 )
        CALL mix_potential (2*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                  & alpha_mix(iter), dr2, tr2/npol, iter, nmix, flmixdpot, convt)
        CALL setmixout(dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag)/2, &
                  & mixin, dvscfin, dbecsum, ndim, 1 )
     ELSE
        CALL mix_potential (2*dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
                  & alpha_mix(iter), dr2, tr2/npol, iter, nmix, flmixdpot, convt)
     ENDIF
     !
     ! NCPP case: dvscfins is a pointer to dvscfin.
     ! USPP case: Interpolate dvscfin from dense to smooth grid
     !            and put the result in dvscfins (needed for the next iteration).
     !
     IF (doublegrid) THEN
        DO is = 1, nspin_mag
           CALL fft_interpolate (dfftp, dvscfin(:,is), dffts, dvscfins(:,is))
        ENDDO
     ENDIF
     !
     ! PAW: compute the response PAW potential 
     ! (see the last term in Eq.(12) in PRB 81, 075123 (2010))
     !
     IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,1)
     !
     ! USPP: With the new change of the HXC response potential dV_HXC 
     ! we compute the integral of the change of this potential and Q.
     ! int3 = \int Q(r) dV_HXC(r) dr
     ! PAW: int3_paw is added to int3 inside of the routine newdq
     !
     IF (okvan) CALL newdq (dvscfin, 1)
     !
     ! Calculate the response occupation matrix
     ! See Eq. (43) in Ref. [1]
     !
     CALL hp_dnsq (lmetq0, iter, convt_chi, dnsscf(:,:,:,:,iq))
     !
     ! Save the response occupation matrix after the first iteration,
     ! which was computed from dpsi corresponding to the perturbing
     ! potential only (dV_HXC=0). This is needed for the calculation of chi0.
     !
     IF ( iter == 1 ) dns0(:,:,:,:,iq) = dnsscf(:,:,:,:,iq)
     !
     ! Compute the average number of iterations 
     ! 
#if defined (__MPI)
     aux_avg(1) = DBLE(ltaver)
     aux_avg(2) = DBLE(lintercall)
     CALL mp_sum ( aux_avg, inter_pool_comm )
     averlt = aux_avg(1) / aux_avg(2)
#else
     averlt = DBLE(ltaver) / DBLE(lintercall)
#endif
     !
     tcpu = get_clock(code)
     !
     WRITE( stdout, '(6x,"Average number of iter. to solve lin. system:",2x,f5.1)') averlt
     WRITE( stdout, '(6x,"Total CPU time :",f8.1,1x,"s")') tcpu
     !
     IF ( check_stop_now() ) CALL hp_stop_smoothly (.FALSE.)
     !
     IF ( iter==niter_max .AND. .NOT.convt_chi) THEN
        WRITE( stdout, '(/6x,"Convergence has not been reached after",1x,i3,1x,"iterations!")') niter_max
        IF ( iter > 1 ) THEN
           WRITE( stdout, '(6x,"The best overall accuracy which was reached :")')
           WRITE( stdout, '(6x,"diff = ",1x,f14.10,1x," iter =",1x,i3)') conv_thr_chi_best, iter_best
        ENDIF
        WRITE( stdout, '(/6x,"Stopping...")')
        CALL hp_stop_smoothly (.TRUE.) 
     ENDIF
     !
     IF (convt_chi) goto 155
     !
  ENDDO  ! loop over the iterations iter
  !
155 CONTINUE
  !
  DEALLOCATE (h_diag)
  DEALLOCATE (aux1)
  DEALLOCATE (aux2)
  DEALLOCATE (dbecsum)
  DEALLOCATE (drhoscf)
  DEALLOCATE (drhoscfh)
  DEALLOCATE (dvscfin)
  DEALLOCATE (dvscfout)
  DEALLOCATE (trace_dns_tot_old)
  IF (doublegrid)       DEALLOCATE (dvscfins)
  IF (ALLOCATED(ldoss)) DEALLOCATE (ldoss)
  IF (ALLOCATED(ldos))  DEALLOCATE (ldos)
  IF (ALLOCATED(becsum1))  DEALLOCATE (becsum1)
  IF (okvan) DEALLOCATE (int3)
  IF (okpaw) THEN
     DEALLOCATE (int3_paw)
     DEALLOCATE (mixin, mixout)
     DEALLOCATE (t)
     DEALLOCATE (tmq)
  ENDIF
  CALL deallocate_bec_type (becp)
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE( tg_dv )
     DEALLOCATE( tg_psic )
  ENDIF
  !
  WRITE( stdout,*) "     "
  WRITE( stdout,*) "     =--------------------------------------------="
  WRITE( stdout,   '(13x,"CONVERGENCE HAS BEEN REACHED")')
  WRITE( stdout,*) "     =--------------------------------------------="
  !
  CALL stop_clock ('hp_solve_linear_system')
  !
  RETURN
  !
END SUBROUTINE hp_solve_linear_system
