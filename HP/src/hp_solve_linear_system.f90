!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
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
  USE klist,                ONLY : lgauss, ltetra, nelec, ngk
  USE gvecs,                ONLY : doublegrid
  USE scf,                  ONLY : rho, v, vrs
  USE fft_base,             ONLY : dfftp, dffts
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE wvfct,                ONLY : nbnd, npwx
  USE uspp,                 ONLY : okvan, nkb, deeq_nc
  USE uspp_param,           ONLY : nhm
  USE becmod,               ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, becp
  USE buffers,              ONLY : save_buffer, get_buffer
  USE noncollin_module,     ONLY : npol, nspin_mag, noncolin, domag
  USE paw_variables,        ONLY : okpaw
  USE paw_onecenter,        ONLY : paw_dpotential
  USE paw_symmetry,         ONLY : paw_dusymmetrize, paw_dumqsymmetrize
  USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE qpoint,               ONLY : nksq, ikks, xq
  USE control_lr,           ONLY : lgamma
  USE units_lr,             ONLY : iuwfc, lrwfc
  USE lrus,                 ONLY : int3, int3_nc, int3_paw, becp1
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE fft_helper_subroutines
  USE fft_interfaces,       ONLY : fft_interpolate
  USE lr_symm_base,         ONLY : irotmq, minus_q, nsymq, rtau
  USE ldaU_lr,              ONLY : dnsscf, vh_u_save, vh_uv_save
  USE ldaU_hp,              ONLY : thresh_init, dns0, trace_dns_tot_old, &
                                   conv_thr_chi_best, iter_best, niter_max, nmix, &
                                   alpha_mix, code, lrdvwfc, iudvwfc, no_metq0
  USE apply_dpot_mod,       ONLY : apply_dpot_allocate, apply_dpot_deallocate
  USE efermi_shift,         ONLY : ef_shift, def
  USE response_kernels,     ONLY : sternheimer_kernel
  USE qpoint_aux,           ONLY : ikmks, ikmkmqs, becpt
  USE lsda_mod,             ONLY : nspin
  USE lr_nc_mag,            ONLY : lr_apply_time_reversal, deeq_nc_save, int3_nc_save
  USE lr_symm_base,         ONLY : lr_npert, upert, upert_mq
  USE ldaU,                 ONLY : lda_plus_u_kind, nsg, v_nsg
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na ! number of the perturbed atom
  INTEGER, INTENT(IN) :: iq ! number of the q point
  !
  REAL(DP),  ALLOCATABLE :: h_diag (:,:) ! diagonal part of the Hamiltonian
  !
  REAL(DP) :: thresh, & ! convergence threshold
              averlt, & ! average number of iterations
              dr2       ! self-consistency error
  !
  REAL(DP) :: dos_ef
  !! density of states at the Fermi level
  !
  REAL(DP), ALLOCATABLE :: becsum1(:,:,:)
  ! 
  COMPLEX(DP), ALLOCATABLE, TARGET :: dvscfin(:,:)
  ! change of the scf potential (input)
  !
  COMPLEX(DP), POINTER :: dvscfins(:,:,:)
  ! change of the scf potential (smooth part only)
  !
  COMPLEX(DP), ALLOCATABLE :: drhoscf  (:,:), &
                              drhoscfh (:,:), &
                              dvscfout (:,:)
  ! change of rho / scf potential (output)
  !
  COMPLEX(DP), ALLOCATABLE :: &
     ldos (:,:),              & ! local density of states at Ef
     ldoss (:,:),             & ! as above, without augmentation charges
     dbecsum (:,:,:,:),       & ! the derivative of becsum
     dbecsum_nc(:,:,:,:,:,:), &
     dbecsum_aux (:,:,:,:),   &
     aux2 (:,:),              & ! auxiliary arrays
     mixin(:), mixout(:)        ! auxiliary arrays for mixing of the response potential
 
  COMPLEX(DP), ALLOCATABLE :: t(:,:,:,:), tmq(:,:,:)
  ! PAW: auxiliary arrays
 
  LOGICAL :: all_conv
  !! True if sternheimer_kernel is converged at all k points
  LOGICAL :: lmetq0,     & ! true if xq=(0,0,0) in a metal
             convt,      & ! not needed for HP 
             convt_chi     ! used instead of convt to control the convergence

  REAL(DP), PARAMETER :: tr2 = 1.D-30 ! threshold parameter

  INTEGER :: iter,       & ! counter on iterations
             ik, ikk,    & ! counter on k points
             ndim,       &
             is,         & ! counter on spin polarizations
             isym,       & ! counter on symmetries
             npw,        & ! number of plane waves at k
             nsolv,      & ! number of linear systems
             isolv,      & ! counter on linear systems    
             ikmk,       & ! index of mk
             nrec          ! the record number for dvpsi
  INTEGER :: nnr 

  REAL(DP) :: tcpu, get_clock ! timing variables
  CHARACTER(LEN=256) :: flmixdpot = 'mixd'
  !
  CALL start_clock ('hp_solve_linear_system')
  !
  WRITE( stdout,*) "     =--------------------------------------------="
  WRITE( stdout, '(13x,"    SOLVE THE LINEAR SYSTEM")')
  WRITE( stdout,*) "     =--------------------------------------------="
  !
  ! Allocate arrays for the SCF density/potential
  !
  ALLOCATE (drhoscf (dffts%nnr, nspin_mag))
  ALLOCATE (drhoscfh(dfftp%nnr, nspin_mag))
  ALLOCATE (dvscfin (dfftp%nnr, nspin_mag))
  ALLOCATE (dvscfout(dfftp%nnr, nspin_mag))
  !
  dvscfin = (0.0_DP, 0.0_DP)
  IF (doublegrid) THEN
     ALLOCATE (dvscfins(dffts%nnr, nspin_mag, 1))
  ELSE
     dvscfins(1:dffts%nnr, 1:nspin_mag, 1:1) => dvscfin
  ENDIF
  nnr = dfftp%nnr
  !$acc enter data create(dvscfins(1:nnr, 1:nspin_mag, 1))
  !
  ! USPP-specific allocations
  !
  IF (okvan) ALLOCATE (int3 ( nhm, nhm, nat, nspin_mag, 1))
  IF (okpaw) ALLOCATE (int3_paw ( nhm, nhm, nat, nspin_mag, 1))
  CALL allocate_bec_type_acc (nkb, nbnd, becp)
  !
  ALLOCATE (dbecsum((nhm*(nhm+1))/2, nat, nspin_mag, 1))
  !
  ! Set symmetry representation in lr_symm_base
  !
  lr_npert = 1
  ALLOCATE(upert(lr_npert, lr_npert, nsymq))
  DO isym = 1, nsymq
     upert(1, 1, isym) = (1.d0, 0.d0)
  ENDDO
  IF (minus_q) THEN
     ALLOCATE(upert_mq(lr_npert, lr_npert))
     upert_mq(1, 1) = (1.d0, 0.d0)
  ENDIF ! minus_q
  !
  nsolv=1
  IF (noncolin.AND.domag) nsolv=2   
  !
  IF (okvan.and.noncolin) ALLOCATE(int3_nc( nhm, nhm, nat, nspin, 1))
  IF (noncolin) ALLOCATE (dbecsum_nc (nhm,nhm, nat , nspin , 1, nsolv))
  !
  IF (noncolin.and.domag.and.okvan) THEN
    ALLOCATE (int3_nc_save( nhm, nhm, nat, nspin_mag, 1, 2))
    ALLOCATE (dbecsum_aux ( (nhm * (nhm + 1))/2 , nat , nspin_mag , 1))
  ENDIF
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
  CALL apply_dpot_allocate()
  ALLOCATE (aux2(npwx*npol, nbnd))
  ALLOCATE (h_diag(npwx*npol, nbnd))
  ALLOCATE (trace_dns_tot_old(nat))
  trace_dns_tot_old(:) = (0.d0, 0.d0)
  !
  convt     = .FALSE.
  convt_chi = .FALSE.
  !
  ! If q=0 for a metal: allocate and compute local DOS and DOS at Ef
  !
  lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
  !
  ! If the user specified in the input that no_metq0=.true. this means
  ! we want to force remove the metallic term at q=0 (e.g. for magnetic insulators
  ! when the smearing is used). We remove the last term in Eq. (22) 
  ! in PRB 103, 045141 (2021) (and also for the response charge density), 
  ! otherwise this term can diverge for magnetic insulators
  ! due to a division by DOS(E_Fermi) which can be vanishing.
  !
  IF (no_metq0) lmetq0 = .FALSE.
  !
  IF (lmetq0) THEN
     ALLOCATE (ldos (dfftp%nnr, nspin_mag))
     ALLOCATE (ldoss(dffts%nnr, nspin_mag))
     ALLOCATE (becsum1 ( (nhm * (nhm + 1))/2, nat, nspin_mag))
     CALL localdos (ldos, ldoss, becsum1, dos_ef)
     IF (.NOT.okpaw) DEALLOCATE (becsum1)
  ENDIF
  !
  ! Compute dV_bare * psi and write to buffer iubar
  !
  DO ik = 1, nksq
     !
     ikk  = ikks(ik)
     npw  = ngk(ikk)
     !
     IF (lsda) current_spin = isk(ikk)
     ! 
     DO isolv = 1, nsolv
        !
        IF (isolv == 1) THEN
           ikmk = ikks(ik)
        ELSE
           ikmk = ikmks(ik)
        ENDIF
        !
        ! Read unperturbed KS wavefuctions psi(k) and psi(k+q)
        !
        IF (nksq > 1 .OR. nsolv == 2) &
           CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
        !
        ! Computes (iter=1) or reads (iter>1) the action of the perturbing
        ! potential on the unperturbed KS wavefunctions: |dvpsi> = dV_pert * |evc>
        ! See Eq. (46) in Ref. [1]
        !
        nrec = ik + (isolv - 1) * nksq
        CALL hp_dvpsi_pert(ik, nrec)
        !
     ENDDO
     !
  ENDDO
  !
  ! The loop of the linear-response calculation
  !
  DO iter = 1, niter_max
     !
     WRITE(stdout,'(/6x,"atom #",i3,3x,"q point #",i4,3x,"iter # ",i3)') na, iq, iter
     !
     drhoscf(:,:)     = (0.d0, 0.d0)
     dvscfout(:,:)    = (0.d0, 0.d0)
     dbecsum(:,:,:,:) = (0.d0, 0.d0)
     !
     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     !
     DO isolv = 1, nsolv
        !
        !  change the sign of the magnetic field if required
        !
        IF (isolv == 2) THEN
           IF (lda_plus_u_kind == 0) THEN
              v%ns_nc (:,:,:,:) = vh_u_save(:,:,:,:,2)
           ELSEIF (lda_plus_u_kind == 2) THEN
              v_nsg(:,:,:,:,:) = vh_uv_save(:,:,:,:,:,2)
           ENDIF
        ENDIF
        !
        ! set threshold for iterative solution of the linear system
        !
        IF ( iter == 1 ) THEN
           ! Starting threshold for iterative solution of the linear system.
           ! A strickt threshold for the first iteration is needed,
           ! because we need dns0 to very high precision.
           thresh = thresh_init * nelec
        ELSE
           ! Threshold for iterative solution of the linear system.
           ! We start with not a strict threshold for iter=2, and then
           ! it decreases with iterations.
           thresh = MIN (1.D-1 * SQRT(dr2), 1.D-2)
         ENDIF
        !
        ! Compute drhoscf, the charge density response to the total potential
        !
         CALL sternheimer_kernel(iter==1, isolv==2, 1, lrdvwfc, iudvwfc, &
            thresh, dvscfins, all_conv, averlt, drhoscf, dbecsum,&
            dbecsum_nc(:,:,:,:,:,isolv), exclude_hubbard=.TRUE.)
        !
        IF ((.NOT. all_conv) .AND. (iter == 1)) THEN
           WRITE(stdout, '(6x, "sternheimer_kernel not converged. Try to increase thresh_init.")')
        ENDIF
        !
        !  reset the original magnetic field if it was changed
        !
        IF (isolv == 2) THEN
           IF (lda_plus_u_kind == 0) THEN
               v%ns_nc (:,:,:,:) = vh_u_save(:,:,:,:,1)
           ELSEIF(lda_plus_u_kind == 2) THEN
               v_nsg (:,:,:,:,:) = vh_uv_save(:,:,:,:,:,1)
           ENDIF
        ENDIF
        !
     ENDDO ! isolv
     !
     IF (nsolv==2) THEN
        drhoscf = drhoscf / 2.0_DP
        dbecsum = dbecsum / 2.0_DP
        dbecsum_nc = dbecsum_nc / 2.0_DP
     ENDIF
     !
     ! USPP: The calculation of dbecsum is distributed across processors (see addusdbec)
     ! Sum over processors the contributions coming from each slice of bands
     !
     IF (noncolin) then
        CALL mp_sum ( dbecsum_nc, intra_pool_comm )
     ELSE
        CALL mp_sum ( dbecsum, intra_pool_comm )
     ENDIF
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
     !  In the noncolinear, spin-orbit case rotate dbecsum
     !
     IF (noncolin.and.okvan) THEN
        CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 1)
        IF (nsolv==2) THEN
           dbecsum_aux=(0.0_DP,0.0_DP)
           CALL set_dbecsum_nc(dbecsum_nc(1,1,1,1,1,2), dbecsum_aux, 1)
           dbecsum(:,:,1,:) = dbecsum(:,:,1,:) + dbecsum_aux(:,:,1,:)
           dbecsum(:,:,2:4,:) = dbecsum(:,:,2:4,:) - dbecsum_aux(:,:,2:4,:)
        ENDIF
     ENDIF
     !
     ! USPP: Compute the total response charge density (standard term + US term)
     !
     IF (okvan) CALL lr_addusddens (drhoscfh, dbecsum)
     !
     call mp_sum ( drhoscf, inter_pool_comm )
     CALL mp_sum ( drhoscfh, inter_pool_comm ) 
     IF (okpaw) CALL mp_sum ( dbecsum, inter_pool_comm )
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
           CALL ef_shift(1, dos_ef, ldos, drhoscfh, dbecsum=dbecsum, becsum1=becsum1)
        ELSE
           CALL ef_shift(1, dos_ef, ldos, drhoscfh)
        ENDIF
        !
        ! Check that def is not too large (it is in Ry). 
        !
        ! Note: This check might be skipped, like in the PHonon code
        IF ( ABS(DBLE(def(1))) < 1.0d-18 .OR. ABS(DBLE(def(1))) > 5.0d0 ) THEN
           !
           IF (noncolin) THEN
              IF (ABS(DBLE(def(1))) > 5.0d0) THEN
                 WRITE( stdout, '(/6x,"WARNING: The Fermi energy shift too big!")')
                 WRITE( stdout, '(6x, "   DOS(E_Fermi) = ",1x,2e12.4)') dos_ef
                 WRITE( stdout, '(6x, "   Fermi_shift  = ",1x,2e12.4)') DBLE(def(1))
                 CALL hp_stop_smoothly (.FALSE.)
              ENDIF
           ELSE
              WRITE( stdout, '(/6x,"WARNING: The Fermi energy shift is zero or too big!")')
              WRITE( stdout, '(6x, "This may happen in two cases:")')
              WRITE( stdout, '(6x, "1. The DOS at the Fermi level is too small:")')
              WRITE( stdout, '(6x, "   DOS(E_Fermi) = ",1x,2e12.4)') dos_ef
              WRITE( stdout, '(6x, "   This means that most likely the system has a gap,")')
              WRITE( stdout, '(6x, "   and hence it should NOT be treated as a metal")')
              WRITE( stdout, '(6x, "   (otherwise numerical instabilities will appear).")')
              WRITE( stdout, '(6x, "2. Numerical instabilities due to too low cutoff")')
              WRITE( stdout, '(6x, "   for hard pseudopotentials.")')
              WRITE( stdout, '(/6x,"Stopping...")')
              WRITE( stdout, '(/6x,"Solution (for magnetic insulators):")')
              WRITE( stdout, '(6x,"Try to use the 2-step scf procedure as in HP/example02")')
              CALL hp_stop_smoothly (.FALSE.)
           ENDIF
        ENDIF    
        !
     ENDIF
     !
     ! Symmetrization of the response charge density.
     !
     CALL hp_psymdvscf (drhoscfh)
     !
     IF ( noncolin.and.domag ) CALL hp_psym_dmag( drhoscfh )
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
     CALL dv_of_drho (dvscfout)
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
           CALL fft_interpolate (dfftp, dvscfin(:,is), dffts, dvscfins(:,is,1))
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
     IF (okvan) then
        CALL newdq (dvscfin, 1)
        IF (noncolin.AND.domag) then
           !
           int3_nc_save(:,:,:,:,:,1)=int3_nc(:,:,:,:,:)
           !
           dvscfin(:,2:4) = -dvscfin(:,2:4)
           IF (okpaw) THEN 
              dbecsum(:,:,2:4,1) = -dbecsum(:,:,2:4,1)
              rho%bec(:,:,2:4) = -rho%bec(:,:,2:4)
           ENDIF
           !
           !   if needed recompute the paw coeffients with the opposite sign of
           !   the magnetic field
           !
           IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,1)
           !
           CALL newdq (dvscfin, 1)
           int3_nc_save(:,:,:,:,:,2) = int3_nc(:,:,:,:,:)
           !
           !  restore the correct sign of the magnetic field.
           !
           dvscfin(:,2:4) = -dvscfin(:,2:4)
           IF (okpaw) THEN 
              dbecsum(:,:,2:4,1) = -dbecsum(:,:,2:4,1)
              rho%bec(:,:,2:4) = -rho%bec(:,:,2:4)
           ENDIF
           !
           !  put into int3_nc the coefficient with +B
           !
           int3_nc(:,:,:,:,:)=int3_nc_save(:,:,:,:,:,1)
        ENDIF
     ENDIF
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
     ! Print the average number of iterations
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
     IF (convt_chi) EXIT
     !
  ENDDO  ! loop over the iterations iter
  !
  CALL apply_dpot_deallocate()
  DEALLOCATE (h_diag)
  DEALLOCATE (aux2)
  DEALLOCATE (dbecsum)
  DEALLOCATE (drhoscf)
  DEALLOCATE (drhoscfh)
  DEALLOCATE (dvscfin)
  DEALLOCATE (dvscfout)
  DEALLOCATE (trace_dns_tot_old)
  DEALLOCATE (upert)
  IF (minus_q) DEALLOCATE(upert_mq)
  !
  IF (ALLOCATED(dbecsum_nc)) DEALLOCATE (dbecsum_nc)
  IF (ALLOCATED(int3_nc)) DEALLOCATE(int3_nc)
  IF (ALLOCATED(int3_nc_save)) DEALLOCATE (int3_nc_save)
  IF (ALLOCATED(dbecsum_aux)) DEALLOCATE (dbecsum_aux)
  !
  !$acc exit data delete(dvscfins)
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
  CALL deallocate_bec_type_acc (becp)
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
