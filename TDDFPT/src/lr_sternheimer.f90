!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE lr_sternheimer
  !
  !    This routine generalizes to finite complex frequencies and 
  !    finite q vectors the routine solve_e of the Quantum ESPRESSO 
  !    distribution. 
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field 
  !    of finite wavevector q and complex frequency omega.
  !    It performs the following tasks:
  !     a) computes the bare potential term  e^{iqr} | psi >
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system at zero
  !        frequency or ccg_many_vectors
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !
CONTAINS

SUBROUTINE one_sternheimer_step(iu, flag)
    !    
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : e2, fpi, rytoev
    USE ions_base,              ONLY : nat
    USE io_global,              ONLY : stdout, ionode
    USE io_files,               ONLY : diropn, nwordwfc, iunwfc
    USE cell_base,              ONLY : tpiba2
    USE fft_interfaces,         ONLY : fwfft
    USE fft_interfaces,         ONLY : fft_interpolate
    USE klist,                  ONLY : lgauss, xk, wk
    USE gvecs,                  ONLY : doublegrid
    USE fft_base,               ONLY : dfftp, dffts
    USE lsda_mod,               ONLY : lsda, nspin, current_spin, isk
    USE spin_orb,               ONLY : domag
    USE wvfct,                  ONLY : nbnd, npwx, g2kin,  et
    USE klist,                  ONLY : ngk, igk_k
    USE check_stop,             ONLY : check_stop_now
    USE buffers,                ONLY : get_buffer, save_buffer
    USE wavefunctions,          ONLY : evc
    USE uspp,                   ONLY : okvan, vkb
    USE uspp_param,             ONLY : nhm
    USE noncollin_module,       ONLY : noncolin, npol, nspin_mag
    USE scf,                    ONLY : rho, v_of_0
    USE gvect,                  ONLY : gg
    USE paw_variables,          ONLY : okpaw
    USE paw_onecenter,          ONLY : paw_dpotential
    USE eqv,                    ONLY : dpsi, dvpsi, evq
    USE units_lr,               ONLY : lrwfc, iuwfc
    USE control_lr,             ONLY : lgamma, alpha_pv, nbnd_occ, &
                                       ext_recover, rec_code, &
                                       lnoloc, convt, tr2_ph, &
                                       alpha_mix, lgamma_gamma, niter_ph, &
                                       flmixdpot, rec_code_read

    USE lrus,                   ONLY : int3_paw, intq, intq_nc
    USE qpoint,                 ONLY : xq, nksq, ikks, ikqs
    USE linear_solvers,         ONLY : ccg_many_vectors
    USE dv_of_drho_lr,          ONLY : dv_of_drho
    USE mp_pools,               ONLY : inter_pool_comm
    USE mp_bands,               ONLY : intra_bgrp_comm
    USE mp_images,              ONLY : root_image, my_image_id
    USE mp,                     ONLY : mp_sum
    USE fft_helper_subroutines, ONLY : fftx_ntgrp
    USE lr_variables,           ONLY : fru, fiu, iundvpsi, iudwf, &
                                       lrdrho, iudrho, n_ipol, lr_verbosity, &
                                       chirr, chirz, chizr, chizz, epsm1, &
                                       current_w, lr1dwf, iu1dwf, itermax!, &
                                       !intq, intq_nc
    USE paw_add_symmetry,       ONLY : paw_deqsymmetrize
    USE wavefunctions,          ONLY : psic
    USE lr_sym_mod,             ONLY : psymeq
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(IN) :: flag   ! if 1 compute the charge-charge and
                                  ! charge magnetization responses
                                  ! if 2 and lsda computes the magnetization
                                  ! magnetization response
    REAL(DP) ::  thresh, anorm, averlt, dr2
    ! thresh: convergence threshold
    ! anorm : the norm of the error
    ! averlt: average number of iterations
    ! dr2   : self-consistency error
    COMPLEX(DP), ALLOCATABLE :: h_diag (:,:)
    COMPLEX(DP), ALLOCATABLE :: h_diag1 (:,:)
    REAL(DP),    ALLOCATABLE :: h_diagr (:,:)
    REAL(DP),    ALLOCATABLE :: h_dia (:,:), s_dia(:,:)
    ! h_diag: diagonal part of the Hamiltonian
    !
    COMPLEX(DP) , ALLOCATABLE, TARGET ::      &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
    COMPLEX(DP) , POINTER ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
    COMPLEX(DP) , ALLOCATABLE ::   &
                   dpsi1(:,:),   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   drhoscfout (:,:), & ! change of the scf charge (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   mixin(:), mixout(:), &  ! auxiliary for paw mixing
                   aux1 (:,:),  ps (:,:), &
                   tg_dv(:,:), &
                   tg_psic(:,:), aux2(:,:), dvpsi1(:,:)
    !
    LOGICAL :: conv_root, exst, all_done_asyn
    ! conv_root: true if linear system is converged
    INTEGER :: kter, iter0, ipol, ibnd, iter, lter, ik, ikk, ikq, &
               ig, is, nrec, ndim, npw, npwq, ios
    ! counters
    INTEGER :: ltaver, lintercall, incr, jpol, v_siz, nwordd0psi
    REAL(DP) :: xqmod2, alpha_pv0
    !
    REAL(DP) :: tcpu, get_clock
    ! timing variables
    !
    COMPLEX(DP) :: w  !frequency
    REAL(DP) :: aa, weight
    LOGICAL :: ldpsi1
    !
    EXTERNAL ch_psi_all, cg_psi
    EXTERNAL ch_psi_all_complex, ccg_psi
    COMPLEX(DP) :: scal_prod
    !

    CALL start_clock ('stern_step')

    w=CMPLX(fru(iu),fiu(iu))
    ldpsi1=ABS(w)>1.D-7
    alpha_pv0=alpha_pv
    alpha_pv=alpha_pv0 + REAL(w)
    !
    ALLOCATE (dvscfin( dfftp%nnr, nspin_mag, 1))
    IF (doublegrid) THEN
       ALLOCATE (dvscfins(dffts%nnr, nspin_mag, 1))
    ELSE
       dvscfins => dvscfin
    ENDIF
    ALLOCATE (dvscfout(dfftp%nnr, nspin_mag, 1))
    ALLOCATE (drhoscfout(dfftp%nnr, nspin_mag))
    IF (okpaw) THEN
       ALLOCATE (mixin(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
       ALLOCATE (mixout(dfftp%nnr*nspin_mag+(nhm*(nhm+1)*nat*nspin_mag)/2) )
    ENDIF
    ALLOCATE (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 1))
    IF (noncolin) ALLOCATE (dbecsum_nc (nhm, nhm, nat, nspin, 1))
    IF (ldpsi1) THEN
       ALLOCATE (dpsi1(npwx*npol,nbnd))
       ALLOCATE (dvpsi1(npwx*npol,nbnd))
       ALLOCATE (h_diag(npwx*npol, nbnd))
       ALLOCATE (h_diag1(npwx*npol, nbnd))
       ALLOCATE (h_dia(npwx,npol))
       ALLOCATE (s_dia(npwx,npol))
    ELSE
       ALLOCATE (h_diagr(npwx*npol, nbnd))
    ENDIF
    ALLOCATE (aux1(dffts%nnr,npol))
    ALLOCATE (aux2(npwx*npol, nbnd))
    IF (okpaw) mixin=(0.0_DP,0.0_DP)
    !

dvpsi =(0.0d0, 0.0d0)

!    IF (rec_code_read == -20.AND.ext_recover) then
!       ! restarting in Electric field calculation
!       IF (okpaw) THEN
!          CALL read_rec(dr2, iter0, 1, dvscfin, dvscfins, dvscfout, dbecsum)
!          CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
!                      mixin, dvscfin, dbecsum, ndim, -1 )
!       ELSE
!          CALL read_rec(dr2, iter0, 1, dvscfin, dvscfins)
!       ENDIF
!    ELSEIF (rec_code_read > -20 .AND. rec_code_read <= -10) then
!       ! restarting in Raman: proceed
!       convt = .true.
!    ELSE
       convt = .false.
       iter0 = 0
!    ENDIF
    !
    incr=1
    IF ( dffts%has_task_groups ) THEN
       !
       v_siz =  dffts%nnr_tg
       ALLOCATE( tg_dv   ( v_siz, nspin_mag ) )
       ALLOCATE( tg_psic( v_siz, npol ) )
       incr = fftx_ntgrp(dffts)
       !
    ENDIF
    !
!    IF ( ionode .AND. fildrho /= ' ') THEN
!       INQUIRE (UNIT = iudrho, OPENED = exst)
!       IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
!       CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
!    ENDIF
    IF (rec_code_read > -20) convt=.TRUE.
    !
    IF (convt) go to 155
    !
    IF ((lgauss.and..not.ldpsi1)) &
            CALL errore ('solve_eq', 'insert a finite frequency', 1)
    !
    IF (lr_verbosity > 5) THEN
       WRITE(stdout,'("<lr_sternheimer_one_step>")')
    ENDIF
    !
    IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))

!    nwordd0psi = 2 * nbnd * npwx * npol * nksq
!    CALL diropn ( iund0psi, 'dvpsi.', nwordd0psi, exst)
!
!    CALL diropn ( iudwf, 'dwf', nwordd0psi, exst)
!    CALL diropn ( iu1dwf, 'mwf', nwordd0psi, exst)
    !
    ! Loop over the iterations
    !
    DO kter = 1, itermax
       !
       FLUSH( stdout)
       iter = kter + iter0

       ltaver = 0
       lintercall = 0
       !
       dvscfout(:,:,:) = (0.0d0, 0.0d0)
       dbecsum(:,:,:,:) = (0.0d0, 0.0d0)
       !
       IF (noncolin) dbecsum_nc = (0.0d0, 0.0d0)
       !
       DO ik = 1, nksq
          !
          ikk  = ikks(ik)
          ikq  = ikqs(ik)
          npw  = ngk(ikk)
          npwq = ngk(ikq)
        if (lsda) current_spin = isk (ikk)
          !
          ! Calculate beta-functions vkb at k+q (Kleinman-Bylander projectors)
          ! The vkb's are needed for the non-local potential in h_psi,
          ! and for the ultrasoft term.
          !
          CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkb)
          !
          ! Read unperturbed wavefuctions evc (wfct at k) 
          ! and evq (wfct at k+q)
          !
          IF (nksq > 1) THEN 
             CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
             CALL get_buffer (evq, nwordwfc, iunwfc, ikq)
          ENDIF
          !
          !  Compute the kinetic energy g2kin: (k+q+G)^2
          !
          CALL g2_kin(ikq)
          !
          ! IF omega non zero
          !
          IF (ldpsi1) THEN
             h_diag=(0.0_DP,0.0_DP)
             h_diag1=(0.0_DP,0.0_DP)
             h_dia=0.0_DP
             s_dia=0.0_DP
             CALL usnldiag( npwq, h_dia, s_dia )
             !

             DO ibnd = 1, nbnd_occ (ikk)
                !
                DO ig = 1, npwq
                   aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                      (et(ibnd,ikk)+w)*s_dia(ig,1)
                   IF (ABS(aa)<1.0_DP) aa=1.0_DP
                   h_diag(ig,ibnd)=CMPLX(1.0d0, 0.d0,kind=DP) / aa
                   aa=g2kin(ig)+v_of_0+h_dia(ig,1)- &
                      (et(ibnd,ikk)-w)*s_dia(ig,1)
                   IF (ABS(aa)<1.0_DP) aa=1.0_DP
                   h_diag1(ig,ibnd)=CMPLX(1.0d0, 0.d0,kind=DP) / aa
                ENDDO
                !
                IF (noncolin) THEN
                   DO ig = 1, npwq
                      aa=g2kin(ig)+v_of_0+h_dia(ig,2)- &
                         (et(ibnd,ikk)+w)*s_dia(ig,2)
                      IF (ABS(aa)<1.0_DP) aa=1.0_DP
                      h_diag(ig+npwx,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa
                      aa=g2kin(ig)+v_of_0+h_dia(ig,2)- &
                         (et(ibnd,ikk)-w)*s_dia(ig,2)
                      IF (ABS(aa)<1.0_DP) aa=1.0_DP
                      h_diag1(ig+npwx,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / aa
                   ENDDO
                ENDIF
             ENDDO
          ELSE
             CALL h_prec (ik, evc, h_diagr)
             !
             DO ibnd = 1, nbnd_occ (ikk)
                !
                DO ig = 1, npwq
                   aa=1.0_DP / h_diagr(ig,ibnd)-et(ibnd,ikk)-REAL(w,KIND=DP)
                   h_diagr(ig,ibnd)=1.d0 /max(1.0d0,aa)
                ENDDO
                IF (noncolin) THEN
                   DO ig = 1, npwq
                      h_diagr(ig+npwx,ibnd)= h_diagr(ig,ibnd)
                   ENDDO
                ENDIF
             ENDDO
          ENDIF
          !
          ! do over polarization
          !
          DO ipol = 1, n_ipol
             !
             nrec = (ipol - 1) * nksq + ik
             !
             IF (iter ==1) THEN
                !
                !  At the first iteration dvbare_q*psi_kpoint is calculated
                !  and written to file
                !
                CALL dveqpsi_us (ik)
                !
                !  with flag=2 the perturbation is a magnetic field along z
                !
                IF (lsda.AND.current_spin==2.AND.flag==2) dvpsi=-dvpsi
!                call save_buffer (dvpsi, lrbar, iubar, nrec)

                CALL save_buffer (dvpsi, nwordwfc, iundvpsi, nrec)
                !
             ELSE
                !
                ! After the first iteration dvbare_q*psi_kpoint is read from file
                !
!                call get_buffer (dvpsi, lrbar, iubar, nrec)
                CALL get_buffer (dvpsi, nwordwfc, iundvpsi, nrec)
                !
                ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
                ! dvscf_q from previous iteration (mix_potential)
                !
                IF ( dffts%has_task_groups ) THEN
                   IF (noncolin) THEN
                      CALL tg_cgather( dffts, dvscfins(:,1,ipol), &
                                                       tg_dv(:,1))
                      IF (domag) THEN
                         DO jpol=2,4
                            CALL tg_cgather( dffts, dvscfins(:,jpol,ipol), &
                                                            tg_dv(:,jpol))
                         ENDDO
                      ENDIF
                   ELSE
                      CALL tg_cgather( dffts, dvscfins(:,current_spin,ipol), &
                                                               tg_dv(:,1))
                   ENDIF
                ENDIF
                !
                aux2=(0.0_DP,0.0_DP)
                DO ibnd = 1, nbnd_occ (ikk), incr
                   IF ( dffts%has_task_groups ) THEN
                      CALL cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, &
                                        nbnd_occ (ikk) )
                      CALL apply_dpot(v_siz, tg_psic, tg_dv, 1)
                      CALL cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, &
                                        nbnd_occ (ikk))
                   ELSE
                      CALL cft_wave (ik, evc (1, ibnd), aux1, +1)
                      CALL apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipol), &
                                                               current_spin)
                      CALL cft_wave (ik, aux2 (1, ibnd), aux1, -1)
                   ENDIF

                ENDDO
                !
                dvpsi=dvpsi+aux2
                !
                CALL adddvscf(ipol,ik)
                !
             ENDIF
             !
             ! Orthogonalize dvpsi to valence states: ps = <evq|dvpsi>
             !
             IF (ldpsi1) THEN
                dvpsi1=dvpsi
                CALL orthogonalize_omega(dvpsi1, evq, ikk, ikq, dpsi, npwq, -w)
             ENDIF
             CALL orthogonalize_omega(dvpsi, evq, ikk, ikq, dpsi, npwq, w)
             ! 
             IF (iter == 1) THEN
                !
                !  At the first iteration dpsi and dvscfin are set to zero,
                !
                dpsi(:,:)=(0.d0,0.d0)
                IF (ldpsi1) dpsi1(:,:)=(0.d0,0.d0)
                dvscfin(:,:,:)=(0.d0,0.d0)
                !
                ! starting threshold for the iterative solution of the linear
                ! system
                !
                thresh = 1.d-2
                IF (lnoloc) thresh = 1.d-5
                !
             ELSE
                !
                ! starting value for  delta_psi is read from iudwf
                !
                nrec = (ipol - 1) * nksq + ik
                CALL get_buffer (dpsi, nwordwfc, iudwf, nrec)
!                call get_buffer (dpsi, lrdwf, iudwf, nrec)
                IF (ldpsi1) CALL get_buffer (dpsi1, nwordwfc, iu1dwf, nrec)
                !
                ! threshold for iterative solution of the linear system
                !
                thresh = min (0.1d0 * sqrt (dr2), thresh)
                !
             ENDIF
             !
             ! iterative solution of the linear system (H-e)*dpsi=dvpsi
             ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
             !
             conv_root = .true.
             !
             current_w=w
             IF (ldpsi1) THEN
                !
                ! Complex or imaginary frequency. Use bicojugate gradient.
                !

                CALL ccgsolve_all (ch_psi_all_complex,ccg_psi,et(1,ikk),dvpsi,dpsi, &
                                    h_diag,npwx,npwq,thresh,ik,lter,conv_root,anorm,&
                                                        nbnd_occ(ikk),npol,current_w)

                !
             ELSE
                !
                ! zero frequency. The standard QE solver
                !
                CALL cgsolve_all (ch_psi_all,cg_psi,et(1,ikk),dvpsi,dpsi, &
                  h_diagr,npwx,npwq,thresh,ik,lter,conv_root,anorm,&
                                                          nbnd_occ(ikk),npol)
                !
             ENDIF
             !
             ltaver = ltaver + lter
             lintercall = lintercall + 1
             IF (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                  &         ' solve_e: root not converged ',es10.3)") ik &
                  &, ibnd, anorm
             !
             ! writes delta_psi on iunit iudwf, k=kpoint,
             !
             nrec = (ipol - 1) * nksq + ik
             CALL save_buffer (dpsi, nwordwfc, iudwf, nrec)
!             call save_buffer(dpsi, lrdwf, iudwf, nrec)
             !

             ! calculates dvscf, sum over k => dvscf_q_ipert
             !
             weight=wk(ikk)
             IF (ldpsi1) THEN
                !
                ! complex frequency, two wavefunctions must be computed
                !
                weight=wk(ikk)/2.0_DP
                !
                ! In this case compute also the wavefunction at frequency -w.
                !
                current_w=-w

                CALL ccgsolve_all (ch_psi_all_complex,ccg_psi,et(1,ikk),dvpsi1,dpsi1, &
                                    h_diag1,npwx,npwq,thresh,ik,lter,conv_root,anorm,&
                                                          nbnd_occ(ikk),npol,current_w)

                ltaver = ltaver + lter
                lintercall = lintercall + 1
                IF (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                  &         ' solve_e: root not converged ',es10.3)") ik &
                  &, ibnd, anorm
                !
                ! writes delta_psi on iunit iudwf, k=kpoint,
                !
                nrec = (ipol - 1) * nksq + ik
                CALL save_buffer(dpsi1, nwordwfc, iu1dwf, nrec)
                !
                ! calculates dvscf, sum over k => dvscf_q_ipert
                !
                CALL DAXPY(npwx*nbnd_occ(ikk)*npol*2, 1.0_DP, dpsi1, 1, dpsi, 1)
                !
             ENDIF
             IF (noncolin) THEN
                CALL incdrhoscf_nc(dvscfout(1,1,ipol), weight, ik, &
                                   dbecsum_nc(1,1,1,1,ipol), dpsi, 1.0d0)
             ELSE

                CALL incdrhoscf (dvscfout(1,current_spin,ipol), weight, &
                           ik, dbecsum(1,1,current_spin,ipol), dpsi)
             ENDIF
             !
          ENDDO   ! on polarizations
          !
!          IF ( with_asyn_images.AND.my_image_id==root_image.AND.ionode ) &
!                             CALL asyn_master(all_done_asyn)
       ENDDO ! on k points
       !

       current_w=w
       !
       !  The calculation of dbecsum is distributed across processors
       !  (see addusdbec) - we sum over processors the contributions
       !  coming from each slice of bands
       !
       IF (noncolin) THEN
          CALL mp_sum ( dbecsum_nc, intra_bgrp_comm )
       ELSE
          CALL mp_sum ( dbecsum, intra_bgrp_comm )
       ENDIF
       !
       IF (doublegrid) then
          DO is=1,nspin_mag
             CALL fft_interpolate (dffts, dvscfout(:,is,1), dfftp, &
                                                   dvscfout(:,is,1))
          ENDDO
       ENDIF
       !
       IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 1)
       !
       CALL addusddenseq (dvscfout, dbecsum)
       !
       !   dvscfout contains the (unsymmetrized) linear charge response
       !   for the three polarizations - symmetrize it
       !
       CALL mp_sum ( dvscfout, inter_pool_comm )
       !
       IF (okpaw) CALL mp_sum ( dbecsum, inter_pool_comm )
       !
       IF (.not. lgamma_gamma) THEN
          CALL psymeq (dvscfout)
       ENDIF
       !
       drhoscfout(:,:)=dvscfout(:,:,1)
       !
       !   save the symmetrized linear charge response to file
       !   calculate the corresponding linear potential response
       !
       IF (lnoloc) THEN
!          CALL dv_of_drho_nlf (dvscfout (1, 1, 1))
       ELSE
          CALL dv_of_drho (dvscfout (1, 1, 1), .FALSE.)
       ENDIF
       !
       !   mix the new potential with the old
       !
       IF (okpaw) THEN
          !
          !  In this case we mix also dbecsum
          !
          CALL setmixout(dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag)/2, &
                         mixout, dvscfout, dbecsum, ndim, -1 )
          CALL mix_potential_eels (2*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, tr2_ph/npol, iter, flmixdpot, &
                         convt)
          CALL setmixout(dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag)/2, &
                         mixin, dvscfin, dbecsum, ndim, 1 )
       ELSE

          CALL mix_potential_eels(2*dfftp%nnr*nspin_mag, dvscfout, dvscfin,  &
               alpha_mix (kter), dr2,  tr2_ph/npol, iter, flmixdpot, convt)

       ENDIF
       !
       IF (doublegrid) then
          DO is=1,nspin_mag
             CALL fft_interpolate (dfftp, dvscfin(:,is,1), dffts, &
                                                 dvscfins(:,is,1))
          ENDDO
       ENDIF
       !
       IF (okpaw) THEN
          IF (noncolin.AND.domag) THEN
!             CALL PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
          ELSE
             !
             ! The presence of c.c. in the formula gives a factor 2.0
             !
             dbecsum=2.0_DP * dbecsum
             IF (.NOT. lgamma_gamma) CALL paw_deqsymmetrize(dbecsum)
             CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,1)
          ENDIF
       ENDIF
       !
       CALL newdq(dvscfin,1)
       !
       CALL mp_sum(ltaver,inter_pool_comm)
       CALL mp_sum(lintercall,inter_pool_comm)
       !
       averlt = DBLE (ltaver) / DBLE (lintercall)
       !
!       tcpu = get_clock ('PHONON')
       tcpu = get_clock ('ccgsolve')
       WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
            &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
       WRITE( stdout, "(5x,' thresh=',es10.3, ' alpha_mix = ',f6.3, &
            &      ' |ddv_scf|^2 = ',es10.3 )") thresh, alpha_mix (kter), dr2
       !
       FLUSH( stdout )
       !
       ! rec_code: state of the calculation
       ! rec_code=-20 Electric Field
       !
       rec_code=-20
!       IF (okpaw) THEN
!          CALL write_rec('solve_e...', 0, dr2, iter, convt, 1, dvscfin, &
!                                                      dvscfout, dbecsum)
!       ELSE
!          CALL write_rec('solve_e...', 0, dr2, iter, convt, 1, dvscfin)
!       ENDIF
       !
!       IF (check_stop_now()) CALL stop_smoothly_ph (.false.)
       !
       IF (convt) goto 155
       !
    ENDDO  ! do over iter_ph
    !
155 CONTINUE
    !
    !  compute here the susceptibility and the inverse of the dielectric
    !  constant
    !
    !  CALL compute_susceptibility(drhoscfout)
    !
!    CLOSE( UNIT = iund0psi)

    DO is=1,nspin_mag
       CALL fwfft ('Rho', drhoscfout(:,is), dfftp)
    ENDDO
    !
    IF (flag==1) THEN
       chirr(iu)=(0.0_DP,0.0_DP)
       chizr(iu)=(0.0_DP,0.0_DP)
       epsm1(iu)=(0.0_DP,0.0_DP)
    ELSE
       chirz(iu)=(0.0_DP,0.0_DP)
       chizz(iu)=(0.0_DP,0.0_DP)
    ENDIF
    !
    xqmod2=(xq(1)**2+xq(2)**2+xq(3)**2)*tpiba2
    !
    IF (ABS(gg(1))<1.d-8) THEN
       IF (flag==1) THEN
          chirr(iu) = drhoscfout(dfftp%nl(1),1)
          IF (lsda) chirr(iu) = chirr(iu) + drhoscfout(dfftp%nl(1),2)
          epsm1(iu) = CMPLX(1.0_DP,0.0_DP)+ chirr(iu)*fpi*e2/xqmod2
          IF (lsda) chizr(iu) = drhoscfout(dfftp%nl(1),1) - &
                                drhoscfout(dfftp%nl(1),2)
       ELSEIF (lsda) THEN
          chizz(iu)=drhoscfout(dfftp%nl(1),1)-drhoscfout(dfftp%nl(1),2)
          chirz(iu)=drhoscfout(dfftp%nl(1),1)+drhoscfout(dfftp%nl(1),2)
       ENDIF
    ENDIF
    !
    IF (flag==1) THEN
       CALL mp_sum(epsm1(iu),intra_bgrp_comm)
       CALL mp_sum(chirr(iu),intra_bgrp_comm)
       CALL mp_sum(chizr(iu),intra_bgrp_comm)
    ELSE
       CALL mp_sum(chizz(iu),intra_bgrp_comm)
       CALL mp_sum(chirz(iu),intra_bgrp_comm)
    ENDIF
    !
    IF (flag==1) THEN
       WRITE(stdout, '(/,6x,"Inverse dielectric constant at &
                          &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
       WRITE(stdout, '(46x,f9.4," +",f9.4," i eV")') current_w * rytoev
       WRITE(stdout,'(/,6x,"epsilon^-1(q,w) =",2f15.6)') epsm1(iu)
       !
       WRITE( stdout, '(/,6x,"Charge-charge susceptibility:")')
       !
       WRITE(stdout,'(/,6x,"chirr(q,w) =",2f15.6)') chirr(iu)
       IF (lsda) THEN
          WRITE(stdout,'(/,6x,"m_z-charge susceptibility:")')
          WRITE(stdout,'(/,6x,"chizr(q,w) =",2f15.6)') chizr(iu)
       ENDIF
       !
    ELSEIF (lsda) THEN
       WRITE( stdout, '(/,6x,"m_z - m_z susceptibility at &
                       &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
       WRITE( stdout, '(43x,f9.4," +",f9.4," i eV")') current_w * rytoev
       WRITE(stdout,'(/,6x,"chizz(q,w) =",2f15.6)') chizz(iu)
       WRITE(stdout,'(/,6x,"chirz(q,w) =",2f15.6)') chirz(iu)
    ENDIF
    !
    deallocate (aux1)
    IF (ldpsi1) THEN
       deallocate (dpsi1)
       deallocate (dvpsi1)
       deallocate (h_diag)
       deallocate (h_diag1)
       deallocate (h_dia)
       deallocate (s_dia)
    ELSE
       deallocate (h_diagr)
    ENDIF
    deallocate (dbecsum)
    deallocate (dvscfout)
    IF (okpaw) THEN
       DEALLOCATE(mixin)
       DEALLOCATE(mixout)
    ENDIF
    deallocate (drhoscfout)
    if (doublegrid) deallocate (dvscfins)
    deallocate (dvscfin)
    if (noncolin) deallocate(dbecsum_nc)
    deallocate(aux2)
    IF ( dffts%has_task_groups ) THEN
       !
       DEALLOCATE( tg_dv  )
       DEALLOCATE( tg_psic)
       !
    ENDIF

!    CLOSE( unit = iund0psi)
!    CLOSE( unit = iudwf)
!    CLOSE( unit = iu1dwf)

    alpha_pv=alpha_pv0
    !
    CALL stop_clock ('stern_step')
    !
    RETURN
    !
END SUBROUTINE one_sternheimer_step

END MODULE lr_sternheimer
