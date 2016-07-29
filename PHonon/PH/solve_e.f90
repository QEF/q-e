!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e
  !-----------------------------------------------------------------------
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field.
  !    It performs the following tasks:
  !     a) computes the bare potential term  x | psi >
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : prefix, diropn
  USE cell_base,             ONLY : tpiba2
  USE klist,                 ONLY : lgauss, xk, wk, ngk, igk_k
  USE gvect,                 ONLY : g
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts, dtgs
  USE fft_parallel,          ONLY : tg_cgather
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,              ONLY : domag
  USE wvfct,                 ONLY : nbnd, npwx, g2kin, et
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer, save_buffer
  USE wavefunctions_module,  ONLY : evc
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : upf, nhm
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag
  USE scf,                   ONLY : rho
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_symmetry,          ONLY : paw_desymmetrize

  USE units_ph,              ONLY : lrdwf, iudwf, lrwfc, iuwfc, lrdrho, &
                                    iudrho
  USE output,                ONLY : fildrho
  USE control_ph,            ONLY : ext_recover, rec_code, &
                                    lnoloc, convt, tr2_ph, nmix_ph, &
                                    alpha_mix, lgamma_gamma, niter_ph, &
                                    flmixdpot, rec_code_read
  USE recover_mod,           ONLY : read_rec, write_rec

  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm, ntask_groups
  USE mp,                    ONLY : mp_sum

  USE lrus,                  ONLY : int3_paw
  USE qpoint,                ONLY : nksq
  USE eqv,                   ONLY : dpsi, dvpsi
  USE control_lr,            ONLY : nbnd_occ, lgamma
  USE dv_of_drho_lr

  implicit none

  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP), allocatable :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian

  complex(DP) , allocatable, target ::      &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  complex(DP) , pointer ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  complex(DP) , allocatable ::   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   dbecsum_nc(:,:,:,:,:), & ! the becsum with dpsi
                   mixin(:), mixout(:), &  ! auxiliary for paw mixing
                   aux1 (:,:),  ps (:,:), &
                   tg_dv(:,:), &
                   tg_psic(:,:), aux2(:,:)

  complex(DP), EXTERNAL :: zdotc      ! the scalar product function

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: npw, npwq
  integer :: kter, iter0, ipol, ibnd, iter, lter, ik, ig, is, nrec, ndim, ios
  ! counters
  integer :: ltaver, lintercall, incr, jpol, v_siz

  real(DP) :: tcpu, get_clock
  ! timing variables

  external ch_psi_all, cg_psi

  call start_clock ('solve_e')
  !
  !  This routine is task group aware
  !
  allocate (dvscfin( dfftp%nnr, nspin_mag, 3))
  if (doublegrid) then
     allocate (dvscfins(dffts%nnr, nspin_mag, 3))
  else
     dvscfins => dvscfin
  endif
  allocate (dvscfout(dfftp%nnr, nspin_mag, 3))
  IF (okpaw) THEN
     ALLOCATE (mixin(dfftp%nnr*nspin_mag*3+(nhm*(nhm+1)*nat*nspin_mag*3)/2) )
     ALLOCATE (mixout(dfftp%nnr*nspin_mag*3+(nhm*(nhm+1)*nat*nspin_mag*3)/2) )
  ENDIF
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 3))
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin, 3))
  allocate (aux1(dffts%nnr,npol))
  allocate (h_diag(npwx*npol, nbnd))
  allocate (aux2(npwx*npol, nbnd))
  IF (okpaw) mixin=(0.0_DP,0.0_DP)

  if (rec_code_read == -20.AND.ext_recover) then
     ! restarting in Electric field calculation
     IF (okpaw) THEN
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins, dvscfout, dbecsum)
        CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixin, dvscfin, dbecsum, ndim, -1 )
     ELSE
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins)
     ENDIF
  else if (rec_code_read > -20 .AND. rec_code_read <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
  endif
  incr=1
  IF ( dtgs%have_task_groups ) THEN
     !
     v_siz =  dtgs%tg_nnr * dtgs%nogrp
     ALLOCATE( tg_dv   ( v_siz, nspin_mag ) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = dtgs%nogrp
     !
  ENDIF
  !
  IF ( ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
  end if
  IF (rec_code_read > -20) convt=.TRUE.
  !
  if (convt) go to 155
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if (lgauss.or..not.lgamma) call errore ('solve_e', &
       'called in the wrong case', 1)

  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph

!     write(6,*) 'kter', kter
     FLUSH( stdout )
     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)
     IF (noncolin) dbecsum_nc=(0.d0,0.d0)

     do ik = 1, nksq
        npw = ngk(ik)
        npwq= npw     ! q=0 always in this routine
        if (lsda) current_spin = isk (ik)
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ik)
        !
        ! compute beta functions and kinetic energy for k-point ik
        ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
        !
        CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
        CALL g2_kin(ik)
        !
        ! compute preconditioning matrix h_diag used by cgsolve_all
        !
        CALL h_prec (ik, evc, h_diag)
        !
        do ipol = 1, 3
           !
           ! computes/reads P_c^+ x psi_kpoint into dvpsi array
           !
           call dvpsi_e (ik, ipol)
           !
           if (iter > 1) then
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
              IF( dtgs%have_task_groups ) THEN
                 IF (noncolin) THEN
                    CALL tg_cgather( dffts, dtgs, dvscfins(:,1,ipol), &
                                                                tg_dv(:,1))
                    IF (domag) THEN
                       DO jpol=2,4
                          CALL tg_cgather( dffts, dtgs, dvscfins(:,jpol,ipol), &
                                                             tg_dv(:,jpol))
                       ENDDO
                    ENDIF
                 ELSE
                    CALL tg_cgather( dffts, dtgs, dvscfins(:,current_spin,ipol), &
                                                             tg_dv(:,1))
                 ENDIF
              ENDIF
              aux2=(0.0_DP,0.0_DP)
              do ibnd = 1, nbnd_occ (ik), incr
                 IF ( dtgs%have_task_groups ) THEN
                    call cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, &
                                      nbnd_occ (ik) )
                    call apply_dpot(v_siz, tg_psic, tg_dv, 1)
                    call cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, &
                                      nbnd_occ (ik))
                 ELSE
                    call cft_wave (ik, evc (1, ibnd), aux1, +1)
                    call apply_dpot(dffts%nnr, aux1, dvscfins(1,1,ipol), current_spin)
                    call cft_wave (ik, aux2 (1, ibnd), aux1, -1)
                 ENDIF
              enddo
              dvpsi=dvpsi+aux2
              !
              call adddvscf(ipol,ik)
              !
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
           !
           CALL orthogonalize(dvpsi, evc, ik, ik, dpsi, npwq, .false.)
           !
           if (iter == 1) then
              !
              !  At the first iteration dpsi and dvscfin are set to zero,
              !
              dpsi(:,:)=(0.d0,0.d0)
              dvscfin(:,:,:)=(0.d0,0.d0)
              !
              ! starting threshold for the iterative solution of the linear
              ! system
              !
              thresh = 1.d-2
              if (lnoloc) thresh = 1.d-5
           else
              ! starting value for  delta_psi is read from iudwf
              !
              nrec = (ipol - 1) * nksq + ik
              call get_buffer (dpsi, lrdwf, iudwf, nrec)
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (0.1d0 * sqrt (dr2), 1.0d-2)
           endif
           !
           ! iterative solution of the linear system (H-e)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
           !

           conv_root = .true.

           call cgsolve_all (ch_psi_all,cg_psi,et(1,ik),dvpsi,dpsi, &
              h_diag,npwx,npw,thresh,ik,lter,conv_root,anorm,nbnd_occ(ik),npol)

           ltaver = ltaver + lter
           lintercall = lintercall + 1
           if (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                &         ' solve_e: root not converged ',es10.3)") ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           nrec = (ipol - 1) * nksq + ik
           call save_buffer(dpsi, lrdwf, iudwf, nrec)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           IF (noncolin) THEN
              call incdrhoscf_nc(dvscfout(1,1,ipol),wk(ik),ik, &
                                 dbecsum_nc(1,1,1,1,ipol), dpsi)
           ELSE
              call incdrhoscf (dvscfout(1,current_spin,ipol), wk(ik), &
                            ik, dbecsum(1,1,current_spin,ipol), dpsi)
           ENDIF
        enddo   ! on polarizations
     enddo      ! on k points
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions
     !  coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        call mp_sum ( dbecsum, intra_bgrp_comm )
     END IF

     if (doublegrid) then
        do is=1,nspin_mag
           do ipol=1,3
              call cinterpolate (dvscfout(1,is,ipol), dvscfout(1,is,ipol), 1)
           enddo
        enddo
     endif
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 3)
     !
     call addusddense (dvscfout, dbecsum)
     !
     !   dvscfout contains the (unsymmetrized) linear charge response
     !   for the three polarizations - symmetrize it
     !
     call mp_sum ( dvscfout, inter_pool_comm )
     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
     if (.not.lgamma_gamma) then
        call psyme (dvscfout)
        IF ( noncolin.and.domag ) CALL psym_dmage(dvscfout)
     endif
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     do ipol=1,3
        if (fildrho.ne.' ') call davcio_drho(dvscfout(1,1,ipol),lrdrho, &
             iudrho,ipol,+1)
        IF (lnoloc) then
           dvscfout(:,:,ipol)=(0.d0,0.d0)
        ELSE
           call dv_of_drho (dvscfout (1, 1, ipol), .false.)
        ENDIF
     enddo
     !
     !   mix the new potential with the old
     !
     IF (okpaw) THEN
     !
     !  In this case we mix also dbecsum
     !
        call setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixout, dvscfout, dbecsum, ndim, -1 )
        call mix_potential (2*3*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, 3*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
        call setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )
     ELSE
        call mix_potential (2*3*dfftp%nnr*nspin_mag, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2, 3 * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
     ENDIF
     if (doublegrid) then
        do is=1,nspin_mag
           do ipol = 1, 3
              call cinterpolate (dvscfin(1,is,ipol),dvscfins(1,is,ipol),-1)
           enddo
        enddo
     endif

     IF (okpaw) THEN
        IF (noncolin) THEN
!           call PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
        ELSE
!
!    The presence of c.c. in the formula gives a factor 2.0
!
           dbecsum=2.0_DP * dbecsum
           IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
           call PAW_dpotential(dbecsum,rho%bec,int3_paw,3)
        ENDIF
     ENDIF

     call newdq(dvscfin,3)

     averlt = DBLE (ltaver) / DBLE (lintercall)

     tcpu = get_clock ('PHONON')
     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / 3
     WRITE( stdout, "(5x,' thresh=',es10.3, ' alpha_mix = ',f6.3, &
          &      ' |ddv_scf|^2 = ',es10.3 )") thresh, alpha_mix (kter), dr2
     !
     FLUSH( stdout )
     !
     ! rec_code: state of the calculation
     ! rec_code=-20 Electric Field
     !
     rec_code=-20
     IF (okpaw) THEN
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 3, dvscfin, &
                                                       dvscfout, dbecsum)
     ELSE
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 3, dvscfin)
     ENDIF

     if (check_stop_now()) call stop_smoothly_ph (.false.)

     if (convt) goto 155

  enddo
155 continue
  deallocate (h_diag)
  deallocate (aux1)
  deallocate (dbecsum)
  deallocate (dvscfout)
  IF (okpaw) THEN
     DEALLOCATE(mixin)
     DEALLOCATE(mixout)
  ENDIF
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  if (noncolin) deallocate(dbecsum_nc)
  deallocate(aux2)
  IF ( dtgs%have_task_groups ) THEN
     !
     DEALLOCATE( tg_dv  )
     DEALLOCATE( tg_psic)
     !
  ENDIF

  call stop_clock ('solve_e')
  return
end subroutine solve_e
