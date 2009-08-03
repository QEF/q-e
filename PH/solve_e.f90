!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
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
  USE ions_base,             ONLY : nat, ntyp => nsp, ityp
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : prefix, iunigk
  USE cell_base,             ONLY : tpiba2
  USE klist,                 ONLY : lgauss, xk, wk
  USE gvect,                 ONLY : nrxx, g
  USE gsmooth,               ONLY : doublegrid, nrxxs
  USE becmod,                ONLY : becp, becp_nc, calbec
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,              ONLY : domag
  USE wvfct,                 ONLY : nbnd, npw, npwx, igk, g2kin,  et
  USE check_stop,            ONLY : check_stop_now
  USE wavefunctions_module,  ONLY : evc
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : upf, nhm
  USE noncollin_module,      ONLY : noncolin, npol
  USE scf,                   ONLY : rho
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential, paw_desymmetrize

  USE eqv,                   ONLY : dpsi, dvpsi, eprec
  USE modes,                 ONLY : max_irr_dim
  USE units_ph,              ONLY : lrdwf, iudwf, lrwfc, iuwfc, lrdrho, &
                                    iudrho
  USE output,                ONLY : fildrho
  USE control_ph,            ONLY : recover, rec_code, &
                                    lnoloc, nbnd_occ, convt, tr2_ph, nmix_ph, &
                                    alpha_mix, lgamma_gamma, niter_ph, &
                                    lgamma, flmixdpot
  USE phus,                  ONLY : int3_paw
  USE qpoint,                ONLY : igkq, npwq, nksq
  USE recover_mod,           ONLY : read_rec, write_rec

  USE mp_global,             ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                    ONLY : mp_sum
  
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
                   aux1 (:,:),  ps (:,:)

  complex(DP), EXTERNAL :: zdotc      ! the scalar product function

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: kter, iter0, ipol, ibnd, jbnd, iter, lter, &
       ik, ig, irr, ir, is, nrec, na, nt, ios
  ! counters
  integer :: ltaver, lintercall

  real(DP) :: tcpu, get_clock
  ! timing variables

  external ch_psi_all, cg_psi

  call start_clock ('solve_e')
  allocate (dvscfin( nrxx, nspin, 3))    
  if (doublegrid) then
     allocate (dvscfins(  nrxxs, nspin, 3))    
  else
     dvscfins => dvscfin
  endif
  allocate (dvscfout( nrxx , nspin, 3))    
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin, 3))    
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin, 3))
  allocate (aux1(nrxxs,npol))    
  allocate (h_diag(npwx*npol, nbnd))    
  if (rec_code == -20.and.recover) then
     ! restarting in Electric field calculation
     CALL read_rec(dr2, iter0, dvscfin, dvscfins, 3)
     convt=.false.
  else if (rec_code > -20 .AND. rec_code <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
  endif
  !
  IF ( ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL DIROPN (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
  end if
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
     CALL flush_unit( stdout )
     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)
     IF (noncolin) dbecsum_nc=(0.d0,0.d0)

     if (nksq.gt.1) rewind (unit = iunigk)
     do ik = 1, nksq
        if (lsda) current_spin = isk (ik)
!        write(6,*) 'current spin', current_spin, ik
        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('solve_e', 'reading igk', abs (ios) )
        endif
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
        npwq = npw
        call init_us_2 (npw, igk, xk (1, ik), vkb)
        !
        ! compute the kinetic energy
        !
        do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ik ) + g (1,igkq (ig)) ) **2 + &
                          (xk (2,ik ) + g (2,igkq (ig)) ) **2 + &
                          (xk (3,ik ) + g (3,igkq (ig)) ) **2 ) * tpiba2
        enddo
        h_diag=0.d0
        do ibnd = 1, nbnd_occ (ik)
           do ig = 1, npw
              h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
           enddo
           IF (noncolin) THEN
              do ig = 1, npw
                 h_diag(ig+npwx,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
              enddo
           END IF
        enddo
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
              do ibnd = 1, nbnd_occ (ik)
                 call cft_wave (evc (1, ibnd), aux1, +1)
                 call apply_dpot(aux1, dvscfins(1,1,ipol), current_spin)
                 call cft_wave (dvpsi (1, ibnd), aux1, -1)
              enddo
              !
              call adddvscf(ipol,ik)
              !
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
           !
           CALL orthogonalize(dvpsi, evc, ik, ik, dpsi)
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
              call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
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
                &         ' solve_e: root not converged ',e10.3)") ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           nrec = (ipol - 1) * nksq + ik
           call davcio (dpsi, lrdwf, iudwf, nrec, + 1)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           IF (noncolin) THEN
              call incdrhoscf_nc(dvscfout(1,1,ipol),wk(ik),ik, &
                                       dbecsum_nc(1,1,1,1,ipol))
           ELSE
              call incdrhoscf (dvscfout(1,current_spin,ipol), wk(ik), &
                            ik, dbecsum(1,1,current_spin,ipol))
           ENDIF
        enddo   ! on polarizations
     enddo      ! on k points
#ifdef __PARA
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions 
     !  coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_pool_comm )
     ELSE
        call mp_sum ( dbecsum, intra_pool_comm )
     END IF
#endif

     if (doublegrid) then
        do is=1,nspin
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
#ifdef __PARA
     call mp_sum ( dvscfout, inter_pool_comm )
#endif
     if (.not.lgamma_gamma) then
#ifdef __PARA
        call psyme (dvscfout)
        IF ( noncolin.and.domag ) CALL psym_dmage(dvscfout)
#else
        call syme (dvscfout)
        IF ( noncolin.and.domag ) CALL sym_dmage(dvscfout)
#endif
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
           call dv_of_drho (0, dvscfout (1, 1, ipol), .false.)
        ENDIF
     enddo
     !
     !   mix the new potential with the old 
     !
     call mix_potential (2 * 3 * nrxx *nspin, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2, 3 * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
     if (doublegrid) then
        do is=1,nspin
           do ipol = 1, 3
              call cinterpolate (dvscfin(1,is,ipol),dvscfins(1,is,ipol),-1)
           enddo
        enddo
     endif

     IF (okpaw) THEN
        IF (noncolin) THEN
!           call PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,max_irr_dim)
        ELSE
!
!    The presence of c.c. in the formula gives a factor 2.0
!
           dbecsum=2.0_DP * dbecsum
           IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
           call PAW_dpotential(dbecsum,rho%bec,int3_paw,3,max_irr_dim)
        ENDIF
     ENDIF

     call newdq(dvscfin,3)

     averlt = DBLE (ltaver) / DBLE (lintercall)
  
     tcpu = get_clock ('PHONON')
     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / 3
     WRITE( stdout, "(5x,' thresh=',e10.3, ' alpha_mix = ',f6.3, &
          &      ' |ddv_scf|^2 = ',e10.3 )") thresh, alpha_mix (kter), dr2
     !
     CALL flush_unit( stdout )
     !
     ! rec_code: state of the calculation
     ! rec_code=-20 Electric Field
     !
     rec_code=-20
     CALL write_rec('solve_e...', irr, dr2, iter, convt, dvscfin, 3)

     if (check_stop_now().and..not.convt) call stop_ph (.false.)

     if (convt) goto 155

  enddo
155 continue
  deallocate (h_diag)
  deallocate (aux1)
  deallocate (dbecsum)
  deallocate (dvscfout)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  if (noncolin) deallocate(dbecsum_nc)


  call stop_clock ('solve_e')
  return
end subroutine solve_e
