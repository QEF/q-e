!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e_fpol ( iw )
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
  USE buffers,               ONLY : get_buffer, save_buffer
  USE check_stop,            ONLY : check_stop_now
  USE wavefunctions_module,  ONLY : evc
  USE cell_base,             ONLY : tpiba2
  USE klist,                 ONLY : ltetra, lgauss, nkstot, wk, xk, ngk, igk_k
  USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
  USE fft_base,              ONLY : dffts, dfftp
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE gvect,                 ONLY : g
  USE gvecs,                 ONLY : doublegrid, nls
  USE becmod,                ONLY : becp, calbec
  USE wvfct,                 ONLY : npwx, nbnd, g2kin, et
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : nhm
  USE control_ph,            ONLY : nmix_ph, tr2_ph, alpha_mix, convt, &
                                    niter_ph, &
                                    rec_code, flmixdpot
  USE output,                ONLY : fildrho
  USE qpoint,                ONLY : nksq
  USE units_ph,              ONLY : lrdwf, iudwf, lrwfc, iuwfc, iudrho, &
                                    lrdrho
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE mp,                    ONLY : mp_sum

  USE eqv,                   ONLY : dpsi, dvpsi
  USE control_lr,            ONLY : nbnd_occ, lgamma
  USE dv_of_drho_lr

  implicit none

  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error

  complex(kind=DP), allocatable :: etc(:,:), h_diag(:,:)
  ! the eigenvalues plus imaginary frequency
  ! the diagonal part of the Hamiltonian which becomes complex now


  complex(DP) , allocatable, target ::      &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  complex(DP) , pointer ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  complex(DP) , allocatable ::   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   auxg (:), aux1 (:),  ps (:,:)

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: npw, npwq
  integer :: kter, iter0, ipol, ibnd, jbnd, iter, lter, &
       ik, ig, irr, ir, is, nrec, ios
  ! counters
  integer :: ltaver, lintercall

  real(DP) :: tcpu
  real(DP) :: eprec1 ! 1.35<ek>, for preconditioning
  real(DP) :: iw     !frequency
  real(dp), external :: ddot, get_clock

  external cch_psi_all, ccg_psi

  if (lsda) call errore ('solve_e_fpol', ' LSDA not implemented', 1)

  call start_clock ('solve_e')
  allocate (dvscfin( dfftp%nnr, nspin, 3))
  if (doublegrid) then
     allocate (dvscfins( dffts%nnr, nspin, 3))
  else
     dvscfins => dvscfin
  endif
  allocate (dvscfout( dfftp%nnr, nspin, 3))
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin, 3))
  allocate (auxg(npwx))
  allocate (aux1(dffts%nnr))
  allocate (ps  (nbnd,nbnd))
  ps (:,:) = (0.d0, 0.d0)
  allocate (h_diag(npwx, nbnd))

  allocate (etc(nbnd, nkstot))
  etc(:,:) = CMPLX( et(:,:), iw ,kind=DP)

  ! restart NOT IMPLEMENTED

  if (rec_code == -20) then
     !read (iunrec) iter0, convt, dr2
     !read (iunrec) dvscfin
     !if (okvan) read (iunrec) int3
     !close (unit = iunrec, status = 'keep')
     !if (doublegrid) then
     !   do is=1,nspin
     !      do ipol=1,3
     !         call cinterpolate (dvscfin(1,is,ipol), dvscfins(1,is,ipol), -1)
     !      enddo
     !   enddo
     !endif
  else if (rec_code > -20 .AND. rec_code <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
  endif
  !
  IF (ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
  end if
  !
  if (convt) go to 155
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if ( (lgauss .or. ltetra) .or. .not.lgamma) call errore ('solve_e_fpol', &
       'called in the wrong case', 1)
  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph

     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)

     do ik = 1, nksq
        if (lsda) current_spin = isk (ik)
        npw = ngk(ik)
        npwq = npw    ! q=0 in ths routine
        !
        ! read unperturbed wavefunctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call get_buffer(evc, lrwfc, iuwfc, ik)
        !
        ! compute beta functions and kinetic energy for k-point ik
        ! needed by h_psi, called by cch_psi_all, called by gmressolve_all
        !
        CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
        CALL g2_kin(ik)
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
                 aux1(:) = (0.d0, 0.d0)
                 do ig = 1, npw
                    aux1 (nls(igk_k(ig,ik)))=evc(ig,ibnd)
                 enddo
                 CALL invfft ('Wave', aux1, dffts)
                 do ir = 1, dffts%nnr
                    aux1(ir)=aux1(ir)*dvscfins(ir,current_spin,ipol)
                 enddo
                 CALL fwfft ('Wave', aux1, dffts)
                 do ig = 1, npwq
                    dvpsi(ig,ibnd)=dvpsi(ig,ibnd)+aux1(nls(igk_k(ig,ik)))
                 enddo
              enddo
              !
              call adddvscf(ipol,ik)
              !
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
           !
           CALL zgemm( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
                (1.d0,0.d0), evc(1,1), npwx, dvpsi(1,1), npwx, (0.d0,0.d0), &
                ps(1,1), nbnd )

           call mp_sum ( ps( :, 1:nbnd_occ(ik) ), intra_bgrp_comm )
           ! dpsi is used as work space to store S|evc>
           !
           CALL calbec (npw, vkb, evc, becp, nbnd_occ(ik) )
           CALL s_psi (npwx, npw, nbnd_occ(ik), evc, dpsi)
           !
           ! |dvpsi> = - (|dvpsi> - S|evc><evc|dvpsi>)
           ! note the change of sign!
           !
           CALL zgemm( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
               (1.d0,0.d0), dpsi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
                dvpsi(1,1), npwx )
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
           else
              ! starting value for  delta_psi is read from iudwf
              !
              nrec = (ipol - 1) * nksq + ik
              call get_buffer(dpsi, lrdwf, iudwf, nrec)
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (0.1d0 * sqrt (dr2), 1.0d-2)
           endif
           !
           ! iterative solution of the linear system (H-e)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           do ibnd = 1, nbnd_occ (ik)
              !
              if ( (abs(iw).lt.0.05) .or. (abs(iw).gt.1.d0) ) then
                 !
                 DO ig = 1, npwq
                    auxg (ig) = g2kin (ig) * evc (ig, ibnd)
                 END DO
                 eprec1 = 1.35_dp*ddot(2*npwq,evc(1,ibnd),1,auxg,1)
                 !
                 do ig = 1, npw
                    h_diag(ig,ibnd)=CMPLX(1.d0, 0.d0,kind=DP) / &
                    CMPLX( max(1.0d0,g2kin(ig)/eprec1)-et(ibnd,ik),-iw ,kind=DP)
                 end do
              else
                 do ig = 1, npw
                    h_diag(ig,ibnd)=CMPLX(1.d0, 0.d0,kind=DP)
                 end do
              endif
              !
           enddo

           conv_root = .true.

           call gmressolve_all (cch_psi_all,ccg_psi,etc(1,ik),dvpsi,dpsi,  &
              h_diag,npwx,npw,thresh,ik,lter,conv_root,anorm,nbnd_occ(ik), 4 )

           ltaver = ltaver + lter
           lintercall = lintercall + 1
           if (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                &         ' solve_e: root not converged ',es10.3)") ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           nrec = (ipol - 1) * nksq + ik
           call save_buffer (dpsi, lrdwf, iudwf, nrec)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           call incdrhoscf (dvscfout(1,current_spin,ipol), wk(ik), &
                            ik, dbecsum(1,1,current_spin,ipol), dpsi)
        enddo   ! on polarizations
     enddo      ! on k points
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions
     !  coming from each slice of bands
     !
     call mp_sum ( dbecsum, intra_bgrp_comm )

     if (doublegrid) then
        do is=1,nspin
           do ipol=1,3
              call cinterpolate (dvscfout(1,is,ipol), dvscfout(1,is,ipol), 1)
           enddo
        enddo
     endif

     call addusddense (dvscfout, dbecsum)
     !
     !   dvscfout contains the (unsymmetrized) linear charge response
     !   for the three polarizations - symmetrize it
     !
     call mp_sum ( dvscfout, inter_pool_comm )
     call psyme (dvscfout)
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     do ipol=1,3
        if (fildrho.ne.' ') call davcio_drho(dvscfout(1,1,ipol),lrdrho, &
             iudrho,ipol,+1)
        call dv_of_drho (dvscfout (1, 1, ipol), .false.)
     enddo
     !
     !   mix the new potential with the old
     !
     call mix_potential (2 * 3 * dfftp%nnr *nspin, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2, 3 * tr2_ph, iter, nmix_ph, flmixdpot, convt)
     if (doublegrid) then
        do is=1,nspin
           do ipol = 1, 3
              call cinterpolate (dvscfin(1,is,ipol),dvscfins(1,is,ipol),-1)
           enddo
        enddo
     endif

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
     ! restart NOT IMPLEMENTED
     !
     !call seqopn (iunrec, 'recover', 'unformatted', exst)
     !
     ! irr: state of the calculation
     ! irr=-20 Electric Field
     !
     !irr = -20
     !
     !write (iunrec) irr
     !
     ! partially calculated results
     !
     !write (iunrec) dyn, dyn00
     !write (iunrec) epsilon, zstareu, zstarue, zstareu0, zstarue0
     !
     ! info on current iteration (iter=0 if potential mixing not available)
     !
     !if (reduce_io) then
     !   write (iunrec) 0, convt, dr2
     !else
     !   write (iunrec) iter, convt, dr2
     !end if
     !write (iunrec) dvscfin
     !if (okvan) write (iunrec) int3

     !close (unit = iunrec, status = 'keep')
     if (check_stop_now()) then
        call stop_smoothly_ph (.false.)
        goto 155
     endif
     if (convt) goto 155
  enddo
155 continue
  deallocate (h_diag)
  deallocate (ps)
  deallocate (aux1)
  deallocate (auxg)
  deallocate (dbecsum)
  deallocate (dvscfout)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  deallocate(etc)

  call stop_clock ('solve_e')
  return
end subroutine solve_e_fpol
