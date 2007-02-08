!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : prefix, iunigk
  use pwcom
  USE check_stop,            ONLY : check_stop_now
  USE wavefunctions_module,  ONLY : evc
  USE kinds,                 ONLY : DP
  USE becmod,                ONLY : becp
  USE uspp_param,            ONLY : nhm
  use phcom
  USE control_flags,         ONLY : reduce_io
  
  implicit none

  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP), allocatable :: h_diag (:,:), eprec(:)
  ! h_diag: diagonal part of the Hamiltonian
  ! eprec : array fo preconditioning

  complex(DP) , allocatable, target ::      &
                   dvscfin (:,:,:)     ! change of the scf potential (input)
  complex(DP) , pointer ::      &
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  complex(DP) , allocatable ::   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   auxg (:), aux1 (:),  ps (:,:)

  complex(DP), EXTERNAL :: ZDOTC      ! the scalar product function

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: kter, iter0, ipol, ibnd, iter, lter, &
       ik, ig, irr, ir, is, nrec, ios
  ! counters
  integer :: ltaver, lintercall

  real(DP) :: tcpu, get_clock
  ! timing variables

  character (len=256) :: flmixdpot
  ! the name of the file with the mixing potential

  external ch_psi_all, cg_psi

  if (lsda) call errore ('solve_e', ' LSDA not implemented', 1)

  call start_clock ('solve_e')
  allocate (dvscfin( nrxx, nspin, 3))    
  if (doublegrid) then
     allocate (dvscfins(  nrxxs, nspin, 3))    
  else
     dvscfins => dvscfin
  endif
  allocate (dvscfout( nrxx , nspin, 3))    
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin, 3))    
  allocate (auxg(npwx))    
  allocate (aux1(nrxxs))    
  allocate (ps  (nbnd,nbnd))    
  ps (:,:) = (0.d0, 0.d0)
  allocate (h_diag(npwx, nbnd))    
  allocate (eprec(nbnd))
  if (irr0 == -20) then
     ! restarting in Electric field calculation
     read (iunrec) iter0, dr2
     read (iunrec) dvscfin
     if (okvan) read (iunrec) int1, int2, int3
     close (unit = iunrec, status = 'keep')
     if (doublegrid) then
        do is=1,nspin
           do ipol=1,3
              call cinterpolate (dvscfin(1,is,ipol), dvscfins(1,is,ipol), -1)
           enddo
        enddo
     endif
     convt = .false.
  else if (irr0 > -20 .AND. irr0 <= -10) then
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
  if (reduce_io) then
     flmixdpot = ' '
  else
     flmixdpot = 'mixd'
  endif
  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph

     iter = kter + iter0
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)

     if (nksq.gt.1) rewind (unit = iunigk)
     do ik = 1, nksq
        if (lsda) current_spin = isk (ik)
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
                    aux1 (nls(igk(ig)))=evc(ig,ibnd)
                 enddo
                 call cft3s (aux1,nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,+2)
                 do ir = 1, nrxxs
                    aux1(ir)=aux1(ir)*dvscfins(ir,current_spin,ipol)
                 enddo
                 call cft3s (aux1,nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
                 do ig = 1, npwq
                    dvpsi(ig,ibnd)=dvpsi(ig,ibnd)+aux1(nls(igkq(ig)))
                 enddo
              enddo
              !
              call adddvscf(ipol,ik)
              !
           endif
           !
           ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
           !
           CALL ZGEMM( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
                (1.d0,0.d0), evc(1,1), npwx, dvpsi(1,1), npwx, (0.d0,0.d0), &
                ps(1,1), nbnd )
#ifdef __PARA
           call reduce (2 * nbnd * nbnd_occ (ik), ps)
#endif
           ! dpsi is used as work space to store S|evc>
           !
           CALL ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, evc)
           CALL s_psi (npwx, npw, nbnd_occ(ik), evc, dpsi)
           !
           ! |dvpsi> = - (|dvpsi> - S|evc><evc|dvpsi>)
           ! note the change of sign!
           !
           CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
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
           do ibnd = 1, nbnd_occ (ik)
              do ig = 1, npw
                 auxg (ig) = g2kin (ig) * evc (ig, ibnd)
              enddo
              eprec (ibnd) = 1.35d0*ZDOTC(npwq,evc(1,ibnd),1,auxg,1)
           enddo
#ifdef __PARA
           call reduce (nbnd_occ (ik), eprec)
#endif
           do ibnd = 1, nbnd_occ (ik)
              do ig = 1, npw
                 h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd))
              enddo
           enddo

           conv_root = .true.

           call cgsolve_all (ch_psi_all,cg_psi,et(1,ik),dvpsi,dpsi, &
              h_diag,npwx,npw,thresh,ik,lter,conv_root,anorm,nbnd_occ(ik),1 )

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
           call incdrhoscf (dvscfout(1,current_spin,ipol), wk(ik), &
                            ik, dbecsum(1,1,current_spin,ipol))
        enddo   ! on polarizations
     enddo      ! on k points
#ifdef __PARA
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions 
     !  coming from each slice of bands
     !
     call reduce (nhm * (nhm + 1) * nat * nspin * 3, dbecsum)
#endif

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
#ifdef __PARA
     call poolreduce (2 * 3 * nrxx *nspin, dvscfout)
#endif
     if (.not.lgamma_gamma) then
#ifdef __PARA
        call psyme (dvscfout)
#else
        call syme (dvscfout)
#endif
     endif
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     do ipol=1,3
        if (fildrho.ne.' ') call davcio_drho(dvscfout(1,1,ipol),lrdrho, &
             iudrho,ipol,+1)
        call dv_of_drho (0, dvscfout (1, 1, ipol), .false.)
     enddo
     IF (lnoloc) dvscfout(:,:,:)=(0.d0,0.d0)
     !
     !   mix the new potential with the old 
     !
     call mix_potential (2 * 3 * nrxx *nspin, dvscfout, dvscfin, alpha_mix ( &
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
     WRITE( stdout, "(5x,' thresh=',e10.3, ' alpha_mix = ',f6.3, &
          &      ' |ddv_scf|^2 = ',e10.3 )") thresh, alpha_mix (kter), dr2
     !
     CALL flush_unit( stdout )
     !
     call seqopn (iunrec, 'recover', 'unformatted', exst)
     !
     ! irr: state of the calculation
     ! irr=-20 Electric Field
     !
     irr = -20
     !
     write (iunrec) irr
     !
     ! partially calculated results
     !
     write (iunrec) dyn, dyn00
     write (iunrec) epsilon, zstareu, zstarue, zstareu0, zstarue0
     !
     ! info on current iteration (iter=0 if potential mixing not available)
     !
     if (reduce_io .or. convt) then
        write (iunrec) 0, dr2
     else
        write (iunrec) iter, dr2
     end if
     write (iunrec) dvscfin
     if (okvan) write (iunrec) int1, int2, int3

     close (unit = iunrec, status = 'keep')
     if (check_stop_now()) then 
        call stop_ph (.false.)
        goto 155
     endif
     if (convt) goto 155

  enddo
155 continue
  deallocate (eprec)
  deallocate (h_diag)
  deallocate (ps)
  deallocate (aux1)
  deallocate (auxg)
  deallocate (dbecsum)
  deallocate (dvscfout)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)

  call stop_clock ('solve_e')
  return
end subroutine solve_e
