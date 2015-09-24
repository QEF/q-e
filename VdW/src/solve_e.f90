!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_e_vdw ( iu )
  !-----------------------------------------------------------------------
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field.
  !    It performs the following tasks:
  !     a) It computes the kinetic energy
  !     b) It adds the term Delta V_{SCF} | psi >
  !     c) It applies P_c^+ to the known part
  !     d) It calls linter to solve the linear system
  !     e) It computes Delta rho, Delta V_{SCF} and symmetrize them
  !
  !
  USE ions_base,             ONLY : nat
  USE cell_base,             ONLY : omega, alat
  USE io_global,             ONLY : stdout
  USE io_files,              ONLY : iunigk, prefix, iunwfc, nwordwfc
  USE fft_base,              ONLY : dffts, dfftp
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE pwcom
  USE scf,                   ONLY : vrs
  USE check_stop,            ONLY : check_stop_now
!  USE wavefunctions_module,  ONLY : evc
  USE kinds,                 ONLY : DP
  USE becmod,                ONLY : becp, calbec
  USE uspp,                  ONLY : okvan, nkb, vkb
  USE uspp_param,            ONLY : nhm
  USE phcom
  USE phus,                  ONLY : becp1
  USE eff_v,                 ONLY : nelecr, veff, et_c, dvext,  evc => evc_veff, &
                                    dpsi_eff
  USE control_vdw
  USE mp_global,            ONLY : intra_pool_comm, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  real(kind=DP) ::  thresh, weight, anorm, averlt, dr2, pola
!  real(kind=DP), allocatable :: h_diag (:,:)
  real(kind=DP), ALLOCATABLE ::  eprec1(:)
  COMPLEX(kind=DP), ALLOCATABLE :: h_diag (:,:)
  ! the diagonal part of the Hamiltonia
  ! the convergence threshold
  ! used for summation over k points
  ! the norm of the error
  ! average number of iterations
  ! cut-off for preconditioning
  ! convergence limit
  ! imag. freq. step

  COMPLEX(kind=DP) , POINTER ::      &
                   dvscfin (:,:,:),  & ! change of the scf potential (input)
                   dvscfins (:,:,:)    ! change of the scf potential (smooth)
  COMPLEX(kind=DP) , ALLOCATABLE ::   &
                   dvscfout (:,:,:), & ! change of the scf potential (output)
                   dbecsum(:,:,:,:), & ! the becsum with dpsi
                   auxg (:,:), aux1 (:), spsi(:), ps (:)

  COMPLEX(kind=DP) :: zdotc      ! the scalar product function

  real(kind=DP) :: iu    ! frequency

  LOGICAL :: conv_root, exst
  ! true if linter is converged
  ! used to open the recover file

  INTEGER :: kter, ipol, ibnd, jbnd, iter, iter0, lter, ltaver, lintercall, &
       ik, ig, irr, ir, is, nrec, nrec1, ios
  ! counter on iterations
  ! counter on perturbations
  ! counter on bands
  ! counter on bands
  ! counter on iterations
  ! counter on iterations of linter
  ! average counter
  ! average number of call to linter
  ! counter on k points
  ! counter on G vectors
  ! the irreducible representation
  ! counter on g vectors
  ! counter on mesh points
  ! the record number
  ! the record number for dpsi
  ! integer variable for I/O control
  ! counter on loop over imag. freq.

  real(kind=DP) :: tcpu, get_clock, &   ! timing variables
                   w1                   ! weight

  CHARACTER (len=256) :: flmixdpot1
  ! the name of the file with the mixing potential
  !
  EXTERNAL ch_psi_all_vdw, pbcg_psi, cg_psi

  IF (lsda) CALL errore ('solve_e', ' LSDA not implemented', 1)

!  CALL init_clocks( .FALSE. )
  CALL start_clock ('solve_e')
  tcpu = get_clock ('VdW')
  !
  ALLOCATE (dvscfin( dfftp%nnr, nspin, 3))
  IF (doublegrid) THEN
     ALLOCATE (dvscfins(  dffts%nnr, nspin, 3))
  ELSE
     dvscfins => dvscfin
  ENDIF
  ALLOCATE (dvscfout( dfftp%nnr , nspin, 3))
  ALLOCATE (dbecsum( nhm*(nhm+1)/2, nat, nspin, 3))
  ALLOCATE (auxg(npwx,1))
  ALLOCATE (aux1(dffts%nnr))
  ALLOCATE (spsi(npwx))
  ALLOCATE (ps  (nbnd))
  ALLOCATE (h_diag(npwx, nbnd))
  ALLOCATE (eprec1(nbnd))
!  if (iter0.ne.0) then
!     if (okvan) read(iunrec) int3
!     read (iunrec) dr2, dvscfin
!     close (unit = iunrec, status = 'keep')
!     if (doublegrid) then
!        do is=1,nspin
!           do ipol=1,3
!              call cinterpolate (dvscfin(1,is,ipol), dvscfins(1,is,ipol), -1)
!           enddo
!        enddo
!     endif
!   else
      iter0 = 0
!  endif
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  IF (degauss/=0.d0.or..not.lgamma) CALL errore ('solve_e', &
       'called in the wrong case', 1)
  !
  ! assign the values of vrs in TFvW
  !
  vrs = veff
  !
  !   The outside loop is over the iterations
  !
  IF (reduce_io) THEN
     flmixdpot1 = ' '
  ELSE
     flmixdpot1 = 'flmixdpot'
  ENDIF
  !
  dr2 = 1.d-6
  !
!  et = 0.d0
  et_c(:,:) = cmplx(et(:,:), iu,kind=dp)
  !
  DO kter = 1, niter_vdw
     iter = kter + iter0
     convt = .true.
     ltaver = 0
     lintercall = 0

     dvscfout(:,:,:)=(0.d0,0.d0)
     dbecsum(:,:,:,:)=(0.d0,0.d0)
!     if (nksq.gt.1) rewind (unit = iunigk)
     DO ik = 1, nksq
        nbnd_occ(ik) = 1
!        if (lsda) current_spin = isk (ik)
!        if (nksq.gt.1) then
!           read (iunigk, err = 100, iostat = ios) npw, igk
!100        call errore ('solve_e', 'reading igk', abs (ios) )
!        endif
        !
        ! generate G-vector index
        !
        CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
!        if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
!        call davcio(evc, nwordwfc, iunwfc, ik, -1)
!        evc = evc_veff
        !
        npwq = npw
        !
        ! we compute the becp terms which are used in the rest of
        !    the code
        !
        IF ( nkb > 0 ) THEN
           CALL init_us_2 (npw, igk, xk (1, ik), vkb)
           CALL calbec( npw, vkb, evc, becp1(ik) )
        ENDIF
        !
        ! compute the kinetic energy
        !
        DO ig = 1, npwq
           g2kin (ig) = ( (xk (1,ik ) + g (1,igkq (ig)) ) **2 + &
                          (xk (2,ik ) + g (2,igkq (ig)) ) **2 + &
                          (xk (3,ik ) + g (3,igkq (ig)) ) **2 ) * tpiba2
        ENDDO
        !
        DO ipol=1,3
           nrec = (ipol - 1) * nksq + ik
           !
           ! computes/reads P_c^+ x psi_kpoint into dvpsi array
           CALL dvpsi_e_vdw (ik, ipol)
           !
           IF (iter==1) THEN
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
           ELSE
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
              DO ibnd = 1, nbnd_occ (ik)
                 aux1(:) = (0.d0, 0.d0)
                 DO ig = 1, npw
                    aux1 (nls(igk(ig)))=evc(ig,ibnd)
                 ENDDO
                 CALL invfft ('Wave', aux1, dffts)
                 DO ir = 1, dffts%nnr
                    aux1(ir)=aux1(ir)*dvscfins(ir,current_spin,ipol)
                 ENDDO
                 CALL fwfft ('Wave', aux1, dffts)
                 DO ig = 1, npwq
                    dvpsi(ig,ibnd)=dvpsi(ig,ibnd)+aux1(nls(igkq(ig)))
                 ENDDO
              ENDDO
              CALL adddvscf(ipol,ik)
              !
              ! starting value for  delta_psi is read from iudwf
              !
              CALL zcopy (npwx, dpsi_eff (1, ipol, 1), 1, dpsi (1, 1), 1)
!              nrec1 = (ipol - 1) * nksq + ik
!              call davcio (dpsi, lrdwf, iudwf, nrec1, - 1)
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (0.1d0 * sqrt (dr2), 1.0d-2)
           ENDIF
           !
           ! Orthogonalize dvpsi
           !
           DO ibnd = 1, nbnd_occ (ik)
              auxg(:,1) = (0.d0, 0.d0)
              DO jbnd = 1, nbnd_occ (ik)
                 ps(jbnd)=-zdotc(npw,evc(1,jbnd),1,dvpsi(1,ibnd),1)
              ENDDO
#ifdef __MPI
              CALL mp_sum( ps, intra_pool_comm )
#endif
              DO jbnd = 1, nbnd_occ (ik)
                 CALL zaxpy (npw, ps (jbnd), evc (1, jbnd), 1, auxg, 1)
              ENDDO
              IF ( nkb > 0 ) CALL calbec (npw, vkb, auxg, becp, 1)
              CALL s_psi (npwx, npw, 1, auxg, spsi)
              CALL daxpy (2*npw, 1.0d0, spsi, 1, dvpsi (1, ibnd), 1)
           ENDDO
           !
           !    Here we change the sign of the known term
           !
           CALL dscal (2*npwx*nbnd, -1.d0, dvpsi, 1)
           !
           ! iterative solution of the linear system (H-e)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           DO ibnd = 1, nbnd_occ (ik)
              DO ig = 1, npw
                 auxg (ig,1) = g2kin (ig) * evc (ig, ibnd)
              ENDDO
              eprec1 (ibnd) = 1.35d0*zdotc(npwq,evc(1,ibnd),1,auxg,1)
           ENDDO
#ifdef __MPI
           CALL mp_sum( eprec1( 1 : nbnd_occ(ik) ), intra_pool_comm )
#endif
           DO ibnd = 1, nbnd_occ (ik)
              DO ig = 1, npw
!                  h_diag(ig,ibnd)=(1.d0, 0.d0)
                 h_diag(ig,ibnd)=(1.d0, 0.d0) / &
                    cmplx(max(1.0d0,g2kin(ig)/eprec1(ibnd))-et(ibnd,ik), -iu,kind=DP)
!                 h_diag(ig,ibnd)=1.d0 / max(1.0d0,g2kin(ig)/eprec1(ibnd))
              ENDDO
           ENDDO
           !
           conv_root = .true.
           !
           CALL gmressolve_all (ch_psi_all_vdw,pbcg_psi,et_c(1,ik),dvpsi,dpsi,  &
              h_diag,npwx,npw,thresh,ik,lter,conv_root,anorm,nbnd_occ(ik), 4 )
           !
           ltaver = ltaver + lter
           lintercall = lintercall + 1
           IF (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4,   &
                &         ' linter: root not converged ',e10.3)") ik &
                &, ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           CALL zcopy (npwx, dpsi (1, 1), 1, dpsi_eff (1, ipol, 1),  1 )
!           nrec1 = (ipol - 1) * nksq + ik
!           call davcio (dpsi, lrdwf, iudwf, nrec1, + 1)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight = wk (ik)
           weight = nelecr
           !
           CALL incdrhoscf_vdw (dvscfout(1,current_spin,ipol), weight, ik, 1)
           !
        ENDDO   ! on perturbation
     ENDDO      ! on k points
     !
     !  keep dpsi for each direction
     !
!  call delta_rho()
!  stop
#ifdef __MPI
     !
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     !
     CALL mp_sum( dbecsum, intra_pool_comm )
#endif

     IF (doublegrid) THEN
        DO is=1,nspin
           DO ipol=1,3
              CALL cinterpolate (dvscfout(1,is,ipol), dvscfout(1,is,ipol), 1)
           ENDDO
        ENDDO
     ENDIF

!     call addusddense (dvscfout, dbecsum)

     !
     !   After the loop over the perturbations we have the change of the pote
     !   for all the modes of this representation. We symmetrize this potenti
     !
#ifdef __MPI
     CALL mp_sum ( dvscfout, inter_pool_comm )
#endif
     fildrho = ' '
     DO ipol=1,3
        IF (fildrho/=' ') CALL davcio_drho(dvscfout(1,1,ipol),lrdrho, &
             iudrho,ipol,+1)
        CALL dv_of_drho_vdw (0, dvscfout (1, 1, ipol), .false.)
     ENDDO
#ifdef __MPI
     CALL psyme (dvscfout)
#else
     CALL syme (dvscfout)
#endif
     !
     !   And we mix with the old potential
     !
     CALL mix_potential (2 * 3 * dfftp%nnr *nspin, dvscfout, dvscfin, al_mix_vdw ( &
          kter), dr2, 3 * tr2_vdw, iter, nmix_vdw, flmixdpot1, convt)
     IF (doublegrid) THEN
        DO is=1,nspin
           DO ipol=1,3
              CALL cinterpolate (dvscfin(1,is,ipol),dvscfins(1,is,ipol),-1)
           ENDDO
        ENDDO
     ENDIF

     CALL newdq(dvscfin,3)

     averlt = dble (ltaver) / dble (lintercall)

     tcpu = get_clock ('VdW')
     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time : ",f7.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / 3
     WRITE( stdout, "(5x,' thresh=',e10.3, ' al_mix_vdw = ',f6.3, &
          &      ' |ddv_scf|^2 = ',e10.3 )") thresh, al_mix_vdw (kter), dr2
#ifdef FLUSH
     FLUSH (6)
#endif

!     call seqopn (iunrec, 'recover', 'unformatted', exst)
!     irr = - 2

!     if (okvan) write (iunrec) int1, int2
!     write (iunrec) dyn, dyn00, epsilon, zstareu, zstarue, zstareu0, &
!          zstarue0
!     if (reduce_io) then
!        write(iunrec) irr, 0, convt, done_irr, comp_irr, ifat
!     else
!        write(iunrec) irr, iter, convt, done_irr, comp_irr, ifat
!        if (okvan) write(iunrec) int3
!        write(iunrec) dr2, dvscfin
!     endif

!     close (unit = iunrec, status = 'keep')
     tcpu = get_clock ('VdW')

     IF (check_stop_now()) THEN
        CALL stop_ph (.false.)
     ENDIF
     IF (convt) GOTO 155

  ENDDO ! of iteration
  !
155 CONTINUE
  !
  DEALLOCATE (eprec1)
  DEALLOCATE (h_diag)
  DEALLOCATE (ps)
  DEALLOCATE (spsi)
  DEALLOCATE (aux1)
  DEALLOCATE (auxg)
  DEALLOCATE (dbecsum)
  DEALLOCATE (dvscfout)
  IF (doublegrid) DEALLOCATE (dvscfins)
  DEALLOCATE (dvscfin)
  !
  CALL stop_clock ('solve_e')
  !
  RETURN
  !
END SUBROUTINE solve_e_vdw
