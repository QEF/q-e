!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE local_dos (iflag, lsign, kpoint, kband, spin_component, &
                      emin, emax, dos)
  !--------------------------------------------------------------------
  !
  !     iflag=0: calculates |psi|^2 for band "kband" at point "kpoint"
  !     iflag=1: calculates the local density of state at e_fermi
  !              (only for metals)
  !     iflag=2: calculates the local density of  electronic entropy
  !              (only for metals with fermi spreading)
  !     iflag=3: calculates the integral of local dos from "emin" to "emax"
  !              (emin, emax in Ry)
  !
  !     lsign:   if true and k=gamma and iflag=0, write |psi|^2 * sign(psi)
  !     spin_component: for iflag=3 and LSDA calculations only
  !                     0 for up+down dos,  1 for up dos, 2 for down dos
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE ener,                 ONLY : ef
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : nl, ngm, g
  USE gvecs,                ONLY : nls, nlsm, doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, wk, xk, &
                                   nkstot, ngk, igk_k
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : sym_rho, sym_rho_init, sym_rho_deallocate
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : lspinorb, fcoef
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE mp_global,            ONLY : me_pool, nproc_pool, my_pool_id, npool
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE becmod,               ONLY : calbec
  USE control_flags,        ONLY : tqr
  USE realus,               ONLY : addusdens_r
  IMPLICIT NONE
  !
  ! input variables
  !
  INTEGER, INTENT(in) :: iflag, kpoint, kband, spin_component
  LOGICAL, INTENT(in) :: lsign
  REAL(DP), INTENT(in) :: emin, emax
  !
  REAL(DP), INTENT(out) :: dos (dfftp%nnr)
  !
  !    local variables
  !
  ! counters for US PPs
  INTEGER :: npw, ikb, jkb, ijkb0, ih, jh, kh, na, ijh, np
  ! counters
  INTEGER :: ir, is, ig, ibnd, ik, irm, isup, isdw, ipol, kkb, is1, is2

  REAL(DP) :: w, w1, modulus, wg_max
  REAL(DP), ALLOCATABLE :: rbecp(:,:), segno(:), maxmod(:)
  COMPLEX(DP), ALLOCATABLE :: becp(:,:),  &
                                   becp_nc(:,:,:), be1(:,:), be2(:,:)
  INTEGER :: who_calculate, iproc
  COMPLEX(DP) :: phase
  REAL(DP), EXTERNAL :: w0gauss, w1gauss
  LOGICAL :: i_am_the_pool
  INTEGER :: which_pool, kpoint_pool
  !
  ! input checks
  !
  IF (noncolin.and. lsign) CALL errore('local_dos','not available',1)
  IF (noncolin.and. gamma_only) CALL errore('local_dos','not available',1)
  !
  IF ( (iflag == 0) .and. ( kband < 1 .or. kband > nbnd ) ) &
       CALL errore ('local_dos', 'wrong band specified', 1)
  IF ( (iflag == 0) .and. ( kpoint < 1 .or. kpoint > nkstot ) ) &
       CALL errore ('local_dos', 'wrong kpoint specified', 1)
  IF (lsign) THEN
     IF (iflag /= 0) CALL errore ('local_dos', 'inconsistent flags', 1)
     IF (sqrt(xk(1,kpoint)**2+xk(2,kpoint)**2+xk(3,kpoint)**2) > 1d-9 )  &
        CALL errore ('local_dos', 'k must be zero', 1)
  ENDIF
  !
  IF (gamma_only) THEN
     ALLOCATE (rbecp(nkb,nbnd))
  ELSE
     IF (noncolin) THEN
        ALLOCATE (becp_nc(nkb,npol,nbnd))
        IF (lspinorb) THEN
          ALLOCATE(be1(nhm,2))
          ALLOCATE(be2(nhm,2))
        ENDIF
     ELSE
        ALLOCATE (becp(nkb,nbnd))
     ENDIF
  ENDIF
  rho%of_r(:,:) = 0.d0
  dos(:) = 0.d0
  becsum(:,:,:) = 0.d0
  IF (lsign) ALLOCATE(segno(dfftp%nnr))
  !
  !   calculate the correct weights
  !
  IF (iflag /= 0.and. iflag /=3 .and. .not.lgauss) CALL errore ('local_dos', &
      'gaussian broadening needed', 1)
  IF (iflag == 2 .and. ngauss /= -99) CALL errore ('local_dos', &
      ' beware: not using Fermi-Dirac function ',  - ngauss)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        IF (iflag == 0) THEN
           wg (ibnd, ik) = 0.d0
        ELSEIF (iflag == 1) THEN
           !    Local density of states at energy emin with broadening emax
           wg(ibnd,ik) = wk(ik) * w0gauss((emin - et(ibnd, ik))/emax, ngauss) / emax
        ELSEIF (iflag == 2) THEN
           wg (ibnd, ik) = - wk (ik) * w1gauss ( (ef - et (ibnd, ik) ) &
                / degauss, ngauss)
        ELSEIF (iflag == 3) THEN
           IF (et (ibnd, ik) <=  emax .and. et (ibnd, ik) >= emin) THEN
              wg (ibnd, ik) = wk (ik)
           ELSE
              wg (ibnd, ik) = 0.d0
           ENDIF
        ELSE
           CALL errore ('local_dos', ' iflag not allowed', abs (iflag) )
        ENDIF
     ENDDO
  ENDDO
  wg_max = MAXVAL(wg(:,:))

  IF ( iflag == 0 .and. npool > 1 ) THEN
     CALL xk_pool( kpoint, nkstot, kpoint_pool,  which_pool )
     IF ( kpoint_pool < 1 .or. kpoint_pool > nks ) &
        CALL errore('local_dos','problems with xk_pool',1)
     i_am_the_pool=(my_pool_id==which_pool)
  ELSE
     i_am_the_pool=.true.
     kpoint_pool=kpoint
  ENDIF

  IF (iflag == 0.and.i_am_the_pool) wg (kband, kpoint_pool) = 1.d0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the density of states
  !
  DO ik = 1, nks
     IF (ik == kpoint_pool .and.i_am_the_pool.or. iflag /= 0) THEN
        IF (lsda) current_spin = isk (ik)
        CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
        npw = ngk(ik)
        CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)

        IF (gamma_only) THEN
           CALL calbec ( npw, vkb, evc, rbecp )
        ELSEIF (noncolin) THEN
           CALL calbec ( npw, vkb, evc, becp_nc )
        ELSE
           CALL calbec ( npw, vkb, evc, becp )
        ENDIF
     !
     !     here we compute the density of states
     !
        DO ibnd = 1, nbnd
         ! Neglect summands with relative weights below machine epsilon
         IF ( wg(ibnd, ik) > epsilon(0.0_DP) * wg_max .and. &
             (ibnd == kband .or. iflag /= 0)) THEN
              IF (noncolin) THEN
                 psic_nc = (0.d0,0.d0)
                 DO ig = 1, npw
                    psic_nc(nls(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
                    psic_nc(nls(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
                 ENDDO
                 DO ipol=1,npol
                    CALL invfft ('Wave', psic_nc(:,ipol), dffts)
                 ENDDO
              ELSE
                 psic(1:dffts%nnr) = (0.d0,0.d0)
                 DO ig = 1, npw
                    psic (nls (igk_k(ig,ik) ) ) = evc (ig, ibnd)
                 ENDDO
                 IF (gamma_only) THEN
                    DO ig = 1, npw
                       psic (nlsm(igk_k (ig,ik) ) ) = conjg(evc (ig, ibnd))
                    ENDDO
                 ENDIF
                 CALL invfft ('Wave', psic, dffts)
              ENDIF
              w1 = wg (ibnd, ik) / omega
!
!  Compute and save the sign of the wavefunction at the gamma point
!
              IF (lsign) THEN
                 IF (gamma_only) THEN
                    !  psi(r) is real by construction
                    segno(1:dffts%nnr) = dble(psic(1:dffts%nnr))
                 ELSE
                    !  determine the phase factor that makes psi(r) real.
                    ALLOCATE(maxmod(nproc_pool))
                    maxmod(me_pool+1)=0.0_DP
                    DO ir = 1, dffts%nnr
                       modulus=abs(psic(ir))
                       IF (modulus > maxmod(me_pool+1)) THEN
                          irm=ir
                          maxmod(me_pool+1)=modulus
                       ENDIF
                    ENDDO
                    who_calculate=1
#if defined(__MPI)
                    CALL mp_sum(maxmod,intra_pool_comm)
                    DO iproc=2,nproc_pool
                       IF (maxmod(iproc)>maxmod(who_calculate)) &
                          who_calculate=iproc
                    ENDDO
#endif
                    IF (maxmod(who_calculate) < 1.d-10) &
                       CALL errore('local_dos','zero wavefunction',1)
                    IF (me_pool+1==who_calculate) &
                          phase = psic(irm)/maxmod(who_calculate)
                    DEALLOCATE(maxmod)
#if defined(__MPI)
                    CALL mp_bcast(phase,who_calculate-1,intra_pool_comm)
#endif
                    segno(1:dffts%nnr) = dble( psic(1:dffts%nnr)*conjg(phase) )
                 ENDIF
                 IF (doublegrid) CALL interpolate (segno, segno, 1)
                 segno(:) = sign( 1.d0, segno(:) )
              ENDIF
              !
              IF (noncolin) THEN
                 DO ipol=1,npol
                    DO ir=1,dffts%nnr
                       rho%of_r(ir,current_spin)=rho%of_r(ir,current_spin)+&
                          w1*(dble(psic_nc(ir,ipol))**2+ &
                             aimag(psic_nc(ir,ipol))**2)
                    ENDDO
                 ENDDO
              ELSE
                 DO ir=1,dffts%nnr
                    rho%of_r(ir,current_spin)=rho%of_r(ir,current_spin) + &
                      w1 * (dble( psic (ir) ) **2 + aimag (psic (ir) ) **2)
                 ENDDO
              ENDIF
        !
        !    If we have a US pseudopotential we compute here the becsum term
        !
              w1 = wg (ibnd, ik)
              ijkb0 = 0
              DO np = 1, ntyp
                IF (upf(np)%tvanp  ) THEN
                  DO na = 1, nat
                    IF (ityp (na) == np) THEN
                      IF (noncolin) THEN
                        IF (upf(np)%has_so) THEN
                          be1=(0.d0,0.d0)
                          be2=(0.d0,0.d0)
                          DO ih = 1, nh(np)
                            ikb = ijkb0 + ih
                            DO kh = 1, nh(np)
                              IF ((nhtol(kh,np)==nhtol(ih,np)).and. &
                                  (nhtoj(kh,np)==nhtoj(ih,np)).and. &
                                  (indv(kh,np)==indv(ih,np))) THEN
                                 kkb=ijkb0 + kh
                                 DO is1=1,2
                                   DO is2=1,2
                                     be1(ih,is1)=be1(ih,is1)+ &
                                           fcoef(ih,kh,is1,is2,np)* &
                                           becp_nc(kkb,is2,ibnd)
                                     be2(ih,is1)=be2(ih,is1)+ &
                                           fcoef(kh,ih,is2,is1,np)* &
                                        conjg(becp_nc(kkb,is2,ibnd))
                                   ENDDO
                                 ENDDO
                              ENDIF
                            ENDDO
                          ENDDO
                        ENDIF
                        ijh = 1
                        DO ih = 1, nh (np)
                          ikb = ijkb0 + ih
                          IF (upf(np)%has_so) THEN
                            becsum(ijh,na,1)=becsum(ijh,na,1)+ w1*    &
                               (be1(ih,1)*be2(ih,1)+be1(ih,2)*be2(ih,2))
                          ELSE
                            becsum(ijh,na,1) = becsum(ijh,na,1)+  &
                             w1*(conjg(becp_nc(ikb,1,ibnd))*      &
                                       becp_nc(ikb,1,ibnd)+       &
                                 conjg(becp_nc(ikb,2,ibnd))*      &
                                       becp_nc(ikb,2,ibnd))
                          ENDIF
                          ijh = ijh + 1
                          DO jh = ih + 1, nh (np)
                            jkb = ijkb0 + jh
                            IF (upf(np)%has_so) THEN
                              becsum(ijh,na,1)=becsum(ijh,na,1) &
                                 + w1*((be1(jh,1)*be2(ih,1)+   &
                                        be1(jh,2)*be2(ih,2))+  &
                                       (be1(ih,1)*be2(jh,1)+   &
                                        be1(ih,2)*be2(jh,2)) )
                            ELSE
                              becsum(ijh,na,1)= becsum(ijh,na,1)+ &
                                   w1*2.d0*dble(conjg(becp_nc(ikb,1,ibnd)) &
                                     *becp_nc(jkb,1,ibnd) + &
                                conjg(becp_nc(ikb,2,ibnd)) &
                                     *becp_nc(jkb,2,ibnd) )
                            ENDIF
                            ijh = ijh + 1
                          ENDDO
                        ENDDO
                      ELSE
                        ijh = 1
                        DO ih = 1, nh (np)
                          ikb = ijkb0 + ih
                          IF (gamma_only) THEN
                              becsum(ijh,na,current_spin) = &
                                    becsum(ijh,na,current_spin) + w1 * &
                                    rbecp(ikb,ibnd)*rbecp(ikb,ibnd)
                          ELSE
                              becsum(ijh,na,current_spin) = &
                                   becsum(ijh,na,current_spin) + w1 * &
                               dble(conjg(becp(ikb,ibnd))*becp(ikb,ibnd))
                          ENDIF
                          ijh = ijh + 1
                          DO jh = ih + 1, nh (np)
                             jkb = ijkb0 + jh
                             IF (gamma_only) THEN
                                becsum(ijh,na,current_spin) = &
                                   becsum(ijh,na,current_spin) + 2.d0*w1 * &
                                   rbecp(ikb,ibnd)*rbecp(jkb,ibnd)
                             ELSE
                                becsum(ijh,na,current_spin) = &
                                  becsum(ijh,na,current_spin) + 2.d0*w1 * &
                                  dble(conjg(becp(ikb,ibnd))*becp(jkb,ibnd))
                             ENDIF
                             ijh = ijh + 1
                          ENDDO
                        ENDDO
                      ENDIF
                      ijkb0 = ijkb0 + nh (np)
                    ENDIF
                  ENDDO
                ELSE
                  DO na = 1, nat
                    IF (ityp (na) == np) ijkb0 = ijkb0 + nh (np)
                  ENDDO
                ENDIF
              ENDDO
           ENDIF
        ENDDO ! loop over bands
    ENDIF 
  ENDDO ! loop over k-points

  IF (gamma_only) THEN
     DEALLOCATE(rbecp)
  ELSE
     IF (noncolin) THEN
        IF (lspinorb) THEN
           DEALLOCATE(be1)
           DEALLOCATE(be2)
        ENDIF
        DEALLOCATE(becp_nc)
     ELSE
        DEALLOCATE(becp)
     ENDIF
  ENDIF
  IF (doublegrid) THEN
     IF (noncolin) THEN
       CALL interpolate(rho%of_r, rho%of_r, 1)
     ELSE
       DO is = 1, nspin
         CALL interpolate(rho%of_r(1, is), rho%of_r(1, is), 1)
       ENDDO
     ENDIF
  ENDIF
  !
  !    Here we add the US contribution to the charge
  !
  IF ( tqr ) THEN
   CALL addusdens_r(rho%of_r(:,:),.false.)
  ELSE
  !
  CALL addusdens(rho%of_r(:,:))
  !
  ENDIF
  !
  IF (nspin == 1 .or. nspin==4) THEN
     is = 1
     dos(:) = rho%of_r (:, is)
  ELSE
     IF ( iflag==3 .and. (spin_component==1 .or. spin_component==2 ) ) THEN
        dos(:) = rho%of_r (:, spin_component)
     ELSE
        isup = 1
        isdw = 2
        dos(:) = rho%of_r (:, isup) + rho%of_r (:, isdw)
     ENDIF
  ENDIF
  IF (lsign) THEN
     dos(:) = dos(:) * segno(:)
     DEALLOCATE(segno)
  ENDIF
#if defined(__MPI)
  CALL mp_sum( dos, inter_pool_comm )
#endif

  IF (iflag == 0 .or. gamma_only) RETURN
  !
  !    symmetrization of the local dos
  !
  CALL sym_rho_init (gamma_only )
  !
  psic(:) = cmplx ( dos(:), 0.0_dp, kind=dp)
  CALL fwfft ('Dense', psic, dfftp)
  rho%of_g(:,1) = psic(nl(:))
  !
  CALL sym_rho (1, rho%of_g)
  !
  psic(:) = (0.0_dp, 0.0_dp)
  psic(nl(:)) = rho%of_g(:,1)
  CALL invfft ('Dense', psic, dfftp)
  dos(:) = dble(psic(:))
  !
  CALL sym_rho_deallocate()
  !
  RETURN

END SUBROUTINE local_dos

!------------------------------------------------------------------------
SUBROUTINE xk_pool( ik, nkstot, ik_pool,  which_pool )
!------------------------------------------------------------------------
!
!  This routine is a simplified version of set_kpoint_vars in
!  xml_io_files. It receives the index ik of a k_point in the complete
!  k point list and return the index within the pool ik_pool, and
!  the number of the pool that has that k point.
!
!
USE mp_global,  ONLY : npool, kunit
!
IMPLICIT NONE

INTEGER, INTENT(in)  :: ik, nkstot
INTEGER, INTENT(out) :: ik_pool, which_pool
!
INTEGER :: nkl, nkr, nkbl
!
!
IF (npool==1) THEN
   which_pool=1
   ik_pool=ik
   RETURN
ENDIF
!
! ... find out number of k points blocks
!
nkbl = nkstot / kunit
!
! ... k points per pool
!
nkl = kunit * ( nkbl / npool )
!
! ... find out the reminder
!
nkr = ( nkstot - nkl * npool ) / kunit
!
! ... calculate the pool and the index within the pool
!
IF (ik<=nkr*(nkl+1)) THEN
   which_pool=(ik-1)/(nkl+1)
   ik_pool=ik-which_pool*(nkl+1)
ELSE
   which_pool=nkr+(ik-nkr*(nkl+1)-1)/nkl
   ik_pool=ik-nkr*(nkl+1)-(which_pool-nkr)*nkl
ENDIF

RETURN
END SUBROUTINE xk_pool
