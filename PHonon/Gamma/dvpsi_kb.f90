!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE dvpsi_kb(ik,nu)
  !----------------------------------------------------------------------
  ! calculates dVion/dtau * psi and stores it in dvpsi
  !
  USE kinds,      ONLY: DP
  USE constants,  ONLY: tpi
  USE atom,       ONLY: rgrid
  USE becmod,     ONLY: calbec
  USE cell_base,  ONLY: omega, tpiba, tpiba2
  USE klist,      ONLY: ngk
  USE ions_base,  ONLY: ntyp => nsp, nat, ityp, tau
  USE uspp_param, ONLY: upf, nh, nhm
  USE uspp,       ONLY: dvan, nkb, vkb
  USE fft_base,   ONLY: dfftp
  USE fft_interfaces, ONLY : invfft
  USE gvect,      ONLY : gstart, nl, nlm, ngl, ngm, g, gg, gl, igtongl
  USE vlocal,     ONLY: vloc
  USE wvfct,      ONLY: nbnd, npwx
  USE wavefunctions_module,  ONLY: evc, psic
  USE cgcom
  !
  IMPLICIT NONE
  INTEGER :: ik, nu
  INTEGER :: npw, ibnd, ir, ih, jkb, ig, na, ng, mu, nt
  COMPLEX(DP), POINTER:: work(:,:), dvloc(:), dvb_cc(:)
  COMPLEX(DP) :: exc
  real(DP), POINTER :: bec1(:,:), bec2(:,:), rhocg(:), dv(:)
  real(DP) :: gu, gtau
  LOGICAL :: has_nlcc
 !
  CALL start_clock('dvpsi_kb')
  !
  has_nlcc=.false.
  rhocg  => auxr
  dv     => auxr
  dvloc  => aux2
  dvb_cc => aux3
  dvloc(:) = (0.d0, 0.d0)
  dvb_cc(:)= (0.d0, 0.d0)
  DO na = 1,nat
     mu = 3*(na-1)
     IF ( u(mu+1,nu)**2+u(mu+2,nu)**2+u(mu+3,nu)**2> 1.0d-12) THEN
        nt=ityp(na)
        IF (upf(nt)%nlcc) CALL drhoc (ngl, gl, omega, tpiba2, rgrid(nt)%mesh,&
                                  rgrid(nt)%r, rgrid(nt)%rab, upf(nt)%rho_atc,&
                                  rhocg )
        has_nlcc = has_nlcc .or. upf(nt)%nlcc
        DO ng = 1,ngm
           gtau = tpi * ( g(1,ng)*tau(1,na) + &
                          g(2,ng)*tau(2,na) + &
                          g(3,ng)*tau(3,na)   )
           gu = tpiba*( g(1,ng)*u(mu+1,nu) + &
                        g(2,ng)*u(mu+2,nu) + &
                        g(3,ng)*u(mu+3,nu)   )
           exc = gu * cmplx(-sin(gtau),-cos(gtau),kind=DP)
           dvloc (nl(ng))=dvloc (nl(ng)) + vloc(igtongl(ng),nt)*exc
           IF (upf(nt)%nlcc) &
              dvb_cc(nl(ng)) = dvb_cc(nl(ng)) + rhocg (igtongl(ng)) * exc
        ENDDO
     ENDIF
  ENDDO
  DO ng = gstart,ngm
     dvloc (nlm(ng))=conjg(dvloc(nl(ng)))
  ENDDO
  !
  !   dVloc/dtau in real space
  !
  CALL invfft ('Dense', dvloc, dfftp)
  DO ir = 1,dfftp%nnr
     dv(ir) =  dble(dvloc(ir))
  ENDDO
  IF (has_nlcc) THEN
     DO ng = gstart,ngm
        dvb_cc (nlm(ng))=conjg(dvb_cc(nl(ng)))
     ENDDO
     CALL invfft ('Dense', dvb_cc, dfftp)
     DO ir = 1,dfftp%nnr
        dv(ir) = dv(ir) +  dble(dvb_cc(ir)) * dmuxc(ir)
     ENDDO
  ENDIF
  !
  !   vloc_psi calculates dVloc/dtau*psi(G)
  !
  npw = ngk(ik)
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  CALL vloc_psi_gamma(npwx, npw, nbnd, evc, dv, dvpsi)
  !
  !   nonlocal (Kleinman-Bylander) contribution.
  !
  jkb=0
  DO nt = 1,ntyp
     ! beware allocations !
     ALLOCATE (work( npwx, nh(nt)))
     ALLOCATE (bec1( nh(nt), nbnd))
     ALLOCATE (bec2( nh(nt), nbnd))
     DO na = 1,nat
        IF (ityp(na) == nt .and. nh(nt) > 0) THEN
           mu =3*(na-1)
           IF ( u(mu+1,nu)**2+u(mu+2,nu)**2+u(mu+3,nu)**2 > 1.0d-12) THEN
              !
              !  first term: sum_l sum_G' [ i V_l(G) V^*_l(G') (G'*u) psi(G')
              !  second term: sum_l sum_G' [-i (G*u) V_l(G) V^*_l(G') psi(G')
              !
              DO ih = 1,nh(nt)
                 DO ig = 1,npw
                    work(ig,ih) = vkb(ig,jkb+ih) * cmplx(0.d0,-1.d0,kind=DP) * &
                                    (tpiba*( g(1,ig)*u(mu+1,nu) +  &
                                             g(2,ig)*u(mu+2,nu) +  &
                                             g(3,ig)*u(mu+3,nu) ) )
                 ENDDO
              ENDDO
              !
              CALL calbec ( npw, work, evc, bec1 )
              CALL calbec ( npw, vkb(:,jkb+1:jkb+nh(nt)), evc, bec2 )
              !
              DO ibnd = 1,nbnd
                 DO ih = 1,nh(nt)
                    bec1(ih,ibnd) = dvan(ih,ih,nt) * bec1(ih,ibnd)
                    bec2(ih,ibnd) = dvan(ih,ih,nt) * bec2(ih,ibnd)
                 ENDDO
              ENDDO
              !
              CALL dgemm ('N', 'N', 2*npw, nbnd, nh(nt), 1.d0, vkb(1,jkb+1), &
                   2*npwx, bec1, max(nh(nt),1), 1.d0, dvpsi, 2*npwx)
              CALL dgemm ('N', 'N', 2*npw, nbnd, nh(nt), 1.d0, work, &
                   2*npwx, bec2, max(nh(nt),1), 1.d0, dvpsi, 2*npwx)
           ENDIF
           jkb = jkb + nh(nt)
        ENDIF
     ENDDO
     DEALLOCATE(work)
     DEALLOCATE(bec2)
     DEALLOCATE(bec1)
  ENDDO
  IF (jkb/=nkb) CALL errore('dvpsi_kb','unexpected error',1)
  !
  CALL stop_clock('dvpsi_kb')
  !
  RETURN
END SUBROUTINE dvpsi_kb
