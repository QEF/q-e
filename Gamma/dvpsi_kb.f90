!
! Copyright (C) 2003-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvpsi_kb(kpoint,nu)
  !----------------------------------------------------------------------
  ! calculates dVion/dtau * psi and stores it in dvpsi
  !
#include "f_defs.h"
  USE kinds,      ONLY: DP
  USE constants,  ONLY: tpi
  USE atom,       ONLY: rgrid
  USE becmod,     ONLY: calbec
  USE cell_base,  ONLY: omega, tpiba, tpiba2
  USE ions_base,  ONLY: ntyp => nsp, nat, ityp, tau
  USE uspp_param, ONLY: upf, nh, nhm
  USE uspp,       ONLY: dvan, nkb, vkb
  USE gvect,      ONLY : gstart, nl, nlm, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, ngm, g, gg, igtongl
  USE vlocal,     ONLY: vloc
  USE wvfct,      ONLY: nbnd, npwx, npw, g2kin, igk
  USE wavefunctions_module,  ONLY: evc, psic
  use cgcom
  !
  implicit none
  integer :: ibnd, ir, ih, jkb, ik, na, nu, ng, mu, nt, kpoint
  complex(DP), pointer:: work(:,:), workcc(:), dvloc(:), dvb_cc(:)
  complex(DP) :: exc
  real(DP), pointer :: bec1(:,:), bec2(:,:), dv(:)
  real(DP) :: gu, gtau
  logical :: has_nlcc
 !
  call start_clock('dvpsi_kb')
  !
  has_nlcc=.false.
  workcc => psic
  dv     => auxr
  dvloc  => aux2
  dvb_cc => aux3
  dvloc(:) = (0.d0, 0.d0)
  dvb_cc(:)= (0.d0, 0.d0)
  do na = 1,nat
     mu = 3*(na-1)
     if ( u(mu+1,nu)**2+u(mu+2,nu)**2+u(mu+3,nu)**2.gt. 1.0d-12) then
        nt=ityp(na)
        if (upf(nt)%nlcc) call drhoc (ngm, gg, omega, tpiba2, rgrid(nt)%mesh,&
                                  rgrid(nt)%dx, rgrid(nt)%r, upf(nt)%rho_atc,&
                                  workcc)
        has_nlcc = has_nlcc .or. upf(nt)%nlcc
        do ng = 1,ngm
           gtau = tpi * ( g(1,ng)*tau(1,na) + &
                          g(2,ng)*tau(2,na) + &
                          g(3,ng)*tau(3,na)   )
           gu = tpiba*( g(1,ng)*u(mu+1,nu) + &
                        g(2,ng)*u(mu+2,nu) + &
                        g(3,ng)*u(mu+3,nu)   )
           exc = gu * CMPLX(-sin(gtau),-cos(gtau))
           dvloc (nl(ng))=dvloc (nl(ng)) + vloc(igtongl(ng),nt)*exc
           if (upf(nt)%nlcc) dvb_cc(nl(ng)) = dvb_cc(nl(ng)) + workcc(ng) * exc
        end do
     end if
  end do
  do ng = gstart,ngm
     dvloc (nlm(ng))=CONJG(dvloc(nl(ng)))
  end do
  !
  !   dVloc/dtau in real space
  !
  call cft3(dvloc, nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  do ir = 1,nrxx
     dv(ir) =  DBLE(dvloc(ir))
  end do
  if (has_nlcc) then
     do ng = gstart,ngm
        dvb_cc (nlm(ng))=CONJG(dvb_cc(nl(ng)))
     end do
     call cft3(dvb_cc,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     do ir = 1,nrxx
        dv(ir) = dv(ir) +  DBLE(dvb_cc(ir)) * dmuxc(ir)
     end do
  end if
  !
  !   vloc_psi calculates dVloc/dtau*psi(G)
  !
  dvpsi(:,:) = (0.d0, 0.d0)
  call vloc_psi(npwx, npw, nbnd, evc, dv, dvpsi)
  !
  !   nonlocal (Kleinman-Bylander) contribution.
  !
  jkb=0
  do nt = 1,ntyp
     ! beware allocations !
     allocate (work( npwx, nh(nt)))    
     allocate (bec1( nh(nt), nbnd))    
     allocate (bec2( nh(nt), nbnd))    
     do na = 1,nat
        if (ityp(na) == nt .and. nh(nt) > 0) then
           mu =3*(na-1)
           if ( u(mu+1,nu)**2+u(mu+2,nu)**2+u(mu+3,nu)**2 > 1.0d-12) then
              !
              !  first term: sum_l sum_G' [ i V_l(G) V^*_l(G') (G'*u) psi(G')
              !  second term: sum_l sum_G' [-i (G*u) V_l(G) V^*_l(G') psi(G')
              !
              do ih = 1,nh(nt)
                 do ik = 1,npw
                    work(ik,ih) = vkb(ik,jkb+ih) * CMPLX(0.d0,-1.d0) * &
                                    (tpiba*( g(1,igk(ik))*u(mu+1,nu) +  &
                                             g(2,igk(ik))*u(mu+2,nu) +  &
                                             g(3,igk(ik))*u(mu+3,nu) ) )
                 end do
              end do
              !
              call calbec ( npw, work, evc, bec1 )
              call calbec ( npw, vkb(:,jkb+1:jkb+nh(nt)), evc, bec2 )
              !
              do ibnd = 1,nbnd
                 do ih = 1,nh(nt)
                    bec1(ih,ibnd) = dvan(ih,ih,nt) * bec1(ih,ibnd)
                    bec2(ih,ibnd) = dvan(ih,ih,nt) * bec2(ih,ibnd)
                 end do
              end do
              !
              call DGEMM ('N', 'N', 2*npw, nbnd, nh(nt), 1.d0, vkb(1,jkb+1), &
                   2*npwx, bec1, max(nh(nt),1), 1.d0, dvpsi, 2*npwx)
              call DGEMM ('N', 'N', 2*npw, nbnd, nh(nt), 1.d0, work, &
                   2*npwx, bec2, max(nh(nt),1), 1.d0, dvpsi, 2*npwx)
           end if
           jkb = jkb + nh(nt)
        end if
     end do
     deallocate(work)
     deallocate(bec2)
     deallocate(bec1)
  end do
  if (jkb.ne.nkb) call errore('dvpsi_kb','unexpected error',1)
  !
  call stop_clock('dvpsi_kb')
  !
  return
end subroutine dvpsi_kb
