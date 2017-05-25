!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE rhod2vkb(dyn0)
  !----------------------------------------------------------------------
  !
  !  calculate the electronic term: <psi|V''|psi>  of the dynamical matrix
  !
  USE kinds, ONLY: dp
  USE constants, ONLY: tpi
  USE ions_base, ONLY : nat, tau, ityp, ntyp => nsp
  USE cell_base, ONLY : tpiba2, tpiba, omega
  USE lsda_mod,  ONLY : current_spin
  USE gvect,  ONLY : ngm, g, igtongl, nl
  USE gvecw,  ONLY: gcutw
  USE wvfct,  ONLY: nbnd, npwx
  USE klist,  ONLY : wk, ngk
  USE scf,    ONLY : rho
  USE vlocal, ONLY: vloc
  USE wavefunctions_module,  ONLY: evc, psic
  USE uspp,   ONLY: nkb, vkb, dvan
  USE uspp_param, ONLY: nh
  USE becmod, ONLY: calbec
  USE cgcom
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE fft_base, ONLY : dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  real(DP) :: dyn0(3*nat,3*nat)
  !
  INTEGER :: npw, i, ih, ibnd, na, nt, nu_i,nu_j,mu_i,mu_j, ir, ng, jkb, ik, &
       ipol, jpol, ijpol
  real(DP) :: weight, fac, gtau
  real(DP), ALLOCATABLE :: dynloc(:,:), dynkb(:,:)
  COMPLEX(DP), ALLOCATABLE :: dvkb(:,:)
  real (DP), ALLOCATABLE :: becp(:,:), becp1(:,:,:), becp2(:,:,:)
  !
  CALL start_clock('rhod2vkb')
  !
  ! contribution from local potential
  !
  ALLOCATE  ( dynloc( 3*nat, nmodes))
  dynloc (:,:) = 0.d0
  DO ir = 1,dfftp%nnr
     psic(ir) = rho%of_r(ir,current_spin)
  ENDDO
  CALL fwfft ('Dense', psic, dfftp)
  DO nu_i = 1,nmodes
     IF (has_equivalent( (nu_i-1)/3+1)==1 ) GOTO 10
     DO na = 1, nat
        mu_i = 3*(na-1)
        IF ( u(mu_i+1,nu_i) == 0.0d0  .and.   &
             u(mu_i+2,nu_i) == 0.0d0  .and.   &
             u(mu_i+3,nu_i) == 0.0d0        ) GOTO 11
        DO nu_j = 1,nmodes
           IF ( u(mu_i+1,nu_j) == 0.0d0  .and. &
                u(mu_i+2,nu_j) == 0.0d0  .and. &
                u(mu_i+3,nu_j) == 0.0d0        ) GOTO 12
           DO ng = 1,ngm
              gtau = tpi * ( g(1,ng)*tau(1,na) + &
                             g(2,ng)*tau(2,na) + &
                             g(3,ng)*tau(3,na)   )
              fac = omega * vloc(igtongl(ng),ityp(na)) * tpiba2 *   &
                        (  dble(psic(nl(ng)))*cos(gtau) -  &
                          aimag(psic(nl(ng)))*sin(gtau)   )
              dynloc(nu_i,nu_j) = dynloc(nu_i,nu_j) + fac *         &
                   ( g(1,ng) * u(mu_i+1,nu_i) +                     &
                     g(2,ng) * u(mu_i+2,nu_i) +                     &
                     g(3,ng) * u(mu_i+3,nu_i)   ) *                 &
                   ( g(1,ng) * u(mu_i+1,nu_j) +                     &
                     g(2,ng) * u(mu_i+2,nu_j) +                     &
                     g(3,ng) * u(mu_i+3,nu_j)   )
           ENDDO
12         CONTINUE
        ENDDO
11      CONTINUE
     ENDDO
10   CONTINUE
  ENDDO
#define GAMMA
#ifdef GAMMA
  dynloc(:,:) = 2.d0 * dynloc(:,:)
#endif
#if defined(__MPI)
  CALL mp_sum( dynloc, intra_pool_comm )
#endif
  !
  !   contribution from nonlocal (Kleinman-Bylander) potential
  !
  ALLOCATE  (dynkb(3*nat,3*nat))
  dynkb=0.d0
  ALLOCATE  ( dvkb( npwx, nkb))
  ALLOCATE  ( becp ( nkb, nbnd))
  ALLOCATE  ( becp1( nkb, nbnd, 3))
  ALLOCATE  ( becp2( nkb, nbnd, 6))
  !
  ik = 1
  npw = ngk(ik)
     ! the sum has four terms which can be reduced to two (note factor 2 in weight):
     !
     ! sum_G sum_G' sum_j sum_l [ psi_j*(G) V_na,l(G)(-iGu_ipol) V^*_na,l(G')( iG'u_jpol) psi_j(G')
     ! sum_G sum_G' sum_j sum_l [ psi_j*(G) V_na,l(G)  V^*_na,l(G') ( iG'u_ipol)( iG'u_jpol) psi_j(G')
     !
     weight = 2.0d0*wk(ik)
     !
     CALL calbec ( npw, vkb, evc, becp )
     !
     !  becp(j,n) = sum_G [ V_n*(G) psi_j(G) ]   n=(na,l)
     !
     DO ipol = 1, 3
        DO jkb = 1, nkb
           DO i = 1,npw
              dvkb(i,jkb) = vkb(i,jkb) * cmplx(0.d0,-tpiba,kind=DP) * g(ipol,i)
           ENDDO
        ENDDO
        !
        CALL calbec ( npw, dvkb, evc, becp1(:,:,ipol) )
        !
        !  becp1(j,n,ipol) = sum_G [ V_n*(G) (iG_ipol) psi_j(G) ]
        !
     ENDDO
     ijpol=0
     DO ipol = 1, 3
        DO jpol = ipol, 3
           DO jkb = 1, nkb
              DO i = 1,npw
                 dvkb(i,jkb) = vkb(i,jkb) * tpiba2 * g(ipol,i)* g(jpol,i)
              ENDDO
           ENDDO
           !
           ijpol=ijpol+1
           CALL calbec ( npw, dvkb, evc, becp2(:,:,ijpol) )
           !
           !  becp2(j,n,ijpol) = sum_G [ V_n*(G) (-iG_ipol) (iG_jpol)  psi_j(G) ]
           !
        ENDDO
     ENDDO
     !
     jkb = 0
     DO nt = 1, ntyp
        DO na =1, nat
           IF (ityp(na)==nt) THEN
              IF (has_equivalent(na)==1 ) GOTO 20
              ijpol = 0
              DO ipol = 1, 3
                 nu_i = 3*(na-1) + ipol
                 DO jpol = ipol, 3
                    nu_j = 3*(na-1) + jpol
                    ijpol = ijpol + 1
                    DO ibnd=1,nbnd
                       DO ih=1,nh(nt)
                          dynkb(nu_i,nu_j) = dynkb(nu_i,nu_j) + &
                               (-becp1(jkb+ih,ibnd,ipol)*becp1(jkb+ih,ibnd,jpol) &
                                +becp2(jkb+ih,ibnd,ijpol)*becp(jkb+ih,ibnd) )    &
                               * dvan(ih,ih,nt) * weight
                       ENDDO
                    ENDDO
                 ENDDO
                 DO jpol = 1, ipol-1
                    nu_j = 3*(na-1) + jpol
                    dynkb(nu_i,nu_j) = dynkb(nu_j,nu_i)
                 ENDDO
              ENDDO
20            CONTINUE
              jkb = jkb + nh(nt)
           ENDIF
        ENDDO
     ENDDO
  !
  DEALLOCATE ( becp2)
  DEALLOCATE ( becp1)
  DEALLOCATE ( becp )
  DEALLOCATE ( dvkb)
  !
  dyn0 (:,:) = 0.d0
  !
  DO nu_i = 1,nmodes
     IF (has_equivalent( (nu_i-1)/3+1)==0 ) THEN
        DO nu_j=1,nmodes
           DO mu_i=1,3*nat
              DO mu_j=1,3*nat
                 dyn0(nu_i,nu_j) = dyn0(nu_i,nu_j) + &
                      dynkb(mu_i,mu_j)*u(mu_i,nu_i)*u(mu_j,nu_j)
              ENDDO
           ENDDO
        ENDDO
        DO nu_j=1,nmodes
           dyn0(nu_i,nu_j) = dyn0(nu_i,nu_j) + dynloc(nu_i,nu_j)
        ENDDO
     ENDIF
  ENDDO
  DEALLOCATE(dynkb)
  DEALLOCATE(dynloc)
  !
  CALL stop_clock('rhod2vkb')
  !
  RETURN
END SUBROUTINE rhod2vkb
