!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE d2ion (nat,ntyp,ityp,zv,tau,alat,omega,                &
     at,bg,g,gg,ngm,nmodes,u,has_equivalent,dyn)
  !-----------------------------------------------------------------------
  !
  !      calculate the ionic contribution to the dynamical matrix
  !      (both real- and reciprocal-space contributions are present)
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : tpi, fpi, e2
  USE io_global,  ONLY : stdout
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  IMPLICIT NONE
  INTEGER :: nat, ntyp, ngm, ityp(nat), nmodes, has_equivalent(nat)
  real(DP)::  tau(3,nat), g(3,ngm), gg(ngm), zv(ntyp), &
       u(3*nat,nmodes), dyn(3*nat,nmodes), at(3,3), bg(3,3), omega, alat
  !
  INTEGER :: nu_i,nu_j, mu_i,mu_j, na,nb, nta,ntb, ng, mxr, nrm, nr, i
  PARAMETER(mxr=50)
  real(DP) :: facg(nat), arg, tpiba2, alpha,  r(3,mxr), r2(mxr), dtau(3), &
              rmax, rr, upperbound, charge, gt2, fac, fnat, df, d2f, ar
  real(DP), EXTERNAL:: qe_erfc
  !
  !
  tpiba2 = (tpi/alat)**2
  !
  charge = 0.d0
  DO na=1, nat
     charge = charge + zv(ityp(na))
  ENDDO
  !
  alpha=0.5d0
  ! appropriate for c60
  WRITE( stdout,'("     d2ion: alpha = ",f6.2)') alpha
  !
  dyn (:,:) = 0.d0
  !
  ! G-space sum here
  !
  DO ng = 1, ngm
     !
     !  for parallel execution: first vector not necessarily G=0
     !
     IF(gg(ng)<1.e-6) GOTO 10
     !
     ! upperbound is a safe upper bound for the error ON THE ENERGY
     !
     upperbound=e2*charge**2*sqrt(2.0d0*alpha/tpi)*             &
          &              qe_erfc(sqrt(tpiba2*gg(ng)/4.d0/alpha))
     IF(upperbound<1.0d-6) GOTO 20
     !
     gt2 = gg(ng)*tpiba2
     fac = -e2*fpi*tpiba2/omega*exp(-gt2/alpha/4.d0)/gt2
     DO na = 1,nat
        nta= ityp(na)
        fnat = 0.0d0
        DO nb= 1,nat
           ntb= ityp(nb)
           arg = tpi*(g(1,ng)*(tau(1,na)-tau(1,nb))+  &
                      g(2,ng)*(tau(2,na)-tau(2,nb))+  &
                      g(3,ng)*(tau(3,na)-tau(3,nb)) )
           facg(nb) = fac*zv(nta)*zv(ntb)*cos(arg)
           fnat = fnat + facg(nb)
        ENDDO
        facg(na) = facg(na) - fnat
        mu_i = 3*(na-1)
        DO nu_i = 1,nmodes
           IF (has_equivalent( (nu_i-1)/3+1 )==1 ) GOTO 15
           arg =  g(1,ng)*u(mu_i+1,nu_i) + &
                  g(2,ng)*u(mu_i+2,nu_i) + &
                  g(3,ng)*u(mu_i+3,nu_i)
           IF (arg==0.0) GOTO 15
           DO nu_j = 1,nmodes
              DO nb= 1,nat
                 mu_j = 3*(nb-1)
                 dyn(nu_i,nu_j) = dyn(nu_i,nu_j) + facg(nb) * arg *  &
                      ( g(1,ng)*u(mu_j+1,nu_j) +   &
                        g(2,ng)*u(mu_j+2,nu_j) +   &
                        g(3,ng)*u(mu_j+3,nu_j)   )
              ENDDO
           ENDDO
15         CONTINUE
        ENDDO
     ENDDO
10   CONTINUE
  ENDDO
  PRINT '("   WARNING: G-sum not converged in d2ion ")'
  PRINT '("   d2ion  : alpha = ",f6.2)', alpha
  !
20 CONTINUE
  !
#define GAMMA
#ifdef GAMMA
  CALL dscal(3*nat*nmodes,2.d0,dyn,1)
#endif
  !
  !   for parallel execution: only node with G=0 calculates R-space term
  !
  IF(gg(1)>1.e-6) GOTO 30
  !
  ! R-space sum here
  !
  rmax=5.0d0/sqrt(alpha)/alat
  !
  ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-12)
  !
  DO na=1, nat
     nta= ityp(na)
     mu_i = 3*(na-1)
     DO nb=1, nat
        IF(nb/=na) THEN
           ntb= ityp(nb)
           mu_j = 3*(nb-1)
           DO i=1,3
              dtau(i)=tau(i,na)-tau(i,nb)
           ENDDO
           !
           ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
           !
           CALL rgen(dtau,rmax,mxr,at,bg,r,r2,nrm)
           DO nr=1, nrm
              rr=sqrt(r2(nr))*alat
              ar = sqrt(alpha)*rr
              d2f = ( 3.d0*qe_erfc(ar) + sqrt(8.d0/tpi)*ar*              &
                     (3.d0+2.d0*ar**2)*exp(-ar**2) ) / rr**5
              df  = (    -qe_erfc(ar) - sqrt(8.d0/tpi)*ar*exp(-ar**2) ) / rr**3
              DO nu_i = 1,nmodes
                 IF (has_equivalent( (nu_i-1)/3+1 )==1 ) GOTO 25
                 arg = r(1,nr)*u(mu_i+1,nu_i) + &
                       r(2,nr)*u(mu_i+2,nu_i) + &
                       r(3,nr)*u(mu_i+3,nu_i)
                 DO nu_j = 1,nmodes
                    dyn(nu_i,nu_j) = dyn(nu_i,nu_j) +            &
                         e2*zv(nta)*zv(ntb) * (d2f*alat *  arg * &
                         alat*( r(1,nr)*u(mu_j+1,nu_j) +         &
                                r(2,nr)*u(mu_j+2,nu_j) +         &
                                r(3,nr)*u(mu_j+3,nu_j)   ) +     &
                         df * ( u(mu_i+1,nu_i)*u(mu_j+1,nu_j) +  &
                                u(mu_i+2,nu_i)*u(mu_j+2,nu_j) +  &
                                u(mu_i+3,nu_i)*u(mu_j+3,nu_j) ) -&
                         d2f*alat * arg *  &
                         alat*( r(1,nr)*u(mu_i+1,nu_j) +         &
                                r(2,nr)*u(mu_i+2,nu_j) +         &
                                r(3,nr)*u(mu_i+3,nu_j)   ) -     &
                         df * ( u(mu_i+1,nu_i)*u(mu_i+1,nu_j) +  &
                                u(mu_i+2,nu_i)*u(mu_i+2,nu_j) +  &
                                u(mu_i+3,nu_i)*u(mu_i+3,nu_j) ) )
                 ENDDO
25               CONTINUE
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
30 CONTINUE
#if defined(__MPI)
  CALL mp_sum( dyn, intra_pool_comm )
#endif
  RETURN
END SUBROUTINE d2ion
