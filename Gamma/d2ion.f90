!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d2ion (nat,ntyp,ityp,zv,tau,alat,omega,                &
     at,bg,g,gg,ngm,nmodes,u,has_equivalent,dyn)
  !-----------------------------------------------------------------------
  !
  !      calculate the ionic contribution to the dynamical matrix
  !      (both real- and reciprocal-space contributions are present)
  !
#include "f_defs.h"
  USE kinds, only : DP
  USE constants, ONLY : tpi, fpi, e2
  USE io_global,  ONLY : stdout
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  implicit none
  integer :: nat, ntyp, ngm, ityp(nat), nmodes, has_equivalent(nat)
  real(DP)::  tau(3,nat), g(3,ngm), gg(ngm), zv(ntyp), &
       u(3*nat,nmodes), dyn(3*nat,nmodes), at(3,3), bg(3,3), omega, alat
  !
  integer :: nu_i,nu_j, mu_i,mu_j, na,nb, nta,ntb, ng, mxr, nrm, nr, i
  parameter(mxr=50)
  real(DP) :: facg(nat), arg, tpiba2, alpha,                       &
       r(3,mxr), r2(mxr), dtau(3), erfc, rmax, rr,                 &
       upperbound, charge, gt2, fac, fnat, df, d2f, ar
  external erfc, rgen
  !
  !
  tpiba2 = (tpi/alat)**2
  !
  charge = 0.d0
  do na=1, nat
     charge = charge + zv(ityp(na))
  end do
  !
  alpha=0.5d0
  ! appropriate for c60
  WRITE( stdout,'("     d2ion: alpha = ",f6.2)') alpha
  !
  dyn (:,:) = 0.d0
  !
  ! G-space sum here
  !
  do ng = 1, ngm
     !
     !  for parallel execution: first vector not necessarily G=0
     !
     if(gg(ng).lt.1.e-6) go to 10
     !
     ! upperbound is a safe upper bound for the error ON THE ENERGY
     !
     upperbound=e2*charge**2*sqrt(2.0d0*alpha/tpi)*             &
          &              erfc(sqrt(tpiba2*gg(ng)/4.d0/alpha))
     if(upperbound.lt.1.0d-6) go to 20
     !
     gt2 = gg(ng)*tpiba2
     fac = -e2*fpi*tpiba2/omega*exp(-gt2/alpha/4.d0)/gt2
     do na = 1,nat
        nta= ityp(na)
        fnat = 0.0d0
        do nb= 1,nat
           ntb= ityp(nb)
           arg = tpi*(g(1,ng)*(tau(1,na)-tau(1,nb))+  &
                      g(2,ng)*(tau(2,na)-tau(2,nb))+  &
                      g(3,ng)*(tau(3,na)-tau(3,nb)) )
           facg(nb) = fac*zv(nta)*zv(ntb)*cos(arg)
           fnat = fnat + facg(nb)
        end do
        facg(na) = facg(na) - fnat
        mu_i = 3*(na-1)
        do nu_i = 1,nmodes
           if (has_equivalent( (nu_i-1)/3+1 ).eq.1 ) go to 15
           arg =  g(1,ng)*u(mu_i+1,nu_i) + &
                  g(2,ng)*u(mu_i+2,nu_i) + &
                  g(3,ng)*u(mu_i+3,nu_i)
           if (arg.eq.0.0) go to 15
           do nu_j = 1,nmodes
              do nb= 1,nat
                 mu_j = 3*(nb-1)
                 dyn(nu_i,nu_j) = dyn(nu_i,nu_j) + facg(nb) * arg *  &
                      ( g(1,ng)*u(mu_j+1,nu_j) +   &
                        g(2,ng)*u(mu_j+2,nu_j) +   &
                        g(3,ng)*u(mu_j+3,nu_j)   )
              end do
           end do
15         continue
        end do
     end do
10   continue
  end do
  print '("   WARNING: G-sum not converged in d2ion ")'
  print '("   d2ion  : alpha = ",f6.2)', alpha
  !
20 continue
  !
#define GAMMA
#ifdef GAMMA
  call DSCAL(3*nat*nmodes,2.d0,dyn,1)
#endif
  !
  !   for parallel execution: only node with G=0 calculates R-space term
  !
  if(gg(1).gt.1.e-6) go to 30
  !
  ! R-space sum here
  !
  rmax=5.0d0/sqrt(alpha)/alat
  !
  ! with this choice terms up to ZiZj*erfc(5) are counted (erfc(5)=2x10^-12)
  !
  do na=1, nat
     nta= ityp(na)
     mu_i = 3*(na-1)
     do nb=1, nat
        if(nb.ne.na) then
           ntb= ityp(nb)
           mu_j = 3*(nb-1)
           do i=1,3
              dtau(i)=tau(i,na)-tau(i,nb)
           end do
           !
           ! generates nearest-neighbors shells r(i)=R(i)-dtau(i)
           !
           call rgen(dtau,rmax,mxr,at,bg,r,r2,nrm)
           do nr=1, nrm
              rr=sqrt(r2(nr))*alat
              ar = sqrt(alpha)*rr
              d2f = ( 3.d0*erfc(ar) + sqrt(8.d0/tpi)*ar*              &
                     (3.d0+2.d0*ar**2)*exp(-ar**2) ) / rr**5
              df  = (    -erfc(ar) - sqrt(8.d0/tpi)*ar*exp(-ar**2) ) / rr**3
              do nu_i = 1,nmodes
                 if (has_equivalent( (nu_i-1)/3+1 ).eq.1 ) go to 25
                 arg = r(1,nr)*u(mu_i+1,nu_i) + &
                       r(2,nr)*u(mu_i+2,nu_i) + &
                       r(3,nr)*u(mu_i+3,nu_i)
                 do nu_j = 1,nmodes
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
                 end do
25               continue
              end do
           end do
        end if
     end do
  end do
  !
30 continue
#ifdef __PARA
  call mp_sum( dyn, intra_pool_comm )
#endif
  return
end subroutine d2ion
