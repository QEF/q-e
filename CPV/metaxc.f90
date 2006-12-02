!
! Copyright (C) 2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE tpssmeta(nnr, nspin,grho,rho,kedtau,etxc)
  !     ===================
  !--------------------------------------------------------------------
  use gvecp, only: ng => ngm
  IMPLICIT NONE
  !
  ! input
  integer nspin , nnr
  real(8)  grho(nnr,3,nspin), rho(nnr,nspin),kedtau(nnr,nspin)
  ! output: excrho: exc * rho ;  E_xc = \int excrho(r) d_r
  ! output: rhor:   contains the exchange-correlation potential
  real(8)  etxc
  REAL(8) :: zeta, rh, grh2
  INTEGER :: k, ipol, is
  REAL(8) :: grho2 (2), sx, sc, v1x, v2x, v3x,v1c, v2c, v3c, &
       v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw ,v2cup(3),v2cdw(3), &
       v3xup, v3xdw,grhoup(3),grhodw(3),&
  segno, arho, atau
  REAL(8), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10
  etxc = 0.d0
  ! calculate the gradient of rho+rho_core in real space
  DO k = 1, nnr
     DO is = 1, nspin
        grho2 (is) = grho(k,1, is)**2 + grho(k,2,is)**2 + grho(k,3, is)**2
     ENDDO
     IF (nspin == 1) THEN
        !
        !    This is the spin-unpolarised case
        !
        arho = ABS (rho (k, 1) )
        segno = SIGN (1.d0, rho (k, 1) )
        atau = kedtau(k,1)
        IF (arho.GT.epsr.AND.grho2 (1) .GT.epsg.AND.ABS(atau).GT.epsr) THEN
           CALL tpsscxc (arho, grho2(1),atau,sx, sc, &
                v1x, v2x,v3x,v1c, v2c,v3c)
           rho (k, 1) =  (v1x + v1c )
           kedtau(k,1)=  (v3x + v3c) *0.5d0
           ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
           DO ipol = 1, 3  
              grho(k,ipol,1) =  (v2x + v2c)*grho (k,ipol,1) 
           ENDDO
           etxc = etxc +  (sx + sc) * segno 
        ELSE  
           DO ipol = 1, 3  
              grho (k, ipol, 1) = 0.d0  
           ENDDO
           kedtau(k,1)=0.d0
        ENDIF
     ELSE
        !
        !    spin-polarised case
        !
        CALL tpsscx_spin(rho (k, 1), rho (k, 2), grho2 (1), grho2 (2), &
             kedtau(k,1),kedtau(k,2),sx, &
             v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
        rh = rho (k, 1) + rho (k, 2)
        IF (rh.GT.epsr) THEN
           zeta = (rho (k, 1) - rho (k, 2) ) / rh
           DO ipol=1,3
              grhoup(ipol)=grho(k,ipol,1)
              grhodw(ipol)=grho(k,ipol,2)
           END DO
           atau=kedtau(k,1)+kedtau(k,2)
           CALL tpsscc_spin(rh,zeta,grhoup,grhodw, &
                atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3c)
        ELSE
           sc = 0.d0  
           v1cup = 0.d0  
           v1cdw = 0.d0  
           v2cup=0.d0
           v2cdw=0.d0
           v3c=0.d0
           !
        ENDIF
        !
        ! first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        rho(k, 1) =  (v1xup + v1cup)
        rho(k, 2) =  (v1xdw + v1cdw) 
        !
        ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        DO ipol = 1, 3  
           grho(k,ipol,1) = (v2xup*grho(k,ipol,1) + v2cup(ipol))
           grho(k,ipol,2) = (v2xdw*grho(k,ipol,2) + v2cdw(ipol)) 
        ENDDO
        kedtau(k,1)=  (v3xup + v3c) *0.5d0
        kedtau(k,2)=  (v3xdw + v3c) *0.5d0
        etxc = etxc +  (sx + sc)
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE tpssmeta

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine tpsscxc(rho,grho,tau, &
     sx,sc,v1x,v2x,v3x,v1c,v2c,v3c)
  !-----------------------------------------------------------------------
  !     tpss metaGGA corrections for exchange and correlation - Hartree a.u.
  !
  !     input:  rho, grho=|\nabla rho|^2
  !     definition:  E_x = \int E_x(rho,grho) dr
  !     output: sx = E_x(rho,grho)
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             sc, v1c, v2c as above for correlation
  !c     v3x= D(E_x)/D(tau)
  !     tau is the kinetic energy density
  !     the same applies to correlation terms
  !
  !     input grho is |\nabla rho|**2
  !
  implicit none  

  real(8) :: rho, grho, tau,sx, sc, v1x, v2x,v3x,v1c,v2c,v3c
  real(8) :: small  
  parameter (small = 1.d-10) 
  ! exchange
  if (rho.le.small) then  
     sx = 0.0d0  
     v1x = 0.0d0  
     v2x = 0.0d0  
     sc = 0.0d0  
     v1c = 0.0d0  
     v2c = 0.0d0  
     v3x = 0.0d0
     v3c=0.d0
     return
  end if
  call metax(rho,grho,tau,sx,v1x,v2x,v3x)
  !
  call metac(rho,grho,tau,sc,v1c,v2c,v3c)
  !
  return  
end subroutine tpsscxc
!c     ==================================================================
subroutine metax(rho,grho2,tau,ex,v1x,v2x,v3x)
  !     ==--------------------------------------------------------------==
  !     ==  meta-GGA exchange potential and energy                            ==
  !  ==                                                              ==
  !  ==--------------------------------------------------------------==
  
  !  NOTA BENE:            E_x(rho,grho)=rho\epsilon_x(rho,grho)
  !  ex = E_x(rho,grho)    NOT \epsilon_x(rho,grho)
  !  v1x= D(E_x)/D(rho)
  !  v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !  v3x= D(E_x)/D( tau )
  !  tau is the kinetic energy density
  !  the same applies to correlation terms
  !  input grho2 is |\nabla rho|^2
  implicit none
  !  INPUT
  real(8) rho,grho2,tau,rs
  !  OUTPUT
  real(8) ex,v1x,v2x,v3x
  !  LOCAL
  real(8) vx_unif,ex_unif
  !  ex_unif:   lda \epsilon_x(rho)
  !  ec_unif:   lda \epsilon_c(rho)
  real(8) small, pi34, third
  parameter (small=1.d-10)
  parameter (pi34 = 0.6203504908994d0, third = 1.d0 / 3.d0)
  !  fx=Fx(p,z)
  !  fxp=d Fx / d p
  !  fxz=d Fx / d z
  real(8) fx,f1x,f2x,f3x
 
  !  ==--------------------------------------------------------------==
  if(abs(tau).lt.small) then
    ex=0.d0
    v1x=0.d0
    v2x=0.d0
    v3x=0.d0
    return
  endif
  rs = pi34/rho**third
  call slater(rs,ex_unif,vx_unif)
  call metaFX(rho,grho2,tau,fx,f1x,f2x,f3x)
  ex =rho*ex_unif
  v1x=vx_unif*fx + ex*f1x
  v2x=ex*f2x
  v3x=ex*f3x
  ex =ex*fx
  
  !  ==--------------------------------------------------------------==
  return
end subroutine metax
!==   ------------------------------------------------------------------
subroutine metac(rho,grho2,tau,ec,v1c,v2c,v3c)
  !==   ------------------------------------------------------------------
  !  meta-GGA correlation energy and potentials
  !==   ------------------------------------------------------------------
  implicit none
  !  INPUT
  real(8) rho, grho2, tau
  !  OUTPUT
  real(8) ec, v1c,v2c,v3c
  !  LOCAL
  real(8) z,z2,tauw,ec_rev,rs
  real(8) d1rev, d2rev, d3rev
  !  d1ec=  D ec_rev / D rho
  !  d2ec=  D ec_rev / D |D rho/ D r| / |\nabla rho|
  !  d3ec=  D ec_rev / D tau
  real(8) cf1,cf2,cf3
  real(8) v1c_pbe, v2c_pbe, ec_pbe
  real(8) v1c_sum, v2c_sum, ec_sum
  real(8) vc_unif,ec_unif
  real(8) dd,cab,cabone
  real(8) rhoup,grhoup,dummy
  real(8) small, pi34,third
  parameter(small=1.0d-10)
  parameter (pi34= 0.75d0/3.141592653589793d+00, &
       third=1.d0/3.d0)
  parameter (dd=2.80d0)  !in unit of Hartree^-1
  parameter (cab=0.53d0, cabone=1.d0+cab)
  !
  if(abs(tau).lt.small) then
     ec=0.d0
     v1c=0.d0
     v2c=0.d0
     v3c=0.d0
     return
  endif
  rhoup=0.5d0*rho
  grhoup=0.5d0*sqrt(grho2)
  if(rhoup.gt.small) then
     call pw_spin((pi34/rhoup)**third,1.d0,ec_unif,vc_unif,dummy)
     if(abs(grhoup).gt.small) then
!1.d0-small to avoid pow_e of 0 in pbec_spin
        call pbec_spin(rhoup,1.d0-small,grhoup**2, & 
             ec_sum,v1c_sum,dummy,v2c_sum)
     else
        ec_sum=0.d0
        v1c_sum=0.d0
        v2c_sum=0.d0
     endif
     ec_sum = ec_sum/rhoup + ec_unif
     v1c_sum = (v1c_sum + vc_unif-ec_sum)/rho !rho, not rhoup
     v2c_sum = v2c_sum/(2.d0*rho)
  else
     ec_sum=0.d0
     v1c_sum=0.d0
     v2c_sum=0.d0
  endif
  !
  rs = (pi34/rho)**third
  call pw (rs, 1, ec_unif, vc_unif)
  !  PBE correlation energy and potential
  !  ec_pbe=rho*H,  not rho*(epsion_c_uinf + H)
  !  v1c_pbe=D (rho*H) /D rho
  !  v2c_pbe= for rho, 2 for 
  call pbec(rho,grho2,ec_pbe,v1c_pbe,v2c_pbe)
  ec_pbe=ec_pbe/rho+ec_unif
  v1c_pbe=(v1c_pbe+vc_unif-ec_pbe)/rho
  v2c_pbe=v2c_pbe/rho
  !
  if(ec_sum .lt. ec_pbe) then
     ec_sum = ec_pbe
     v1c_sum= v1c_pbe
     v2c_sum= v2c_pbe
  endif
  !
  tauw=0.1250d0*grho2/rho
  z=tauw/tau
  z2=z*z
  !  
  ec_rev=ec_pbe*(1+cab*z2)-cabone*z2*ec_sum
  d1rev = v1c_pbe + (cab*v1c_pbe-cabone*v1c_sum)*z2  &
       -(ec_pbe*cab - ec_sum*cabone)*2.d0*z2/rho
  d2rev = v2c_pbe + (cab*v2c_pbe-cabone*v2c_sum)*z2  &
       +(ec_pbe*cab - ec_sum*cabone)*4.d0*z2/grho2
  d3rev = -(ec_pbe*cab - ec_sum*cabone)*2.d0*z2/tau
  !
  cf1=1.d0+dd*ec_rev*z2*z
  cf2=rho*(1.d0+2.d0*z2*z*dd*ec_rev)
  cf3=ec_rev*ec_rev*3.d0*dd*z2*z
  v1c=ec_rev*cf1 + cf2*d1rev-cf3
  !
  cf3=cf3*rho
  v2c=cf2*d2rev + cf3*2.d0/grho2
  v3c=cf2*d3rev - cf3/tau

  ec=rho*ec_rev*(1.d0+dd*ec_rev*z2*z)  !-rho*ec_unif
  v1c=v1c !-vc_unif
  !  ==--------------------------------------------------------------==
  return
end subroutine metac
!  ==================================================================
subroutine metaFX(rho,grho2,tau,fx,f1x,f2x,f3x)
  !  ==================================================================
  implicit none
  !  INPUT
  !  charge density, square of gradient of rho, and kinetic energy density
  real(8) rho, grho2, tau
  !  OUTPUT
  !  fx = Fx(p,z)
  !  f1x=D (Fx) / D rho
  !  f2x=D (Fx) / D ( D rho/D r_alpha) /|nabla rho|
  !  f3x=D (Fx) / D tau
  real(8) fx, f1x, f2x, f3x
  !  LOCAL
  real(8) x, p, z, qb, al, localdp, dz
  real(8) dfdx, dxdp, dxdz, dqbdp, daldp, dqbdz, daldz
  real(8) fxp, fxz  ! fxp =D fx /D p
  real(8) tauw, tau_unif
  !  work variables
  real(8) xf1,xf2
  real(8) xfac1, xfac2, xfac3,xfac4,xfac5,xfac6,xfac7,z2
  !
  real(8) pi, THRD, ee, cc, kk, bb,miu,fac1,small
  parameter(pi=3.141592653589793d0)
  parameter(THRD=0.3333333333333333d0)
  parameter(ee=1.537d0)
  parameter(cc=1.59096d0)
  parameter(kk=0.804d0)
  parameter(bb=0.40d0)
  parameter(miu=0.21951d0)
  parameter(fac1=9.57078000062731d0) !fac1=(3*pi^2)^(2/3)
  parameter(small=1.00d-6)
  !==------------------------------------------------------------- 
  tauw=0.125d0*grho2/rho
  z=tauw/tau
  
  p=sqrt(grho2)/rho**THRD/rho
  p=p*p/(fac1*4.0d0)
  tau_unif=0.3d0*fac1*rho**(5.0d0/3.0d0)
  al=(tau-tauw)/tau_unif
  al=abs(al)  !make sure al is always .gt. 0.d0
  qb=0.45d0*(al-1.d0)/sqrt(1.d0+bb*al*(al-1.d0))
  qb=qb+2.0d0*THRD*p
  
  !  calculate x(p,z) and fx
  z2=z*z
  xf1=10.0d0/81.0d0
  xfac1=xf1+cc*z2/(1+z2)**2.d0
  xfac2=146.0d0/2025.0d0
  xfac3=sqrt(0.5d0*(0.36d0*z2+p*p))
  xfac4=xf1*xf1/kk
  xfac5=2.d0*sqrt(ee)*xf1*0.36d0
  xfac6=xfac1*p+xfac2*qb**2.0d0-73.0d0/405.0d0*qb*xfac3
  xfac6=xfac6+xfac4*p**2.0d0+xfac5*z2+ee*miu*p**3.d0
  xfac7=(1+sqrt(ee)*p)
  x=xfac6/(xfac7*xfac7)
  !  fx=kk-kk/(1.d0+x/kk)
  fx=1.d0 + kk-kk/(1.d0+x/kk)
  
  !  calculate the derivatives of fx w.r.t p and z
  dfdx=(kk/(kk+x))**2.d0
  daldp=5.0d0*THRD*(tau/tauw-1.0d0)
  !   daldz=-0.50d0*THRD*
  !  * (tau/(2.0d0*fac1*rho**THRD*0.1250d0*sqrt(grho2)))**2.0d0 
  daldz=-5.0d0*THRD*p/z2
  dqbdz=0.45d0*(0.50*bb*(al-1.d0)+1.d0)
  dqbdz=dqbdz/(1.d0+bb*al*(al-1.d0))**1.5d0
  
  dqbdp=dqbdz*daldp+2.d0*THRD
  dqbdz=dqbdz*daldz
  !  calculate d x /d p
  xf1=73.0d0/405.0d0/xfac3*0.50d0*qb
  xf2=2.0d0*xfac2*qb-73.0d0/405.d0*xfac3
  
  dxdp=-xf1*p
  dxdp=dxdp+xfac1+xf2*dqbdp
  dxdp=dxdp+2.0d0*xfac4*p
  dxdp=dxdp+3.0d0*ee*miu*p*p
  dxdp=dxdp/(xfac7*xfac7)-2.0d0*x*sqrt(ee)/xfac7
  !  d x/ dz
  dxdz=-xf1*0.36d0*z
  xfac1=cc*2.0d0*z*(1-z2)/(1+z2)**3.0d0
  dxdz=dxdz+xfac1*p+xf2*dqbdz
  dxdz=dxdz+xfac5*2.0d0*z
  dxdz=dxdz/(xfac7*xfac7)
  
  fxp=dfdx*dxdp
  fxz=dfdx*dxdz
  !  calculate f1x
  localdp=-8.0d0*THRD*p/rho  ! D p /D rho
  dz=-z/rho             ! D z /D rho
  f1x=fxp*localdp+fxz*dz
  !  f2x
  localdp=2.0d0/(fac1*4.0d0*rho**(8.0d0/3.0d0))
  dz=2.0d0*0.125d0/(rho*tau)
  f2x=fxp*localdp + fxz*dz
  !   f3x
  localdp=0.d0
  dz=-z/tau
  f3x=fxz*dz
  
  !  

  !==---------------------------------------------------------------
  return
end subroutine metaFX

!-----------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine tpsscx_spin(rhoup,rhodw,grhoup2,grhodw2,tauup,taudw,sx,&
     v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
  !-----------------------------------------------------------------------
  !     TPSS metaGGA for exchange - Hartree a.u.
  !
  implicit none  
  !
  !     dummy arguments
  !
  real(8) :: rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, v1xdw, &
       v2xup, v2xdw
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !
  real(8):: tauup,taudw, &! up and down kinetic energy density 
       v3xup,v3xdw         ! derivatives of exchange wr. tau
  real(8) :: small  
  parameter (small = 1.d-10)  
  real(8) :: rho, sxup, sxdw  
  !
  ! exchange
  rho = rhoup + rhodw  
  if (rhoup.gt.small.and.sqrt(abs(grhoup2)).gt.small &
       .and. abs(tauup).gt.small) then
     call metax(2.d0*rhoup,4.d0*grhoup2, &
          2.0d0*tauup,sxup,v1xup,v2xup,v3xup)
  else
     sxup=0.d0
     v1xup=0.d0
     v2xup=0.d0
     v3xup=0.d0
  endif
  if (rhodw.gt.small.and.sqrt(abs(grhodw2)).gt.small &
       .and. abs(taudw).gt.small) then
     call metax(2.d0*rhodw,4.d0*grhodw2,  &
          2.0d0*taudw,sxdw,v1xdw,v2xdw,v3xdw)
  else
     sxdw=0.d0
     v1xdw=0.d0
     v2xdw=0.d0
     v3xdw=0.d0
  endif
  sx=0.5d0*(sxup+sxdw)
  v2xup=2.d0*v2xup
  v2xdw=2.d0*v2xdw
  !
  return  
end subroutine tpsscx_spin
!
!-----------------------------------------------------------------------
subroutine tpsscc_spin(rho,zeta,grhoup,grhodw, &
     atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3c)
!-----------------------------------------------------------------------
!     tpss metaGGA for correlations - Hartree a.u.
!
  implicit none  
!
!     dummy arguments
!
  real(8) :: rho, zeta, grhoup(3),grhodw(3), sc, v1cup, v1cdw, v2c  
  ! the total charge
  ! the magnetization
  ! the gradient of the charge
  ! exchange and correlation energies
  ! derivatives of correlation wr. rho
  ! derivatives of correlation wr. grho
  real(8) :: atau,v2cup(3),v2cdw(3),v3c,grho_vec(3),grho !grho=grho2
  real(8) :: small  
  integer :: ipol
  parameter (small = 1.d-10)  
  !
  !
! vector
  grho_vec=grhoup+grhodw
  grho=0.d0
  do ipol=1,3
     grho = grho + grho_vec(ipol)**2
  end do
!
!
  if (rho.le.small.or.abs (zeta) .gt.1.d0.or.sqrt (abs (grho) ) &
       .le.small.or.abs(atau).lt.small) then
     sc = 0.0d0  
     v1cup = 0.0d0  
     v1cdw = 0.0d0  
     v2cup(:) = 0.0d0  
     v2cdw(:) = 0.0d0
     v3c = 0.d0
  else
     call metac_spin(rho,zeta,grhoup,grhodw, &
          atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3c)
  end if
  !
  return  
end subroutine tpsscc_spin
!

!---------------------------------------------------------------
subroutine metac_spin(rho,zeta,grhoup,grhodw, &
     tau,sc,v1up,v1dw,v2up,v2dw,v3)
!---------------------------------------------------------------
  implicit none
!     input
  real(8) :: rho, zeta,grhoup(3),grhodw(3), tau
!     output
  real(8) :: sc, v1up, v1dw, v2up(3), v2dw(3), v3
!     local
  real(8) :: rhoup, rhodw,tauw,grhovec(3),grho2,grho,&
       grhoup2,grhodw2
  !grhovec   vector gradient of rho
  !grho    mod of gradient of rho
  real(8) :: ec_u, vcup_u, vcdw_u 
  real(8) :: ec_pbe, v1up_pbe, v1dw_pbe,v2up_pbe(3),v2dw_pbe(3)
  real(8) :: ecup_0, v1up_0, v2up_0(3),v2_tmp
  real(8) :: ecdw_0, v1dw_0, v2dw_0(3)
  real(8) :: ec_rev, cab, aa, bb, aa2
  real(8) :: z2,z,ca0,dca0da,dcabda,dcabdb
  real(8) :: term(3),term1,term2,term3
  real(8) :: drev1up, drev1dw,drev2up(3),drev2dw(3),drev3
  real(8) :: sum, dsum1up, dsum1dw,dsum2up(3),dsum2dw(3)
  real(8) :: dcab1up, dcab1dw,dcab2up(3),dcab2dw(3)
  real(8) :: db1up,   db1dw,  db2up(3),  db2dw(3)
  real(8) :: da1up,   da1dw
  real(8) :: ecup_til,ecdw_til
  real(8) :: v1up_uptil, v1up_dwtil, v2up_uptil(3),v2up_dwtil(3)
  real(8) :: v1dw_uptil, v1dw_dwtil, v2dw_uptil(3),v2dw_dwtil(3)
  real(8) :: small, pi34, p43, third, fac
  parameter(small=1.0d-10, &
       fac=3.09366772628013593097d0**2)
!     fac = (3*PI**2)**(2/3)
  parameter (pi34= 0.75d0 / 3.141592653589793d+00, &
       p43=4.0d0/3.d0,third=1.d0/3.d0)
  integer:: ipol
!-----------
  rhoup=(1+zeta)*0.50d0*rho
  rhodw=(1-zeta)*0.50d0*rho
  grho2=0.d0
  grhoup2=0.d0
  grhodw2=0.d0
  do ipol=1,3
     grhovec(ipol)=grhoup(ipol)+grhodw(ipol)
     grho2=grho2+grhovec(ipol)**2
     grhoup2=grhoup2+grhoup(ipol)**2
     grhodw2=grhodw2+grhodw(ipol)**2
  end do 
  grho=sqrt(grho2)
!     
  if(rho.gt.small) then
     v2_tmp=0.d0
     call pw_spin((pi34/rho)**third,zeta,ec_u,vcup_u,vcdw_u)
     if((abs(grho).gt.small) .and. (zeta .le. 1.d0)) then
        call pbec_spin(rho,zeta,grho2, &
             ec_pbe,v1up_pbe,v1dw_pbe,v2_tmp)
     else
        ec_pbe=0.d0
        v1up_pbe=0.d0
        v1dw_pbe=0.d0
        v2up_pbe=0.d0
     endif
     ec_pbe = ec_pbe/rho+ec_u
!     v1xx_pbe = D_epsilon_c/ D_rho_xx   :xx= up, dw
     v1up_pbe = (v1up_pbe+vcup_u-ec_pbe)/rho
     v1dw_pbe = (v1dw_pbe+vcdw_u-ec_pbe)/rho
!     v2xx_pbe = (D_Ec / D grho)/rho = (D_Ec/ D |grho| /|grho|)*grho/rho
     v2up_pbe = v2_tmp/rho*grhovec
!    v2dw === v2up for PBE
     v2dw_pbe = v2up_pbe
  else
     ec_pbe=0.d0
     v1up_pbe=0.d0
     v1dw_pbe=0.d0
     v2up_pbe=0.d0
     v2dw_pbe=0.d0
  endif
!     ec_pbe(rhoup,0,grhoup,0)
  if(rhoup.gt.small) then
     v2_tmp=0.d0
     call pw_spin((pi34/rhoup)**third,1.d0,ec_u,vcup_u,vcdw_u)
     if(sqrt(grhoup2).gt.small) then
        call pbec_spin(rhoup,1.d0-small,grhoup2,&
             ecup_0,v1up_0,v1dw_0,v2_tmp)
     else
        ecup_0=0.d0
        v1up_0=0.d0
        v2up_0=0.d0
     endif
     ecup_0 = ecup_0/rhoup + ec_u
     v1up_0 = (v1up_0 + vcup_u-ecup_0)/rhoup
     v2up_0 = v2_tmp/rhoup*grhoup
  else
     ecup_0 = 0.d0
     v1up_0 = 0.d0
     v2up_0=0.d0
  endif
!
  if(ecup_0.gt.ec_pbe) then
     ecup_til = ecup_0
     v1up_uptil=v1up_0
     v2up_uptil=v2up_0
     v1up_dwtil=0.d0
     v2up_dwtil=0.d0
  else
     ecup_til = ec_pbe
     v1up_uptil= v1up_pbe
     v1up_dwtil= v1dw_pbe
     v2up_uptil= v2up_pbe
     v2up_dwtil= v2up_pbe
  endif
!     ec_pbe(rhodw,0,grhodw,0)
!     zeta = 1.d0
  if(rhodw.gt.small) then
     v2_tmp=0.d0
     call pw_spin((pi34/rhodw)**third,-1.d0,ec_u,vcup_u,vcdw_u)
     if(sqrt(grhodw2).gt.small) then
        call pbec_spin(rhodw,-1.d0+small,grhodw2,&
             ecdw_0,v1up_0,v1dw_0,v2_tmp)
     else
        ecdw_0=0.d0
        v1dw_0=0.d0
        v2dw_0=0.d0
     endif
     ecdw_0 = ecdw_0/rhodw + ec_u
     v1dw_0 = (v1dw_0 + vcdw_u-ecdw_0)/rhodw
     v2dw_0 = v2_tmp/rhodw*grhodw
  else
     ecdw_0 = 0.d0
     v1dw_0 = 0.d0
     v2dw_0 = 0.d0
  endif
!     
  if(ecdw_0.gt.ec_pbe) then
     ecdw_til = ecdw_0
     v1dw_dwtil=v1dw_0
     v2dw_dwtil=v2dw_0
     v1dw_uptil=0.d0
     v2dw_uptil=0.d0
  else
     ecdw_til = ec_pbe
     v1dw_dwtil= v1dw_pbe
     v2dw_dwtil= v2dw_pbe
     v1dw_uptil= v1up_pbe
     v2dw_uptil= v2dw_pbe
  endif
!cccccccccccccccccccccccccccccccccccccccccc-------checked
  sum=(rhoup*ecup_til+rhodw*ecdw_til)/rho
  dsum1up=(ecup_til-ecdw_til)*rhodw/rho**2 &
       +  (rhoup*v1up_uptil + rhodw*v1dw_uptil)/rho
  dsum1dw=(ecdw_til-ecup_til)*rhoup/rho**2 &
       +  (rhodw*v1dw_dwtil + rhoup*v1up_dwtil)/rho
! vector  
  dsum2up=(rhoup*v2up_uptil + rhodw*v2dw_uptil)/rho
  dsum2dw=(rhodw*v2dw_dwtil + rhoup*v2up_dwtil)/rho
!ccccccccccccccccccccccccccccccccccccccccc---------checked
  aa=zeta
!  bb=(rho*(grhoup-grhodw) - (rhoup-rhodw)*grho)**2 &
!       /(4.d0*fac*rho**(14.d0/3.d0))
  bb=0.d0
  do ipol=1,3
     term(ipol)= rhodw*grhoup(ipol)-rhoup*grhodw(ipol)
     bb=bb+ term(ipol)**2
  end do
!vector
  term=term/(fac*rho**(14.d0/3.d0))   
  bb=bb/(fac*rho**(14.d0/3.d0))
!      bb=(rhodw*grhoup-rhoup*grhodw)**2/fac*rho**(-14.0d0/3.d0)
  aa2=aa*aa
  ca0=0.53d0+aa2*(0.87d0+aa2*(0.50d0+aa2*2.26d0))
  dca0da = aa*(1.74d0+aa2*(2.0d0+aa2*13.56d0))
  if(abs(aa).le.1.d0-small) then
     term3 =(1.d0+aa)**(-p43) + (1.d0-aa)**(-p43)
     term1=(1.d0+bb*0.50*term3)
     term2=(1.d0+aa)**(-7.d0/3.d0) + (1.d0-aa)**(-7.d0/3.d0)
     cab =ca0/term1**4
     dcabda = (dca0da/ca0 + 8.d0/3.d0*bb*term2/term1)*cab
     dcabdb = -2.d0*cab*term3/term1
  else
     cab=0.d0
     dcabda=0.d0
     dcabdb=0.d0
  endif
  da1up=2.d0*rhodw/rho**2
  da1dw=-2.d0*rhoup/rho**2
  db1up=-2.d0*(grhodw(1)*term(1)+grhodw(2)*term(2)+grhodw(3)*term(3)) &
       -14.d0/3.d0*bb/rho
  db1dw= 2.d0*(grhoup(1)*term(1)+grhoup(2)*term(2)+grhoup(3)*term(3)) &
       -14.d0/3.d0*bb/rho
  !vector, not scalar
  db2up= term*rhodw*2.d0
  db2dw=-term*rhoup*2.d0
!
  dcab1up = dcabda*da1up + dcabdb*db1up
  dcab1dw = dcabda*da1dw + dcabdb*db1dw
  !vector, not scalar
  dcab2up = dcabdb*db2up
  dcab2dw = dcabdb*db2dw
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc------checked
  tauw=0.1250d0*grho2/rho
  z=tauw/tau
  z2=z*z
!     
  term1=1.d0+cab*z2
  term2=(1.d0+cab)*z2
  ec_rev = ec_pbe*term1-term2*sum
!     
  drev1up=v1up_pbe*term1 + &
       ec_pbe*(z2*dcab1up - 2.d0*cab*z2/rho) &
       + (2.d0*term2/rho - z2*dcab1up)*sum &
       - term2*dsum1up 
!     
  drev1dw=v1dw_pbe*term1 + &
       ec_pbe*(z2*dcab1dw - 2.d0*cab*z2/rho)  &
       + (2.d0*term2/rho - z2*dcab1dw)*sum  &
       - term2*dsum1dw
!
! vector, not scalar
  drev2up=v2up_pbe*term1 + &
       ec_pbe*(z2*dcab2up+0.5d0*cab*z/(rho*tau)*grhovec)& 
       - (term2*4.d0/grho2*grhovec + z2*dcab2up)*sum &
       - term2*dsum2up
  drev2dw=v2dw_pbe*term1 +  &
       ec_pbe*(z2*dcab2dw+0.5d0*cab*z/(rho*tau)*grhovec) &
       - (term2*4.d0/grho2*grhovec + z2*dcab2dw)*sum  &
       - term2*dsum2dw
!
  drev3 = ((1.d0+cab)*sum-ec_pbe*cab)*2.d0*z2/tau
!ccccccccccccccccccccccccccccccccccccccccccccccccccc----checked
  term1=ec_rev*(1.d0+2.8d0*ec_rev*z2*z)
  term2=(1.d0+5.6d0*ec_rev*z2*z)*rho
  term3=-8.4d0*ec_rev*ec_rev*z2*z
!
  v1up = term1 + term2*drev1up + term3
  v1dw = term1 + term2*drev1dw + term3
!
  term3=term3*rho
  v3   = term2*drev3 + term3/tau
!
  term3=-2.d0*term3/grho2    !grho/|grho|^2 = 1/grho
  v2up = term2*drev2up + term3*grhovec
  v2dw = term2*drev2dw + term3*grhovec
  !
  !
  !  call pw_spin((pi34/rho)**third,zeta,ec_u,vcup_u,vcdw_u)
  sc=rho*ec_rev*(1.d0+2.8d0*ec_rev*z2*z) !-rho*ec_u
  !  v1up=v1up-vcup_u
  !  v1dw=v1dw-vcdw_u
  
  return
end subroutine metac_spin
