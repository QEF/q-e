!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------------
!
!                        META-GGA FUNCTIONALS
!
!  Available functionals : 
!           - TPSS (Tao, Perdew, Staroverov & Scuseria)
!           - TB09 (via libxc)
!           - M06L
!
!=========================================================================
!   
!-------------------------------------------------------------------------
!
!                             TPSS
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine tpsscxc( rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c )
  !-----------------------------------------------------------------------
  !     TPSS metaGGA corrections for exchange and correlation - Hartree a.u.
  !
  !       input:  rho, grho=|\nabla rho|^2, tau = kinetic energy density
  !               definition:  E_x = \int E_x(rho,grho) dr
  !               output: sx = E_x(rho,grho)
  !                       v1x= D(E_x)/D(rho)
  !                       v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !                       sc, v1c, v2c as above for correlation
  !                       v3x= D(E_x)/D(tau)
  !
  USE kinds,            ONLY : DP
#if defined(__LIBXC)
  use xc_f90_types_m
  use xc_f90_lib_m
#endif
  implicit none  
  real(DP), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c
#if defined(__LIBXC)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  integer :: size = 1
  integer :: func_id = 202  !
  real(dp) :: lapl_rho, vlapl_rho ! not used in TPSS
 
  lapl_rho = grho

  ! exchange  
  func_id = 202
  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)    
  call xc_f90_mgga_exc_vxc(xc_func, size, rho, grho, lapl_rho, 0.5_dp*tau,&
                               sx, v1x, v2x, vlapl_rho, v3x)  
  call xc_f90_func_end(xc_func)	
  
  sx  = sx * rho
  v2x = v2x*2.0_dp
  v3x = 0.5_dp*v3x

  ! correlation
  func_id = 231  ! Perdew, Tao, Staroverov & Scuseria correlation  
  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)   
  call xc_f90_mgga_exc_vxc(xc_func,size , rho, grho, lapl_rho, 0.5_dp*tau,&
                               sc, v1c, v2c, vlapl_rho, v3c)  
  call xc_f90_func_end(xc_func)

  sc  = sc * rho
  v2c = v2c*2.0_dp
  v3c = 0.5_dp*v3c       

#else
  real(DP), parameter :: small = 1.E-10_DP

  if (rho.le.small) then  
     sx = 0.0_DP  
     v1x = 0.0_DP  
     v2x = 0.0_DP 
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
     v3x = 0.0_DP
     v3c=0.0_DP
     return
  end if
  ! exchange
  call metax(rho,grho,tau,sx,v1x,v2x,v3x)
  ! correlation
  call metac(rho,grho,tau,sc,v1c,v2c,v3c)
  !
#endif
  !
  return  
end subroutine tpsscxc
!-------------------------------------------------------------------------
subroutine metax(rho,grho2,tau,ex,v1x,v2x,v3x)
  !    --------------------------------------------------------------==
  !     ==  TPSS meta-GGA exchange potential and energy                    
  !  ==                                                              ==
  !  ==--------------------------------------------------------------==
  
  USE kinds,            ONLY : DP

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
  real(DP) :: rho,grho2,tau,rs
  !  OUTPUT
  real(DP) :: ex,v1x,v2x,v3x
  !  LOCAL
  real(DP) :: vx_unif,ex_unif
  !  ex_unif:   lda \epsilon_x(rho)
  !  ec_unif:   lda \epsilon_c(rho)
  real(DP) :: small, pi34, third
  parameter (small=1.E-10_DP)
  parameter (pi34 = 0.6203504908994_DP, third = 1.0_DP / 3.0_DP)
  !  fx=Fx(p,z)
  !  fxp=d Fx / d p
  !  fxz=d Fx / d z
  real(DP) fx,f1x,f2x,f3x
 
  !  ==--------------------------------------------------------------==
  if(abs(tau).lt.small) then
    ex=0.0_DP
    v1x=0.0_DP
    v2x=0.0_DP
    v3x=0.0_DP
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
  !  TPSS meta-GGA correlation energy and potentials
  !==   ------------------------------------------------------------------
  USE kinds,            ONLY : DP
  implicit none
  !  INPUT
  real(DP) :: rho, grho2, tau
  !  OUTPUT
  real(DP) :: ec, v1c,v2c,v3c
  !  LOCAL
  real(DP) :: z,z2,tauw,ec_rev,rs
  real(DP) :: d1rev, d2rev, d3rev
  !  d1ec=  D ec_rev / D rho
  !  d2ec=  D ec_rev / D |D rho/ D r| / |\nabla rho|
  !  d3ec=  D ec_rev / D tau
  real(DP) :: cf1,cf2,cf3
  real(DP) :: v1c_pbe, v2c_pbe, ec_pbe
  real(DP) :: v1c_sum, v2c_sum, ec_sum
  real(DP) :: vc_unif,ec_unif
  real(DP) :: dd,cab,cabone
  real(DP) :: rhoup,grhoup,dummy
  real(DP) :: small, pi34,third
  parameter(small=1.0E-10_DP)
  parameter (pi34= 0.75_DP/3.141592653589793_DP, &
       third=1.0_DP/3.0_DP)
  parameter (dd=2.80_DP)  !in unit of Hartree^-1
  parameter (cab=0.53_DP, cabone=1.0_DP+cab)
  !
   if(abs(tau).lt.small) then
     ec=0.0_DP
     v1c=0.0_DP
     v2c=0.0_DP
     v3c=0.0_DP
     return
  endif
  rhoup=0.5_DP*rho
  grhoup=0.5_DP*SQRT(grho2)
  if(rhoup.gt.small) then
     call pw_spin((pi34/rhoup)**third,1.0_DP,ec_unif,vc_unif,dummy)
     if(abs(grhoup).gt.small) then
!1.0_DP-small to avoid pow_e of 0 in pbec_spin
        call pbec_spin(rhoup,1.0_DP-small,grhoup**2,1,& 
             ec_sum,v1c_sum,dummy,v2c_sum)
     else
        ec_sum=0.0_DP
        v1c_sum=0.0_DP
        v2c_sum=0.0_DP
     endif
     ec_sum = ec_sum/rhoup + ec_unif
     v1c_sum = (v1c_sum + vc_unif-ec_sum)/rho !rho, not rhoup
     v2c_sum = v2c_sum/(2.0_DP*rho)
  else
     ec_sum=0.0_DP
     v1c_sum=0.0_DP
     v2c_sum=0.0_DP
  endif
  !
  rs = (pi34/rho)**third
  call pw (rs, 1, ec_unif, vc_unif)
  !  PBE correlation energy and potential
  !  ec_pbe=rho*H,  not rho*(epsion_c_uinf + H)
  !  v1c_pbe=D (rho*H) /D rho
  !  v2c_pbe= for rho, 2 for 
  call pbec(rho,grho2,1,ec_pbe,v1c_pbe,v2c_pbe)
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
  tauw=0.1250_DP*grho2/rho
  z=tauw/tau
  z2=z*z
  !  
  ec_rev=ec_pbe*(1+cab*z2)-cabone*z2*ec_sum
  d1rev = v1c_pbe + (cab*v1c_pbe-cabone*v1c_sum)*z2  &
       -(ec_pbe*cab - ec_sum*cabone)*2.0_DP*z2/rho
  d2rev = v2c_pbe + (cab*v2c_pbe-cabone*v2c_sum)*z2  &
       +(ec_pbe*cab - ec_sum*cabone)*4.0_DP*z2/grho2
  d3rev = -(ec_pbe*cab - ec_sum*cabone)*2.0_DP*z2/tau
  !
  cf1=1.0_DP+dd*ec_rev*z2*z
  cf2=rho*(1.0_DP+2.0_DP*z2*z*dd*ec_rev)
  cf3=ec_rev*ec_rev*3.0_DP*dd*z2*z
  v1c=ec_rev*cf1 + cf2*d1rev-cf3
  !
  cf3=cf3*rho
  v2c=cf2*d2rev + cf3*2.0_DP/grho2
  v3c=cf2*d3rev - cf3/tau

  ec=rho*ec_rev*(1.0_DP+dd*ec_rev*z2*z)  !-rho*ec_unif
  v1c=v1c !-vc_unif
  !  ==--------------------------------------------------------------==
  return
end subroutine metac
!-------------------------------------------------------------------------
subroutine metaFX(rho,grho2,tau,fx,f1x,f2x,f3x)
  !-------------------------------------------------------------------------
  USE kinds,            ONLY : DP
  implicit none
  !  INPUT
  !  charge density, square of gradient of rho, and kinetic energy density
  real(DP) rho, grho2, tau
  !  OUTPUT
  !  fx = Fx(p,z)
  !  f1x=D (Fx) / D rho
  !  f2x=D (Fx) / D ( D rho/D r_alpha) /|nabla rho|
  !  f3x=D (Fx) / D tau
  real(DP) fx, f1x, f2x, f3x
  !  LOCAL
  real(DP) x, p, z, qb, al, localdp, dz
  real(DP) dfdx, dxdp, dxdz, dqbdp, daldp, dqbdz, daldz
  real(DP) fxp, fxz  ! fxp =D fx /D p
  real(DP) tauw, tau_unif
  !  work variables
  real(DP) xf1,xf2
  real(DP) xfac1, xfac2, xfac3,xfac4,xfac5,xfac6,xfac7,z2
  !
  real(DP) pi, THRD, ee, cc, kk, bb,miu,fac1,small
  parameter(pi=3.141592653589793_DP)
  parameter(THRD=0.3333333333333333_DP)
  parameter(ee=1.537_DP)
  parameter(cc=1.59096_DP)
  parameter(kk=0.804_DP)
  parameter(bb=0.40_DP)
  parameter(miu=0.21951_DP)
  parameter(fac1=9.57078000062731_DP) !fac1=(3*pi^2)^(2/3)
  parameter(small=1.0E-6_DP)
  !==------------------------------------------------------------- 
  tauw=0.125_DP*grho2/rho
  z=tauw/tau
  
  p=sqrt(grho2)/rho**THRD/rho
  p=p*p/(fac1*4.0_DP)
  tau_unif=0.3_DP*fac1*rho**(5.0_DP/3.0_DP)
  al=(tau-tauw)/tau_unif
  al=abs(al)  !make sure al is always .gt. 0.0_DP
  qb=0.45_DP*(al-1.0_DP)/sqrt(1.0_DP+bb*al*(al-1.0_DP))
  qb=qb+2.0_DP*THRD*p
  
  !  calculate x(p,z) and fx
  z2=z*z
  xf1=10.0_DP/81.0_DP
  xfac1=xf1+cc*z2/(1+z2)**2.0_DP
  xfac2=146.0_DP/2025.0_DP
  xfac3=sqrt(0.5_DP*(0.36_DP*z2+p*p))
  xfac4=xf1*xf1/kk
  xfac5=2.0_DP*sqrt(ee)*xf1*0.36_DP
  xfac6=xfac1*p+xfac2*qb**2.0_DP-73.0_DP/405.0_DP*qb*xfac3
  xfac6=xfac6+xfac4*p**2.0_DP+xfac5*z2+ee*miu*p**3.0_DP
  xfac7=(1+sqrt(ee)*p)
  x=xfac6/(xfac7*xfac7)
  !  fx=kk-kk/(1.0_DP+x/kk)
  fx=1.0_DP + kk-kk/(1.0_DP+x/kk)
  
  !  calculate the derivatives of fx w.r.t p and z
  dfdx=(kk/(kk+x))**2.0_DP
  daldp=5.0_DP*THRD*(tau/tauw-1.0_DP)
  !   daldz=-0.50_DP*THRD*
  !  * (tau/(2.0_DP*fac1*rho**THRD*0.1250_DP*sqrt(grho2)))**2.0_DP
  daldz=-5.0_DP*THRD*p/z2
  dqbdz=0.45_DP*(0.50_DP*bb*(al-1.0_DP)+1.0_DP)
  dqbdz=dqbdz/(1.0_DP+bb*al*(al-1.0_DP))**1.5_DP
  
  dqbdp=dqbdz*daldp+2.0_DP*THRD
  dqbdz=dqbdz*daldz
  !  calculate d x /d p
  xf1=73.0_DP/405.0_DP/xfac3*0.50_DP*qb
  xf2=2.0_DP*xfac2*qb-73.0_DP/405.0_DP*xfac3
  
  dxdp=-xf1*p
  dxdp=dxdp+xfac1+xf2*dqbdp
  dxdp=dxdp+2.0_DP*xfac4*p
  dxdp=dxdp+3.0_DP*ee*miu*p*p
  dxdp=dxdp/(xfac7*xfac7)-2.0_DP*x*sqrt(ee)/xfac7
  !  d x/ dz
  dxdz=-xf1*0.36_DP*z
  xfac1=cc*2.0_DP*z*(1-z2)/(1+z2)**3.0_DP
  dxdz=dxdz+xfac1*p+xf2*dqbdz
  dxdz=dxdz+xfac5*2.0_DP*z
  dxdz=dxdz/(xfac7*xfac7)
  
  fxp=dfdx*dxdp
  fxz=dfdx*dxdz
  !  calculate f1x
  localdp=-8.0_DP*THRD*p/rho  ! D p /D rho
  dz=-z/rho             ! D z /D rho
  f1x=fxp*localdp+fxz*dz
  !  f2x
  localdp=2.0_DP/(fac1*4.0_DP*rho**(8.0_DP/3.0_DP))
  dz=2.0_DP*0.125_DP/(rho*tau)
  f2x=fxp*localdp + fxz*dz
  !   f3x
  localdp=0.0_DP
  dz=-z/tau
  f3x=fxz*dz
  

  return
end subroutine metaFX

!-------------------------------------------------------------------
subroutine tpsscx_spin(rhoup,rhodw,grhoup2,grhodw2,tauup,taudw,sx,&
     v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
  !-----------------------------------------------------------------
  !     TPSS metaGGA for exchange - Hartree a.u.
  !
  USE kinds,            ONLY : DP
  implicit none  
  !
  !     dummy arguments
  !
  real(DP) :: rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, v1xdw, &
       v2xup, v2xdw
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !
  real(DP):: tauup,taudw, &! up and down kinetic energy density 
       v3xup,v3xdw         ! derivatives of exchange wr. tau
  real(DP) :: small  
  parameter (small = 1.E-10_DP)  
  real(DP) :: rho, sxup, sxdw  
  !
  ! exchange
  rho = rhoup + rhodw  
  if (rhoup.gt.small.and.sqrt(abs(grhoup2)).gt.small &
       .and. abs(tauup).gt.small) then
     call metax(2.0_DP*rhoup,4.0_DP*grhoup2, &
          2.0_DP*tauup,sxup,v1xup,v2xup,v3xup)
  else
     sxup=0.0_DP
     v1xup=0.0_DP
     v2xup=0.0_DP
     v3xup=0.0_DP
  endif
  if (rhodw.gt.small.and.sqrt(abs(grhodw2)).gt.small &
       .and. abs(taudw).gt.small) then
     call metax(2.0_DP*rhodw,4.0_DP*grhodw2,  &
          2.0_DP*taudw,sxdw,v1xdw,v2xdw,v3xdw)
  else
     sxdw=0.0_DP
     v1xdw=0.0_DP
     v2xdw=0.0_DP
     v3xdw=0.0_DP
  endif
  sx=0.5_DP*(sxup+sxdw)
  v2xup=2.0_DP*v2xup
  v2xdw=2.0_DP*v2xdw
  !
  return  
end subroutine tpsscx_spin
!
!-----------------------------------------------------------------------
subroutine tpsscc_spin(rho,zeta,grhoup,grhodw, &
     atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3cup, v3cdw)
!-----------------------------------------------------------------------
!     tpss metaGGA for correlations - Hartree a.u.
!
  USE kinds,            ONLY : DP
  implicit none  
!
!     dummy arguments
!
  real(DP) :: rho, zeta, grhoup(3),grhodw(3), sc, v1cup, v1cdw, v3c  
  ! the total charge
  ! the magnetization
  ! the gradient of the charge
  ! exchange and correlation energies
  ! derivatives of correlation wr. rho
  ! derivatives of correlation wr. grho
  
  real(DP), dimension(3) :: v2cup, v2cdw, grho_vec
  real(DP) :: atau,v3cup, v3cdw, grho !grho=grho2
  real(DP) :: small  
  integer :: ipol
  parameter (small = 1.E-10_DP)  
  !
  !
! vector
  grho_vec=grhoup+grhodw
  grho=0.0_DP
  do ipol=1,3
     grho    = grho    + grho_vec(ipol)**2
  end do
!
!
  if (rho.le.small.or.abs (zeta) .gt.1.0_DP.or.sqrt (abs (grho) ) &
       .le.small.or.abs(atau).lt.small) then
     
     sc    = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
 
     v2cup(:) = 0.0_DP
     v2cdw(:) = 0.0_DP

     v3cup = 0.0_DP
     v3cdw = 0.0_DP

     v3c = 0.0_DP
     
  else
     call metac_spin(rho,zeta,grhoup,grhodw, &
          atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3c)
  end if
  !
  !
  v3cup = v3c
  v3cdw = v3c
  !
  return  
end subroutine tpsscc_spin
!

!---------------------------------------------------------------
subroutine metac_spin(rho,zeta,grhoup,grhodw, &
     tau,sc,v1up,v1dw,v2up,v2dw,v3)
!---------------------------------------------------------------
  USE kinds,            ONLY : DP
  implicit none
!     input
  real(DP) :: rho, zeta,grhoup(3),grhodw(3), tau
!     output
  real(DP) :: sc, v1up, v1dw, v2up(3), v2dw(3), v3
!     local
  real(DP) :: rhoup, rhodw,tauw,grhovec(3),grho2,grho,&
       grhoup2,grhodw2
  !grhovec   vector gradient of rho
  !grho    mod of gradient of rho
  real(DP) :: ec_u, vcup_u, vcdw_u 
  real(DP) :: ec_pbe, v1up_pbe, v1dw_pbe,v2up_pbe(3),v2dw_pbe(3)
  real(DP) :: ecup_0, v1up_0, v2up_0(3),v2_tmp
  real(DP) :: ecdw_0, v1dw_0, v2dw_0(3)
  real(DP) :: ec_rev, cab, aa, bb, aa2
  real(DP) :: z2,z,ca0,dca0da,dcabda,dcabdb
  real(DP) :: term(3),term1,term2,term3
  real(DP) :: drev1up, drev1dw,drev2up(3),drev2dw(3),drev3
  real(DP) :: sum, dsum1up, dsum1dw,dsum2up(3),dsum2dw(3)
  real(DP) :: dcab1up, dcab1dw,dcab2up(3),dcab2dw(3)
  real(DP) :: db1up,   db1dw,  db2up(3),  db2dw(3)
  real(DP) :: da1up,   da1dw
  real(DP) :: ecup_til,ecdw_til
  real(DP) :: v1up_uptil, v1up_dwtil, v2up_uptil(3),v2up_dwtil(3)
  real(DP) :: v1dw_uptil, v1dw_dwtil, v2dw_uptil(3),v2dw_dwtil(3)
  real(DP) :: small, pi34, p43, third, fac
  parameter(small=1.0E-10_DP, &
       fac=3.09366772628013593097_DP**2)
!     fac = (3*PI**2)**(2/3)
  parameter (pi34= 0.75_DP / 3.141592653589793_DP, &
       p43=4.0_DP/3.0_DP,third=1.0_DP/3.0_DP)
  integer:: ipol
!-----------
  rhoup=(1+zeta)*0.5_DP*rho
  rhodw=(1-zeta)*0.5_DP*rho
  grho2=0.0_DP
  grhoup2=0.0_DP
  grhodw2=0.0_DP
  do ipol=1,3
     grhovec(ipol)=grhoup(ipol)+grhodw(ipol)
     grho2=grho2+grhovec(ipol)**2
     grhoup2=grhoup2+grhoup(ipol)**2
     grhodw2=grhodw2+grhodw(ipol)**2
  end do 
  grho=sqrt(grho2)
!     
  if(rho.gt.small) then
     v2_tmp=0.0_DP
     call pw_spin((pi34/rho)**third,zeta,ec_u,vcup_u,vcdw_u)
     if((abs(grho).gt.small) .and. (zeta .le. 1.0_DP)) then
        call pbec_spin(rho,zeta,grho2,1,&
             ec_pbe,v1up_pbe,v1dw_pbe,v2_tmp)
     else
        ec_pbe=0.0_DP
        v1up_pbe=0.0_DP
        v1dw_pbe=0.0_DP
        v2up_pbe=0.0_DP
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
     ec_pbe=0.0_DP
     v1up_pbe=0.0_DP
     v1dw_pbe=0.0_DP
     v2up_pbe=0.0_DP
     v2dw_pbe=0.0_DP
  endif
!     ec_pbe(rhoup,0,grhoup,0)
  if(rhoup.gt.small) then
     v2_tmp=0.0_DP
     call pw_spin((pi34/rhoup)**third,1.0_DP,ec_u,vcup_u,vcdw_u)
     if(sqrt(grhoup2).gt.small) then
        call pbec_spin(rhoup,1.0_DP-small,grhoup2,1,&
             ecup_0,v1up_0,v1dw_0,v2_tmp)
     else
        ecup_0=0.0_DP
        v1up_0=0.0_DP
        v2up_0=0.0_DP
     endif
     ecup_0 = ecup_0/rhoup + ec_u
     v1up_0 = (v1up_0 + vcup_u-ecup_0)/rhoup
     v2up_0 = v2_tmp/rhoup*grhoup
  else
     ecup_0 = 0.0_DP
     v1up_0 = 0.0_DP
     v2up_0 = 0.0_DP
  endif
!
  if(ecup_0.gt.ec_pbe) then
     ecup_til = ecup_0
     v1up_uptil=v1up_0
     v2up_uptil=v2up_0
     v1up_dwtil=0.0_DP
     v2up_dwtil=0.0_DP
  else
     ecup_til = ec_pbe
     v1up_uptil= v1up_pbe
     v1up_dwtil= v1dw_pbe
     v2up_uptil= v2up_pbe
     v2up_dwtil= v2up_pbe
  endif
!     ec_pbe(rhodw,0,grhodw,0)
!     zeta = 1.0_DP
  if(rhodw.gt.small) then
     v2_tmp=0.0_DP
     call pw_spin((pi34/rhodw)**third,-1.0_DP,ec_u,vcup_u,vcdw_u)
     if(sqrt(grhodw2).gt.small) then
        call pbec_spin(rhodw,-1.0_DP+small,grhodw2,1,&
             ecdw_0,v1up_0,v1dw_0,v2_tmp)
     else
        ecdw_0=0.0_DP
        v1dw_0=0.0_DP
        v2dw_0=0.0_DP
     endif
     ecdw_0 = ecdw_0/rhodw + ec_u
     v1dw_0 = (v1dw_0 + vcdw_u-ecdw_0)/rhodw
     v2dw_0 = v2_tmp/rhodw*grhodw
  else
     ecdw_0 = 0.0_DP
     v1dw_0 = 0.0_DP
     v2dw_0 = 0.0_DP
  endif
!     
  if(ecdw_0.gt.ec_pbe) then
     ecdw_til = ecdw_0
     v1dw_dwtil=v1dw_0
     v2dw_dwtil=v2dw_0
     v1dw_uptil=0.0_DP
     v2dw_uptil=0.0_DP
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
!       /(4.0_DP*fac*rho**(14.0_DP/3.0_DP))
  bb=0.0_DP
  do ipol=1,3
     term(ipol)= rhodw*grhoup(ipol)-rhoup*grhodw(ipol)
     bb=bb+ term(ipol)**2
  end do
!vector
  term=term/(fac*rho**(14.0_DP/3.0_DP))   
  bb=bb/(fac*rho**(14.0_DP/3.0_DP))
!      bb=(rhodw*grhoup-rhoup*grhodw)**2/fac*rho**(-14.0_DP/3.0_DP)
  aa2=aa*aa
  ca0=0.53_DP+aa2*(0.87_DP+aa2*(0.50_DP+aa2*2.26_DP))
  dca0da = aa*(1.74_DP+aa2*(2.0_DP+aa2*13.56_DP))
  if(abs(aa).le.1.0_DP-small) then
     term3 =(1.0_DP+aa)**(-p43) + (1.0_DP-aa)**(-p43)
     term1=(1.0_DP+bb*0.50_DP*term3)
     term2=(1.0_DP+aa)**(-7.0_DP/3.0_DP) + (1.0_DP-aa)**(-7.0_DP/3.0_DP)
     cab =ca0/term1**4
     dcabda = (dca0da/ca0 + 8.0_DP/3.0_DP*bb*term2/term1)*cab
     dcabdb = -2.0_DP*cab*term3/term1
  else
     cab=0.0_DP
     dcabda=0.0_DP
     dcabdb=0.0_DP
  endif
  da1up=2.0_DP*rhodw/rho**2
  da1dw=-2.0_DP*rhoup/rho**2
  db1up=-2.0_DP*(grhodw(1)*term(1)+grhodw(2)*term(2)+grhodw(3)*term(3)) &
       -14.0_DP/3.0_DP*bb/rho
  db1dw= 2.0_DP*(grhoup(1)*term(1)+grhoup(2)*term(2)+grhoup(3)*term(3)) &
       -14.0_DP/3.0_DP*bb/rho
  !vector, not scalar
  db2up= term*rhodw*2.0_DP
  db2dw=-term*rhoup*2.0_DP
!
  dcab1up = dcabda*da1up + dcabdb*db1up
  dcab1dw = dcabda*da1dw + dcabdb*db1dw
  !vector, not scalar
  dcab2up = dcabdb*db2up
  dcab2dw = dcabdb*db2dw
!cccccccccccccccccccccccccccccccccccccccccccccccccccccc------checked
  tauw=0.1250_DP*grho2/rho
  z=tauw/tau
  z2=z*z
!     
  term1=1.0_DP+cab*z2
  term2=(1.0_DP+cab)*z2
  ec_rev = ec_pbe*term1-term2*sum
!     
  drev1up=v1up_pbe*term1 + &
       ec_pbe*(z2*dcab1up - 2.0_DP*cab*z2/rho) &
       + (2.0_DP*term2/rho - z2*dcab1up)*sum &
       - term2*dsum1up 
!     
  drev1dw=v1dw_pbe*term1 + &
       ec_pbe*(z2*dcab1dw - 2.0_DP*cab*z2/rho)  &
       + (2.0_DP*term2/rho - z2*dcab1dw)*sum  &
       - term2*dsum1dw
!
! vector, not scalar
  drev2up=v2up_pbe*term1 + &
       ec_pbe*(z2*dcab2up+0.5_DP*cab*z/(rho*tau)*grhovec)& 
       - (term2*4.0_DP/grho2*grhovec + z2*dcab2up)*sum &
       - term2*dsum2up
  drev2dw=v2dw_pbe*term1 +  &
       ec_pbe*(z2*dcab2dw+0.5_DP*cab*z/(rho*tau)*grhovec) &
       - (term2*4.0_DP/grho2*grhovec + z2*dcab2dw)*sum  &
       - term2*dsum2dw
!
  drev3 = ((1.0_DP+cab)*sum-ec_pbe*cab)*2.0_DP*z2/tau
!ccccccccccccccccccccccccccccccccccccccccccccccccccc----checked
  term1=ec_rev*(1.0_DP+2.8_DP*ec_rev*z2*z)
  term2=(1.0_DP+5.6_DP*ec_rev*z2*z)*rho
  term3=-8.4_DP*ec_rev*ec_rev*z2*z
!
  v1up = term1 + term2*drev1up + term3
  v1dw = term1 + term2*drev1dw + term3
!
  term3=term3*rho
  v3   = term2*drev3 + term3/tau
!
  term3=-2.0_DP*term3/grho2    !grho/|grho|^2 = 1/grho
  v2up = term2*drev2up + term3*grhovec
  v2dw = term2*drev2dw + term3*grhovec
  !
  !
  !  call pw_spin((pi34/rho)**third,zeta,ec_u,vcup_u,vcdw_u)
  sc=rho*ec_rev*(1.0_DP+2.8_DP*ec_rev*z2*z) !-rho*ec_u
  !  v1up=v1up-vcup_u
  !  v1dw=v1dw-vcdw_u
  
  return
end subroutine metac_spin
!
!-------------------------------------------------------------------------
!
!                          END TPSSS
!-------------------------------------------------------------------------
!
!=========================================================================
!
!-------------------------------------------------------------------------
!
!                             M06L
!
!
!           input:  - rho
!                   - grho2=|\nabla rho|^2
!                   - tau = the kinetic energy density
!                           It is defined as summ_i( |nabla phi_i|**2 ) 
!
!           definition:  E_x = \int ex dr
!
!           output:     ex (rho, grho, tau)
!                       v1x= D(E_x)/D(rho)
!                       v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
!                            ( v2x = 1/|grho| * dsx / d|grho| = 2 *  dsx / dgrho2 )
!                       v3x= D(E_x)/D(tau)
!   
!                       ec, v1c, v2c, v3c as above for correlation
!
!-------------------------------------------------------------------------
!
subroutine m06lxc (rho, grho2, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c)

  !-----------------------------------------------------------------------
  !    
  !
  USE kinds,        ONLY : dp

  implicit none  

  real(dp), intent(in)  :: rho, grho2, tau
  real(dp), intent(out) :: ex, ec, v1x, v2x,v3x,v1c,v2c,v3c
  !
  real(dp)              :: rhoa, rhob, grho2a, grho2b, taua, taub, v1cb, v2cb, v3cb
  real(dp), parameter   :: zero = 0.0_dp, two = 2.0_dp, four = 4.0_dp
  !
  !
  rhoa   = rho   / two   ! one component only
  rhob   = rhoa
  !
  grho2a = grho2 / four
  grho2b = grho2a
  !
  taua   = tau * two * 0.5_dp ! Taua, which is Tau_sigma is half Tau
  taub   = taua               ! Tau is defined as summ_i( |nabla phi_i|**2 ) 
                              ! in the M06L routine
  !
  call m06lx (rhoa, grho2a, taua, ex, v1x, v2x, v3x)
  !
  ex  = two * ex  ! Add the two components up + dw
  !
  v2x = 0.5_dp * v2x
  !
  call m06lc (rhoa, rhob, grho2a, grho2b, taua, taub, ec, v1c, v2c, v3c, &
  &           v1cb, v2cb, v3cb)
  !
  !
  v2c = 0.5_dp * v2c
  !
end subroutine m06lxc

!-------------------------------------------------------------------------
!
subroutine m06lxc_spin (rhoup, rhodw, grhoup2, grhodw2, tauup, taudw,      &
          &             ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
          &             v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw)

  !-----------------------------------------------------------------------
  !    
  !
  USE kinds,        ONLY : dp

  implicit none  

  real(dp), intent(in)  :: rhoup, rhodw, grhoup2, grhodw2, tauup, taudw
  real(dp), intent(out) :: ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
   &                       v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw
  !
  real(dp)              :: exup, exdw, taua, taub
  real(dp), parameter   :: zero = 0.0_dp, two = 2.0_dp
  !
  !
  !
  taua   = tauup * two       ! Tau is defined as summ_i( |nabla phi_i|**2 ) 
  taub   = taudw * two       ! in the rest of the routine
  !
  call m06lx (rhoup, grhoup2, taua, exup, v1xup, v2xup, v3xup)
  call m06lx (rhodw, grhodw2, taub, exdw, v1xdw, v2xdw, v3xdw)
  !
  ex = exup + exdw
  !
  !
  call m06lc (rhoup, rhodw, grhoup2, grhodw2, taua, taub,     &
  &           ec, v1cup, v2cup, v3cup, v1cdw, v2cdw, v3cdw)
  !
  !
  !
end subroutine m06lxc_spin


!===============================  M06L exchange ==========================


subroutine m06lx (rho, grho2, tau, ex, v1x, v2x, v3x)

!_________________________________________________________________________
  
  use kinds, ONLY  : dp
  use constants, ONLY : pi
  
  implicit none

  real(dp), intent(in)   :: rho, grho2, tau 
  real(dp), intent(out)  :: ex, v1x, v2x, v3x
  
  real(dp)               :: v1x_unif,ex_unif, ex_pbe, &
  &                         sx_pbe, v1x_pbe, v2x_pbe
  !
  !  ex_unif:   lda \epsilon_x(rho)
  ! v2x = 1/|grho| * dsx / d|grho| = 2 *  dsx / dgrho2
  !
  real(dp), parameter     :: zero = 0._dp, one = 1.0_dp, two=2.0_dp, three = 3.0_dp,  &
  &                          four = 4.0_dp, five = 5.0_dp, six = 6.0_dp,              &
  &                          eight = 8.0_dp,                                          &
  &                          f12 = one/two, f13 = one/three, f23 = two/three,         &
  &                          f53 = five/three, f83 = eight/three, f43 = four/three,   &
  &                          pi34 = pi*three/four, pi2 = pi*pi,                       &
  &                          small=1.d-10
  
  real(dp)                :: d0, d1, d2, d3, d4, d5, CF, CT, CX, alpha
  real(dp), dimension(0:11)   &
  &                       :: at
  integer                 :: i
  !
  ! 
  ! VSXC98 variables (LDA part)
  !
  real(dp)                :: xs, xs2, grho, rhom83, rho13, rho43, zs, gh
  real(dp)                :: hg, dhg_dxs2, dhg_dzs
  real(dp)                :: dxs2_drho, dxs2_dgrho2, dzs_drho, dzs_dtau
  real(dp)                :: ex_vs98, v1x_vs98, v2x_vs98, v3x_vs98, v2x_vs98_g   
  !
  ! GGA and MGGA variables
  !  
  real(dp)                :: tau_unif, ts, ws, fws, dfws, dfws_drho, dfws_dtau,     &
  &                          dws_dts, dts_dtau, dts_drho
  !
  ! _________________________________________________________________________________________ 
    
  ! set parameters
  
  at(0)  =    3.987756d-01
  at(1)  =    2.548219d-01
  at(2)  =    3.923994d-01 
  at(3)  =   -2.103655d+00
  at(4)  =   -6.302147d+00
  at(5)  =    1.097615d+01
  at(6)  =    3.097273d+01
  at(7)  =   -2.318489d+01
  at(8)  =   -5.673480d+01
  at(9)  =    2.160364d+01
  at(10) =    3.421814d+01
  at(11) =   -9.049762d+00
  
  d0     =    6.012244d-01
  d1     =    4.748822d-03
  d2     =   -8.635108d-03
  d3     =   -9.308062d-06
  d4     =    4.482811d-05
  d5     =    zero
  
  alpha = 1.86726d-03
!___________________________________________________

  if (rho < small .and. tau < small) then
    ex = zero
       v1x = zero
       v2x = zero
       v3x = zero
       return
       
  end if
  
! _________VSXC98 functional (LDA part)_____________
!
! set variables 
  
  
  CF  =   (three/five) * (six*pi2)**f23
  CT  =   CF / two
  CX  =  -(three/two) * (three/(four*pi))**f13  ! Cx LSDA
  
!  if (rho >= small .and. grho>=small) then

  grho = sqrt(grho2) 
  rho43 = rho**f43
  rho13 = rho**f13
  rhom83 = one/rho**f83
  xs  = grho / rho43
  xs2 = xs * xs  
  zs  = tau/rho**f53 - CF
  gh  = one + alpha * (xs2 + zs)

   
  if (gh >= small) then
    call gvt4 (xs2, zs, d0, d1, d2, d3, d4, d5, alpha, hg, dhg_dxs2, dhg_dzs)
  else 
    hg        = zero
    dhg_dxs2  = zero
    dhg_dzs   = zero
  end if
  
  dxs2_drho   = -f83*xs2/rho
  dxs2_dgrho2 =  rhom83
  dzs_drho    = -f53*tau*rhom83
  dzs_dtau    =  one/rho**f53
  
  ex_unif     =  CX * rho43
  ex_vs98     =  ex_unif * hg
  v1x_vs98    =  CX * ( f43 * hg * rho**f13 ) +         &
&                ex_unif * ( dhg_dxs2*dxs2_drho + dhg_dzs*dzs_drho ) 
  v2x_vs98    =  two * ex_unif * dhg_dxs2 * dxs2_dgrho2
  v3x_vs98    =  ex_unif * dhg_dzs * dzs_dtau
  
!____________________mo6lx functional____________________________

  tau_unif  = CF * rho**f53  ! Tau is define as summ_i( |nabla phi_i|**2 )
  ts        = tau_unif / tau
  ws        = (ts - one)/(ts + one)
  
  fws       = zero
  dfws      = zero
  
  do i = 0, 11
    fws  =  fws + at(i)*ws**i
    dfws = dfws + i*at(i)*ws**(i-1)
  end do
  
  dws_dts   = two/((ts+1)**2)
  dts_drho  = ( (six*pi*pi*rho)**f23 )/tau
  dts_dtau  = -ts/tau
  dfws_drho = dfws*dws_dts*dts_drho
  dfws_dtau = dfws*dws_dts*dts_dtau

  call pbex_m06l (two*rho, four*grho2, sx_pbe, v1x_pbe, v2x_pbe)
  
  v1x_unif = f43 * CX * rho13
  
  sx_pbe   = f12 * sx_pbe 
  v1x_pbe  = v1x_pbe  + v1x_unif 
  v2x_pbe  = two * v2x_pbe
  

  ex_pbe   = sx_pbe + ex_unif 
  
!________energy and potential_____________________________  
   
  ex  = ex_vs98  + ex_pbe*fws  
                        
  v1x = v1x_vs98 + v1x_pbe*fws + ex_pbe*dfws_drho 
  v2x = v2x_vs98 + v2x_pbe*fws               
  v3X = v3x_vs98 + ex_pbe*dfws_dtau    
               
!__________________________________________________________

end subroutine m06lx

!__________________________________________________________

subroutine pbex_m06l (rho, grho2, sx, v1x, v2x)

  !---------------------------------------------------------------
  !
  ! PBE exchange (without Slater exchange):
  ! J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  !
  ! v2x = 1/|grho| * dsx / d|grho| = 2 *  dsx / dgrho2
  !
  USE kinds
  USE constants, ONLY : pi

  implicit none
  
  real(dp)      :: rho, grho2, sx, v1x, v2x
  ! input: charge and squared gradient
  ! output: energy
  ! output: potential
  integer       :: iflag
  
  ! local variables
  real(dp)      :: grho, rho43, xs, xs2, dxs2_drho, dxs2_dgrho2
  real(dp)      :: CX, denom, C1, C2, ex, Fx, dFx_dxs2, dex_drho
  
  real(dp), parameter ::     mu=0.21951_dp, ka=0.804_dp, one = 1.0_dp, two=2.0_dp, three = 3.0_dp,   &
  &                          four = 4.0_dp, six = 6.0_dp, eight = 8.0_dp,              &
  &                          f13 = one/three, f23 = two/three, f43 = four/three,       &
  &                          f34=three/four, f83 = eight/three
  
!_____________________________________________________________________

  CX    =  f34 * (three/pi)**f13            ! Cx LDA
  denom =  four * (three*pi**two)**f23
  C1    =  mu / denom
  C2    =  mu / (ka * denom)

  grho  = sqrt(grho2) 
  rho43 = rho**f43
  xs    = grho / rho43
  xs2   = xs * xs  
  
  dxs2_drho  = -f83 * xs2 / rho
  dxs2_dgrho2 = one /rho**f83
  
  ex        = - CX * rho43
  dex_drho  = - f43 * CX * rho**f13
  
  Fx        =   C1*xs2 / (one + C2*xs2)
  dFx_dxs2  =   C1 / (one + C2*xs2)**2
  !
  !   Energy
  !
  sx = Fx * ex
  !
  !   Potential
  !
  v1x = dex_drho * Fx  +  ex * dFx_dxs2 * dxs2_drho
  v2x = two * ex * dFx_dxs2* dxs2_dgrho2
  !
  !
end subroutine pbex_m06l


!===============================  M06L correlation ==========================
!
!-------------------------------------------------------------------------
!
subroutine m06lc (rhoa, rhob, grho2a, grho2b, taua, taub, ec, v1c_up, v2c_up, v3c_up,   &
&                 v1c_dw, v2c_dw, v3c_dw)
!-------------------------------------------------------------------------
!
  
  use kinds, only     : dp
  use constants, only : pi
  implicit none
!-------------------------------------------------------------------------
  
  real(dp), intent(in)    :: rhoa, rhob, grho2a, grho2b, taua, taub
  real(dp), intent(out)   :: ec, v1c_up, v2c_up, v3c_up, v1c_dw, v2c_dw, v3c_dw
  !
  real(dp), parameter     :: zero = 0._dp, one = 1.0_dp, two=2.0_dp, three = 3.0_dp,    &
  &                          four = 4.0_dp, five = 5.0_dp, six = 6.0_dp,                &
  &                          eight = 8.0_dp,                                            &
  &                          f12 = one/two, f13 = one/three, f23 = two/three,           &
  &                          f53 = five/three, f83 = eight/three, f43 = four/three,     &
  &                          pi34 = three/(four*pi), pi2 = pi*pi, f35 = three/five,     &
  &                          small=1.d-10
  !
  ! parameters of the MO6Lc functional
  !
  real(dp), dimension(0:4):: cs, cab
  ! 
  real(dp)                :: ds0, ds1, ds2, ds3, ds4, ds5, CF, alpha, Ds,                &
  &                          dab0, dab1, dab2, dab3, dab4, dab5, gama_ab, gama_s,        &
  &                          alpha_s, alpha_ab 
  !
  ! functions and variables
  !
  real(dp)                :: ec_pw_a, ec_pw_b, ec_pw_ab, vc_pw_a, vc_pw_b, vv,           &
  &                          vc_pw_ab, vc_pw_up, vc_pw_dw, Ecaa, Ecbb, Ecab,             &
  &                          Ec_UEG_ab, Ec_UEG_aa, Ec_UEG_bb, decab_drhoa, decab_drhob,  &
  &                          v1_ab_up, v1_ab_dw, v2_ab_up, v2_ab_dw, v3_ab_up, v3_ab_dw, &
  &                          v1_aa_up, v2_aa_up, v3_aa_up, v1_bb_dw, v2_bb_dw, v3_bb_dw            
  !  
  real(dp)                :: xsa, xs2a, rsa, grhoa, xsb, xs2b, grhob, rsb, zsa, zsb,     &
  &                          xs2ab, zsab, zeta, rho, rs,                                 &
  &                          dxs2a_drhoa, dxs2b_drhob, dxs2a_dgrhoa2, dxs2b_dgrhob2,     &
  &                          dzsa_drhoa, dzsb_drhob, dzsa_dtaua, dzsb_dtaub
  !
  real(dp)                :: hga, dhga_dxs2a, dhga_dzsa, hgb, dhgb_dxs2b, dhgb_dzsb,     &
  &                          hgab, dhgab_dxs2ab, dhgab_dzsab,                            &
  &                          Dsa, Dsb, dDsa_dxs2a, dDsa_dzsa, dDsb_dxs2b, dDsb_dzsb,     &
  &                          gsa, gsb, gsab, dgsa_dxs2a, dgsb_dxs2b, dgsab_dxs2ab, num

  integer                 :: ifunc

!_____________________________________________________________________________________
 
  dab0     =   3.957626d-01
  dab1     =  -5.614546d-01
  dab2     =   1.403963d-02
  dab3     =   9.831442d-04
  dab4     =  -3.577176d-03
  dab5     =   zero

  cab(0)   =   6.042374d-01
  cab(1)   =   1.776783d+02
  cab(2)   =  -2.513252d+02
  cab(3)   =   7.635173d+01
  cab(4)   =  -1.255699d+01   

  gama_ab  =   0.0031_dp
  alpha_ab =   0.00304966_dp

  
  ds0      =   4.650534d-01
  ds1      =   1.617589d-01
  ds2      =   1.833657d-01
  ds3      =   4.692100d-04
  ds4      =  -4.990573d-03
  ds5      =   zero

  cs(0)    =   5.349466d-01
  cs(1)    =   5.396620d-01
  cs(2)    =  -3.161217d+01
  cs(3)    =   5.149592d+01
  cs(4)    =  -2.919613d+01

  gama_s   =   0.06_dp
  alpha_s  =   0.00515088_dp
  
  CF  =   f35 * (six*pi2)**f23
  
  ifunc = 1     ! iflag=1  J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
 
!______________Ecaa_____________________________________________________
   
  if (rhoa < small .and. taua < small ) then
  
    Ecaa     = zero
    v1_aa_up = zero
    v2_aa_up = zero
    v3_aa_up = zero
    
  else
  
    rsa   = (pi34/rhoa)**f13
    grhoa = sqrt(grho2a)
    xsa   = grhoa / rhoa**f43
    xs2a  = xsa * xsa
    zsa   = taua/rhoa**f53 - CF
    
    dxs2a_drhoa   = -f83*xs2a/rhoa
    dxs2a_dgrhoa2 =  one/(rhoa**f83)
  
    dzsa_drhoa   = -f53*taua/(rhoa**f83)
    dzsa_dtaua   =  one/rhoa**f53
    
    Dsa        = one - xs2a/(four * (zsa + CF))
    dDsa_dxs2a = - one/(four * (zsa + CF))
    dDsa_dzsa  = xs2a/(four * (zsa + CF)**2)
  
    ec_pw_a = zero
    vc_pw_a = zero
  
    call pw_spin (rsa, one, ec_pw_a, vc_pw_a, vv)
    call gvt4  (xs2a, zsa, ds0, ds1, ds2, ds3, ds4, ds5, alpha_s, hga, dhga_dxs2a, dhga_dzsa)
    call gfunc (cs, gama_s, xs2a, gsa, dgsa_dxs2a)
  
    Ec_UEG_aa  = rhoa*ec_pw_a
    num        = (dgsa_dxs2a + dhga_dxs2a)*Dsa + (gsa + hga)*dDsa_dxs2a 
    !
    !
    Ecaa       = Ec_UEG_aa * (gsa + hga) * Dsa
   
    v1_aa_up   = vc_pw_a * (gsa + hga) * Dsa                            &
    &          + Ec_UEG_aa * num * dxs2a_drhoa                          &
    &          + Ec_UEG_aa * (dhga_dzsa*Dsa + (gsa + hga)*dDsa_dzsa) * dzsa_drhoa

    v2_aa_up   = two * Ec_UEG_aa * num * dxs2a_dgrhoa2

    v3_aa_up   = Ec_UEG_aa * (dhga_dzsa*Dsa + (gsa + hga)*dDsa_dzsa) * dzsa_dtaua
    !
  end if
!
!______________Ecbb_____________________________________________________
  
  if (rhob < small .and. taub < small) then
  
    Ecbb     = zero
    v1_bb_dw = zero
    v2_bb_dw = zero
    v3_bb_dw = zero
    
  else
  
    rsb   = (pi34/rhob)**f13
    grhob = sqrt(grho2b)
    xsb   = grhob / rhob**f43
    xs2b  = xsb * xsb
    zsb   = taub/rhob**f53 - CF
  
    dxs2b_drhob   = -f83*xs2b/rhob
    dxs2b_dgrhob2 =  one /rhob**f83
  
    dzsb_drhob   = -f53*taub/(rhob**f83)
    dzsb_dtaub   =  one/rhob**f53
  
    Dsb        = one - xs2b/(four * (zsb + CF))
    dDsb_dxs2b = - one/(four * (zsb + CF))
    dDsb_dzsb  =  xs2b/(four * (zsb + CF)**2)
  
    call pw_spin (rsb, one, ec_pw_b, vc_pw_b, vv)
    call gvt4  (xs2b, zsb, ds0, ds1, ds2, ds3, ds4, ds5, alpha_s, hgb, dhgb_dxs2b, dhgb_dzsb)
    call gfunc (cs, gama_s, xs2b, gsb, dgsb_dxs2b)
  
    Ec_UEG_bb  = rhob*ec_pw_b
    num        = (dgsb_dxs2b + dhgb_dxs2b)*Dsb + (gsb + hgb)*dDsb_dxs2b
    !
    !  
    Ecbb       = Ec_UEG_bb * (gsb + hgb) * Dsb

    v1_bb_dw   = vc_pw_b * (gsb + hgb) * Dsb                             &
    &          + Ec_UEG_bb * num * dxs2b_drhob                           &
    &          + Ec_UEG_bb * (dhgb_dzsb*Dsb + (gsb + hgb)*dDsb_dzsb)*dzsb_drhob
  
    v2_bb_dw   = two * Ec_UEG_bb * num * dxs2b_dgrhob2
  
    v3_bb_dw   = Ec_UEG_bb * (dhgb_dzsb*Dsb + (gsb + hgb)*dDsb_dzsb)*dzsb_dtaub
    !
  end if
!
!________________Ecab____________________________________________

  if (rhoa < small .and. rhob < small) then
  
    Ecab     = zero
    v1_ab_up = zero
    v1_ab_dw = zero
    v2_ab_up = zero
    v2_ab_dw = zero
    v3_ab_up = zero
    v3_ab_dw = zero

  else

    xs2ab = xs2a + xs2b
    zsab  = zsa + zsb
    rho   = rhoa + rhob
    zeta  = (rhoa - rhob)/rho
    rs    = (pi34/rho)**f13
  
    call gvt4 (xs2ab, zsab, dab0, dab1, dab2, dab3, dab4, dab5, alpha_ab, hgab, dhgab_dxs2ab, dhgab_dzsab)
    call pw_spin (rs, zeta, ec_pw_ab, vc_pw_up, vc_pw_dw)
    call gfunc (cab, gama_ab, xs2ab, gsab, dgsab_dxs2ab)
    
    decab_drhoa  =  vc_pw_up - vc_pw_a
    decab_drhob  =  vc_pw_dw - vc_pw_b
  
    Ec_UEG_ab    = ec_pw_ab*rho - ec_pw_a*rhoa - ec_pw_b*rhob
    !
    !
    Ecab         = Ec_UEG_ab * (gsab + hgab)
 
    v1_ab_up = decab_drhoa * (gsab + hgab)                               &
    &        + Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2a_drhoa   &
    &        + Ec_UEG_ab * dhgab_dzsab * dzsa_drhoa
  
    v1_ab_dw = decab_drhob * (gsab + hgab)                               &
    &        + Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2b_drhob   &
    &        + Ec_UEG_ab * dhgab_dzsab * dzsb_drhob
  
    v2_ab_up = two * Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2a_dgrhoa2
    v2_ab_dw = two * Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2b_dgrhob2
  
    v3_ab_up = Ec_UEG_ab * dhgab_dzsab * dzsa_dtaua
    v3_ab_dw = Ec_UEG_ab * dhgab_dzsab * dzsb_dtaub
    !
  end if
!
!___________________ec and vc_____________________________________________

  ec     = Ecaa + Ecbb + Ecab 
 
  v1c_up = v1_aa_up + v1_ab_up
  v2c_up = v2_aa_up + v2_ab_up
  v3c_up = v3_aa_up + v3_ab_up
  
  v1c_dw = v1_bb_dw + v1_ab_dw
  v2c_dw = v2_bb_dw + v2_ab_dw
  v3c_dw = v3_bb_dw + v3_ab_dw
  
!__________________________________________________________________________
  
  contains
!__________________________________________________________________________  
  
  subroutine gfunc (cspin, gama, xspin, gs, dgs_dx)
  
    implicit none
    real(dp), dimension (0:4), intent(in) :: cspin
    real(dp), intent(in)                  :: xspin, gama
    real(dp), intent(out)                 :: gs, dgs_dx
    !
    real(dp)               :: de, d2, x1, x2, x3, x4
    real(dp), parameter    :: one=1.0d0, two=2.0d0, three=3.0d0, four=4.0d0
    
    !__________________
    
    de  = one/(one + gama*xspin)
    d2 = de**2
    x1 = gama * xspin * de
    x2 = x1**2
    x3 = x1**3
    x4 = x1**4
    
    gs     = cspin(0) + cspin(1)*x1 + cspin(2)*x2 + cspin(3)*x3 + cspin(4)*x4
    dgs_dx = gama*d2* (cspin(1) + two*cspin(2)*x1 + three*cspin(3)*x2 + four*cspin(4)*x3)
       
  end subroutine gfunc

 !___________________________________________________________________
  
end subroutine m06lc

!___________________________________________________________________


subroutine gvt4 (x, z, a, b, c, d, e, f, alpha, hg, dh_dx, dh_dz)
  
  use kinds, only  : dp
  
  implicit none
  
  real(dp), intent(in)   :: X, z, a, b, c, d, e, f, alpha
  real(dp), intent(out) :: hg, dh_dx, dh_dz
  
  real(dp)              :: gamma, gamma2, gamma3
  real(dp), parameter   :: one=1.0_dp, two=2.0_dp, three=3.0_dp
  
  gamma  = one + alpha*(x+z)
  gamma2 = gamma*gamma
  gamma3 = gamma2*gamma
  
  hg = a/gamma + (b*x + c*z)/gamma2 + (d*x*x + e*x*z + f*z*z)/gamma3
  
  dh_dx = ( -alpha*a + b + (two*x*(d - alpha*b) + z*(e - two*alpha*c))/ gamma  &
&         - three*alpha*(d*x*x + e*x*z + f*z*z)/gamma2  )/gamma2

  dh_dz = ( -alpha*a + c + (two*z*(f - alpha*c) + x*(e -two*alpha*b))/ gamma  &
&         - three*alpha*(d*x*x + e*x*z + f*z*z)/gamma2  )/gamma2
  
  return
  
end subroutine gvt4

!-------------------------------------------------------------------------
!
!                          END M06L
!
!=========================================================================
!
!                          TB09
!-----------------------------------------------------------------------

subroutine tb09cxc(rho, grho, tau, sx, sc, v1x, v2x,v3x,v1c, v2c,v3c)

  USE kinds,            ONLY : DP
#if defined(__LIBXC)
  use xc_f90_types_m
  use xc_f90_lib_m
#endif
  implicit none  
  real(DP), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c
#if defined(__LIBXC)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  integer :: size = 1
  integer :: func_id = 208  !Tran & Blaha correction to Becke & Johnson
  real(dp) :: lapl_rho, vlapl_rho ! not used in TB09
  
  lapl_rho = grho
  
  ! ---------------------------- Exchange
  
  func_id = 208  !Tran & Blaha correction to Becke & Johnson  -- TB09 
 
  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)    
  call xc_f90_mgga_vxc(xc_func, size, rho, grho, lapl_rho, 0.5_dp*tau, &
      v1x, v2x, vlapl_rho, v3x)  
  call xc_f90_func_end(xc_func)	
  
  sx = 0.0d0		       
  v2x = v2x*2.0_dp
  v3x = v3x*0.5_dp

  ! ---------------------------- Correlation  
  
  func_id = 231  ! Perdew, Tao, Staroverov & Scuseria correlation  -- TPSS
  
  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)   
  call xc_f90_mgga_exc_vxc(xc_func, size, rho, grho, lapl_rho, 0.5_dp*tau, &
     sc, v1c, v2c, vlapl_rho, v3c)  
  call xc_f90_func_end(xc_func)

  sc = sc * rho			       
  v2c = v2c*2.0_dp
  v3c = v3c*0.5_dp
#else
  call errore('tb09','need libxc',1)
#endif
end subroutine tb09cxc

!c     ==================================================================




!======================================================================
! SCAN meta-GGA
!======================================================================
subroutine SCANcxc(rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c)
  USE kinds,            ONLY : DP
  USE funct,            ONLY : libxc_major, libxc_minor, libxc_micro, get_libxc_version
#if defined(__LIBXC)
  use xc_f90_types_m
  use xc_f90_lib_m
#endif
  implicit none  
  real(DP), intent(in) :: rho, grho, tau
  real(dp), intent(out):: sx, sc, v1x, v2x, v3x, v1c, v2c, v3c

#if defined(__LIBXC)
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  integer, SAVE :: major=0, minor=0, micro=0
  integer :: size = 1
  integer :: func_id
  real(dp) :: lapl_rho, vlapl_rho ! not used?

  if (libxc_major == 0) call get_libxc_version
  if (libxc_major < 3 .or. (libxc_major == 3 .and. libxc_minor /= -1)) & 
      call errore('SCAN meta-GGA','please, recompile with LibXC trunk (i.e. >3.0.0))',1)
 
  lapl_rho = grho

  ! exchange
  func_id = 263  ! XC_MGGA_X_SCAN
 
  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)
  call xc_f90_mgga_exc_vxc(xc_func, size, rho, grho, lapl_rho, 0.5_dp*tau,&
                           sx, v1x, v2x, vlapl_rho, v3x)  
  call xc_f90_func_end(xc_func)
  sx = 0.0d0
  v2x = v2x*2.0_dp
  v3x = v3x*0.5_dp

  ! correlation  
  func_id = 267  ! XC_MGGA_C_SCAN
  call xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)   
  call xc_f90_mgga_exc_vxc(xc_func, size, rho, grho, lapl_rho, 0.5_dp*tau, &
                           sc, v1c, v2c, vlapl_rho, v3c)  
  call xc_f90_func_end(xc_func)

  sc = sc * rho
  v2c = v2c*2.0_dp
  v3c = v3c*0.5_dp

#else
  call errore('SCAN meta-GGA','please, recompile with LibXC trunk (i.e. >3.0.0))',1)
#endif

end subroutine SCANcxc

