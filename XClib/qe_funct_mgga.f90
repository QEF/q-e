!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
MODULE metagga
!-------------------------------------------------------------------------
!! MetaGGA functionals. Available functionals:
!
!! * TPSS (Tao, Perdew, Staroverov & Scuseria)
!! * M06L
!
USE exch_lda, ONLY: slater
USE corr_lda, ONLY: pw, pw_spin   
USE corr_gga, ONLY: pbec, pbec_spin
!
 CONTAINS
!                             TPSS
!
!-------------------------------------------------------------------------
SUBROUTINE tpsscxc( rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c )
  !-----------------------------------------------------------------------
  !! TPSS metaGGA corrections for exchange and correlation - Hartree a.u.  
  !! Definition:  \(E_x = \int E_x(\text{rho},\text{grho}) dr\)
  !
  USE kind_l,            ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the charge density
  REAL(DP), INTENT(IN) :: grho
  !! grho = |\nabla rho|^2
  REAL(DP), INTENT(IN) :: tau
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: sx
  !! sx = E_x(rho,grho)
  REAL(DP), INTENT(OUT) :: sc
  !! sc = E_c(rho,grho)
  REAL(DP), INTENT(OUT) :: v1x
  !! v1x = D(E_x)/D(rho)
  REAL(DP), INTENT(OUT) :: v2x
  !! v2x = D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3x
  !! v3x = D(E_x)/D(tau)
  REAL(DP), INTENT(OUT) :: v1c
  !! v1c = D(E_c)/D(rho)
  REAL(DP), INTENT(OUT) :: v2c
  !! v2c = D(E_c)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3c
  !! v3c = D(E_c)/D(tau)
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  IF (rho <= small) THEN
     sx  = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
     sc  = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
     v3x = 0.0_DP
     v3c = 0.0_DP
     RETURN
  ENDIF
  !
  ! exchange
  CALL metax( rho, grho, tau, sx, v1x, v2x, v3x )
  ! correlation
  CALL metac( rho, grho, tau, sc, v1c, v2c, v3c )
  !
  RETURN
  !
END SUBROUTINE tpsscxc
!
!
!-------------------------------------------------------------------------
SUBROUTINE metax( rho, grho2, tau, ex, v1x, v2x, v3x )
  !--------------------------------------------------------------------
  !! TPSS meta-GGA exchange potential and energy.
  !
  !! NOTE: E_x(rho,grho) = rho\epsilon_x(rho,grho) ;
  !! ex = E_x(rho,grho)    NOT \epsilon_x(rho,grho) ;
  !! v1x= D(E_x)/D(rho) ;
  !! v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho| ;
  !! v3x= D(E_x)/D( tau ) ;
  !! tau is the kinetic energy density ;
  !! the same applies to correlation terms ;
  !! input grho2 is |\nabla rho|^2 .
  !
  USE kind_l,       ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the charge density
  REAL(DP), INTENT(IN) :: grho2
  !! grho2 = |\nabla rho|^2
  REAL(DP), INTENT(IN) :: tau
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: ex
  !! ex = E_x(rho,grho)
  REAL(DP), INTENT(OUT) :: v1x
  !! v1x = D(E_x)/D(rho)
  REAL(DP), INTENT(OUT) :: v2x
  !! v2x = D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3x
  !! v3x = D(E_x)/D(tau)
  !
  ! ... local variables
  !
  REAL(DP) :: rs, vx_unif, ex_unif
  !  ex_unif:   lda \epsilon_x(rho)
  !  ec_unif:   lda \epsilon_c(rho)
  REAL(DP), PARAMETER :: small=1.E-10_DP
  REAL(DP), PARAMETER :: pi34=0.6203504908994_DP, third=1.0_DP/3.0_DP
  !  fx=Fx(p,z)
  !  fxp=d Fx / d p
  !  fxz=d Fx / d z
  REAL(DP) :: fx, f1x, f2x, f3x
  !
  IF (ABS(tau) < small) THEN
    ex  = 0.0_DP
    v1x = 0.0_DP
    v2x = 0.0_DP
    v3x = 0.0_DP
    RETURN
  ENDIF
  !
  rs = pi34/rho**third
  CALL slater( rs, ex_unif, vx_unif )
  CALL metaFX( rho, grho2, tau, fx, f1x, f2x, f3x )
  !
  ex = rho*ex_unif
  v1x = vx_unif*fx + ex*f1x
  v2x = ex*f2x
  v3x = ex*f3x
  ex  = ex*fx
  !
  RETURN
  !
END SUBROUTINE metax
!
!
!------------------------------------------------------------------
SUBROUTINE metac( rho, grho2, tau, ec, v1c, v2c, v3c )
  !--------------------------------------------------------------
  !! TPSS meta-GGA correlation energy and potentials.
  !
  USE kind_l,  ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the charge density
  REAL(DP), INTENT(IN) :: grho2
  !! grho2 = |\nabla rho|^2
  REAL(DP), INTENT(IN) :: tau
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: ec
  !! ec = E_c(rho,grho)
  REAL(DP), INTENT(OUT) :: v1c
  !! v1c = D(E_c)/D(rho)
  REAL(DP), INTENT(OUT) :: v2c
  !! v2c = D(E_c)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3c
  !! v3c = D(E_c)/D(tau)
  !
  ! ... local variables
  !
  REAL(DP) :: z, z2, tauw, ec_rev
  REAL(DP) :: d1rev, d2rev, d3rev
  !  d1ec=  D ec_rev / D rho
  !  d2ec=  D ec_rev / D |D rho/ D r| / |\nabla rho|
  !  d3ec=  D ec_rev / D tau
  REAL(DP) :: cf1, cf2, cf3
  REAL(DP) :: v1c_pbe, v2c_pbe, ec_pbe
  REAL(DP) :: v1c_sum(2), v2c_sum, ec_sum
  !
  REAL(DP) :: rs
  REAL(DP) :: ec_unif, vc_unif
  REAL(DP) :: vc_unif_s(2)
  !
  REAL(DP) :: rhoup, grhoup
  !
  REAL(DP), PARAMETER :: small=1.0E-10_DP
  REAL(DP), PARAMETER :: pi34=0.75_DP/3.141592653589793_DP, &
                         third=1.0_DP/3.0_DP
  REAL(DP), PARAMETER :: dd=2.80_DP  !in unit of Hartree^-1
  REAL(DP), PARAMETER :: cab=0.53_DP, cabone=1.0_DP+cab
  !
  IF (ABS(tau) < small) THEN
     ec  = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
     v3c = 0.0_DP
     RETURN
  ENDIF
  !
  rhoup = 0.5_DP*rho
  grhoup = 0.5_DP*SQRT(grho2)
  !
  IF (rhoup > small) THEN
     !
     rs = (pi34/rhoup)**third
     CALL pw_spin( rs, 1.0_DP, ec_unif, vc_unif_s(1), vc_unif_s(2) )
     !
     IF (ABS(grhoup) > small) THEN
        ! zeta=1.0_DP-small to avoid pow_e of 0 in pbec_spin
        CALL pbec_spin( rhoup, 1.0_DP-small, grhoup**2, 1, &
                        ec_sum, v1c_sum(1), v1c_sum(2), v2c_sum )
     ELSE
        ec_sum  = 0.0_DP
        v1c_sum = 0.0_DP
        v2c_sum = 0.0_DP
     ENDIF
     ec_sum = ec_sum/rhoup + ec_unif
     v1c_sum(1) = (v1c_sum(1) + vc_unif_s(1)-ec_sum)/rho !rho, not rhoup
     v2c_sum = v2c_sum/(2.0_DP*rho)
  ELSE
     ec_sum  = 0.0_DP
     v1c_sum = 0.0_DP
     v2c_sum = 0.0_DP
  ENDIF
  !
  rs = (pi34/rho)**third
  CALL pw( rs, 1, ec_unif, vc_unif )
  !
  !  ... PBE correlation energy and potential:
  !  ec_pbe=rho*H,  not rho*(epsion_c_uinf + H)
  !  v1c_pbe=D (rho*H) /D rho
  !  v2c_pbe= for rho, 2 for
  CALL pbec( rho, grho2, 1, ec_pbe, v1c_pbe, v2c_pbe )
  !
  ec_pbe  = ec_pbe/rho+ec_unif
  v1c_pbe = (v1c_pbe+vc_unif-ec_pbe)/rho
  v2c_pbe = v2c_pbe/rho
  !
  IF (ec_sum < ec_pbe) THEN
     ec_sum = ec_pbe
     v1c_sum(1) = v1c_pbe
     v2c_sum = v2c_pbe
  ENDIF
  !
  tauw = 0.1250_DP * grho2/rho
  z = tauw/tau
  z2 = z*z
  !
  ec_rev = ec_pbe*(1+cab*z2)-cabone*z2*ec_sum
  !
  d1rev  = v1c_pbe + (cab*v1c_pbe-cabone*v1c_sum(1))*z2  &
                   - (ec_pbe*cab - ec_sum*cabone)*2.0_DP*z2/rho
  d2rev  = v2c_pbe + (cab*v2c_pbe-cabone*v2c_sum)*z2     &
                   + (ec_pbe*cab - ec_sum*cabone)*4.0_DP*z2/grho2
  d3rev  = -(ec_pbe*cab - ec_sum*cabone)*2.0_DP*z2/tau
  !
  cf1 = 1.0_DP + dd*ec_rev*z2*z
  cf2 = rho*(1.0_DP+2.0_DP*z2*z*dd*ec_rev)
  cf3 = ec_rev*ec_rev*3.0_DP*dd*z2*z
  v1c = ec_rev*cf1 + cf2*d1rev-cf3
  !
  cf3 = cf3*rho
  v2c = cf2*d2rev + cf3*2.0_DP/grho2
  v3c = cf2*d3rev - cf3/tau
  ec = rho*ec_rev*(1.0_DP+dd*ec_rev*z2*z)  !-rho*ec_unif(1)
  !
  RETURN
  !
END SUBROUTINE metac
!
!
!-------------------------------------------------------------------------
SUBROUTINE metaFX( rho, grho2, tau, fx, f1x, f2x, f3x )
  !-------------------------------------------------------------------------
  !! FX calculation.
  !
  USE kind_l,           ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! charge density
  REAL(DP), INTENT(IN) :: grho2
  !! square of gradient of rho
  REAL(DP), INTENT(IN) :: tau
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: fx
  !! fx = Fx(p,z)
  REAL(DP), INTENT(OUT) :: f1x
  !! f1x=D (Fx) / D rho
  REAL(DP), INTENT(OUT) :: f2x
  !! f2x=D (Fx) / D ( D rho/D r_alpha) /|nabla rho|
  REAL(DP), INTENT(OUT) :: f3x
  !! f3x=D (Fx) / D tau
  !
  ! ... local variables
  !
  REAL(DP) x, p, z, qb, al, localdp, dz
  REAL(DP) dfdx, dxdp, dxdz, dqbdp, daldp, dqbdz, daldz
  REAL(DP) fxp, fxz  ! fxp =D fx /D p
  REAL(DP) tauw, tau_unif
  !  work variables
  REAL(DP) xf1, xf2
  REAL(DP) xfac1, xfac2, xfac3, xfac4, xfac5, xfac6, xfac7, z2
  !
  REAL(DP), PARAMETER :: pi=3.141592653589793_DP
  REAL(DP), PARAMETER :: THRD=0.3333333333333333_DP
  REAL(DP), PARAMETER :: ee=1.537_DP
  REAL(DP), PARAMETER :: cc=1.59096_DP
  REAL(DP), PARAMETER :: kk=0.804_DP
  REAL(DP), PARAMETER :: bb=0.40_DP
  REAL(DP), PARAMETER :: miu=0.21951_DP
  REAL(DP), PARAMETER :: fac1=9.57078000062731_DP !fac1=(3*pi^2)^(2/3)
  REAL(DP), PARAMETER :: small=1.0E-6_DP
  !
  tauw = 0.125_DP*grho2/rho
  z = tauw/tau
  !
  p = SQRT(grho2)/rho**THRD/rho
  p = p*p/(fac1*4.0_DP)
  tau_unif = 0.3_DP*fac1*rho**(5.0_DP/3.0_DP)
  al = (tau-tauw)/tau_unif
  al = ABS(al)  !make sure al is always .gt. 0.0_DP
  qb = 0.45_DP*(al-1.0_DP)/SQRT(1.0_DP+bb*al*(al-1.0_DP))
  qb = qb+2.0_DP*THRD*p
  !
  !  calculate x(p,z) and fx
  z2 = z*z
  xf1 = 10.0_DP/81.0_DP
  xfac1 = xf1+cc*z2/(1+z2)**2.0_DP
  xfac2 = 146.0_DP/2025.0_DP
  xfac3 = SQRT(0.5_DP*(0.36_DP*z2+p*p))
  xfac4 = xf1*xf1/kk
  xfac5 = 2.0_DP*SQRT(ee)*xf1*0.36_DP
  xfac6 = xfac1*p+xfac2*qb**2.0_DP-73.0_DP/405.0_DP*qb*xfac3
  xfac6 = xfac6+xfac4*p**2.0_DP+xfac5*z2+ee*miu*p**3.0_DP
  xfac7 = (1+SQRT(ee)*p)
  x = xfac6/(xfac7*xfac7)
  !  fx=kk-kk/(1.0_DP+x/kk)
  fx = 1.0_DP + kk-kk/(1.0_DP+x/kk)
  !
  !  calculate the derivatives of fx w.r.t p and z
  dfdx = (kk/(kk+x))**2.0_DP
  daldp = 5.0_DP*THRD*(tau/tauw-1.0_DP)
  !
  !   daldz=-0.50_DP*THRD*
  !  * (tau/(2.0_DP*fac1*rho**THRD*0.1250_DP*sqrt(grho2)))**2.0_DP
  daldz = -5.0_DP*THRD*p/z2
  dqbdz = 0.45_DP*(0.50_DP*bb*(al-1.0_DP)+1.0_DP)
  dqbdz = dqbdz/(1.0_DP+bb*al*(al-1.0_DP))**1.5_DP
  !
  dqbdp = dqbdz*daldp+2.0_DP*THRD
  dqbdz = dqbdz*daldz
  !
  !  calculate dx/dp
  xf1 = 73.0_DP/405.0_DP/xfac3*0.50_DP*qb
  xf2 = 2.0_DP*xfac2*qb-73.0_DP/405.0_DP*xfac3
  !
  dxdp = -xf1*p
  dxdp = dxdp+xfac1+xf2*dqbdp
  dxdp = dxdp+2.0_DP*xfac4*p
  dxdp = dxdp+3.0_DP*ee*miu*p*p
  dxdp = dxdp/(xfac7*xfac7)-2.0_DP*x*SQRT(ee)/xfac7
  !
  !  dx/dz
  dxdz = -xf1*0.36_DP*z
  xfac1 = cc*2.0_DP*z*(1-z2)/(1+z2)**3.0_DP
  dxdz = dxdz+xfac1*p+xf2*dqbdz
  dxdz = dxdz+xfac5*2.0_DP*z
  dxdz = dxdz/(xfac7*xfac7)
  !
  fxp = dfdx*dxdp
  fxz = dfdx*dxdz
  !
  !  calculate f1x
  localdp = -8.0_DP*THRD*p/rho  ! D p /D rho
  dz = -z/rho                   ! D z /D rho
  f1x = fxp*localdp+fxz*dz
  !
  !  f2x
  localdp = 2.0_DP/(fac1*4.0_DP*rho**(8.0_DP/3.0_DP))
  dz = 2.0_DP*0.125_DP/(rho*tau)
  f2x = fxp*localdp + fxz*dz
  !
  !  f3x
  localdp = 0.0_DP
  dz = -z/tau
  f3x = fxz*dz
  !
  RETURN
  !
END SUBROUTINE metaFX
!
!
!---------------------------------------------------------------------------
SUBROUTINE tpsscx_spin( rhoup, rhodw, grhoup2, grhodw2, tauup, taudw, sx, &
                        v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw )
  !-----------------------------------------------------------------------
  !! TPSS metaGGA for exchange - Hartree a.u.
  !
  USE kind_l,            ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rhoup
  !! up charge
  REAL(DP), INTENT(IN) :: rhodw
  !! down charge
  REAL(DP), INTENT(IN) :: grhoup2
  !! up gradient of the charge
  REAL(DP), INTENT(IN) :: grhodw2
  !! down gradient of the charge
  REAL(DP), INTENT(IN) :: tauup
  !! up kinetic energy density
  REAL(DP), INTENT(IN) :: taudw
  !! down kinetic energy density
  REAL(DP), INTENT(OUT) :: sx
  !! exchange energy
  REAL(DP), INTENT(OUT) :: v1xup
  !! derivatives of exchange wr. rho - up
  REAL(DP), INTENT(OUT) :: v1xdw
  !! derivatives of exchange wr. rho - down
  REAL(DP), INTENT(OUT) :: v2xup
  !! derivatives of exchange wr. grho - up
  REAL(DP), INTENT(OUT) :: v2xdw
  !! derivatives of exchange wr. grho - down
  REAL(DP), INTENT(OUT) :: v3xup
  !! derivatives of exchange wr. tau - up
  REAL(DP), INTENT(OUT) :: v3xdw
  !! derivatives of exchange wr. tau - down
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  REAL(DP) :: rho, sxup, sxdw
  !
  ! exchange
  rho = rhoup + rhodw
  IF (rhoup>small .AND. SQRT(ABS(grhoup2))>small &
       .AND. ABS(tauup) > small) THEN
     CALL metax( 2.0_DP*rhoup, 4.0_DP*grhoup2, &
                 2.0_DP*tauup, sxup, v1xup, v2xup, v3xup )
  ELSE
     sxup = 0.0_DP
     v1xup = 0.0_DP
     v2xup = 0.0_DP
     v3xup = 0.0_DP
  ENDIF
  !
  IF (rhodw > small .AND. SQRT(ABS(grhodw2)) > small &
       .AND. ABS(taudw) > small) THEN
     CALL metax( 2.0_DP*rhodw, 4.0_DP*grhodw2, &
                 2.0_DP*taudw, sxdw, v1xdw, v2xdw, v3xdw )
  ELSE
     sxdw = 0.0_DP
     v1xdw = 0.0_DP
     v2xdw = 0.0_DP
     v3xdw = 0.0_DP
  ENDIF
  !
  sx = 0.5_DP*(sxup+sxdw)
  v2xup = 2.0_DP*v2xup
  v2xdw = 2.0_DP*v2xdw
  !
  RETURN
  !
END SUBROUTINE tpsscx_spin
!
!
!---------------------------------------------------------------------------
SUBROUTINE tpsscc_spin( rho, zeta, grhoup, grhodw, &
                        atau, sc, v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw )
  !--------------------------------------------------------------------------
  !! TPSS metaGGA for correlations - Hartree a.u.
  !
  USE kind_l,       ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the total charge
  REAL(DP), INTENT(IN) :: zeta
  !! the magnetization
  REAL(DP), INTENT(IN) :: atau
  !! the total kinetic energy density
  REAL(DP), INTENT(IN) :: grhoup(3)
  !! the gradient of the charge - up
  REAL(DP), INTENT(IN) :: grhodw(3)
  !! the gradient of the charge - down
  REAL(DP), INTENT(OUT) :: sc
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1cup
  !! derivatives of correlation wr. rho - up
  REAL(DP), INTENT(OUT) :: v1cdw
  !! derivatives of correlation wr. rho - down
  REAL(DP), INTENT(OUT) :: v2cup(3)
  !! derivatives of correlation wr. grho - up
  REAL(DP), INTENT(OUT) :: v2cdw(3)
  !! derivatives of correlation wr. grho - down
  REAL(DP), INTENT(OUT) :: v3cup
  !! derivatives of correlation wr. tau - up
  REAL(DP), INTENT(OUT) :: v3cdw
  !! derivatives of correlation wr. tau - down
  !
  ! ... local variables
  !
  REAL(DP) :: grho_vec(3)
  REAL(DP) :: v3c, grho !grho=grho2
  INTEGER :: ipol
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  ! vector
  grho_vec = grhoup + grhodw
  grho = 0.0_DP
  !
  DO ipol = 1, 3
     grho = grho + grho_vec(ipol)**2
  ENDDO
  !
  IF (rho <= small .OR. ABS(zeta) > 1.0_DP .OR. SQRT(ABS(grho)) <= small &
       .OR. ABS(atau) < small )  THEN
     !
     sc = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
     !
     v2cup(:) = 0.0_DP
     v2cdw(:) = 0.0_DP
     !
     v3cup = 0.0_DP
     v3cdw = 0.0_DP
     !
     v3c = 0.0_DP
  ELSE
     CALL metac_spin( rho, zeta, grhoup, grhodw, &
                        atau, sc, v1cup, v1cdw, v2cup, v2cdw, v3c )
  ENDIF
  !
  v3cup = v3c
  v3cdw = v3c
  !
  RETURN
  !
END SUBROUTINE tpsscc_spin
!
!
!---------------------------------------------------------------
SUBROUTINE metac_spin( rho, zeta, grhoup, grhodw, &
                       tau, sc, v1up, v1dw, v2up, v2dw, v3 )
  !---------------------------------------------------------------
  !! TPSS meta-GGA correlation energy and potentials - polarized case.
  !
  USE kind_l,    ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the total charge
  REAL(DP), INTENT(IN) :: zeta
  !! the magnetization
  REAL(DP), INTENT(IN) :: grhoup(3)
  !! the gradient of the charge - up
  REAL(DP), INTENT(IN) :: grhodw(3)
  !! the gradient of the charge - down
  REAL(DP), INTENT(IN) :: tau
  !! the total kinetic energy density
  REAL(DP), INTENT(OUT) :: sc
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1up
  !! derivatives of correlation wr. rho - up
  REAL(DP), INTENT(OUT) :: v1dw
  !! derivatives of correlation wr. rho - down
  REAL(DP), INTENT(OUT) :: v2up(3)
  !! derivatives of correlation wr. grho - up
  REAL(DP), INTENT(OUT) :: v2dw(3)
  !! derivatives of correlation wr. grho - down
  REAL(DP), INTENT(OUT) :: v3
  !! derivatives of correlation wr. tau
  !
  ! ... local variables
  !
  REAL(DP) :: rhoup, rhodw, tauw, grhovec(3), grho2, grho, &
              grhoup2, grhodw2
  !
  REAL(DP) :: rs, ec_u, vc_u(2), v1_0v(2), v1_pbe(2)
  !
  !grhovec   vector gradient of rho
  !grho    mod of gradient of rho
  !REAL(DP) :: ec_u, vcup_u, vcdw_u
  REAL(DP) :: ec_pbe, v1up_pbe, v1dw_pbe, v2up_pbe(3), v2dw_pbe(3)
  REAL(DP) :: ecup_0, v1up_0, v2up_0(3), v2_tmp
  REAL(DP) :: ecdw_0, v1dw_0, v2dw_0(3)
  REAL(DP) :: ec_rev, cab, aa, bb, aa2
  !
  REAL(DP) :: z2, z, ca0, dca0da, dcabda, dcabdb
  REAL(DP) :: term(3), term1, term2, term3
  !
  REAL(DP) :: drev1up, drev1dw, drev2up(3), drev2dw(3), drev3
  REAL(DP) :: sum, dsum1up, dsum1dw, dsum2up(3), dsum2dw(3)
  !
  REAL(DP) :: dcab1up, dcab1dw, dcab2up(3), dcab2dw(3)
  REAL(DP) :: db1up,   db1dw,   db2up(3),   db2dw(3)
  REAL(DP) :: da1up,   da1dw
  REAL(DP) :: ecup_til, ecdw_til
  !
  REAL(DP) :: v1up_uptil, v1up_dwtil, v2up_uptil(3), v2up_dwtil(3)
  REAL(DP) :: v1dw_uptil, v1dw_dwtil, v2dw_uptil(3), v2dw_dwtil(3)
  !
  REAL(DP), PARAMETER :: small=1.0E-10_DP, &
                         fac=3.09366772628013593097_DP**2
  !                      fac = (3*PI**2)**(2/3)
  REAL(DP), PARAMETER :: pi34=0.75_DP/3.141592653589793_DP, &
                         p43=4.0_DP/3.0_DP, third=1.0_DP/3.0_DP
  INTEGER :: ipol
  !
  rhoup = (1+zeta)*0.5_DP*rho
  rhodw = (1-zeta)*0.5_DP*rho
  grho2   = 0.0_DP
  grhoup2 = 0.0_DP
  grhodw2 = 0.0_DP
  !
  DO ipol = 1, 3
     grhovec(ipol) = grhoup(ipol) + grhodw(ipol)
     grho2   = grho2   + grhovec(ipol)**2
     grhoup2 = grhoup2 + grhoup(ipol)**2
     grhodw2 = grhodw2 + grhodw(ipol)**2
  ENDDO
  !
  grho = SQRT(grho2)
  !
  IF (rho > small) THEN
     !
     v2_tmp = 0.0_DP
     rs = (pi34/rho)**third
     CALL pw_spin( rs, zeta, ec_u, vc_u(1), vc_u(2) )
     !
     IF ((ABS(grho) > small) .AND. (zeta <= 1.0_DP)) THEN
        CALL pbec_spin( rho, zeta, grho2, 1, ec_pbe, v1_pbe(1), v1_pbe(2), v2_tmp )
        v1up_pbe=v1_pbe(1)
        v1dw_pbe=v1_pbe(2)
     ELSE
        ec_pbe   = 0.0_DP
        v1up_pbe = 0.0_DP
        v1dw_pbe = 0.0_DP
        v2up_pbe = 0.0_DP
     ENDIF
     !
     ec_pbe = ec_pbe/rho+ec_u
     ! v1xx_pbe = D_epsilon_c/ D_rho_xx   :xx= up, dw
     v1up_pbe = (v1up_pbe+vc_u(1)-ec_pbe)/rho
     v1dw_pbe = (v1dw_pbe+vc_u(2)-ec_pbe)/rho
     ! v2xx_pbe = (D_Ec / D grho)/rho = (D_Ec/ D|grho|/|grho|)*grho/rho
     v2up_pbe = v2_tmp/rho*grhovec
     ! v2dw === v2up for PBE
     v2dw_pbe = v2up_pbe
  ELSE
     ec_pbe   = 0.0_DP
     v1up_pbe = 0.0_DP
     v1dw_pbe = 0.0_DP
     v2up_pbe = 0.0_DP
     v2dw_pbe = 0.0_DP
  ENDIF
  !
  ! ec_pbe(rhoup,0,grhoup,0)
  IF (rhoup > small) THEN
     v2_tmp = 0.0_DP
     !
     rs = (pi34/rhoup)**third
     CALL pw_spin( rs, 1.0_DP, ec_u, vc_u(1), vc_u(2) )
     !
     IF (SQRT(grhoup2) > small) THEN
        CALL pbec_spin( rhoup, 1.0_DP-small, grhoup2, 1, &
                        ecup_0, v1_0v(1), v1_0v(2), v2_tmp )
        v1up_0 = v1_0v(1)
        v1dw_0 = v1_0v(2)
     ELSE
        ecup_0 = 0.0_DP
        v1up_0 = 0.0_DP
        v2up_0 = 0.0_DP
     ENDIF
     !
     ecup_0 = ecup_0/rhoup + ec_u
     v1up_0 = (v1up_0 + vc_u(1)-ecup_0)/rhoup
     v2up_0 = v2_tmp/rhoup*grhoup
  ELSE
     ecup_0 = 0.0_DP
     v1up_0 = 0.0_DP
     v2up_0 = 0.0_DP
  ENDIF
  !
  IF (ecup_0 > ec_pbe) THEN
     ecup_til = ecup_0
     v1up_uptil = v1up_0
     v2up_uptil = v2up_0
     v1up_dwtil = 0.0_DP
     v2up_dwtil = 0.0_DP
  ELSE
     ecup_til = ec_pbe
     v1up_uptil = v1up_pbe
     v1up_dwtil = v1dw_pbe
     v2up_uptil = v2up_pbe
     v2up_dwtil = v2up_pbe
  ENDIF
  ! ec_pbe(rhodw,0,grhodw,0)
  ! zeta = 1.0_DP
  IF (rhodw > small) THEN
     v2_tmp = 0.0_DP
     !
     rs = (pi34/rhodw)**third
     CALL pw_spin( rs, -1.0_DP, ec_u, vc_u(1), vc_u(2) )
     !
     IF (SQRT(grhodw2) > small) THEN
        !
        CALL pbec_spin( rhodw, -1.0_DP+small, grhodw2, 1, &
                        ecdw_0, v1_0v(1), v1_0v(2), v2_tmp )
        !
        v1up_0 = v1_0v(1)
        v1dw_0 = v1_0v(2)
        !
     ELSE
        ecdw_0 = 0.0_DP
        v1dw_0 = 0.0_DP
        v2dw_0 = 0.0_DP
     ENDIF
     !
     ecdw_0 = ecdw_0/rhodw + ec_u
     v1dw_0 = (v1dw_0 + vc_u(2)-ecdw_0)/rhodw
     v2dw_0 = v2_tmp/rhodw*grhodw
  ELSE
     ecdw_0 = 0.0_DP
     v1dw_0 = 0.0_DP
     v2dw_0 = 0.0_DP
  ENDIF
  !
  IF (ecdw_0 > ec_pbe) THEN
     ecdw_til = ecdw_0
     v1dw_dwtil = v1dw_0
     v2dw_dwtil = v2dw_0
     v1dw_uptil = 0.0_DP
     v2dw_uptil = 0.0_DP
  ELSE
     ecdw_til = ec_pbe
     v1dw_dwtil = v1dw_pbe
     v2dw_dwtil = v2dw_pbe
     v1dw_uptil = v1up_pbe
     v2dw_uptil = v2dw_pbe
  ENDIF
  !cccccccccccccccccccccccccccccccccccccccccc-------checked
  sum = (rhoup*ecup_til+rhodw*ecdw_til)/rho
  dsum1up = (ecup_til-ecdw_til)*rhodw/rho**2  &
            + (rhoup*v1up_uptil + rhodw*v1dw_uptil)/rho
  dsum1dw = (ecdw_til-ecup_til)*rhoup/rho**2  &
            + (rhodw*v1dw_dwtil + rhoup*v1up_dwtil)/rho
  ! vector
  dsum2up = (rhoup*v2up_uptil + rhodw*v2dw_uptil)/rho
  dsum2dw = (rhodw*v2dw_dwtil + rhoup*v2up_dwtil)/rho
  !ccccccccccccccccccccccccccccccccccccccccc---------checked
  aa = zeta
  !  bb = (rho*(grhoup-grhodw) - (rhoup-rhodw)*grho)**2 &
  !       /(4.0_DP*fac*rho**(14.0_DP/3.0_DP))
  bb = 0.0_DP
  DO ipol = 1, 3
     term(ipol) = rhodw*grhoup(ipol)-rhoup*grhodw(ipol)
     bb = bb + term(ipol)**2
  ENDDO
  ! vector
  term = term/(fac*rho**(14.0_DP/3.0_DP))
  bb = bb/(fac*rho**(14.0_DP/3.0_DP))
  ! bb = (rhodw*grhoup-rhoup*grhodw)**2/fac*rho**(-14.0_DP/3.0_DP)
  aa2 = aa*aa
  ca0 = 0.53_DP+aa2*(0.87_DP+aa2*(0.50_DP+aa2*2.26_DP))
  dca0da = aa*(1.74_DP+aa2*(2.0_DP+aa2*13.56_DP))
  !
  IF (ABS(aa) <= 1.0_DP-small) THEN
     term3 = (1.0_DP+aa)**(-p43) + (1.0_DP-aa)**(-p43)
     term1 = (1.0_DP+bb*0.50_DP*term3)
     term2 = (1.0_DP+aa)**(-7.0_DP/3.0_DP) + (1.0_DP-aa)**(-7.0_DP/3.0_DP)
     cab = ca0/term1**4
     dcabda = (dca0da/ca0 + 8.0_DP/3.0_DP*bb*term2/term1)*cab
     dcabdb = -2.0_DP*cab*term3/term1
  ELSE
     cab = 0.0_DP
     dcabda = 0.0_DP
     dcabdb = 0.0_DP
  ENDIF
  !
  da1up = 2.0_DP*rhodw/rho**2
  da1dw = -2.0_DP*rhoup/rho**2
  db1up = -2.0_DP*(grhodw(1)*term(1)+grhodw(2)*term(2)+grhodw(3)*term(3)) &
          -14.0_DP/3.0_DP*bb/rho
  db1dw =  2.0_DP*(grhoup(1)*term(1)+grhoup(2)*term(2)+grhoup(3)*term(3)) &
          -14.0_DP/3.0_DP*bb/rho
  !vector, not scalar
  db2up =  term*rhodw*2.0_DP
  db2dw = -term*rhoup*2.0_DP
  !
  dcab1up = dcabda*da1up + dcabdb*db1up
  dcab1dw = dcabda*da1dw + dcabdb*db1dw
  !vector, not scalar
  dcab2up = dcabdb*db2up
  dcab2dw = dcabdb*db2dw
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccc------checked
  tauw = 0.1250_DP*grho2/rho
  z = tauw/tau
  z2 = z*z
  !
  term1 = 1.0_DP+cab*z2
  term2 = (1.0_DP+cab)*z2
  ec_rev = ec_pbe*term1-term2*sum
  !
  drev1up = v1up_pbe*term1 &
            + ec_pbe*(z2*dcab1up - 2.0_DP*cab*z2/rho) &
            + (2.0_DP*term2/rho - z2*dcab1up)*sum &
            - term2*dsum1up
  !
  drev1dw = v1dw_pbe*term1 &
            + ec_pbe*(z2*dcab1dw - 2.0_DP*cab*z2/rho)  &
            + (2.0_DP*term2/rho - z2*dcab1dw)*sum  &
            - term2*dsum1dw
  !
  ! vector, not scalar
  drev2up = v2up_pbe*term1 &
            + ec_pbe*(z2*dcab2up+0.5_DP*cab*z/(rho*tau)*grhovec) &
            - (term2*4.0_DP/grho2*grhovec + z2*dcab2up)*sum &
            - term2*dsum2up
  !
  drev2dw = v2dw_pbe*term1 &
            + ec_pbe*(z2*dcab2dw+0.5_DP*cab*z/(rho*tau)*grhovec) &
            - (term2*4.0_DP/grho2*grhovec + z2*dcab2dw)*sum  &
            - term2*dsum2dw
  !
  drev3 = ((1.0_DP+cab)*sum-ec_pbe*cab)*2.0_DP*z2/tau
  !ccccccccccccccccccccccccccccccccccccccccccccccccccc----checked
  term1 = ec_rev*(1.0_DP+2.8_DP*ec_rev*z2*z)
  term2 = (1.0_DP+5.6_DP*ec_rev*z2*z)*rho
  term3 = -8.4_DP*ec_rev*ec_rev*z2*z
  !
  v1up = term1 + term2*drev1up + term3
  v1dw = term1 + term2*drev1dw + term3
  !
  term3 = term3*rho
  v3 = term2*drev3 + term3/tau
  !
  term3 = -2.0_DP*term3/grho2    !grho/|grho|^2 = 1/grho
  v2up = term2*drev2up + term3*grhovec
  v2dw = term2*drev2dw + term3*grhovec
  !
  !  call pw_spin((pi34/rho)**third,zeta,ec_u,vcup_u,vcdw_u)
  sc = rho*ec_rev*(1.0_DP+2.8_DP*ec_rev*z2*z) !-rho*ec_u
  !  v1up=v1up-vcup_u
  !  v1dw=v1dw-vcdw_u
  !
  RETURN
  !
END SUBROUTINE metac_spin
!
!                          END TPSSS
!=========================================================================
!
!                           --- M06L ---
!
!           input:  - rho
!                   - grho2=|\nabla rho|^2
!                   - tau = the input kinetic energy density
!                           It is defined as 0.5summ_i( |nabla phi_i|**2 )
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
SUBROUTINE m06lxc( rho, grho2, tau, ex, ec, v1x, v2x, v3x, v1c, v2c, v3c )
  !-----------------------------------------------------------------------
  !! Wrapper to M06L exchange+correlation routines (unpolarized).
  !
  USE kind_l,        ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the total charge
  REAL(DP), INTENT(IN) :: grho2
  !! square of gradient of rho
  REAL(DP), INTENT(IN) :: tau
  !! the total kinetic energy density
  REAL(DP), INTENT(OUT) :: ex
  !! exchange energy
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1x
  !! derivatives of exchange wr. rho
  REAL(DP), INTENT(OUT) :: v2x
  !! derivatives of exchange wr. grho
  REAL(DP), INTENT(OUT) :: v3x
  !! derivatives of exchange wr. tau
  REAL(DP), INTENT(OUT) :: v1c
  !! derivatives of correlation wr. rho
  REAL(DP), INTENT(OUT) :: v2c
  !! derivatives of correlation wr. grho
  REAL(DP), INTENT(OUT) :: v3c
  !! derivatives of correlation wr. tau
  !
  ! ... local variables
  !
  REAL(DP) :: rhoa, rhob, grho2a, grho2b, taua, taub, v1cb, v2cb, v3cb
  REAL(DP), PARAMETER :: zero = 0.0_dp, two = 2.0_dp, four = 4.0_dp
  !
  !
  rhoa = rho / two   ! one component only
  rhob = rhoa
  !
  grho2a = grho2 / four
  grho2b = grho2a
  !
  taua = tau * two * 0.5_dp ! Taua, which is Tau_sigma is half Tau
  taub = taua               ! Tau is defined as summ_i( |nabla phi_i|**2 )
                            ! in the following M06L routines
  !
  CALL m06lx( rhoa, grho2a, taua, ex, v1x, v2x, v3x )
  !
  ex = two * ex  ! Add the two components up + dw
  !
  v2x = 0.5_dp * v2x
  v3x = 2.0_dp * v3x
  !
  CALL m06lc( rhoa, rhob, grho2a, grho2b, taua, taub, ec, v1c, v2c, v3c, &
              v1cb, v2cb, v3cb )
  !
  v2c = 0.5_dp * v2c
  v3c = 2.0_dp * v3c
  !
END SUBROUTINE m06lxc
!
!
!-------------------------------------------------------------------------
SUBROUTINE m06lxc_spin( rhoup, rhodw, grhoup2, grhodw2, tauup, taudw,      &
                        ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
                        v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw )
  !-----------------------------------------------------------------------
  !! Wrapper to M06L exchange+correlation routines (polarized).
  !
  USE kind_l,        ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN)  :: rhoup, rhodw, grhoup2, grhodw2, tauup, taudw
  REAL(DP), INTENT(OUT) :: ex, ec, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw,  &
                           v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw
  !
  REAL(DP) :: exup, exdw, taua, taub
  REAL(DP), PARAMETER :: zero = 0.0_dp, two = 2.0_dp
  !
  taua = tauup * two      ! Tau is defined as summ_i( |nabla phi_i|**2 )
  taub = taudw * two      ! in the rest of the routine
  !
  CALL m06lx( rhoup, grhoup2, taua, exup, v1xup, v2xup, v3xup )
  CALL m06lx( rhodw, grhodw2, taub, exdw, v1xdw, v2xdw, v3xdw )
  !
  ex = exup + exdw
  v3xup = 2.0_dp * v3xup
  v3xdw = 2.0_dp * v3xdw
  !
  CALL m06lc( rhoup, rhodw, grhoup2, grhodw2, taua, taub, &
              ec, v1cup, v2cup, v3cup, v1cdw, v2cdw, v3cdw )
  v3cup = 2.0_dp * v3cup
  v3cdw = 2.0_dp * v3cdw
  !
END SUBROUTINE m06lxc_spin
! !
!
!===============================  M06L exchange ==========================
!
!-------------------------------------------------------------------------------
SUBROUTINE m06lx( rho, grho2, tau, ex, v1x, v2x, v3x )
  !---------------------------------------------------------------------------
  !! M06L exchange.
  !
  USE kind_l,       ONLY : DP
  USE constants_l,   ONLY : pi
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the total charge
  REAL(DP), INTENT(IN) :: grho2
  !! square of gradient of rho
  REAL(DP), INTENT(IN) :: tau
  !! the total kinetic energy density
  REAL(DP), INTENT(OUT) :: ex
  !! exchange energy
  REAL(DP), INTENT(OUT) :: v1x
  !! derivatives of exchange wr. rho
  REAL(DP), INTENT(OUT) :: v2x
  !! derivatives of exchange wr. grho
  REAL(DP), INTENT(OUT) :: v3x
  !! derivatives of exchange wr. tau
  !
  ! ... local variables
  !
  REAL(DP) :: v1x_unif, ex_unif, ex_pbe, sx_pbe, v1x_pbe, v2x_pbe
  ! ex_unif:   lda \epsilon_x(rho)
  ! v2x = 1/|grho| * dsx / d|grho| = 2 *  dsx / dgrho2
  !
  REAL(DP), PARAMETER :: zero = 0._dp,  one  = 1.0_dp, two = 2.0_dp, three = 3.0_dp,  &
                         four = 4.0_dp, five = 5.0_dp, six = 6.0_dp, eight = 8.0_dp,  &
                         f12 = one/two,    f13 = one/three,   f23 = two/three,  &
                         f53 = five/three, f83 = eight/three, f43 = four/three, &
                         pi34 = pi*three/four, pi2 = pi*pi, small = 1.d-10
  !
  REAL(DP) :: d0, d1, d2, d3, d4, d5, CF, CT, CX, alpha
  REAL(DP), DIMENSION(0:11) :: at
  INTEGER :: i
  !
  ! VSXC98 variables (LDA part)
  !
  REAL(DP) :: xs, xs2, grho, rhom83, rho13, rho43, zs, gh
  REAL(DP) :: hg, dhg_dxs2, dhg_dzs
  REAL(DP) :: dxs2_drho, dxs2_dgrho2, dzs_drho, dzs_dtau
  REAL(DP) :: ex_vs98, v1x_vs98, v2x_vs98, v3x_vs98
  !
  ! GGA and MGGA variables
  !
  REAL(DP) :: tau_unif, ts, ws, fws, dfws, dfws_drho, dfws_dtau, &
              dws_dts, dts_dtau, dts_drho
  !
  ! set parameters
  !
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
  !
  d0     =    6.012244d-01
  d1     =    4.748822d-03
  d2     =   -8.635108d-03
  d3     =   -9.308062d-06
  d4     =    4.482811d-05
  d5     =    zero
  !
  alpha = 1.86726d-03
  !
  IF (rho < small .OR. tau < small) THEN  !.AND.?
     ex = zero
     v1x = zero
     v2x = zero
     v3x = zero
     RETURN
  ENDIF
  !
  ! ... VSXC98 functional (LDA part)
  !
  ! set variables
  !
  CF =  (three/five) * (six*pi2)**f23
  CT =  CF / two
  CX = -(three/two) * (three/(four*pi))**f13  ! Cx LSDA
  !
  !  IF (rho >= small .AND. grho>=small) THEN
  !
  grho = SQRT(grho2)
  rho43 = rho**f43
  rho13 = rho**f13
  rhom83 = one / rho**f83
  !
  xs  = grho / rho43
  xs2 = xs * xs
  zs  = tau / rho**f53 - CF
  gh  = one + alpha * (xs2 + zs)
  !
  IF (gh >= small) THEN
    CALL gvt4( xs2, zs, d0, d1, d2, d3, d4, d5, alpha, hg, dhg_dxs2, dhg_dzs )
  ELSE
    hg = zero
    dhg_dxs2 = zero
    dhg_dzs  = zero
  ENDIF
  !
  dxs2_drho   = -f83*xs2/rho
  dxs2_dgrho2 =  rhom83
  dzs_drho = -f53*tau*rhom83
  dzs_dtau =  one/rho**f53
  !
  ex_unif  =  CX * rho43
  ex_vs98  =  ex_unif * hg
  v1x_vs98 =  CX * ( f43 * hg * rho**f13 ) + &
              ex_unif * ( dhg_dxs2*dxs2_drho + dhg_dzs*dzs_drho )
  v2x_vs98 =  two * ex_unif * dhg_dxs2 * dxs2_dgrho2
  v3x_vs98 =  ex_unif * dhg_dzs * dzs_dtau
  !
  ! ... mo6lx functional
  !
  tau_unif = CF * rho**f53  ! Tau is define as summ_i( |nabla phi_i|**2 )
  ts = tau_unif / tau
  ws = (ts - one)/(ts + one)
  !
  fws  = zero
  dfws = zero
  !
  DO i = 0, 11
    fws  =  fws + at(i)*ws**i
    dfws = dfws + i*at(i)*ws**(i-1)
  ENDDO
  !
  dws_dts  = two/((ts+1)**2)
  dts_drho = ( (six*pi*pi*rho)**f23 )/tau
  dts_dtau = -ts/tau
  dfws_drho = dfws*dws_dts*dts_drho
  dfws_dtau = dfws*dws_dts*dts_dtau
  !
  CALL pbex_m06l( two*rho, four*grho2, sx_pbe, v1x_pbe, v2x_pbe )
  !
  v1x_unif = f43 * CX * rho13
  !
  sx_pbe  = f12 * sx_pbe
  v1x_pbe = v1x_pbe  + v1x_unif
  v2x_pbe = two * v2x_pbe
  !
  ex_pbe  = sx_pbe + ex_unif
  !
  ! ... energy and potential
  !
  ex = ex_vs98  + ex_pbe*fws
  !
  v1x = v1x_vs98 + v1x_pbe*fws + ex_pbe*dfws_drho
  v2x = v2x_vs98 + v2x_pbe*fws
  v3X = v3x_vs98 + ex_pbe*dfws_dtau
  !
END SUBROUTINE m06lx
!
!
!-------------------------------------------------------------------
SUBROUTINE pbex_m06l( rho, grho2, sx, v1x, v2x )
  !---------------------------------------------------------------
  !! PBE exchange (without Slater exchange):
  !! J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  !
  !! v2x = 1/|grho| * dsx / d|grho| = 2 *  dsx / dgrho2
  !
  USE kind_l
  USE constants_l,   ONLY : pi
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! charge density
  REAL(DP), INTENT(IN) :: grho2
  !! squared gradient
  REAL(DP), INTENT(OUT) :: sx
  !! energy
  REAL(DP), INTENT(OUT) :: v1x
  !! potential w.r. rho
  REAL(DP), INTENT(OUT) :: v2x
  !! potential w.r. grho
  !
  ! ... local variables
  !
  REAL(DP) :: grho, rho43, xs, xs2, dxs2_drho, dxs2_dgrho2
  REAL(DP) :: CX, denom, C1, C2, ex, Fx, dFx_dxs2, dex_drho
  !
  REAL(DP), PARAMETER :: mu = 0.21951_dp, ka = 0.804_dp, &
                         one  = 1.0_dp, two = 2.0_dp, three = 3.0_dp, &
                         four = 4.0_dp, six = 6.0_dp, eight = 8.0_dp, &
                         f13 = one/three,  f23 = two/three,  f43 = four/three, &
                         f34 = three/four, f83 = eight/three
  !
  CX    =  f34 * (three/pi)**f13            ! Cx LDA
  denom =  four * (three*pi**two)**f23
  C1    =  mu / denom
  C2    =  mu / (ka * denom)
  !
  grho  = SQRT(grho2)
  rho43 = rho**f43
  xs    = grho / rho43
  xs2   = xs * xs
  !
  dxs2_drho = -f83 * xs2 / rho
  dxs2_dgrho2 = one /rho**f83
  !
  ex = - CX * rho43
  dex_drho = - f43 * CX * rho**f13
  !
  Fx = C1*xs2 / (one + C2*xs2)
  dFx_dxs2 = C1 / (one + C2*xs2)**2
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
END SUBROUTINE pbex_m06l
!
!
!===============================  M06L correlation ==========================
!
!---------------------------------------------------------------------------------
SUBROUTINE m06lc( rhoa, rhob, grho2a, grho2b, taua, taub, ec, v1c_up, v2c_up, &
                  v3c_up, v1c_dw, v2c_dw, v3c_dw )
  !-------------------------------------------------------------------------------
  !! M06L correlation.
  !
  USE kind_l,        ONLY : DP
  USE constants_l,    ONLY : pi
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rhoa
  !! charge density up
  REAL(DP), INTENT(IN) :: rhob
  !! charge density down
  REAL(DP), INTENT(IN) :: grho2a
  !! squared gradient up
  REAL(DP), INTENT(IN) :: grho2b
  !! squared gradient down
  REAL(DP), INTENT(IN) :: taua
  !! kinetic energy density up
  REAL(DP), INTENT(IN) :: taub
  !! kinetic energyt density down
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1c_up
  !! corr. potential w.r. rho - up
  REAL(DP), INTENT(OUT) :: v2c_up
  !! corr. potential w.r. rho - up
  REAL(DP), INTENT(OUT) :: v3c_up
  !! corr. potential w.r. rho - up
  REAL(DP), INTENT(OUT) :: v1c_dw
  !! corr. potential w.r. rho - down
  REAL(DP), INTENT(OUT) :: v2c_dw
  !! corr. potential w.r. grho - down
  REAL(DP), INTENT(OUT) :: v3c_dw
  !! corr. potential w.r. tau - down
  !
  ! ... local variables
  !
  REAL(DP) :: rs, zeta
  REAL(DP) :: vc_v(2)
  !
  REAL(DP), PARAMETER :: zero = 0._dp,  one  = 1.0_dp, two = 2.0_dp, three = 3.0_dp, &
                         four = 4.0_dp, five = 5.0_dp, six = 6.0_dp, eight = 8.0_dp, &
                         f12 = one/two, f13 = one/three, f23 = two/three,            &
                         f53 = five/three, f83 = eight/three, f43 = four/three,      &
                         pi34 = three/(four*pi), pi2 = pi*pi, f35 = three/five,      &
                         small = 1.d-10
  !
  ! parameters of the MO6Lc functional
  !
  REAL(DP), DIMENSION(0:4):: cs, cab
  !
  REAL(DP) :: ds0, ds1, ds2, ds3, ds4, ds5, CF,         &
              dab0, dab1, dab2, dab3, dab4, dab5, gama_ab, gama_s, &
              alpha_s, alpha_ab
  !
  ! functions and variables
  !
  REAL(DP) :: ec_pw_a, ec_pw_b, ec_pw_ab
  !
  REAL(DP) :: vv, vc_pw_a, vc_pw_b, vc_pw_up, vc_pw_dw, Ecaa, Ecbb, Ecab, &
              Ec_UEG_ab, Ec_UEG_aa, Ec_UEG_bb, decab_drhoa, decab_drhob,            &
              v1_ab_up, v1_ab_dw, v2_ab_up, v2_ab_dw, v3_ab_up, v3_ab_dw,           &
              v1_aa_up, v2_aa_up, v3_aa_up, v1_bb_dw, v2_bb_dw, v3_bb_dw
  !
  REAL(DP) :: xsa, xs2a, rsa, grhoa, xsb, xs2b, grhob, rsb, zsa, zsb, &
              xs2ab, zsab, rho,                                       &
              dxs2a_drhoa, dxs2b_drhob, dxs2a_dgrhoa2, dxs2b_dgrhob2, &
              dzsa_drhoa, dzsb_drhob, dzsa_dtaua, dzsb_dtaub
  !
  REAL(DP) :: hga, dhga_dxs2a, dhga_dzsa, hgb, dhgb_dxs2b, dhgb_dzsb,   &
              hgab, dhgab_dxs2ab, dhgab_dzsab,                          &
              Dsa, Dsb, dDsa_dxs2a, dDsa_dzsa, dDsb_dxs2b, dDsb_dzsb,   &
              gsa, gsb, gsab, dgsa_dxs2a, dgsb_dxs2b, dgsab_dxs2ab, num
  !
  INTEGER :: ifunc
  !
  !
  dab0 =  3.957626d-01
  dab1 = -5.614546d-01
  dab2 =  1.403963d-02
  dab3 =  9.831442d-04
  dab4 = -3.577176d-03
  dab5 =  zero
  !
  cab(0) =  6.042374d-01
  cab(1) =  1.776783d+02
  cab(2) = -2.513252d+02
  cab(3) =  7.635173d+01
  cab(4) = -1.255699d+01
  !
  gama_ab  = 0.0031_dp
  alpha_ab = 0.00304966_dp
  !
  ds0 =  4.650534d-01
  ds1 =  1.617589d-01
  ds2 =  1.833657d-01
  ds3 =  4.692100d-04
  ds4 = -4.990573d-03
  ds5 =  zero
  !
  cs(0) =  5.349466d-01
  cs(1) =  5.396620d-01
  cs(2) = -3.161217d+01
  cs(3) =  5.149592d+01
  cs(4) = -2.919613d+01
  !
  gama_s  = 0.06_dp
  alpha_s = 0.00515088_dp
  !
  CF = f35 * (six*pi2)**f23
  !
  ifunc = 1     ! iflag=1  J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  ! ... Ecaa
  !
  IF (rhoa < small .OR. taua < small) THEN  !.AND.?
    !
    Ecaa = zero
    v1_aa_up = zero
    v2_aa_up = zero
    v3_aa_up = zero
    ec_pw_a = zero
    vc_pw_a = zero
    xs2a = zero
    zsa = zero
    dxs2a_drhoa   = zero
    dxs2a_dgrhoa2 = zero
    dzsa_drhoa    = zero
    dzsa_dtaua    = zero
    !
  ELSE
    !
    rsa   = (pi34/rhoa)**f13
    grhoa = SQRT(grho2a)
    xsa   = grhoa / rhoa**f43
    xs2a  = xsa * xsa
    zsa   = taua/rhoa**f53 - CF
    !
    dxs2a_drhoa   = -f83*xs2a/rhoa
    dxs2a_dgrhoa2 =  one/(rhoa**f83)
    !
    dzsa_drhoa   = -f53*taua/(rhoa**f83)
    dzsa_dtaua   =  one/rhoa**f53
    !
    Dsa        = one - xs2a/(four * (zsa + CF))
    dDsa_dxs2a = - one/(four * (zsa + CF))
    dDsa_dzsa  = xs2a/(four * (zsa + CF)**2)
    !
    ec_pw_a = zero
    vc_pw_a = zero
    !
    rs   = rsa
    zeta = one
    CALL pw_spin( rs, zeta, ec_pw_a, vc_v(1), vc_v(2) )
    vc_pw_a = vc_v(1)
    vv      = vc_v(2)
    !
    CALL gvt4( xs2a, zsa, ds0, ds1, ds2, ds3, ds4, ds5, alpha_s, hga, dhga_dxs2a, dhga_dzsa )
    CALL gfunc( cs, gama_s, xs2a, gsa, dgsa_dxs2a )
    !
    Ec_UEG_aa = rhoa*ec_pw_a
    num = (dgsa_dxs2a + dhga_dxs2a)*Dsa + (gsa + hga)*dDsa_dxs2a
    !
    !
    Ecaa = Ec_UEG_aa * (gsa + hga) * Dsa
    !
    v1_aa_up = vc_pw_a * (gsa + hga) * Dsa +                          &
               Ec_UEG_aa * num * dxs2a_drhoa +                        &
               Ec_UEG_aa * (dhga_dzsa*Dsa + (gsa + hga)*dDsa_dzsa) * dzsa_drhoa
    !
    v2_aa_up = two * Ec_UEG_aa * num * dxs2a_dgrhoa2
    !
    v3_aa_up = Ec_UEG_aa * (dhga_dzsa*Dsa + (gsa + hga)*dDsa_dzsa) * dzsa_dtaua
    !
  ENDIF
  !
  ! ... Ecbb
  !
  IF (rhob < small .OR. taub < small) THEN
    !
    Ecbb = zero
    v1_bb_dw = zero
    v2_bb_dw = zero
    v3_bb_dw = zero
    ec_pw_b = zero
    vc_pw_b = zero
    xs2b = zero
    zsb = zero
    dxs2b_drhob   = zero
    dxs2b_dgrhob2 = zero
    dzsb_drhob    = zero
    dzsb_dtaub    = zero
    !
  ELSE
    !
    rsb   = (pi34/rhob)**f13
    grhob = SQRT(grho2b)
    xsb   = grhob / rhob**f43
    xs2b  = xsb * xsb
    zsb   = taub/rhob**f53 - CF

    dxs2b_drhob   = -f83*xs2b/rhob
    dxs2b_dgrhob2 =  one /rhob**f83

    dzsb_drhob = -f53*taub/(rhob**f83)
    dzsb_dtaub =  one/rhob**f53

    Dsb        = one - xs2b/(four * (zsb + CF))
    dDsb_dxs2b = - one/(four * (zsb + CF))
    dDsb_dzsb  =  xs2b/(four * (zsb + CF)**2)
    !
    zeta = one
    rs   = rsb
    CALL pw_spin( rs, zeta, ec_pw_b, vc_v(1), vc_v(2) )
    vc_pw_b = vc_v(1)
    vv      = vc_v(2)
    !
    CALL gvt4( xs2b, zsb, ds0, ds1, ds2, ds3, ds4, ds5, alpha_s, hgb, dhgb_dxs2b, dhgb_dzsb )
    CALL gfunc( cs, gama_s, xs2b, gsb, dgsb_dxs2b )
    !
    Ec_UEG_bb = rhob*ec_pw_b
    num = (dgsb_dxs2b + dhgb_dxs2b)*Dsb + (gsb + hgb)*dDsb_dxs2b
    !
    Ecbb = Ec_UEG_bb * (gsb + hgb) * Dsb
    !
    v1_bb_dw = vc_pw_b * (gsb + hgb) * Dsb +    &
               Ec_UEG_bb * num * dxs2b_drhob +  &
               Ec_UEG_bb * (dhgb_dzsb*Dsb + (gsb + hgb)*dDsb_dzsb)*dzsb_drhob
    !
    v2_bb_dw = two * Ec_UEG_bb * num * dxs2b_dgrhob2
    !
    v3_bb_dw = Ec_UEG_bb * (dhgb_dzsb*Dsb + (gsb + hgb)*dDsb_dzsb)*dzsb_dtaub
    !
  ENDIF
  !
  ! ... Ecab
  !
  IF (rhoa < small .AND. rhob < small) THEN
    !
    Ecab = zero
    v1_ab_up = zero
    v1_ab_dw = zero
    v2_ab_up = zero
    v2_ab_dw = zero
    v3_ab_up = zero
    v3_ab_dw = zero
    !
  ELSE
    !
    xs2ab = xs2a + xs2b
    zsab  = zsa + zsb
    rho   = rhoa + rhob
    zeta  = (rhoa - rhob)/rho
    rs    = (pi34/rho)**f13
    !
    CALL gvt4( xs2ab, zsab, dab0, dab1, dab2, dab3, dab4, dab5, alpha_ab, hgab, &
               dhgab_dxs2ab, dhgab_dzsab )
    !
    CALL pw_spin( rs, zeta, ec_pw_ab, vc_v(1), vc_v(2) )
    vc_pw_up=vc_v(1) ; vc_pw_dw=vc_v(2)
    !
    CALL gfunc( cab, gama_ab, xs2ab, gsab, dgsab_dxs2ab )
    !
    decab_drhoa = vc_pw_up - vc_pw_a
    decab_drhob = vc_pw_dw - vc_pw_b
    !
    Ec_UEG_ab = ec_pw_ab*rho - ec_pw_a*rhoa - ec_pw_b*rhob
    !
    Ecab = Ec_UEG_ab * (gsab + hgab)
    !
    v1_ab_up = decab_drhoa * (gsab + hgab) +                             &
               Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2a_drhoa + &
               Ec_UEG_ab * dhgab_dzsab * dzsa_drhoa
    !
    v1_ab_dw = decab_drhob * (gsab + hgab) +                             &
               Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2b_drhob + &
               Ec_UEG_ab * dhgab_dzsab * dzsb_drhob
    !
    v2_ab_up = two * Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2a_dgrhoa2
    v2_ab_dw = two * Ec_UEG_ab * (dgsab_dxs2ab + dhgab_dxs2ab) * dxs2b_dgrhob2

    v3_ab_up = Ec_UEG_ab * dhgab_dzsab * dzsa_dtaua
    v3_ab_dw = Ec_UEG_ab * dhgab_dzsab * dzsb_dtaub
    !
  ENDIF
  !
  ! ... ec and vc
  !
  ec = Ecaa + Ecbb + Ecab
  !
  v1c_up = v1_aa_up + v1_ab_up
  v2c_up = v2_aa_up + v2_ab_up
  v3c_up = v3_aa_up + v3_ab_up
  !
  v1c_dw = v1_bb_dw + v1_ab_dw
  v2c_dw = v2_bb_dw + v2_ab_dw
  v3c_dw = v3_bb_dw + v3_ab_dw
  !
END SUBROUTINE m06lc
  !
  !
  !-------------------------------------------------------------------
SUBROUTINE gfunc( cspin, gama, xspin, gs, dgs_dx )
    !-----------------------------------------------------------------
    !
    USE kind_l,   ONLY : DP
    !
    IMPLICIT NONE
    !
    !$acc routine seq
    !
    REAL(DP), INTENT(IN)  :: cspin(0:4)
    REAL(DP), INTENT(IN)  :: xspin, gama
    REAL(DP), INTENT(OUT) :: gs, dgs_dx
    !
    REAL(DP) :: de, d2, x1, x2, x3, x4
    REAL(DP), PARAMETER :: one=1.0d0, two=2.0d0, &
                           three=3.0d0, four=4.0d0
    !
    de = one/(one + gama*xspin)
    d2 = de**2
    x1 = gama * xspin * de
    x2 = x1**2
    x3 = x1**3
    x4 = x1**4
    !
    gs = cspin(0) + cspin(1)*x1 + cspin(2)*x2 + cspin(3)*x3 + cspin(4)*x4
    dgs_dx = gama*d2*(cspin(1) + two*cspin(2)*x1 + three*cspin(3)*x2 + four*cspin(4)*x3)
    !
  END SUBROUTINE gfunc
  !
!
!
!-------------------------------------------------------------------------
SUBROUTINE gvt4( x, z, a, b, c, d, e, f, alpha, hg, dh_dx, dh_dz )
  !----------------------------------------------------------------------
  !
  USE kind_l,    ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: X, z, a, b, c, d, e, f, alpha
  REAL(DP), INTENT(OUT) :: hg, dh_dx, dh_dz
  !
  REAL(DP) :: gamma, gamma2, gamma3
  REAL(DP), PARAMETER :: one=1.0_dp, two=2.0_dp, three=3.0_dp
  !
  gamma  = one + alpha*(x+z)
  gamma2 = gamma*gamma
  gamma3 = gamma2*gamma
  !
  hg = a/gamma + (b*x + c*z)/gamma2 + (d*x*x + e*x*z + f*z*z)/gamma3
  !
  dh_dx = ( -alpha*a + b + (two*x*(d - alpha*b) + z*(e - two*alpha*c))/ gamma &
          - three*alpha*(d*x*x + e*x*z + f*z*z)/gamma2  )/gamma2
  !
  dh_dz = ( -alpha*a + c + (two*z*(f - alpha*c) + x*(e -two*alpha*b))/ gamma  &
          - three*alpha*(d*x*x + e*x*z + f*z*z)/gamma2  )/gamma2
  !
  RETURN
  !
END SUBROUTINE gvt4
!
!                          END M06L
!=========================================================================
END MODULE
