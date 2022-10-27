!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------------
MODULE corr_gga
  !-------------------------------------------------------------------------------
  !! GGA correlation functionals.
  !
  USE corr_lda, ONLY: pw, pw_spin
  !
CONTAINS
!
!
!-----------------------------------------------------------------------
SUBROUTINE perdew86( rho, grho, sc, v1c, v2c )
  !-----------------------------------------------------------------------
  !! Perdew gradient correction on correlation: PRB 33, 8822 (1986).
  !
  USE kind_l, ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sc, v1c, v2c
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: p1=0.023266_DP, p2=7.389d-6, p3=8.723_DP, &
                         p4=0.472_DP
  REAL(DP), PARAMETER :: pc1=0.001667_DP, pc2=0.002568_DP, pci=pc1 + pc2
  REAL(DP), PARAMETER :: third=1._DP/3._DP, pi34=0.6203504908994_DP
                                          ! pi34=(3/4pi)^(1/3)
  REAL(DP) :: rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs
  REAL(DP) :: dcna, dcnb, dcn, phi, ephi
  !
  rho13 = rho**third
  rho43 = rho13**4
  rs  = pi34 / rho13
  rs2 = rs * rs
  rs3 = rs * rs2
  !
  cna = pc2 + p1 * rs + p2 * rs2
  cnb = 1._DP + p3 * rs + p4 * rs2 + 1.d4 * p2 * rs3
  cn = pc1 + cna / cnb
  !
  drs  = - third * pi34 / rho43
  dcna = (p1 + 2._DP * p2 * rs) * drs
  dcnb = (p3 + 2._DP * p4 * rs + 3.d4 * p2 * rs2) * drs
  dcn  = dcna / cnb - cna / (cnb * cnb) * dcnb
  !
  phi = 0.192_DP * pci / cn * SQRT(grho) * rho**(-7._DP/6._DP)
  ! SdG: in the original paper 1.745*0.11=0.19195 is used
  ephi = EXP( -phi )
  !
  sc  = grho / rho43 * cn * ephi
  v1c = sc * ( (1._DP+phi) * dcn / cn - ((4._DP/3._DP)-(7._DP/ &
                  6._DP)*phi) / rho )
  v2c = cn * ephi / rho43 * (2._DP - phi)
  !
  RETURN
  !
END SUBROUTINE perdew86
!
!
!-----------------------------------------------------------------------
SUBROUTINE ggac( rho, grho, sc, v1c, v2c )
  !-----------------------------------------------------------------------
  !! Perdew-Wang GGA (PW91) correlation part
  !
  USE kind_l,    ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sc, v1c, v2c
  !
  ! ... local variables
  !
  REAL(DP) :: rs, ec, vc
  !
  REAL(DP), PARAMETER :: al=0.09_DP,  pa=0.023266_DP, pb=7.389d-6, &
                         pc=8.723_DP, pd=0.472_DP,                 &
                         cx=-0.001667_DP, cxc0=0.002568_DP, cc0=-cx+cxc0
  !
  REAL(DP), PARAMETER :: third=1._DP/3._DP, pi34=0.6203504908994_DP, &
                         nu=15.755920349483144_DP, be=nu*cc0,        &
                         xkf=1.919158292677513_DP, xks=1.128379167095513_DP
                         ! pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
                         ! xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  REAL(DP) :: kf, ks, rs1, rs2, rs3, t, expe, af, bf, y, xy, qy, s1
  REAL(DP) :: h0, dh0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, &
              dh1, ddh1
  !
  rs = pi34 / rho**third
  !
  CALL pw( rs, 1, ec, vc )
  !
  rs1 = rs
  rs2 = rs1 * rs1
  rs3 = rs1 * rs2
  !
  kf = xkf / rs1
  ks = xks * SQRT(kf)
  t = SQRT(grho) / (2._DP * ks * rho)
  !
  expe = EXP( - 2._DP * al * ec / (be * be) )
  af = 2._DP * al / be * (1._DP / (expe-1._DP) )
  bf = expe * (vc - ec)
  !
  y = af * t * t
  xy = (1._DP + y) / (1._DP + y + y * y)
  qy = y * y * (2._DP + y) / (1._DP + y + y * y)**2
  s1 = 1._DP + 2._DP * al / be * t * t * xy
  !
  h0 = be * be / (2._DP * al) * LOG(s1)
  dh0 = be * t * t / s1 * ( - 7._DP / 3._DP * xy - qy * (af * bf / &
                            be-7._DP / 3._DP) )
  ddh0 = be / (2._DP * ks * ks * rho) * (xy - qy) / s1
  !
  ee = - 100._DP * (ks / kf * t)**2
  !
  cna = cxc0 + pa * rs1 + pb * rs2
  dcna = pa * rs1 + 2._DP * pb * rs2
  cnb = 1._DP + pc * rs1 + pd * rs2 + 1.d4 * pb * rs3
  dcnb = pc * rs1 + 2._DP * pd * rs2 + 3.d4 * pb * rs3
  cn = cna / cnb - cx
  dcn = dcna / cnb - cna * dcnb / (cnb * cnb)
  !
  h1 = nu * (cn - cc0 - 3._DP / 7._DP * cx) * t * t * EXP(ee)
  dh1 = - third * ( h1 * (7._DP + 8._DP * ee) + nu * t * t * EXP(ee) &
                    * dcn )
  ddh1 = 2._DP * h1 * (1._DP + ee) * rho / grho
  !
  sc = rho * (h0 + h1)
  v1c = h0 + h1 + dh0 + dh1
  v2c = ddh0 + ddh1
  !
  RETURN
  !
END SUBROUTINE ggac
!
!
!-----------------------------------------------------------------------
SUBROUTINE glyp( rho, grho, sc, v1c, v2c )
  !-----------------------------------------------------------------------
  !! Lee Yang Parr: gradient correction part.
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sc, v1c, v2c
  !
  ! ... local varibles
  !
  REAL(DP), PARAMETER :: a=0.04918_DP, b=0.132_DP, c=0.2533_DP, &
                         d=0.349_DP
  REAL(DP) :: rhom13, rhom43, rhom53, om, xl, ff, dom, dxl
  !
  rhom13 = rho**(-1._DP/3._DP)
  om = EXP(-c*rhom13) / (1._DP+d*rhom13)
  xl = 1._DP + (7._DP/3._DP) * ( c*rhom13 + d * rhom13 / (1._DP + &
                                 d * rhom13) )
  ff = a * b * grho / 24._DP
  rhom53 = rhom13**5
  !
  sc = ff * rhom53 * om * xl
  !
  dom = - om * (c + d+c * d * rhom13) / (1.d0 + d * rhom13)
  dxl = (7.d0 / 3.d0) * (c + d+2.d0 * c * d * rhom13 + c * d * d * &
       rhom13**2) / (1.d0 + d * rhom13) **2
  rhom43 = rhom13**4
  !
  v1c = - ff * rhom43 / 3.d0 * ( 5.d0 * rhom43 * om * xl + rhom53 * &
                                    dom * xl + rhom53 * om * dxl )
  v2c = 2.d0 * sc / grho
  !
  RETURN
  !
END SUBROUTINE glyp
!
!
!---------------------------------------------------------------
SUBROUTINE pbec( rho, grho, iflag, sc, v1c, v2c )
  !---------------------------------------------------------------
  !! PBE correlation (without LDA part)
  !
  !! * iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  !! * iflag=2: J.P.Perdew et al., PRL 100, 136406 (2008).
  !! * iflag=3: L. Chiodo et al, PRL 108, 126402 (2012)  (PBEQ2D)
  !
  USE kind_l,    ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  INTEGER,  INTENT(IN) :: iflag
  REAL(DP), INTENT(IN) :: rho, grho
  ! input: charge and squared gradient
  REAL(DP), INTENT(OUT) :: sc, v1c, v2c
  ! output: energy, potential
  REAL(DP), PARAMETER :: ga = 0.0310906908696548950_DP
  REAL(DP) :: be(3)
  !             pbe           pbesol   pbeq2d
  DATA be / 0.06672455060314922_DP, 0.046_DP, 0.06672455060314922_DP/
  REAL(DP), PARAMETER :: third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0
  REAL(DP), PARAMETER :: xkf = 1.919158292677513d0, xks = 1.128379167095513d0
  ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  !
  REAL(DP) :: rs, ec, vc
  !
  REAL(DP) :: kf, ks, t, expe, af, bf, y, xy, qy
  REAL(DP) :: s1, h0, dh0, ddh0, sc2D, v1c2D, v2c2D
  !
  rs = pi34 / rho**third
  !
  CALL pw( rs, 1, ec, vc )
  !
  kf = xkf / rs
  ks = xks * SQRT(kf)
  t = SQRT(grho) / (2._DP * ks * rho)
  expe = EXP( - ec / ga )
  af = be(iflag) / ga * (1._DP / (expe-1._DP))
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1._DP + y) / (1._DP + y + y * y)
  qy = y * y * (2._DP + y) / (1._DP + y + y * y)**2
  s1 = 1._DP + be(iflag) / ga * t * t * xy
  h0 = ga * LOG(s1)
  dh0 = be(iflag) * t * t / s1 * ( - 7._DP / 3._DP * xy - qy * (af * bf / &
        be(iflag)-7._DP / 3._DP) )
  ddh0 = be(iflag) / (2._DP * ks * ks * rho) * (xy - qy) / s1
  !
  sc  = rho * h0
  v1c = h0 + dh0
  v2c = ddh0
  ! q2D
  IF (iflag == 3) THEN
     CALL cpbe2d( rho, grho, sc2D, v1c2D, v2c2D )
     sc  = sc  + sc2D
     v1c = v1c + v1c2D
     v2c = v2c + v2c2D
  ENDIF
  !
  RETURN
  !
END SUBROUTINE pbec


! ===========> SPIN <===========
!
!-----------------------------------------------------------------------
SUBROUTINE perdew86_spin( rho, zeta, grho, sc, v1c_up, v1c_dw, v2c )
  !---------------------------------------------------------------------
  !! Perdew gradient correction on correlation: PRB 33, 8822 (1986)
  !! spin-polarized case.
  !
  USE kind_l,    ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the total charge density
  REAL(DP), INTENT(IN) :: zeta
  !! the magnetization
  REAL(DP), INTENT(IN) :: grho
  !! the gradient of the charge squared
  REAL(DP), INTENT(OUT) :: sc
  !! correlation energies
  REAL(DP), INTENT(OUT) :: v1c_up, v1c_dw !v1c(2)
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT) :: v2c
  !! derivatives of correlation wr. grho
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: p1=0.023266_DP, p2=7.389D-6, p3=8.723_DP, &
                         p4=0.472_DP
  REAL(DP), PARAMETER :: pc1=0.001667_DP, pc2 = 0.002568_DP, pci=pc1+pc2
  REAL(DP), PARAMETER :: third=1._DP/3._DP, pi34=0.6203504908994_DP
  !                                         pi34=(3/4pi)^(1/3)
  !
  REAL(DP) :: rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs
  REAL(DP) :: dcna, dcnb, dcn, phi, ephi, dd, ddd
  !
  rho13 = rho**third
  rho43 = rho13**4
  rs  = pi34 / rho13
  rs2 = rs * rs
  rs3 = rs * rs2
  !
  cna = pc2 + p1 * rs + p2 * rs2
  cnb = 1._DP + p3 * rs + p4 * rs2 + 1.D4 * p2 * rs3
  cn  = pc1 + cna / cnb
  !
  drs  = -third * pi34 / rho43
  dcna = (p1 + 2._DP * p2 * rs) * drs
  dcnb = (p3 + 2._DP * p4 * rs + 3.D4 * p2 * rs2) * drs
  dcn  = dcna / cnb - cna / (cnb * cnb) * dcnb
  phi  = 0.192_DP * pci / cn * SQRT(grho) * rho**(-7._DP/6._DP)
  !SdG: in the original paper 1.745*0.11=0.19195 is used
  !
  dd  = (2._DP)**third * SQRT( ( (1.0_DP + zeta) * 0.5_DP)**(5.0_DP/ &
       3._DP) + ( (1.0_DP - zeta) * 0.5d0)**(5._DP/3._DP) )
  !
  ddd = (2._DP)**(-4.0_DP/3._DP) * 5._DP * ( ((1._DP + zeta) &
       * 0.5_DP)**(2._DP/3._DP) - ((1._DP - zeta)*0.5d0)**(2._DP/ &
       3._DP))/(3._DP*dd)
  !
  ephi = EXP(-phi)
  sc = grho / rho43 * cn * ephi / dd
  !
  v1c_up = sc * ( (1._DP + phi) * dcn / cn - ( (4._DP / 3._DP) - &
       (7._DP/6._DP)*phi)/rho) - sc * ddd / dd * (1._DP - zeta)/rho
  v1c_dw = sc * ( (1._DP + phi) * dcn / cn - ( (4._DP / 3._DP) - &
       (7._DP/6._DP)*phi)/rho) + sc * ddd / dd * (1._DP + zeta)/rho
  v2c = cn * ephi / rho43 * (2._DP - phi) / dd
  !
  RETURN
  !
END SUBROUTINE perdew86_spin
!
!
!-----------------------------------------------------------------------
SUBROUTINE ggac_spin( rho, zeta, grho, sc, v1c_up, v1c_dw, v2c )
  !---------------------------------------------------------------------
  !! Perdew-Wang GGA (PW91) correlation part - spin-polarized
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! the total charge density
  REAL(DP), INTENT(IN) :: zeta
  !! the magnetization
  REAL(DP), INTENT(IN) :: grho
  !! the gradient of the charge squared
  REAL(DP), INTENT(OUT) :: sc
  !! correlation energies
  REAL(DP), INTENT(OUT) :: v1c_up, v1c_dw
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT) :: v2c
  !! derivatives of correlation wr. grho
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: al=0.09_DP,  pa=0.023266_DP, pb=7.389D-6, &
                         pc=8.723_DP, pd=0.472_DP
  REAL(DP), PARAMETER :: cx=-0.001667_DP, cxc0=0.002568_DP, cc0=-cx+cxc0
  REAL(DP), PARAMETER :: third=1._DP/3._DP, pi34=0.6203504908994_DP
  !                                         pi34=(3/4pi)^(1/3)
  REAL(DP), PARAMETER :: nu=15.755920349483144_DP, be=nu*cc0
  !                      nu=(16/pi)*(3 pi^2)^(1/3)
  REAL(DP), PARAMETER :: xkf=1.919158292677513_DP, xks=1.128379167095513_DP
  !                      xkf=(9 pi/4)^(1/3)      , xks=sqrt(4/pi)
  !
  REAL(DP) :: rs, ec, vc_up, vc_dn
  REAL(DP) :: kf, ks, rs2, rs3, t, expe, af, y, xy, qy, s1, h0, ddh0, ee, &
              cn, dcn, cna, dcna, cnb, dcnb, h1, dh1, ddh1, fz, fz2, fz3, &
              fz4, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, dh0zdw, dh1zup, &
              dh1zdw
  !
  rs = pi34 / rho**third
  !
  CALL pw_spin( rs, zeta, ec, vc_up, vc_dn )
  !
  rs2 = rs * rs
  rs3 = rs * rs2
  kf = xkf / rs
  ks = xks * SQRT(kf)
  !
  fz = 0.5_DP * ( (1._DP+zeta)**(2._DP/3._DP) + (1._DP-zeta)**(2._DP/3._DP) )
  fz2 = fz  * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1._DP+zeta)**(-1._DP/3._DP) - (1._DP-zeta)**(-1._DP/3._DP) )/3._DP
  !
  t  = SQRT(grho) / (2._DP * fz * ks * rho)
  expe = EXP( - 2._DP * al * ec / (fz3 * be * be) )
  af = 2._DP * al / be * (1._DP / (expe-1._DP) )
  bfup = expe * (vc_up - ec) / fz3
  bfdw = expe * (vc_dn - ec) / fz3
  !
  y = af * t * t
  xy = (1._DP + y) / (1._DP + y + y * y)
  qy = y * y * (2._DP + y) / (1._DP + y + y * y)**2
  s1 = 1._DP + 2._DP * al / be * t * t * xy
  !
  h0 = fz3 * be * be / (2._DP * al) * LOG(s1)
  dh0up  = be * t * t * fz3 / s1 * ( - 7._DP / 3._DP * xy - qy *   &
           (af * bfup / be-7._DP / 3._DP) )
  dh0dw  = be * t * t * fz3 / s1 * ( - 7._DP / 3._DP * xy - qy *   &
           (af * bfdw / be-7._DP / 3._DP) )
  dh0zup = (3._DP * h0 / fz - be * t * t * fz2 / s1 * (2._DP *xy - &
           qy * (3._DP * af * expe * ec / fz3 / be+2._DP) ) ) *    &
           dfz * (1._DP - zeta )
  dh0zdw = -(3._DP * h0 / fz - be * t * t * fz3 / s1 * (2._DP*xy - &
           qy * (3._DP * af * expe * ec / fz3 / be + 2._DP) ) ) *  &
           dfz * (1._DP + zeta )
  ddh0   = be * fz / (2._DP * ks * ks * rho) * (xy - qy) / s1
  !
  ee = -100._DP * fz4 * (ks / kf * t)**2
  cna  = cxc0 + pa * rs + pb * rs2
  dcna = pa * rs + 2._DP * pb * rs2
  cnb  = 1._DP + pc * rs + pd * rs2 + 1.D4 * pb * rs3
  dcnb = pc * rs + 2._DP * pd * rs2 + 3.D4 * pb * rs3
  cn  = cna / cnb - cx
  dcn = dcna / cnb - cna * dcnb / (cnb * cnb)
  h1  = nu * (cn - cc0 - 3._DP / 7._DP * cx) * fz3 * t * t * EXP(ee)
  dh1 = - third * (h1 * (7._DP + 8._DP * ee) + fz3 * nu * t * t * &
        EXP(ee) * dcn)
  !
  ddh1 = 2._DP * h1 * (1._DP + ee) * rho / grho
  dh1zup =  (1._DP - zeta) * dfz * h1 * (1._DP + 2._DP * ee / fz)
  dh1zdw = -(1._DP + zeta) * dfz * h1 * (1._DP + 2._DP * ee / fz)
  !
  sc     = rho * (h0 + h1)
  v1c_up = h0 + h1 + dh0up + dh1 + dh0zup + dh1zup
  v1c_dw = h0 + h1 + dh0up + dh1 + dh0zdw + dh1zdw
  v2c    = ddh0 + ddh1
  !
  !
  RETURN
  !
END SUBROUTINE ggac_spin
!
!
!-------------------------------------------------------------------
SUBROUTINE pbec_spin( rho, zeta, grho, iflag, sc, v1c_up, v1c_dw, v2c )
  !-----------------------------------------------------------------
  !! PBE correlation (without LDA part) - spin-polarized.
  !
  !! * iflag = 1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996);
  !! * iflag = 2: J.P.Perdew et al., PRL 100, 136406 (2008).
  !
  USE kind_l, ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  INTEGER, INTENT(IN) :: iflag
  !! see main comments
  REAL(DP), INTENT(IN) :: rho
  !! the total charge density
  REAL(DP), INTENT(IN) :: zeta
  !! the magnetization
  REAL(DP), INTENT(IN) :: grho
  !! the gradient of the charge squared
  REAL(DP), INTENT(OUT) :: sc
  !! correlation energies
  REAL(DP), INTENT(OUT) :: v1c_up, v1c_dw
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT) :: v2c
  !! derivatives of correlation wr. grho
  !
  ! ... local variables
  !
  REAL(DP) :: rs, ec, vc_up, vc_dn
  !
  REAL(DP), PARAMETER :: ga=0.031091_DP
  REAL(DP) :: be(2)
  DATA be / 0.06672455060314922_DP, 0.046_DP /
  REAL(DP), PARAMETER :: third=1.D0/3.D0, pi34=0.6203504908994_DP
  !                                         pi34=(3/4pi)^(1/3)
  REAL(DP), PARAMETER :: xkf=1.919158292677513_DP, xks=1.128379167095513_DP
  !                      xkf=(9 pi/4)^(1/3)      , xks=sqrt(4/pi)
  REAL(DP) :: kf, ks, t, expe, af, y, xy, qy, s1, h0, ddh0
  REAL(DP) :: fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, &
              dh0zup, dh0zdw
  !
  !
  rs = pi34 / rho**third
  !
  CALL pw_spin( rs, zeta, ec, vc_up, vc_dn )
  !
  kf = xkf / rs
  ks = xks * SQRT(kf)
  !
  fz = 0.5_DP*( (1._DP+zeta)**(2._DP/3._DP) + (1._DP-zeta)**(2._DP/3._DP) )
  fz2 = fz * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1._DP+zeta)**(-1._DP/3._DP) - (1._DP - zeta)**(-1._DP/3._DP) ) &
         / 3._DP
  !
  t  = SQRT(grho) / (2._DP * fz * ks * rho)
  expe = EXP( - ec / (fz3 * ga) )
  af   = be(iflag) / ga * (1._DP / (expe-1._DP) )
  bfup = expe * (vc_up - ec) / fz3
  bfdw = expe * (vc_dn - ec) / fz3
  !
  y  = af * t * t
  xy = (1._DP + y) / (1._DP + y + y * y)
  qy = y * y * (2._DP + y) / (1._DP + y + y * y)**2
  s1 = 1._DP + be(iflag) / ga * t * t * xy
  !
  h0 = fz3 * ga * LOG(s1)
  !
  dh0up = be(iflag) * t * t * fz3 / s1 * ( -7._DP/3._DP * xy - qy * &
          (af * bfup / be(iflag)-7._DP/3._DP) )
  !
  dh0dw = be(iflag) * t * t * fz3 / s1 * ( -7._DP/3._DP * xy - qy * &
          (af * bfdw / be(iflag)-7._DP/3._DP) )
  !
  dh0zup =   (3._DP * h0 / fz - be(iflag) * t * t * fz2 / s1 *  &
             (2._DP * xy - qy * (3._DP * af * expe * ec / fz3 / &
             be(iflag)+2._DP) ) ) * dfz * (1._DP - zeta)
  !
  dh0zdw = - (3._DP * h0 / fz - be(iflag) * t * t * fz2 / s1 *  &
             (2._DP * xy - qy * (3._DP * af * expe * ec / fz3 / &
             be(iflag)+2._DP) ) ) * dfz * (1._DP + zeta)
  !
  ddh0 = be(iflag) * fz / (2._DP * ks * ks * rho) * (xy - qy) / s1
  !
  sc     = rho * h0
  v1c_up = h0 + dh0up + dh0zup
  v1c_dw = h0 + dh0dw + dh0zdw
  v2c    = ddh0
  !
  !
  RETURN
  !
END SUBROUTINE pbec_spin
!
!
!------------------------------------------------------------------------
SUBROUTINE lsd_glyp( rho_in_up, rho_in_dw, grho_up, grho_dw, grho_ud, sc, &
                     v1c_up, v1c_dw, v2c_up, v2c_dw, v2c_ud )
  !----------------------------------------------------------------------
  !! Lee, Yang, Parr: gradient correction part.
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho_in_up, rho_in_dw
  !! the total charge density
  REAL(DP), INTENT(IN) :: grho_up, grho_dw
  !! the gradient of the charge squared
  REAL(DP), INTENT(IN) :: grho_ud
  !! gradient off-diagonal term up-down
  REAL(DP), INTENT(OUT) :: sc
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1c_up, v1c_dw
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT) :: v2c_up, v2c_dw
  !! derivatives of correlation wr. grho
  REAL(DP), INTENT(OUT) :: v2c_ud
  !! derivative of correlation wr. grho, off-diag. term
  !
  ! ... local variables
  !
  REAL(DP) :: ra, rb, rho, grhoaa, grhoab, grhobb
  REAL(DP) :: rm3, dr, or, dor, der, dder
  REAL(DP) :: dlaa, dlab, dlbb, dlaaa, dlaab, dlaba, dlabb, dlbba, dlbbb
  REAL(DP), PARAMETER :: a=0.04918_DP, b=0.132_DP, c=0.2533_DP, d=0.349_DP
  !
  ra = rho_in_up
  rb = rho_in_dw
  rho = ra + rb
  rm3 = rho**(-1._DP/3._DP)
  !
  dr = ( 1._DP + d*rm3 )
  or = EXP(-c*rm3) / dr * rm3**11._DP
  dor = -1._DP/3._DP * rm3**4 * or * (11._DP/rm3-c-d/dr)
  !
  der  = c*rm3 + d*rm3/dr
  dder = 1._DP/3._DP * (d*d*rm3**5/dr/dr - der/rho)
  !
  dlaa = -a*b*or*( ra*rb/9._DP*(1._DP-3*der-(der-11._DP)*ra/rho)-rb*rb )
  dlab = -a*b*or*( ra*rb/9._DP*(47._DP-7._DP*der)-4._DP/3._DP*rho*rho  )
  dlbb = -a*b*or*( ra*rb/9._DP*(1._DP-3*der-(der-11._DP)*rb/rho)-ra*ra )
  !
  dlaaa = dor/or*dlaa-a*b*or*(rb/9._DP*(1._DP-3*der-(der-11._DP)*ra/rho)-     &
          ra*rb/9._DP*((3._DP+ra/rho)*dder+(der-11._DP)*rb/rho/rho))
  dlaab = dor/or*dlaa-a*b*or*(ra/9._DP*(1._DP-3._DP*der-(der-11._DP)*ra/rho)- &
          ra*rb/9._DP*((3._DP+ra/rho)*dder-(der-11._DP)*ra/rho/rho)-2._DP*rb)
  dlaba = dor/or*dlab-a*b*or*(rb/9._DP*(47._DP-7._DP*der)-7._DP/9._DP*ra*rb*dder- &
          8._DP/3._DP*rho)
  dlabb = dor/or*dlab-a*b*or*(ra/9._DP*(47._DP-7._DP*der)-7._DP/9._DP*ra*rb*dder- &
          8._DP/3._DP*rho)
  dlbba = dor/or*dlbb-a*b*or*(rb/9._DP*(1._DP-3._DP*der-(der-11._DP)*rb/rho)- &
          ra*rb/9._DP*((3._DP+rb/rho)*dder-(der-11._DP)*rb/rho/rho)-2._DP*ra)
  dlbbb = dor/or*dlbb-a*b*or*(ra/9._DP*(1._DP-3*der-(der-11._DP)*rb/rho)-     &
          ra*rb/9._DP*((3._DP+rb/rho)*dder+(der-11._DP)*ra/rho/rho))
  !
  grhoaa = grho_up
  grhobb = grho_dw
  grhoab = grho_ud
  !
  sc     = dlaa *grhoaa + dlab *grhoab + dlbb *grhobb
  v1c_up = dlaaa*grhoaa + dlaba*grhoab + dlbba*grhobb
  v1c_dw = dlaab*grhoaa + dlabb*grhoab + dlbbb*grhobb
  v2c_up = 2._DP*dlaa
  v2c_dw = 2._DP*dlbb
  v2c_ud = dlab
  !
  !
  RETURN
  !
END SUBROUTINE lsd_glyp
!
!
!---------------------------------------------------------------
SUBROUTINE cpbe2d( rho, grho, sc, v1c, v2c )
  !---------------------------------------------------------------
  !! 2D correction (last term of Eq. 5, PRL 108, 126402 (2012))
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN)  :: rho, grho
  REAL(DP), INTENT(OUT) :: sc, v1c, v2c
  !
  ! ... local variables
  !
  REAL(8), PARAMETER :: pi=3.14159265358979323846d0
  REAL(DP), PARAMETER :: ex1=0.333333333333333333_DP, ex2=1.166666666666667_DP
  REAL(DP), PARAMETER :: ex3=ex2+1.0_DP
  REAL(DP) :: fac1, fac2, zeta, phi, gr, rs, drsdn, akf, aks, t, dtdn, dtdgr
  REAL(DP) :: p, a, g, alpha1, beta1,beta2,beta3,beta4, dgdrs, epsc, depscdrs
  REAL(DP) :: c, gamma1, beta, aa, cg, h, dhdaa, dhdt, dhdrs
  REAL(DP) :: epscpbe, depscpbedrs, depscpbedt, a0,a1,a2, b0,b1,b2, c0,c1,c2
  REAL(DP) :: e0,e1,e2, f0,f1,f2, g0,g1,g2, h0,h1,h2, d0,d1,d2, ff, dffdt
  REAL(DP) :: rs3d, rs2d, drs2ddrs3d, eps2d, deps2ddrs2, depsGGAdrs, depsGGAdt
  REAL(DP) :: drs2ddt, rs2, ec, decdn, decdgr, daadepsc
  !
  fac1 = (3.d0*pi*pi)**ex1
  fac2 = SQRT(4.d0*fac1/pi)
  !
  zeta = 0.d0
  phi  = 1.d0
  !
  gr = SQRT(grho)
  !
  rs = (3.d0/4.d0/pi/rho)**ex1
  drsdn = -DBLE(3 ** (0.1D1 / 0.3D1)) * DBLE(2 ** (0.1D1 / 0.3D1)) * &
  0.3141592654D1 ** (-0.1D1 / 0.3D1) * (0.1D1 / rho) ** (-0.2D1 /    &
  0.3D1) / rho ** 2 / 0.6D1
  !
  akf  = (3.d0*pi*pi*rho)**(1.d0/3.d0)
  aks  = DSQRT(4.d0*akf/pi)
  t    = gr/2.d0 / phi / aks / rho
  dtdn = -7.d0/6.d0 * gr/2.d0 / phi/DSQRT(4.d0/pi)/   &
         ((3.d0*pi*pi)**(1.d0/6.d0)) / (rho**(13.d0/6.d0))
  dtdgr = 1.d0/2.d0/phi/aks/rho
  !
  ! for the LDA correlation
  p = 1.d0
  A = 0.0310906908696548950_DP
  alpha1 = 0.21370d0
  beta1  = 7.5957d0
  beta2  = 3.5876d0
  beta3  = 1.6382d0
  beta4  = 0.49294d0
  G = -0.2D1 * A * DBLE(1 + alpha1 * rs) * LOG(0.1D1 + 0.1D1 / A / ( &
  beta1 * SQRT(DBLE(rs)) + DBLE(beta2 * rs) + DBLE(beta3 * rs ** (   &
  0.3D1 / 0.2D1)) + DBLE(beta4 * rs ** (p + 1))) / 0.2D1)
  !
  dGdrs = -0.2D1 * A * alpha1 * LOG(0.1D1 + 0.1D1 / A / (beta1 * SQRT(rs) &
   + beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + beta4 * rs **           &
  (p + 1)) / 0.2D1) + (0.1D1 + alpha1 * rs) / (beta1 * SQRT(rs) +         &
  beta2 * rs + beta3 * rs ** (0.3D1 / 0.2D1) + beta4 * rs ** (p + 1))     &
  ** 2 * (beta1 * rs ** (-0.1D1 / 0.2D1) / 0.2D1 + beta2 + 0.3D1 /        &
  0.2D1 * beta3 * SQRT(rs) + beta4 * rs ** (p + 1) * DBLE(p + 1) /        &
  rs) / (0.1D1 + 0.1D1 / A / (beta1 * SQRT(rs) + beta2 * rs + beta3 *     &
  rs ** (0.3D1 / 0.2D1) + beta4 * rs ** (p + 1)) / 0.2D1)
  !
  epsc = G
  depscdrs = dGdrs
  !
  ! PBE
  c = 1.d0
  gamma1 = 0.0310906908696548950_dp
  beta = 0.06672455060314922_dp
  !
  AA = beta / gamma1 / (EXP(-epsc / gamma1 / phi ** 3) - 0.1D1)
  cg = beta / gamma1 ** 2 / (EXP(-epsc/ gamma1 / phi ** 3) - 0.1D1) &
       ** 2 / phi ** 3 * EXP(-epsc / gamma1 / phi ** 3)
  dAAdepsc = cg
  !
  IF (t <= 10.d0) THEN
     !
     H = DBLE(gamma1) * phi ** 3 * LOG(DBLE(1 + beta / gamma1 * t ** 2   &
         * (1 + AA * t ** 2) / (1 + c * AA * t ** 2 + AA ** 2 * t ** 4)))
     !
     dHdAA = gamma1 * phi ** 3 * (beta / gamma1 * t ** 4 / (1 + c * AA   &
      * t ** 2 + AA ** 2 * t ** 4) - beta / gamma1 * t ** 2 * (1 + AA *  &
      t ** 2) / (1 + c * AA * t ** 2 + AA ** 2 * t ** 4) ** 2 * (c * t **&
      2 + 2 * AA * t ** 4)) / (1 + beta / gamma1 * t ** 2 * (1 + AA *    &
      t ** 2) / (1 + c * AA * t ** 2 + AA ** 2 * t ** 4))
     !
     dHdt = gamma1 * phi ** 3 * (2 * beta / gamma1 * t * (1 + AA * t **  &
      2) / (1 + c * AA * t ** 2 + AA ** 2 * t ** 4) + 2 * beta / gamma1  &
      * t ** 3 * AA / (1 + c * AA * t ** 2 + AA ** 2 * t ** 4) - beta /  &
      gamma1 * t ** 2 * (1 + AA * t ** 2) / (1 + c * AA * t ** 2 + AA ** &
      2 * t ** 4) ** 2 * (2 * c * AA * t + 4 * AA ** 2 * t ** 3)) / (1   &
      + beta / gamma1 * t ** 2 * (1 + AA * t ** 2) / (1 + c * AA * t **  &
      2 + AA ** 2 * t ** 4))
     !
  ELSE
     !
     H = gamma1 * (phi**3) * DLOG(1.d0+(beta/gamma1)*(1.d0/AA))
     !
     dHdAA = gamma1 * (phi**3) * 1.d0/(1.d0+(beta/gamma1)*(1.d0/AA))* &
             (beta/gamma1) * (-1.d0/AA/AA)
     !
     dHdt = 0.d0
     !
  ENDIF
  !
  dHdrs = dHdAA*dAAdepsc*depscdrs
  !
  epscPBE     = epsc+H
  depscPBEdrs = depscdrs+dHdrs
  depscPBEdt  = dHdt
  !
  ! START THE 2D CORRECTION
  !
  beta = 1.3386d0
  a0 = -0.1925d0
  a1 =  0.117331d0
  a2 =  0.0234188d0
  b0 =  0.0863136d0
  b1 = -0.03394d0
  b2 = -0.037093d0
  c0 =  0.057234d0
  c1 = -0.00766765d0
  c2 =  0.0163618d0
  e0 =  1.0022d0
  e1 =  0.4133d0
  e2 =  1.424301d0
  f0 = -0.02069d0
  f1 =  0.d0
  f2 =  0.d0
  g0 =  0.340d0
  g1 =  0.0668467d0
  g2 =  0.d0
  h0 =  0.01747d0
  h1 =  0.0007799d0
  h2 =  1.163099d0
  d0 = -a0*h0
  d1 = -a1*h1
  d2 = -a2*h2
  !
  ff = t ** 4 * (1 + t ** 2) / (1000000 + t ** 6)
  dffdt = 4 * t ** 3 * (1 + t ** 2) / (1000000 + t ** 6) + 2 * t ** &
  5 / (1000000 + t ** 6) - 6 * t ** 9 * (1 + t ** 2) / (1000000 + t &
  ** 6) ** 2
!
  rs3d=rs
  rs2d = 0.4552100000D0 * DBLE(3 ** (0.7D1 / 0.12D2)) * DBLE(4 ** ( &
  0.5D1 / 0.12D2)) * (0.1D1 / pi) ** (-0.5D1 / 0.12D2) * rs3d ** (  &
  0.5D1 / 0.4D1) * SQRT(t)

  cg = 0.5690125000D0 * DBLE(3 ** (0.7D1 / 0.12D2)) * DBLE(4 ** (       &
  0.5D1 / 0.12D2)) * (0.1D1 / pi) ** (-0.5D1 / 0.12D2) * rs3d ** (0.1D1 &
   / 0.4D1) * SQRT(t)
  drs2ddrs3d=cg

  cg = 0.2276050000D0 * DBLE(3 ** (0.7D1 / 0.12D2)) * DBLE(4 ** (       &
  0.5D1 / 0.12D2)) * DBLE((1 / pi) ** (-0.5D1 / 0.12D2)) * DBLE(rs3d ** &
   (0.5D1 / 0.4D1)) * DBLE(t ** (-0.1D1 / 0.2D1))
  drs2ddt=cg
  rs2=rs2d
  !
  eps2d = (EXP(-beta * rs2) - 0.1D1) * (-0.2D1 / 0.3D1 * SQRT(0.2D1)   &
   * DBLE((1 + zeta) ** (0.3D1 / 0.2D1) + (1 - zeta) ** (0.3D1 /       &
  0.2D1)) / pi / rs2 + 0.4D1 / 0.3D1 * (0.1D1 + 0.3D1 / 0.8D1 * DBLE(  &
  zeta ** 2) + 0.3D1 / 0.128D3 * DBLE(zeta ** 4)) * SQRT(0.2D1) / pi / &
   rs2) + a0 + (b0 * rs2 + c0 * rs2 ** 2 + d0 * rs2 ** 3) * LOG(0.1D1  &
   + 0.1D1 / (e0 * rs2 + f0 * rs2 ** (0.3D1 / 0.2D1) + g0 * rs2 **     &
  2 + h0 * rs2 ** 3)) + (a1 + (b1 * rs2 + c1 * rs2 ** 2 + d1 * rs2 **  &
   3) * LOG(0.1D1 + 0.1D1 / (e1 * rs2 + f1 * rs2 ** (0.3D1 / 0.2D1)    &
   + g1 * rs2 ** 2 + h1 * rs2 ** 3))) * DBLE(zeta ** 2) + (a2 + (b2    &
  * rs2 + c2 * rs2 ** 2 + d2 * rs2 ** 3) * LOG(0.1D1 + 0.1D1 / (e2 *   &
   rs2 + f2 * rs2 ** (0.3D1 / 0.2D1) + g2 * rs2 ** 2 + h2 * rs2 ** 3   &
  ))) * DBLE(zeta ** 4)
  !
  cg = -beta * EXP(-beta * rs2) * (-0.2D1 / 0.3D1 * SQRT(0.2D1) *      &
  DBLE((1 + zeta) ** (0.3D1 / 0.2D1) + (1 - zeta) ** (0.3D1 / 0.2D1))  &
  / pi / rs2 + 0.4D1 / 0.3D1 * (0.1D1 + 0.3D1 / 0.8D1 * DBLE(zeta **   &
   2) + 0.3D1 / 0.128D3 * DBLE(zeta ** 4)) * SQRT(0.2D1) / pi / rs2)   &
   + (EXP(-beta * rs2) - 0.1D1) * (0.2D1 / 0.3D1 * SQRT(0.2D1) * DBLE  &
  ((1 + zeta) ** (0.3D1 / 0.2D1) + (1 - zeta) ** (0.3D1 / 0.2D1)) /    &
   pi / rs2 ** 2 - 0.4D1 / 0.3D1 * (0.1D1 + 0.3D1 / 0.8D1 * DBLE(zeta  &
   ** 2) + 0.3D1 / 0.128D3 * DBLE(zeta ** 4)) * SQRT(0.2D1) / pi /     &
  rs2 ** 2) + (b0 + 0.2D1 * c0 * rs2 + 0.3D1 * d0 * rs2 ** 2) * LOG(   &
  0.1D1 + 0.1D1 / (e0 * rs2 + f0 * rs2 ** (0.3D1 / 0.2D1) + g0 * rs2   &
   ** 2 + h0 * rs2 ** 3)) - (b0 * rs2 + c0 * rs2 ** 2 + d0 * rs2 **    &
  3) / (e0 * rs2 + f0 * rs2 ** (0.3D1 / 0.2D1) + g0 * rs2 ** 2 + h0    &
  * rs2 ** 3) ** 2 * (e0 + 0.3D1 / 0.2D1 * f0 * SQRT(rs2) + 0.2D1 *    &
  g0 * rs2 + 0.3D1 * h0 * rs2 ** 2) / (0.1D1 + 0.1D1 / (e0 * rs2 + f0  &
   * rs2 ** (0.3D1 / 0.2D1) + g0 * rs2 ** 2 + h0 * rs2 ** 3)) + ((     &
  b1 + 0.2D1 * c1 * rs2 + 0.3D1 * d1 * rs2 ** 2) * LOG(0.1D1 + 0.1D1   &
  / (e1 * rs2 + f1 * rs2 ** (0.3D1 / 0.2D1) + g1 * rs2 ** 2 + h1 *     &
  rs2 ** 3)) - (b1 * rs2 + c1 * rs2 ** 2 + d1 * rs2 ** 3) / (e1 * rs2  &
   + f1 * rs2 ** (0.3D1 / 0.2D1) + g1 * rs2 ** 2 + h1 * rs2 ** 3) **   &
   2 * (e1 + 0.3D1 / 0.2D1 * f1 * SQRT(rs2) + 0.2D1 * g1 * rs2 +       &
  0.3D1 * h1 * rs2 ** 2) / (0.1D1 + 0.1D1 / (e1 * rs2 + f1 * rs2 ** (  &
  0.3D1 / 0.2D1) + g1 * rs2 ** 2 + h1 * rs2 ** 3))) * DBLE(zeta ** 2)  &
  + ((b2 + 0.2D1 * c2 * rs2 + 0.3D1 * d2 * rs2 ** 2) * LOG(0.1D1 +     &
  0.1D1 / (e2 * rs2 + f2 * rs2 ** (0.3D1 / 0.2D1) + g2 * rs2 ** 2 + h2 &
   * rs2 ** 3)) - (b2 * rs2 + c2 * rs2 ** 2 + d2 * rs2 ** 3) / (e2     &
  * rs2 + f2 * rs2 ** (0.3D1 / 0.2D1) + g2 * rs2 ** 2 + h2 * rs2 **    &
  3) ** 2 * (e2 + 0.3D1 / 0.2D1 * f2 * SQRT(rs2) + 0.2D1 * g2 * rs2    &
  + 0.3D1 * h2 * rs2 ** 2) / (0.1D1 + 0.1D1 / (e2 * rs2 + f2 * rs2 **  &
   (0.3D1 / 0.2D1) + g2 * rs2 ** 2 + h2 * rs2 ** 3))) * DBLE(zeta **   &
   4)
  !
  deps2ddrs2=cg
  !
  ! GGA-2D
  !
  depsGGAdrs = ff*(-depscPBEdrs+deps2ddrs2*drs2ddrs3d)
  depsGGAdt  = dffdt*(-epscPBE+eps2d)+ff*      &
               (-depscPBEdt+deps2ddrs2*drs2ddt)
  !
  ec = rho*(ff*(-epscPBE+eps2d))
  !
  decdn = ff*(-epscPBE+eps2d)+rho*depsGGAdrs*drsdn+ &
          rho*depsGGAdt*dtdn
  !
  decdgr = rho*depsGGAdt*dtdgr
  !
  sc = ec
  v1c = decdn
  v2c = decdgr/gr
  !
  RETURN
  !
END SUBROUTINE cpbe2d
! !
END MODULE

