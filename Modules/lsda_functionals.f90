!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine becke88_spin (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! Becke exchange: A.D. Becke, PRA 38, 3098 (1988) - Spin polarized case
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP) :: rho, grho, sx, v1x, v2x
  ! input: charge
  ! input: gradient
  ! output: the up and down energies
  ! output: first part of the potential
  ! output: the second part of the potential
  !
  real(DP) :: beta, third
  parameter (beta = 0.0042d0, third = 1.d0 / 3.d0)
  real(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee
  !
  rho13 = rho**third
  rho43 = rho13**4
  xs = sqrt (grho) / rho43
  xs2 = xs * xs
  sa2b8 = sqrt (1.0d0 + xs2)
  shm1 = log (xs + sa2b8)
  dd = 1.0d0 + 6.0d0 * beta * xs * shm1
  dd2 = dd * dd
  ee = 6.0d0 * beta * xs2 / sa2b8 - 1.d0
  sx = grho / rho43 * ( - beta / dd)
  v1x = - (4.d0 / 3.d0) * xs2 * beta * rho13 * ee / dd2
  v2x = beta * (ee-dd) / (rho43 * dd2)
  !
  return
end subroutine becke88_spin
!
!-----------------------------------------------------------------------
subroutine perdew86_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  ! Perdew gradient correction on correlation: PRB 33, 8822 (1986)
  ! spin-polarized case
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
  real(DP) :: p1, p2, p3, p4, pc1, pc2, pci
  parameter (p1 = 0.023266d0, p2 = 7.389d-6, p3 = 8.723d0, p4 = &
       0.472d0)
  parameter (pc1 = 0.001667d0, pc2 = 0.002568d0, pci = pc1 + pc2)
  real(DP) :: third, pi34
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  ! pi34=(3/4pi)^(1/3)
  !
  real(DP) :: rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs
  real(DP) :: dcna, dcnb, dcn, phi, ephi, dd, ddd
  !
  rho13 = rho**third
  rho43 = rho13**4
  rs = pi34 / rho13
  rs2 = rs * rs
  rs3 = rs * rs2
  cna = pc2 + p1 * rs + p2 * rs2
  cnb = 1.d0 + p3 * rs + p4 * rs2 + 1.d4 * p2 * rs3
  cn = pc1 + cna / cnb
  drs = - third * pi34 / rho43
  dcna = (p1 + 2.d0 * p2 * rs) * drs
  dcnb = (p3 + 2.d0 * p4 * rs + 3.d4 * p2 * rs2) * drs
  dcn = dcna / cnb - cna / (cnb * cnb) * dcnb
  phi = 0.192d0 * pci / cn * sqrt (grho) * rho** ( - 7.d0 / 6.d0)
  !SdG: in the original paper 1.745*0.11=0.19195 is used
  dd = (2.d0) **third * sqrt ( ( (1.d0 + zeta) * 0.5d0) ** (5.d0 / &
       3.d0) + ( (1.d0 - zeta) * 0.5d0) ** (5.d0 / 3.d0) )
  ddd = (2.d0) ** ( - 4.d0 / 3.d0) * 5.d0 * ( ( (1.d0 + zeta) &
       * 0.5d0) ** (2.d0 / 3.d0) - ( (1.d0 - zeta) * 0.5d0) ** (2.d0 / &
       3.d0) ) / (3.d0 * dd)
  ephi = exp ( - phi)
  sc = grho / rho43 * cn * ephi / dd
  v1cup = sc * ( (1.d0 + phi) * dcn / cn - ( (4.d0 / 3.d0) - &
       (7.d0 / 6.d0) * phi) / rho) - sc * ddd / dd * (1.d0 - zeta) &
       / rho
  v1cdw = sc * ( (1.d0 + phi) * dcn / cn - ( (4.d0 / 3.d0) - &
       (7.d0 / 6.d0) * phi) / rho) + sc * ddd / dd * (1.d0 + zeta) &
       / rho
  v2c = cn * ephi / rho43 * (2.d0 - phi) / dd
  !
  return
end subroutine perdew86_spin
!
!-----------------------------------------------------------------------
subroutine ggac_spin(rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91) correlation part - spin-polarized
  !
  USE kinds, ONLY: DP
  implicit none
  real(DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
  !
  real(dp) :: rs, ec
  real(dp) :: vc(2)
  !
  real(DP) :: al, pa, pb, pc, pd, cx, cxc0, cc0
  parameter (al = 0.09d0, pa = 0.023266d0, pb = 7.389d-6, pc = &
       8.723d0, pd = 0.472d0)
  parameter (cx = - 0.001667d0, cxc0 = 0.002568d0, cc0 = - cx + &
       cxc0)
  real(DP) :: third, pi34, nu, be, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (nu = 15.755920349483144d0, be = nu * cc0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
  ! xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(DP) :: kf, ks, rs2, rs3, vcup, vcdw, t, expe, af, y, &
       xy, qy, s1, h0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, dh1, &
       ddh1, fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, &
       dh0zdw, dh1zup, dh1zdw
  !
  rs = pi34 / rho**third
  rs2 = rs * rs
  rs3 = rs * rs2
  !
  call pw_spin( rs, zeta, ec, vc )
  !
  kf = xkf / rs
  ks = xks * sqrt(kf)
  fz = 0.5d0 * ( (1.d0 + zeta) ** (2.d0 / 3.d0) + (1.d0 - zeta) ** ( &
       2.d0 / 3.d0) )
  fz2 = fz * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1.d0 + zeta) ** ( - 1.d0 / 3.d0) - (1.d0 - zeta) ** ( - &
       1.d0 / 3.d0) ) / 3.d0
  t = sqrt(grho) / (2.d0 * fz * ks * rho)
  expe = exp( - 2.d0 * al * ec / (fz3 * be * be) )
  af = 2.d0 * al / be * (1.d0 / (expe-1.d0) )
  bfup = expe * (vc(1) - ec) / fz3
  bfdw = expe * (vc(2) - ec) / fz3
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y)**2
  s1 = 1.d0 + 2.d0 * al / be * t * t * xy
  h0 = fz3 * be * be / (2.d0 * al) * log(s1)
  dh0up = be * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfup / be-7.d0 / 3.d0) )
  dh0dw = be * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfdw / be-7.d0 / 3.d0) )
  dh0zup = (3.d0 * h0 / fz - be * t * t * fz2 / s1 * (2.d0 * xy - &
       qy * (3.d0 * af * expe * ec / fz3 / be+2.d0) ) ) * dfz * (1.d0 - &
       zeta)
  dh0zdw = - (3.d0 * h0 / fz - be * t * t * fz3 / s1 * (2.d0 * xy - &
       qy * (3.d0 * af * expe * ec / fz3 / be+2.d0) ) ) * dfz * (1.d0 + &
       zeta)
  ddh0 = be * fz / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  ee = - 100.d0 * fz4 * (ks / kf * t) **2
  cna = cxc0 + pa * rs + pb * rs2
  dcna = pa * rs + 2.d0 * pb * rs2
  cnb = 1.d0 + pc * rs + pd * rs2 + 1.d4 * pb * rs3
  dcnb = pc * rs + 2.d0 * pd * rs2 + 3.d4 * pb * rs3
  cn = cna / cnb - cx
  dcn = dcna / cnb - cna * dcnb / (cnb * cnb)
  h1 = nu * (cn - cc0 - 3.d0 / 7.d0 * cx) * fz3 * t * t * exp(ee)
  dh1 = - third * (h1 * (7.d0 + 8.d0 * ee) + fz3 * nu * t * t * exp &
       (ee) * dcn)
  ddh1 = 2.d0 * h1 * (1.d0 + ee) * rho / grho
  dh1zup = (1.d0 - zeta) * dfz * h1 * (1.d0 + 2.d0 * ee / fz)
  dh1zdw = - (1.d0 + zeta) * dfz * h1 * (1.d0 + 2.d0 * ee / fz)
  sc = rho * (h0 + h1)
  v1cup = h0 + h1 + dh0up + dh1 + dh0zup + dh1zup
  v1cdw = h0 + h1 + dh0up + dh1 + dh0zdw + dh1zdw
  v2c = ddh0 + ddh1
  return
end subroutine ggac_spin
!
!---------------------------------------------------------------
subroutine pbec_spin (rho, zeta, grho, iflag, sc, v1cup, v1cdw, v2c)
  !---------------------------------------------------------------
  !
  ! PBE correlation (without LDA part) - spin-polarized
  ! iflag = 1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  ! iflag = 2: J.P.Perdew et al., PRL 100, 136406 (2008)
  !
  USE kinds, ONLY : DP
  implicit none
  integer, intent(in) :: iflag
  real(DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
  !
  real(dp) :: rs, ec
  real(dp) :: vc(2)
  !
  real(DP) :: ga, be(2)
  parameter (ga = 0.031091d0)
  data be / 0.06672455060314922_dp,  0.046_dp /
  real(DP) :: third, pi34, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(DP) :: kf, ks, t, expe, af, y, xy, qy, &
       s1, h0, ddh0
  real(DP) :: fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, &
       dh0zup, dh0zdw
  !
  rs = pi34 / rho**third
  !
  call pw_spin( rs, zeta, ec, vc )
  !
  kf = xkf / rs
  ks = xks * sqrt(kf)
  fz = 0.5d0 * ( (1.d0 + zeta) ** (2.d0 / 3.d0) + (1.d0 - zeta) ** ( &
       2.d0 / 3.d0) )
  fz2 = fz * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1.d0 + zeta) ** ( - 1.d0 / 3.d0) - (1.d0 - zeta) ** ( - &
       1.d0 / 3.d0) ) / 3.d0
  t = sqrt(grho) / (2.d0 * fz * ks * rho)
  expe = exp( - ec / (fz3 * ga) )
  af = be(iflag) / ga * (1.d0 / (expe-1.d0) )
  bfup = expe * (vc(1) - ec) / fz3
  bfdw = expe * (vc(2) - ec) / fz3
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + be(iflag) / ga * t * t * xy
  h0 = fz3 * ga * log(s1)
  dh0up = be(iflag) * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfup / be(iflag)-7.d0 / 3.d0) )
  dh0dw = be(iflag) * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfdw / be(iflag)-7.d0 / 3.d0) )
  dh0zup = (3.d0 * h0 / fz - be(iflag) * t * t * fz2 / s1 * (2.d0 * xy - &
  qy * (3.d0 * af * expe * ec / fz3 / be(iflag)+2.d0) ) ) * dfz * (1.d0 - zeta)
  dh0zdw = - (3.d0 * h0 / fz - be(iflag) * t * t * fz2 / s1 * (2.d0 * xy - &
  qy * (3.d0 * af * expe * ec / fz3 / be(iflag)+2.d0) ) ) * dfz * (1.d0 + zeta)
  !
  ddh0 = be(iflag) * fz / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  sc = rho * h0
  v1cup = h0 + dh0up + dh0zup
  v1cdw = h0 + dh0dw + dh0zdw
  v2c = ddh0
  return
end subroutine pbec_spin
!
