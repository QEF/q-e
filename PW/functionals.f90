!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine xc (rho, ex, ec, vx, vc)
  !-----------------------------------------------------------------------
  !     lda exchange and correlation functionals - Hartree a.u.
  !
  !     exchange   :  Slater (alpha=2/3)
  !     correlation:  Ceperley-Alder (Perdew-Zunger parameters)
  !                   Vosko-Wilk-Nusair
  !                   Lee-Yang-Parr
  !                   Perdew-Wang
  !                   Wigner
  !                   Hedin-Lundqvist
  !                   Ortiz-Ballone (Perdew-Zunger formula)
  !                   Ortiz-Ballone (Perdew-Wang formula)
  !                   Gunnarsson-Lundqvist
  !
  !     input : rho=rho(r)
  !     definitions: E_x = \int E_x(rho) dr, E_x(rho) = rho\epsilon_c(rho)
  !                  same for correlation
  !     output: ex = \epsilon_x(rho) ( NOT E_x(rho) )
  !             vx = dE_x(rho)/drho  ( NOT d\epsilon_x(rho)/drho )
  !             ec, vc as above for correlation
  !
  use parameters
  use funct
  implicit none

  real(kind=DP) :: rho, ec, vc, ex, vx
  !
  real(kind=DP) :: small, pi34, third
  parameter (small = 1.d-10)
  parameter (pi34 = 0.6203504908994d0, third = 1.d0 / 3.d0)
  ! pi34=(3/4pi)^(1/3)
  real(kind=DP) :: rs
  !
  if (rho.le.small) then
     ec = 0.0d0
     vc = 0.0d0
     ex = 0.0d0
     vx = 0.0d0
     return
  else
     rs = pi34 / rho**third
     ! rs as in the theory of metals: rs=(3/(4pi rho))^(1/3)
  endif
  !..exchange
  if (iexch.eq.1) then
     call slater (rs, ex, vx)
  else
     ex = 0.0d0
     vx = 0.0d0
  endif
  !..correlation
  if (icorr.eq.1) then
     call pz (rs, 1, ec, vc)
  elseif (icorr.eq.2) then
     call vwn (rs, ec, vc)
  elseif (icorr.eq.3) then
     call lyp (rs, ec, vc)
  elseif (icorr.eq.4) then
     call pw (rs, 1, ec, vc)
  elseif (icorr.eq.5) then
     call wigner (rs, ec, vc)
  elseif (icorr.eq.6) then
     call hl (rs, ec, vc)
  elseif (icorr.eq.7) then
     call pz (rs, 2, ec, vc)
  elseif (icorr.eq.8) then
     call pw (rs, 2, ec, vc)
  elseif (icorr.eq.9) then
     call gl (rs, ec, vc)
  else
     ec = 0.0d0
     vc = 0.0d0
  endif
  !
  return
end subroutine xc
!
!-----------------------------------------------------------------------
subroutine gcxc (rho, grho, sx, sc, v1x, v2x, v1c, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     exchange  :  Becke88
  !                  GGA (Generalized Gradient Approximation), PW91
  !                  PBE
  !                  revPBE
  !     correlation: Perdew86
  !                  GGA (PW91)
  !                  Lee-Yang-Parr
  !                  PBE
  !
  !     input:  rho, grho=|\nabla rho|^2
  !     definition:  E_x = \int E_x(rho,grho) dr
  !     output: sx = E_x(rho,grho)
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             sc, v1c, v2c as above for correlation
  !
  use funct
  use parameters
  implicit none

  real(kind=DP) :: rho, grho, sx, sc, v1x, v2x, v1c, v2c
  real(kind=DP) :: small

  parameter (small = 1.d-10)
  ! exchange
  if (rho.le.small) then
     sx = 0.0d0
     v1x = 0.0d0
     v2x = 0.0d0
  elseif (igcx.eq.1) then
     call becke88 (rho, grho, sx, v1x, v2x)
  elseif (igcx.eq.2) then
     call ggax (rho, grho, sx, v1x, v2x)
  elseif (igcx.eq.3) then
     call pbex (rho, grho, 1, sx, v1x, v2x)
  elseif (igcx.eq.4) then
     call pbex (rho, grho, 2, sx, v1x, v2x)
  else
     sx = 0.0d0
     v1x = 0.0d0
     v2x = 0.0d0
  endif
  ! correlation
  if (rho.le.small) then
     sc = 0.0d0
     v1c = 0.0d0
     v2c = 0.0d0
  elseif (igcc.eq.1) then
     call perdew86 (rho, grho, sc, v1c, v2c)
  elseif (igcc.eq.2) then
     call ggac (rho, grho, sc, v1c, v2c)
  elseif (igcc.eq.3) then
     call glyp (rho, grho, sc, v1c, v2c)
  elseif (igcc.eq.4) then
     call pbec (rho, grho, sc, v1c, v2c)
  else
     sc = 0.0d0
     v1c = 0.0d0
     v2c = 0.0d0
  endif
  !
  return
end subroutine gcxc
!
!-----------------------------------------------------------------------
subroutine slater (rs, ex, vx)
  !-----------------------------------------------------------------------
  !        Slater exchange with alpha=2/3
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, ex, vx
  real(kind=DP) :: f, alpha
  parameter (f = - 0.687247939924714d0)
  ! f = -9/8*(3/2pi)^(2/3)
  parameter (alpha = 2.0d0 / 3.0d0)
  !
  ex = f * alpha / rs
  vx = 4.d0 / 3.d0 * f * alpha / rs
  !
  return
end subroutine slater
!
!-----------------------------------------------------------------------
subroutine pz (rs, iflag, ec, vc)
  !-----------------------------------------------------------------------
  !     LDA parameterization form Monte Carlo data
  !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, ec, vc
  integer :: iflag
  !
  real(kind=DP) :: a (2), b (2), c (2), d (2), gc (2), b1 (2), b2 (2)
  real(kind=DP) :: lnrs, rs12, ox, dox
  !
  data a / 0.0311d0, 0.031091d0 /, b / - 0.048d0, - 0.046644d0 /, &
       c / 0.0020d0, 0.00419d0 /, d / - 0.0116d0, - 0.00983d0 /
  data gc / - 0.1423d0, - 0.103756d0 /, b1 / 1.0529d0, 0.56371d0 /, &
       b2 / 0.3334d0, 0.27358d0 /
  !
  if (rs.lt.1.0d0) then
     ! high density formula
     lnrs = log (rs)
     ec = a (iflag) * lnrs + b (iflag) + c (iflag) * rs * lnrs + d ( &
          iflag) * rs
     vc = a (iflag) * lnrs + (b (iflag) - a (iflag) / 3.d0) + 2.d0 / &
          3.d0 * c (iflag) * rs * lnrs + (2.d0 * d (iflag) - c (iflag) ) &
          / 3.d0 * rs
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     ox = 1.d0 + b1 (iflag) * rs12 + b2 (iflag) * rs
     dox = 1.d0 + 7.d0 / 6.d0 * b1 (iflag) * rs12 + 4.d0 / 3.d0 * &
          b2 (iflag) * rs
     ec = gc (iflag) / ox
     vc = ec * dox / ox
  endif
  !
  return
end subroutine pz
!
!-----------------------------------------------------------------------
subroutine vwn (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, ec, vc
  real(kind=DP) :: a, b, c, x0
  parameter (a = 0.0310907, b = 3.72744, c = 12.9352, x0 = - &
       0.10498)
  real(kind=DP) :: q, f1, f2, f3, rs12, fx, qx, tx, tt
  !
  q = sqrt (4.d0 * c - b * b)
  f1 = 2.d0 * b / q
  f2 = b * x0 / (x0 * x0 + b * x0 + c)
  f3 = 2.d0 * (2.d0 * x0 + b) / q
  rs12 = sqrt (rs)
  fx = rs + b * rs12 + c
  qx = atan (q / (2.d0 * rs12 + b) )
  ec = a * (log (rs / fx) + f1 * qx - f2 * (log ( (rs12 - x0) **2 / &
       fx) + f3 * qx) )
  tx = 2.d0 * rs12 + b
  tt = tx * tx + q * q
  vc = ec - rs12 * a / 6.d0 * (2.d0 / rs12 - tx / fx - 4.d0 * b / &
       tt - f2 * (2.d0 / (rs12 - x0) - tx / fx - 4.d0 * (2.d0 * x0 + b) &
       / tt) )
  !
  return
end subroutine vwn
!-----------------------------------------------------------------------
subroutine lyp (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     C. Lee, W. Yang, and R.G. Parr, PRB 37, 785 (1988)
  !     LDA part only
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, ec, vc
  real(kind=DP) :: a, b, c, d, pi43
  parameter (a = 0.04918d0, b = 0.132d0 * 2.87123400018819108d0)
  ! pi43 = (4pi/3)^(1/3)
  parameter (pi43 = 1.61199195401647d0, c = 0.2533d0 * pi43, d = &
       0.349d0 * pi43)
  real(kind=DP) :: ecrs, ox
  !
  ecrs = b * exp ( - c * rs)
  ox = 1.d0 / (1.d0 + d * rs)
  ec = - a * ox * (1.d0 + ecrs)
  vc = ec - rs / 3.d0 * a * ox * (d * ox + ecrs * (d * ox + c) )
  !
  return
end subroutine lyp
!
!-----------------------------------------------------------------------
subroutine pw (rs, iflag, ec, vc)
  !-----------------------------------------------------------------------
  !     iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, ec, vc
  integer :: iflag
  !
  real(kind=DP) :: a, b1, b2, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, b1 = 7.5957d0, b2 = 3.5876d0, c0 = a, &
       c1 = 0.046644d0, c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, &
       d1 = 1.4408d0)
  real(kind=DP) :: lnrs, rs12, rs32, rs2, om, dom, olog
  real(kind=DP) :: a1 (2), b3 (2), b4 (2)
  data a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, - 0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /
  !
  ! high- and low-density formulae implemented but not used in PW case
  ! (reason: inconsistencies in PBE/PW91 functionals)
  !
  if (rs.lt.1d0.and.iflag.eq.2) then
     ! high density formula
     lnrs = log (rs)
     ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs
     vc = c0 * lnrs - (c1 + c0 / 3.d0) + 2.d0 / 3.d0 * c2 * rs * &
          lnrs - (2.d0 * c3 + c2) / 3.d0 * rs
  elseif (rs.gt.100.d0.and.iflag.eq.2) then
     ! low density formula
     ec = - d0 / rs + d1 / rs**1.5d0
     vc = - 4.d0 / 3.d0 * d0 / rs + 1.5d0 * d1 / rs**1.5d0
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     rs32 = rs * rs12
     rs2 = rs**2
     om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 (iflag) * rs32 + b4 ( &
          iflag) * rs2)
     dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 ( &
          iflag) * rs32 + 2.d0 * b4 (iflag) * rs2)
     olog = log (1.d0 + 1.0d0 / om)
     ec = - 2.d0 * a * (1.d0 + a1 (iflag) * rs) * olog
     vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 (iflag) * rs) &
          * olog - 2.d0 / 3.d0 * a * (1.d0 + a1 (iflag) * rs) * dom / &
          (om * (om + 1.d0) )
  endif
  !
  return
end subroutine pw
!
!-----------------------------------------------------------------------
subroutine wigner (rs, ec, vc)
  !-----------------------------------------------------------------------
  !        Wigner correlation
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, ec, vc
  real(kind=DP) :: pi34, rho13
  parameter (pi34 = 0.6203504908994d0)
  ! pi34=(3/4pi)^(1/3), rho13=rho^(1/3)
  !
  rho13 = pi34 / rs
  vc = - rho13 * ( (0.943656d0 + 8.8963d0 * rho13) / (1.d0 + &
       12.57d0 * rho13) **2)
  ec = - 0.738d0 * rho13 * (0.959d0 / (1.d0 + 12.57d0 * rho13) )
  !
  return
end subroutine wigner
!
!-----------------------------------------------------------------------
subroutine hl (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     L. Hedin and  B.I. Lundqvist,  J. Phys. C 4, 2064 (1971)
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, ec, vc
  real(kind=DP) :: a, x
  !
  a = log (1.0d0 + 21.d0 / rs)
  x = rs / 21.0d0
  ec = a + (x**3 * a - x * x) + x / 2.d0 - 1.0d0 / 3.0d0
  ec = - 0.0225d0 * ec
  vc = - 0.0225d0 * a
  !
  return
end subroutine hl
!
!-----------------------------------------------------------------------
subroutine gl (rs, ec, vc)
  !-----------------------------------------------------------------------
  !  O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976)
  !
  use parameters
  implicit none
  real(kind=DP) :: rs, vc, ec
  real(kind=DP) :: c, r, x
  parameter (c = 0.0333, r = 11.4)
  ! c=0.0203, r=15.9 for the paramagnetic case
  !
  x = rs / r
  vc = - c * log (1.d0 + 1.d0 / x)
  ec = - c * ( (1.d0 + x**3) * log (1.d0 + 1.d0 / x) - 1.0d0 / &
       3.0d0 + x * (0.5d0 - x) )
  !
  return
end subroutine gl
!
!-----------------------------------------------------------------------
subroutine becke88 (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! Becke exchange: A.D. Becke, PRA 38, 3098 (1988)
  ! only gradient-corrected part, no Slater term included
  !
  use parameters
  implicit none
  real(kind=DP) :: rho, grho, sx, v1x, v2x
  real(kind=DP) :: beta, third, two13
  parameter (beta = 0.0042d0)
  parameter (third = 1.d0 / 3.d0, two13 = 1.259921049894873d0)
  ! two13 = 2^(1/3)
  real(kind=DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee
  !
  rho13 = rho**third
  rho43 = rho13**4
  xs = two13 * sqrt (grho) / rho43
  xs2 = xs * xs
  sa2b8 = sqrt (1.0d0 + xs2)
  shm1 = log (xs + sa2b8)
  dd = 1.0d0 + 6.0d0 * beta * xs * shm1
  dd2 = dd * dd
  ee = 6.0d0 * beta * xs2 / sa2b8 - 1.d0
  sx = two13 * grho / rho43 * ( - beta / dd)
  v1x = - (4.d0 / 3.d0) / two13 * xs2 * beta * rho13 * ee / dd2
  v2x = two13 * beta * (ee-dd) / (rho43 * dd2)
  !
  return
end subroutine becke88
!
!-----------------------------------------------------------------------
subroutine ggax (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91), exchange part:
  ! J.P. Perdew et al.,PRB 46, 6671 (1992)
  !
  use parameters
  implicit none
  real(kind=DP) :: rho, grho, sx, v1x, v2x
  real(kind=DP) :: f1, f2, f3, f4, f5
  parameter (f1 = 0.19645d0, f2 = 7.7956d0, f3 = 0.2743d0, f4 = &
       0.1508d0, f5 = 0.004d0)
  real(kind=DP) :: fp1, fp2
  parameter (fp1 = - 0.019292021296426d0, fp2 = 0.161620459673995d0)
  ! fp1 = -3/(16 pi)*(3 pi^2)^(-1/3)
  ! fp2 = (1/2)(3 pi^2)**(-1/3)
  real(kind=DP) :: rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das, &
       dbs, dls
  !
  rhom43 = rho** ( - 4.d0 / 3.d0)
  s = fp2 * sqrt (grho) * rhom43
  s2 = s * s
  s3 = s2 * s
  s4 = s2 * s2
  exps = f4 * exp ( - 100.d0 * s2)
  as = f3 - exps - f5 * s2
  sa2b8 = sqrt (1.0d0 + f2 * f2 * s2)
  shm1 = log (f2 * s + sa2b8)
  bs = 1.d0 + f1 * s * shm1 + f5 * s4
  das = (200.d0 * exps - 2.d0 * f5) * s
  dbs = f1 * (shm1 + f2 * s / sa2b8) + 4.d0 * f5 * s3
  dls = (das / as - dbs / bs)
  sx = fp1 * grho * rhom43 * as / bs
  v1x = - 4.d0 / 3.d0 * sx / rho * (1.d0 + s * dls)
  v2x = fp1 * rhom43 * as / bs * (2.d0 + s * dls)
  !
  return
end subroutine ggax
!
!-----------------------------------------------------------------------
subroutine perdew86 (rho, grho, sc, v1c, v2c)
  !-----------------------------------------------------------------------
  ! Perdew gradient correction on correlation: PRB 33, 8822 (1986)
  !
  use parameters
  implicit none
  real(kind=DP) :: rho, grho, sc, v1c, v2c
  real(kind=DP) :: p1, p2, p3, p4, pc1, pc2, pci
  parameter (p1 = 0.023266d0, p2 = 7.389d-6, p3 = 8.723d0, p4 = &
       0.472d0)
  parameter (pc1 = 0.001667d0, pc2 = 0.002568d0, pci = pc1 + pc2)
  real(kind=DP) :: third, pi34
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  ! pi34=(3/4pi)^(1/3)
  real(kind=DP) :: rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs
  real(kind=DP) :: dcna, dcnb, dcn, phi, ephi
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
  ! SdG: in the original paper 1.745*0.11=0.19195 is used
  ephi = exp ( - phi)
  sc = grho / rho43 * cn * ephi
  v1c = sc * ( (1.d0 + phi) * dcn / cn - ( (4.d0 / 3.d0) - (7.d0 / &
       6.d0) * phi) / rho)
  v2c = cn * ephi / rho43 * (2.d0 - phi)
  !
  return
end subroutine perdew86
!
!-----------------------------------------------------------------------
subroutine glyp (rho, grho, sc, v1c, v2c)
  !-----------------------------------------------------------------------
  ! Lee Yang Parr: gradient correction part
  !
  use parameters
  implicit none
  real(kind=DP) :: rho, grho, sc, v1c, v2c
  real(kind=DP) :: a, b, c, d
  parameter (a = 0.04918d0, b = 0.132d0, c = 0.2533d0, d = 0.349d0)
  real(kind=DP) :: rhom13, rhom43, rhom53, om, xl, ff, dom, dxl
  !
  rhom13 = rho** ( - 1.d0 / 3.d0)
  om = exp ( - c * rhom13) / (1.d0 + d * rhom13)
  xl = 1.d0 + (7.d0 / 3.d0) * (c * rhom13 + d * rhom13 / (1.d0 + d * &
       rhom13) )
  ff = a * b * grho / 24.d0
  rhom53 = rhom13**5
  sc = ff * rhom53 * om * xl
  dom = - om * (c + d+c * d * rhom13) / (1.d0 + d * rhom13)
  dxl = (7.d0 / 3.d0) * (c + d+2.d0 * c * d * rhom13 + c * d * d * &
       rhom13**2) / (1.d0 + d * rhom13) **2
  rhom43 = rhom13**4
  v1c = - ff * rhom43 / 3.d0 * (5.d0 * rhom43 * om * xl + rhom53 * &
       dom * xl + rhom53 * om * dxl)
  v2c = 2.d0 * sc / grho
  !
  return
end subroutine glyp
!
!-----------------------------------------------------------------------
subroutine ggac (rho, grho, sc, v1c, v2c)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91) correlation part
  !
  use parameters
  implicit none
  real(kind=DP) :: rho, grho, sc, v1c, v2c
  real(kind=DP) :: al, pa, pb, pc, pd, cx, cxc0, cc0
  parameter (al = 0.09d0, pa = 0.023266d0, pb = 7.389d-6, pc = &
       8.723d0, pd = 0.472d0)
  parameter (cx = - 0.001667d0, cxc0 = 0.002568d0, cc0 = - cx + &
       cxc0)
  real(kind=DP) :: third, pi34, nu, be, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (nu = 15.755920349483144d0, be = nu * cc0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
  ! xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(kind=DP) :: kf, ks, rs, rs2, rs3, ec, vc, t, expe, af, bf, y, xy, &
       qy, s1
  real(kind=DP) :: h0, dh0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, &
       dh1, ddh1
  !
  rs = pi34 / rho**third
  rs2 = rs * rs
  rs3 = rs * rs2
  call pw (rs, 1, ec, vc)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  t = sqrt (grho) / (2.d0 * ks * rho)
  expe = exp ( - 2.d0 * al * ec / (be * be) )
  af = 2.d0 * al / be * (1.d0 / (expe-1.d0) )
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + 2.d0 * al / be * t * t * xy
  h0 = be * be / (2.d0 * al) * log (s1)
  dh0 = be * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
       be-7.d0 / 3.d0) )
  ddh0 = be / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  ee = - 100.d0 * (ks / kf * t) **2
  cna = cxc0 + pa * rs + pb * rs2
  dcna = pa * rs + 2.d0 * pb * rs2
  cnb = 1.d0 + pc * rs + pd * rs2 + 1.d4 * pb * rs3
  dcnb = pc * rs + 2.d0 * pd * rs2 + 3.d4 * pb * rs3
  cn = cna / cnb - cx
  dcn = dcna / cnb - cna * dcnb / (cnb * cnb)
  h1 = nu * (cn - cc0 - 3.d0 / 7.d0 * cx) * t * t * exp (ee)
  dh1 = - third * (h1 * (7.d0 + 8.d0 * ee) + nu * t * t * exp (ee) &
       * dcn)
  ddh1 = 2.d0 * h1 * (1.d0 + ee) * rho / grho
  sc = rho * (h0 + h1)
  v1c = h0 + h1 + dh0 + dh1
  v2c = ddh0 + ddh1
  !
  return
end subroutine ggac
!
!---------------------------------------------------------------
subroutine pbex (rho, grho, iflag, sx, v1x, v2x)
  !---------------------------------------------------------------
  !
  ! PBE exchange (without Slater exchange):
  ! iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  ! iflag=2  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)
  !
  use parameters
  implicit none
  real(kind=DP) :: rho, grho, sx, v1x, v2x
  ! input: charge and squared gradient
  ! output: energy
  ! output: potential
  integer :: iflag
  ! local variables
  real(kind=DP) :: kf, agrho, s1, s2, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  real(kind=DP) :: dxunif, dfx, f1, f2, f3, dfx1
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  real(kind=DP) :: pi, third, c1, c2, c5
  parameter (pi = 3.14159265358979d0, third = 1.d0 / 3.d0, c1 = &
       0.75d0 / pi, c2 = 3.093667726280136d0, c5 = 4.d0 * third)
  ! parameters of the functional
  real(kind=DP) :: k (2), mu
  data k / 0.804d0, 1.2450D0 /, mu / 0.21951d0 /
  !
  agrho = sqrt (grho)
  kf = c2 * rho**third
  dsg = 0.5d0 / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  !
  !   Energy
  !
  f1 = s2 * mu / k (iflag)
  f2 = 1.d0 + f1
  f3 = k (iflag) / f2
  fx = k (iflag) - f3
  exunif = - c1 * kf
  sx = exunif * fx
  !
  !   Potential
  !
  dxunif = exunif * third
  dfx1 = f2 * f2
  dfx = 2.d0 * mu * s1 / dfx1
  v1x = sx + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho

  sx = sx * rho
  return
end subroutine pbex
!
!---------------------------------------------------------------
subroutine pbec (rho, grho, sc, v1c, v2c)
  !---------------------------------------------------------------
  !
  ! PBE correlation (without LDA part)
  ! J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  !
  use parameters
  implicit none
  real(kind=DP) :: rho, grho, sc, v1c, v2c
  real(kind=DP) :: ga, be
  parameter (ga = 0.031091d0, be = 0.066725d0)
  real(kind=DP) :: third, pi34, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(kind=DP) :: kf, ks, rs, ec, vc, t, expe, af, bf, y, xy, qy
  real(kind=DP) :: s1, h0, dh0, ddh0
  !
  rs = pi34 / rho**third
  call pw (rs, 1, ec, vc)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  t = sqrt (grho) / (2.d0 * ks * rho)
  expe = exp ( - ec / ga)
  af = be / ga * (1.d0 / (expe-1.d0) )
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + be / ga * t * t * xy
  h0 = ga * log (s1)
  dh0 = be * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
       be-7.d0 / 3.d0) )
  ddh0 = be / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  sc = rho * h0
  v1c = h0 + dh0
  v2c = ddh0
  !
  return
end subroutine pbec
