!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine xc_spin (rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     lsd exchange and correlation functionals - Hartree a.u.
  !
  !     exchange  :  Slater (alpha=2/3)
  !     correlation: Ceperley & Alder (Perdew-Zunger parameters)
  !                  Perdew & Wang
  !
  !     input : rho = rhoup(r)+rhodw(r)
  !             zeta=(rhoup(r)-rhodw(r))/rho
  !
  USE kinds
  use funct
  implicit none

  real(kind=DP) :: rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw
  !
  real(kind=DP), parameter :: small= 1.d-10, third = 1.d0/3.d0, &
       pi34= 0.6203504908994d0 ! pi34=(3/4pi)^(1/3)
  real(kind=DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0d0
     vcup = 0.0d0
     vcdw = 0.0d0
     ex = 0.0d0
     vxup = 0.0d0
     vxdw = 0.0d0
     return
  else
     rs = pi34 / rho**third
  endif
  !..exchange
  if (iexch == 1) then
     call slater_spin (rho, zeta, ex, vxup, vxdw)
  elseif (iexch == 2) then
     call slater1_spin (rho, zeta, ex, vxup, vxdw)
  ELSEIF (iexch == 3) THEN
     call slater_rxc_spin ( rho, zeta, ex, vxup, vxdw )
  else
     ex = 0.0d0
     vxup = 0.0d0
     vxdw = 0.0d0
  endif
  !..correlation
  if (icorr == 0) then
     ec = 0.0d0
     vcup = 0.0d0
     vcdw = 0.0d0
  elseif (icorr == 1) then
     call pz_spin (rs, zeta, ec, vcup, vcdw)
  elseif (icorr == 4) then
     call pw_spin (rs, zeta, ec, vcup, vcdw)
  else
     call errore ('lsda_functional', 'not implemented', icorr)
  endif
  !
  return
end subroutine xc_spin
!
!-----------------------------------------------------------------------
subroutine gcx_spin (rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, &
     v1xdw, v2xup, v2xdw)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange - Hartree a.u.
  !     Implemented:  Becke88, GGA (PW91), PBE, revPBE
  !
  use funct
  USE kinds
  implicit none
  !
  !     dummy arguments
  !
  real(kind=DP) :: rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, v1xdw, &
       v2xup, v2xdw
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !
  real(kind=DP), parameter :: small = 1.d-10
  real(kind=DP) :: rho, sxup, sxdw
  integer :: iflag
  !
  !
  ! exchange
  rho = rhoup + rhodw
  if (rho <= small .or. igcx == 0) then
     sx = 0.0d0
     v1xup = 0.0d0
     v2xup = 0.0d0
     v1xdw = 0.0d0
     v2xdw = 0.0d0
  elseif (igcx == 1) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.d0
        v1xup = 0.d0
        v2xup = 0.d0
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.d0
        v1xdw = 0.d0
        v2xdw = 0.d0
     endif
     sx = sxup + sxdw
  elseif (igcx == 2) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call ggax (2.d0 * rhoup, 4.d0 * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.d0
        v1xup = 0.d0
        v2xup = 0.d0
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call ggax (2.d0 * rhodw, 4.d0 * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.d0
        v1xdw = 0.d0
        v2xdw = 0.d0
     endif
     sx = 0.5d0 * (sxup + sxdw)
     v2xup = 2.d0 * v2xup
     v2xdw = 2.d0 * v2xdw
  elseif (igcx == 3 .or. igcx == 4) then
     ! igcx=3: PBE  igcx=4: revised PBE
     if (igcx == 3) then
        iflag = 1
     else
        iflag = 2
     endif
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call pbex (2.d0 * rhoup, 4.d0 * grhoup2, iflag, sxup, v1xup, v2xup)
     else
        sxup = 0.d0
        v1xup = 0.d0
        v2xup = 0.d0
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call pbex (2.d0 * rhodw, 4.d0 * grhodw2, iflag, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.d0
        v1xdw = 0.d0
        v2xdw = 0.d0
     endif
     sx = 0.5d0 * (sxup + sxdw)
     v2xup = 2.d0 * v2xup
     v2xdw = 2.d0 * v2xdw
  else
     call errore ('gcx_spin', 'not implemented', igcx)
  endif
  !
  return
end subroutine gcx_spin
!
!-----------------------------------------------------------------------
subroutine gcc_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for correlations - Hartree a.u.
  !     Implemented:  Perdew86, GGA (PW91), PBE
  !
  use funct
  USE kinds
  implicit none
  !
  !     dummy arguments
  !
  real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
  ! the total charge
  ! the magnetization
  ! the gradient of the charge squared
  ! exchange and correlation energies
  ! derivatives of correlation wr. rho
  ! derivatives of correlation wr. grho

  real(kind=DP), parameter :: small = 1.d-10, epsr=1.d-6
  !
  !
  if ( abs(zeta) > 1.d0 ) then
     sc = 0.0d0
     v1cup = 0.0d0
     v1cdw = 0.0d0
     v2c = 0.0d0
     return
  else
     !
     ! ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
     zeta = SIGN( MIN( ABS( zeta ), ( 1.D0 - epsr ) ) , zeta )
  endif

  if (rho <= small .or. sqrt(abs(grho)) <= small) then
     sc = 0.0d0
     v1cup = 0.0d0
     v1cdw = 0.0d0
     v2c = 0.0d0
  elseif (igcc == 1) then
     call perdew86_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 2) then
     call ggac_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 3 .or. igcc > 4) then
     call errore ('lsda_functionals', 'not implemented', igcc)
  elseif (igcc == 4) then
        call pbec_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  else
     sc = 0.0d0
     v1cup = 0.0d0
     v1cdw = 0.0d0
     v2c = 0.0d0
  endif
  !
  return
end subroutine gcc_spin
!
!-----------------------------------------------------------------------
subroutine pz_polarized (rs, ec, vc)
  !-----------------------------------------------------------------------
  !     J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     spin-polarized energy and potential
  !
  USE kinds
  implicit none
  real(kind=DP) :: rs, ec, vc
  real(kind=DP) :: a, b, c, d, gc, b1, b2
  parameter (a = 0.01555d0, b = - 0.0269d0, c = 0.0007d0, d = &
       - 0.0048d0, gc = - 0.0843d0, b1 = 1.3981d0, b2 = 0.2611d0)
  real(kind=DP) :: lnrs, rs12, ox, dox
  LOGICAL       :: xc_rel
  REAL(KIND=DP), PARAMETER :: xcprefact = 0.022575584, pi34 = 0.6203504908994d0 
  REAL(KIND=DP) :: betha, etha, csi, prefact
  !
  if (rs.lt.1.0d0) then
     ! high density formula
     lnrs = log (rs)
     ec = a * lnrs + b + c * rs * lnrs + d * rs
     vc = a * lnrs + (b - a / 3.d0) + 2.d0 / 3.d0 * c * rs * lnrs + &
          (2.d0 * d-c) / 3.d0 * rs
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     ox = 1.d0 + b1 * rs12 + b2 * rs
     dox = 1.d0 + 7.d0 / 6.d0 * b1 * rs12 + 4.d0 / 3.d0 * b2 * rs
     ec = gc / ox
     vc = ec * dox / ox
  endif
  !
!  IF ( lxc_rel ) THEN
!     betha = prefact * pi34 / rs
!     etha = DSQRT( 1 + betha**2 )
!     csi = betha + etha
!     prefact = 1.0D0 - (3.0D0/2.0D0) * ( (betha*etha - log(csi))/betha**2 )**2
!     ec = ec * prefact
!     vc = vc * prefact
!  ENDIF
  return
end subroutine pz_polarized
!
!-----------------------------------------------------------------------
subroutine pz_spin (rs, zeta, ec, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds
  implicit none
  real(kind=DP) :: rs, zeta, ec, vcup, vcdw
  !
  real(kind=DP) :: ecu, vcu, ecp, vcp, fz, dfz
  real(kind=DP) :: p43, third
  parameter (p43 = 4.0d0 / 3.d0, third = 1.d0 / 3.d0)
  !
  ! unpolarized part (Perdew-Zunger formula)
  call pz (rs, 1, ecu, vcu)
  ! polarization contribution
  call pz_polarized (rs, ecp, vcp)
  !
  fz = ( (1.0d0 + zeta) **p43 + (1.d0 - zeta) **p43 - 2.d0) / &
       (2.d0**p43 - 2.d0)
  dfz = p43 * ( (1.0d0 + zeta) **third- (1.d0 - zeta) **third) &
       / (2.d0**p43 - 2.d0)
  !
  ec = ecu + fz * (ecp - ecu)
  vcup = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (1.d0 - zeta)
  vcdw = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * ( - 1.d0 - &
       zeta)
  !
  return
end subroutine pz_spin
!
!-----------------------------------------------------------------------
subroutine pw_spin (rs, zeta, ec, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds
  implicit none
  real(kind=DP) :: rs, zeta, ec, vcup, vcdw
  ! xc parameters, unpolarised
  real(kind=DP) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, a1 = 0.21370d0, b1 = 7.5957d0, b2 = &
       3.5876d0, b3 = 1.6382d0, b4 = 0.49294d0, c0 = a, c1 = 0.046644d0, &
       c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, d1 = 1.4408d0)
  ! xc parameters, polarised
  real(kind=DP) :: ap, a1p, b1p, b2p, b3p, b4p, c0p, c1p, c2p, c3p, d0p, &
       d1p
  parameter (ap = 0.015545d0, a1p = 0.20548d0, b1p = 14.1189d0, b2p &
       = 6.1977d0, b3p = 3.3662d0, b4p = 0.62517d0, c0p = ap, c1p = &
       0.025599d0, c2p = 0.00319d0, c3p = 0.00384d0, d0p = 0.3287d0, d1p &
       = 1.7697d0)
  ! xc parameters, antiferro
  real(kind=DP) :: aa, a1a, b1a, b2a, b3a, b4a, c0a, c1a, c2a, c3a, d0a, &
       d1a
  parameter (aa = 0.016887d0, a1a = 0.11125d0, b1a = 10.357d0, b2a = &
       3.6231d0, b3a = 0.88026d0, b4a = 0.49671d0, c0a = aa, c1a = &
       0.035475d0, c2a = 0.00188d0, c3a = 0.00521d0, d0a = 0.2240d0, d1a &
       = 0.3969d0)
  real(kind=DP) :: fz0
  parameter (fz0 = 1.709921d0)
  real(kind=DP) :: rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz
  real(kind=DP) :: om, dom, olog, epwc, vpwc
  real(kind=DP) :: omp, domp, ologp, epwcp, vpwcp
  real(kind=DP) :: oma, doma, ologa, alpha, vpwca
  !
  !     if(rs.lt.0.5d0) then
  ! high density formula (not implemented)
  !
  !     else if(rs.gt.100.d0) then
  ! low density formula  (not implemented)
  !
  !     else
  ! interpolation formula
  zeta2 = zeta * zeta
  zeta3 = zeta2 * zeta
  zeta4 = zeta3 * zeta
  rs12 = sqrt (rs)
  rs32 = rs * rs12
  rs2 = rs**2
  ! unpolarised
  om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
  dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 * rs32 &
       + 2.d0 * b4 * rs2)
  olog = log (1.d0 + 1.0d0 / om)
  epwc = - 2.d0 * a * (1.d0 + a1 * rs) * olog
  vpwc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 * rs) * olog - 2.d0 / &
       3.d0 * a * (1.d0 + a1 * rs) * dom / (om * (om + 1.d0) )
  ! polarized
  omp = 2.d0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
  domp = 2.d0 * ap * (0.5d0 * b1p * rs12 + b2p * rs + 1.5d0 * b3p * &
       rs32 + 2.d0 * b4p * rs2)
  ologp = log (1.d0 + 1.0d0 / omp)
  epwcp = - 2.d0 * ap * (1.d0 + a1p * rs) * ologp
  vpwcp = - 2.d0 * ap * (1.d0 + 2.d0 / 3.d0 * a1p * rs) * ologp - &
       2.d0 / 3.d0 * ap * (1.d0 + a1p * rs) * domp / (omp * (omp + 1.d0) &
       )
  ! antiferro
  oma = 2.d0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
  doma = 2.d0 * aa * (0.5d0 * b1a * rs12 + b2a * rs + 1.5d0 * b3a * &
       rs32 + 2.d0 * b4a * rs2)
  ologa = log (1.d0 + 1.0d0 / oma)
  alpha = 2.d0 * aa * (1.d0 + a1a * rs) * ologa
  vpwca = + 2.d0 * aa * (1.d0 + 2.d0 / 3.d0 * a1a * rs) * ologa + &
       2.d0 / 3.d0 * aa * (1.d0 + a1a * rs) * doma / (oma * (oma + 1.d0) &
       )
  !
  fz = ( (1.d0 + zeta) ** (4.d0 / 3.d0) + (1.d0 - zeta) ** (4.d0 / &
       3.d0) - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2.d0)
  dfz = ( (1.d0 + zeta) ** (1.d0 / 3.d0) - (1.d0 - zeta) ** (1.d0 / &
       3.d0) ) * 4.d0 / (3.d0 * (2.d0** (4.d0 / 3.d0) - 2.d0) )
  !
  ec = epwc + alpha * fz * (1.d0 - zeta4) / fz0 + (epwcp - epwc) &
       * fz * zeta4
  !
  vcup = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
       * fz * zeta4 + (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
       zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
       * (1.d0 - zeta)

  vcdw = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
       * fz * zeta4 - (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
       zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
       * (1.d0 + zeta)
  !      endif
  !
  return
end subroutine pw_spin
!
!-----------------------------------------------------------------------
subroutine becke88_spin (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! Becke exchange: A.D. Becke, PRA 38, 3098 (1988) - Spin polarized case
  !
  USE kinds
  implicit none
  real(kind=DP) :: rho, grho, sx, v1x, v2x
  ! input: charge
  ! input: gradient
  ! output: the up and down energies
  ! output: first part of the potential
  ! output: the second part of the potential
  !
  real(kind=DP) :: beta, third
  parameter (beta = 0.0042d0, third = 1.d0 / 3.d0)
  real(kind=DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee
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
  USE kinds
  implicit none
  real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
  real(kind=DP) :: p1, p2, p3, p4, pc1, pc2, pci
  parameter (p1 = 0.023266d0, p2 = 7.389d-6, p3 = 8.723d0, p4 = &
       0.472d0)
  parameter (pc1 = 0.001667d0, pc2 = 0.002568d0, pci = pc1 + pc2)
  real(kind=DP) :: third, pi34
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  ! pi34=(3/4pi)^(1/3)
  !
  real(kind=DP) :: rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs
  real(kind=DP) :: dcna, dcnb, dcn, phi, ephi, dd, ddd
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
subroutine ggac_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91) correlation part - spin-polarized
  !
  USE kinds
  implicit none
  real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
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
  real(kind=DP) :: kf, ks, rs, rs2, rs3, ec, vcup, vcdw, t, expe, af, y, &
       xy, qy, s1, h0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, dh1, &
       ddh1, fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, &
       dh0zdw, dh1zup, dh1zdw
  !
  rs = pi34 / rho**third
  rs2 = rs * rs
  rs3 = rs * rs2
  call pw_spin (rs, zeta, ec, vcup, vcdw)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  fz = 0.5d0 * ( (1.d0 + zeta) ** (2.d0 / 3.d0) + (1.d0 - zeta) ** ( &
       2.d0 / 3.d0) )
  fz2 = fz * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1.d0 + zeta) ** ( - 1.d0 / 3.d0) - (1.d0 - zeta) ** ( - &
       1.d0 / 3.d0) ) / 3.d0
  t = sqrt (grho) / (2.d0 * fz * ks * rho)
  expe = exp ( - 2.d0 * al * ec / (fz3 * be * be) )
  af = 2.d0 * al / be * (1.d0 / (expe-1.d0) )
  bfup = expe * (vcup - ec) / fz3
  bfdw = expe * (vcdw - ec) / fz3
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + 2.d0 * al / be * t * t * xy
  h0 = fz3 * be * be / (2.d0 * al) * log (s1)
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
  h1 = nu * (cn - cc0 - 3.d0 / 7.d0 * cx) * fz3 * t * t * exp (ee)
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
subroutine pbec_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !---------------------------------------------------------------
  !
  ! PBE correlation (without LDA part) - spin-polarized
  ! J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  !
  USE kinds
  implicit none
  real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
  real(kind=DP) :: ga, be
  parameter (ga = 0.031091d0, be = 0.066725d0)
  real(kind=DP) :: third, pi34, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(kind=DP) :: kf, ks, rs, ec, vcup, vcdw, t, expe, af, y, xy, qy, &
       s1, h0, ddh0
  real(kind=DP) :: fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, &
       dh0zup, dh0zdw
  !
  rs = pi34 / rho**third
  call pw_spin (rs, zeta, ec, vcup, vcdw)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  fz = 0.5d0 * ( (1.d0 + zeta) ** (2.d0 / 3.d0) + (1.d0 - zeta) ** ( &
       2.d0 / 3.d0) )
  fz2 = fz * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1.d0 + zeta) ** ( - 1.d0 / 3.d0) - (1.d0 - zeta) ** ( - &
       1.d0 / 3.d0) ) / 3.d0
  t = sqrt (grho) / (2.d0 * fz * ks * rho)
  expe = exp ( - ec / (fz3 * ga) )
  af = be / ga * (1.d0 / (expe-1.d0) )
  bfup = expe * (vcup - ec) / fz3
  bfdw = expe * (vcdw - ec) / fz3
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + be / ga * t * t * xy
  h0 = fz3 * ga * log (s1)
  dh0up = be * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfup / be-7.d0 / 3.d0) )
  dh0dw = be * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfdw / be-7.d0 / 3.d0) )
  dh0zup = (3.d0 * h0 / fz - be * t * t * fz2 / s1 * (2.d0 * xy - &
       qy * (3.d0 * af * expe * ec / fz3 / be+2.d0) ) ) * dfz * (1.d0 - &
       zeta)
  dh0zdw = - (3.d0 * h0 / fz - be * t * t * fz2 / s1 * (2.d0 * xy - &
       qy * (3.d0 * af * expe * ec / fz3 / be+2.d0) ) ) * dfz * (1.d0 + &
       zeta)

  ddh0 = be * fz / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  sc = rho * h0
  v1cup = h0 + dh0up + dh0zup
  v1cdw = h0 + dh0dw + dh0zdw
  v2c = ddh0
  return
end subroutine pbec_spin
!
!-----------------------------------------------------------------------
subroutine slater_spin (rho, zeta, ex, vxup, vxdw)
  !-----------------------------------------------------------------------
  !     Slater exchange with alpha=2/3, spin-polarized case
  !
  USE kinds
  implicit none
  real(kind=DP) :: rho, zeta, ex, vxup, vxdw
  real(kind=DP) :: f, alpha, third, p43
  parameter (f = - 1.10783814957303361d0, alpha = 2.0d0 / 3.0d0)
  ! f = -9/8*(3/pi)^(1/3)
  parameter (third = 1.d0 / 3.d0, p43 = 4.d0 / 3.d0)
  real(kind=DP) :: exup, exdw, rho13
  !
  rho13 = ( (1.d0 + zeta) * rho) **third
  exup = f * alpha * rho13
  vxup = p43 * f * alpha * rho13
  rho13 = ( (1.d0 - zeta) * rho) **third
  exdw = f * alpha * rho13
  vxdw = p43 * f * alpha * rho13
  ex = 0.5d0 * ( (1.d0 + zeta) * exup + (1.d0 - zeta) * exdw)
  !
  return
end subroutine slater_spin

!-----------------------------------------------------------------------
SUBROUTINE slater_rxc_spin ( rho, Z, ex, vxup, vxdw )
  !-----------------------------------------------------------------------
  !     Slater exchange with alpha=2/3, relativistic exchange case
  !
  USE kinds
  IMPLICIT none
  real (kind=DP):: rho, ex, vxup, vxdw
  !
  real(kind=DP), PARAMETER :: ZERO=0.D0, ONE=1.D0, PFIVE=.5D0, &
       OPF=1.5D0, C014=0.014D0, pi = 3.14159265358979d0
  real (kind=DP):: rs, trd, ftrd, tftm, a0, alp, z, fz, fzp, vxp, exp, &
       beta, sb, alb, vxf, exf

  TRD = ONE/3
  FTRD = 4*TRD
  TFTM = 2**FTRD-2
  A0 = (4/(9*PI))**TRD
  
  !      X-alpha parameter:
  ALP = 2 * TRD
  
  IF (rho <=  ZERO) THEN
     EX = ZERO
     vxup  = ZERO
     vxdw  = ZERO
     RETURN
  ELSE
     FZ = ((1+Z)**FTRD+(1-Z)**FTRD-2)/TFTM
     FZP = FTRD*((1+Z)**TRD-(1-Z)**TRD)/TFTM
  ENDIF
  RS = (3 / (4*PI*rho) )**TRD
  VXP = -3*ALP/(2*PI*A0*RS)
  EXP = 3*VXP/4
  
  BETA = C014/RS
  SB = SQRT(1+BETA*BETA)
  ALB = LOG(BETA+SB)
  VXP = VXP * (-PFIVE + OPF * ALB / (BETA*SB))
  EXP = EXP * (ONE-OPF*((BETA*SB-ALB)/BETA**2)**2)
  
  VXF = 2**TRD*VXP
  EXF = 2**TRD*EXP
  vxup  = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP)
  vxdw  = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP)
  EX    = EXP + FZ*(EXF-EXP)
       
END SUBROUTINE slater_rxc_spin


!-----------------------------------------------------------------------
subroutine slater1_spin (rho, zeta, ex, vxup, vxdw)
  !-----------------------------------------------------------------------
  !     Slater exchange with alpha=2/3, spin-polarized case
  !
  use kinds, only: dp
  implicit none
  real(kind=DP) :: rho, zeta, ex, vxup, vxdw
  real(kind=DP), parameter :: f = - 1.10783814957303361d0, alpha = 1.0d0, &
       third = 1.d0 / 3.d0, p43 = 4.d0 / 3.d0
  ! f = -9/8*(3/pi)^(1/3)
  real(kind=DP) :: exup, exdw, rho13
  !
  rho13 = ( (1.d0 + zeta) * rho) **third
  exup = f * alpha * rho13
  vxup = p43 * f * alpha * rho13
  rho13 = ( (1.d0 - zeta) * rho) **third
  exdw = f * alpha * rho13
  vxdw = p43 * f * alpha * rho13
  ex = 0.5d0 * ( (1.d0 + zeta) * exup + (1.d0 - zeta) * exdw)
  !
  return
end subroutine slater1_spin
