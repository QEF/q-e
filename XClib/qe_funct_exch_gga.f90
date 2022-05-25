!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
MODULE exch_gga
!------------------------------------------------------------------------
!! GGA exchange functionals
!
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE becke88( rho, grho, sx, v1x, v2x )
  !-----------------------------------------------------------------------
  !! Becke exchange: A.D. Becke, PRA 38, 3098 (1988)
  !! only gradient-corrected part, no Slater term included
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee
  REAL(DP), PARAMETER :: beta=0.0042_DP
  REAL(DP), PARAMETER :: third=1._DP/3._DP, two13=1.259921049894873_DP
                                          ! two13= 2^(1/3)
  !
  rho13 = rho**third
  rho43 = rho13**4
  !
  xs = two13 * SQRT(grho)/rho43
  xs2 = xs * xs
  !
  sa2b8 = SQRT(1.0_DP + xs2)
  shm1 = LOG(xs + sa2b8)
  !
  dd = 1.0_DP + 6.0_DP * beta * xs * shm1
  dd2 = dd * dd
  !
  ee = 6.0_DP * beta * xs2 / sa2b8 - 1._DP
  sx = two13 * grho / rho43 * ( - beta / dd)
  !
  v1x = - (4._DP/3._DP) / two13 * xs2 * beta * rho13 * ee / dd2
  v2x = two13 * beta * (ee-dd) / (rho43 * dd2)
  !
  RETURN
  !
END SUBROUTINE becke88
!
!
!-----------------------------------------------------------------------
SUBROUTINE ggax( rho, grho, sx, v1x, v2x )
  !-----------------------------------------------------------------------
  !! Perdew-Wang GGA (PW91), exchange part:
  !! J.P. Perdew et al.,PRB 46, 6671 (1992)
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das, &
              dbs, dls
  REAL(DP), PARAMETER :: f1=0.19645_DP, f2=7.7956_DP, f3=0.2743_DP, &
                         f4=0.1508_DP,  f5=0.004_DP
  REAL(DP), PARAMETER :: fp1=-0.019292021296426_DP, fp2=0.161620459673995_DP
                       ! fp1= -3/(16 pi)*(3 pi^2)^(-1/3)
                       ! fp2= (1/2)(3 pi^2)**(-1/3)
  !
  rhom43 = rho**(-4.d0/3.d0)
  s  = fp2 * SQRT(grho) * rhom43
  s2 = s * s
  s3 = s2 * s
  s4 = s2 * s2
  !
  exps  = f4 * EXP( - 100.d0 * s2)
  as    = f3 - exps - f5 * s2
  sa2b8 = SQRT(1.0d0 + f2 * f2 * s2)
  shm1  = LOG(f2 * s + sa2b8)
  bs    = 1.d0 + f1 * s * shm1 + f5 * s4
  !
  das = (200.d0 * exps - 2.d0 * f5) * s
  dbs = f1 * (shm1 + f2 * s / sa2b8) + 4.d0 * f5 * s3
  dls = (das / as - dbs / bs)
  !
  sx  = fp1 * grho * rhom43 * as / bs
  v1x = - 4.d0 / 3.d0 * sx / rho * (1.d0 + s * dls)
  v2x = fp1 * rhom43 * as / bs * (2.d0 + s * dls)
  !
  RETURN
  !
END SUBROUTINE ggax
!
!
!---------------------------------------------------------------
SUBROUTINE pbex( rho, grho, iflag, sx, v1x, v2x )
  !---------------------------------------------------------------
  !! PBE exchange (without Slater exchange):
  !! iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  !! iflag=2  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)
  !! iflag=3  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
  !! iflag=4  PBEQ2D: L. Chiodo et al., PRL 108, 126402 (2012)
  !! iflag=5  optB88: Klimes et al., J. Phys. Cond. Matter, 22, 022201 (2010)
  !! iflag=6  optB86b: Klimes et al., Phys. Rev. B 83, 195131 (2011)
  !! iflag=7  ev: Engel and Vosko, PRB 47, 13164 (1991)
  !! iflag=8  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
  !! iflag=9  W31X: D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)
  !
  USE kind_l,      ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  INTEGER, INTENT(IN) :: iflag
  REAL(DP), INTENT(IN) :: rho, grho
  ! input: charge and squared gradient
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  ! output: energy, potential
  !
  ! ... local variables
  !
  REAL(DP) :: kf, agrho, s1, s2, ds, dsg, exunif, fx, sx_s
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  ! auxiliary variable for energy calculation
  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1
  REAL(DP) :: p, amu, ab, c, dfxdp, dfxds, s, ak
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0
  REAL(DP), PARAMETER :: third=1._DP/3._DP, c1=0.75_DP/pi,        &
                         c2=3.093667726280136_DP, c5=4._DP*third, &
                         c6=c2*2.51984210_DP, c7=0.8_DP
                         ! (3pi^2)^(1/3)*2^(4/3)
  ! parameters of the functional
  REAL(DP) :: k(9), mu(9), ev(6)
  !         pbe         revpbe       pbesol     pbeq2d     optB88   optB86b
  !         ev          rpbe         W31x
  DATA k  / 0.804_DP,   1.2450_DP,   0.804_DP , 0.804_DP,  1.2_DP,  0.0_DP,       &
            0.000_DP,   0.8040_DP,   1.10_DP /,                                   &
       mu / 0.2195149727645171_DP, 0.2195149727645171_DP, 0.12345679012345679_DP, &
            0.12345679012345679_DP, 0.22_DP, 0.1234_DP, 0.000_DP,                 &
            0.2195149727645171_DP, 0.12345679012345679_DP /,                      &
       ev / 1.647127_DP, 0.980118_DP, 0.017399_DP, 1.523671_DP, 0.367229_DP,      &
            0.011282_DP /  ! a and b parameters of Engel and Vosko
  !
  SELECT CASE( iflag )
  CASE( 4 )
     !
     agrho = SQRT(grho)
     kf = c2 * rho**third
     dsg = 0.5_DP / kf
     s1 = agrho * dsg / rho
     p = s1*s1
     s = s1
     ak = 0.804_DP
     amu = 10._DP/81._DP
     ab = 0.5217_DP
     c = 2._DP
     fx =  ak - ak / (1.0_DP + amu * p / ak)  + p**2 * (1.0_DP + p)/       &
            (10**c + p**3) * ( -1.0_DP - ak + ak / (1.0_DP + amu * p / ak) &
           + ab * p ** (-0.1D1/ 0.4D1) )
     !
     exunif = - c1 * kf
     sx_s = exunif * fx
     !
     dxunif = exunif * third
     !
     dfxdp = DBLE(1 / (1 + amu * p / ak) ** 2 * amu) + DBLE(2 * p * (1   &
     + p) / (10 ** c + p ** 3) * (-1 - ak + ak / (1 + amu * p / ak) + ab &
     * p ** (-0.1d1 / 0.4D1))) + DBLE(p ** 2 / (10 ** c + p ** 3) * (    &
     -1 - ak + ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) -  &
     DBLE(3 * p ** 4 * (1 + p) / (10 ** c + p ** 3) ** 2 * (-1 - ak +    &
     ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) + DBLE(p **  &
     2) * DBLE(1 + p) / DBLE(10 ** c + p ** 3) * (-DBLE(1 / (1 + amu *   &
     p / ak) ** 2 * amu) - DBLE(ab * p ** (-0.5d1 / 0.4D1)) / 0.4D1)
     !
     dfxds = dfxdp*2._DP*s
     dfx = dfxds
     ds = - c5 * s1
     !
     v1x = sx_s + dxunif * fx + exunif * dfx * ds
     v2x = exunif * dfx * dsg / agrho
     sx  = sx_s * rho
     !
  CASE( 5, 9 )
     !
     agrho = SQRT(grho)
     kf = c2 * rho**third
     dsg = 0.5_DP / kf
     s1 = agrho * dsg / rho
     ab = mu(iflag) / k(iflag)
     p = s1*c6
     c = LOG(p + SQRT(p*p+1)) ! asinh(p)
     dfx1 = 1 + ab*s1*c
     fx = mu(iflag)*s1*s1/dfx1
     !
     exunif = - c1 * kf
     sx_s = exunif * fx
     !
     dxunif = exunif * third
     !
     dfx = 2*fx/s1-fx/dfx1*(ab*c+ab*s1/SQRT(p*p+1)*c6)
     ds  = - c5 * s1
     !
     v1x = sx_s + dxunif * fx + exunif * dfx * ds
     v2x = exunif * dfx * dsg / agrho
     sx  = sx_s * rho
     !
  CASE( 6 )
     !
     agrho = SQRT(grho)
     kf = c2 * rho**third
     dsg = 0.5_DP / kf
     s1 = agrho * dsg / rho
     p = mu(iflag)*s1*s1
     fx =  p / ( 1._DP + p )**c7
     !
     exunif = - c1 * kf
     sx_s = exunif * fx
     !
     dxunif = exunif * third
     !
     dfx = 2*mu(iflag)*s1*fx*(1+(1-c7)*p)/(p*(1+p))
     ds = - c5 * s1
     !
     v1x = sx_s + dxunif * fx + exunif * dfx * ds
     v2x = exunif * dfx * dsg / agrho
     sx  = sx_s * rho
     !
  CASE( 7 )
     !
     agrho = SQRT(grho)
     kf = c2 * rho**third
     dsg = 0.5_DP / kf
     s1 = agrho * dsg / rho
     s2 = s1 * s1
     s = s2*s2
     f1 =  1._DP + ev(1)*s2 + ev(2)*s + ev(3)*s*s2
     f2 =  1._DP + ev(4)*s2 + ev(5)*s + ev(6)*s*s2
     fx = f1 / f2 - 1._DP
     !
     exunif = - c1 * kf
     sx_s = exunif * fx
     !
     dxunif = exunif * third
     ds = - c5 * s1
     !
     dfx  =  ev(1) + 2*ev(2)*s2 + 3*ev(3)*s
     dfx1 =  ev(4) + 2*ev(5)*s2 + 3*ev(6)*s
     dfx  = 2 * s1 * ( dfx - f1*dfx1/f2 ) / f2
     !
     v1x = sx_s + dxunif * fx + exunif * dfx * ds
     v2x = exunif * dfx * dsg / agrho
     sx  = sx_s * rho
     !
  CASE(8)
     !
     agrho = SQRT(grho)
     kf = c2 * rho**third
     dsg = 0.5_DP / kf
     s1 = agrho * dsg / rho
     s2 = s1 * s1
     f1 = exp( - mu(iflag) * s2 / k(iflag) )
     f2 = 1._DP - f1
     fx = k(iflag) * f2
     !
     exunif = - c1 * kf
     sx_s = exunif * fx
     !
     dxunif = exunif * third
     ds = - c5 * s1
     !
     dfx = 2._DP * mu(iflag) * s1 * exp( - mu(iflag) * s2 / k(iflag) )
     !
     v1x = sx_s + dxunif * fx + exunif * dfx * ds
     v2x = exunif * dfx * dsg / agrho
     sx  = sx_s * rho
     !
  CASE DEFAULT
     !
     agrho = SQRT(grho)
     kf = c2 * rho**third
     dsg = 0.5_DP / kf
     s1 = agrho * dsg / rho
     s2 = s1 * s1
     f1 = s2 * mu(iflag) / k(iflag)
     f2 = 1._DP + f1
     f3 = k(iflag) / f2
     fx = k(iflag) - f3
     !
     exunif = - c1 * kf
     sx_s = exunif * fx
     !
     dxunif = exunif * third
     ds = - c5 * s1
     !
     dfx1 = f2 * f2
     dfx = 2._DP * mu(iflag) * s1 / dfx1
     !
     v1x = sx_s + dxunif * fx + exunif * dfx * ds
     v2x = exunif * dfx * dsg / agrho
     sx  = sx_s * rho
     !
  END SELECT
  !
  !
  RETURN
  !
END SUBROUTINE pbex
!
!
!----------------------------------------------------------------------------
SUBROUTINE hcth( rho, grho, sx, v1x, v2x )
  !--------------------------------------------------------------------------
  !! HCTH/120, JCP 109, p. 6264 (1998)
  !! Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
  !! Present release: Mauro Boero, Tsukuba, 11/05/2004
  !
  !! * rhoa = rhob = 0.5 * rho
  !! * grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
  !! * sx  : total exchange correlation energy at point r
  !! * v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)
  !! * v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0
  REAL(DP), PARAMETER :: o3 = 1.0d0/3.0d0, o34 = 4.0d0/3.0d0, fr83 = 8.d0/3.d0
  REAL(DP) :: cg0(6), cg1(6), caa(6), cab(6), cx(6)
  REAL(DP) :: r3q2, r3pi, gr, rho_o3, rho_o34, xa, xa2, ra, rab,        &
       dra_drho, drab_drho, g, dg, era1, dera1_dra, erab0, derab0_drab, &
       ex, dex_drho, uaa, uab, ux, ffaa, ffab,  dffaa_drho, dffab_drho, &
       denaa, denab, denx, f83rho, bygr, gaa, gab, gx, taa, tab, txx,   &
       dgaa_drho, dgab_drho, dgx_drho, dgaa_dgr, dgab_dgr, dgx_dgr
  !
  r3q2 = 2.d0**(-o3)
  r3pi = (3.d0/pi)**o3
  ! ... coefficients for pw correlation
  cg0(1) = 0.031091d0
  cg0(2) = 0.213700d0
  cg0(3) = 7.595700d0
  cg0(4) = 3.587600d0
  cg0(5) = 1.638200d0
  cg0(6) = 0.492940d0
  cg1(1) = 0.015545d0
  cg1(2) = 0.205480d0
  cg1(3) =14.118900d0
  cg1(4) = 6.197700d0
  cg1(5) = 3.366200d0
  cg1(6) = 0.625170d0
  ! ... hcth-19-4 ...
  caa(1) =  0.489508d+00
  caa(2) = -0.260699d+00
  caa(3) =  0.432917d+00
  caa(4) = -0.199247d+01
  caa(5) =  0.248531d+01
  caa(6) =  0.200000d+00
  cab(1) =  0.514730d+00
  cab(2) =  0.692982d+01
  cab(3) = -0.247073d+02
  cab(4) =  0.231098d+02
  cab(5) = -0.113234d+02
  cab(6) =  0.006000d+00
  cx(1)  =  0.109163d+01
  cx(2)  = -0.747215d+00
  cx(3)  =  0.507833d+01
  cx(4)  = -0.410746d+01
  cx(5)  =  0.117173d+01
  cx(6)  =  0.004000d+00
  !  ... ... ... ... ...
  !
  gr = DSQRT(grho)
  rho_o3  = rho**(o3)
  rho_o34 = rho**(o34)
  xa = 1.25992105d0*gr/rho_o34
  xa2 = xa*xa
  ra = 0.781592642d0/rho_o3
  rab = r3q2*ra
  dra_drho = -0.260530881d0/rho_o34
  drab_drho = r3q2*dra_drho
  CALL pwcorr( ra, cg1, g, dg )
  era1 = g
  dera1_dra = dg
  CALL pwcorr( rab, cg0, g, dg )
  erab0 = g
  derab0_drab = dg
  ex = -0.75d0*r3pi*rho_o34
  dex_drho = -r3pi*rho_o3
  uaa = caa(6)*xa2
  uaa = uaa/(1.0d0+uaa)
  uab = cab(6)*xa2
  uab = uab/(1.0d0+uab)
  ux = cx(6)*xa2
  ux = ux/(1.0d0+ux)
  ffaa = rho*era1
  ffab = rho*erab0-ffaa
  dffaa_drho = era1 + rho*dera1_dra*dra_drho
  dffab_drho = erab0 + rho*derab0_drab*drab_drho - dffaa_drho
  ! mb-> i-loop removed
  denaa = 1.d0 / (1.0d0+caa(6)*xa2)
  denab = 1.d0 / (1.0d0+cab(6)*xa2)
  denx  = 1.d0 / (1.0d0+cx(6)*xa2)
  f83rho = fr83 / rho
  bygr = 2.0d0/gr
  gaa = caa(1)+uaa*(caa(2)+uaa*(caa(3)+uaa*(caa(4)+uaa*caa(5))))
  gab = cab(1)+uab*(cab(2)+uab*(cab(3)+uab*(cab(4)+uab*cab(5))))
  gx  = cx(1)+ux*(cx(2)+ux*(cx(3)+ux*(cx(4)+ux*cx(5))))
  taa = denaa*uaa*(caa(2)+uaa*(2.d0*caa(3)+uaa &
        *(3.d0*caa(4)+uaa*4.d0*caa(5))))
  tab = denab*uab*(cab(2)+uab*(2.d0*cab(3)+uab &
        *(3.d0*cab(4)+uab*4.d0*cab(5))))
  txx = denx*ux*(cx(2)+ux*(2.d0*cx(3)+ux &
        *(3.d0*cx(4)+ux*4.d0*cx(5))))
  dgaa_drho = -f83rho*taa
  dgab_drho = -f83rho*tab
  dgx_drho  = -f83rho*txx
  dgaa_dgr  =  bygr*taa
  dgab_dgr  =  bygr*tab
  dgx_dgr   =  bygr*txx
  ! mb
  sx  = ex*gx + ffaa*gaa + ffab*gab
  v1x = dex_drho*gx + ex*dgx_drho          &
             + dffaa_drho*gaa + ffaa*dgaa_drho &
             + dffab_drho*gab + ffab*dgab_drho
  v2x = (ex*dgx_dgr + ffaa*dgaa_dgr + ffab*dgab_dgr) / gr
  !
  RETURN
  !
END SUBROUTINE hcth
    !
    !-------------------------------------------------------
    SUBROUTINE pwcorr( r, c, g, dg )
      !-----------------------------------------------------
      !
      USE kind_l,   ONLY: DP
      !
      IMPLICIT NONE
      !
      !$acc routine seq
      !
      REAL(DP), INTENT(IN)  :: r, c(6)
      REAL(DP), INTENT(OUT) :: g, dg
      !
      ! ... local variables
      !
      REAL(DP) :: r12, r32, r2, rb, drb, sb
      !
      r12 = DSQRT(r)
      r32 = r*r12
      r2  = r*r
      rb  = c(3)*r12 + c(4)*r + c(5)*r32 + c(6)*r2
      sb  = 1.0d0 + 1.0d0/(2.0d0*c(1)*rb)
      g   = -2.0d0 * c(1) * (1.0d0+c(2)*r) * DLOG(sb)
      drb = c(3)/(2.0d0*r12) + c(4) + 1.5d0*c(5)*r12 + 2.0d0*c(6)*r
      dg  = (1.0d0+c(2)*r)*drb/(rb*rb*sb) - 2.0d0*c(1)*c(2)*DLOG(sb)
      !
      RETURN
      !
    END SUBROUTINE pwcorr
!
!
!-----------------------------------------------------------------------------
SUBROUTINE optx( rho, grho, sx, v1x, v2x )
  !---------------------------------------------------------------------------
  !! OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
  !! Present release: Mauro Boero, Tsukuba, 10/9/2002
  !
  !! rhoa = rhob = 0.5 * rho in LDA implementation
  !! grho is the SQUARE of the gradient of rho! --> gr=sqrt(grho)
  !! sx  : total exchange correlation energy at point r
  !! v1x : d(sx)/drho
  !! v2x : 1/gr*d(sx)/d(gr)
  !
  USE kind_l,   ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: small=1.D-30, smal2=1.D-10
  ! ... coefficients and exponents
  REAL(DP), PARAMETER :: o43=4.0d0/3.0d0, two13=1.259921049894873D0,      &
       two53=3.174802103936399D0, gam=0.006D0, a1cx=0.9784571170284421D0, &
       a2=1.43169D0
  REAL(DP) :: gr, rho43, xa, gamx2, uden, uu
  !
  ! ... OPTX in compact form
  !
  gr = MAX(grho,smal2)
  rho43 = rho**o43
  xa = two13*DSQRT(gr)/rho43
  gamx2 = gam*xa*xa
  uden = 1.d+00/(1.d+00+gamx2)
  uu = a2*gamx2*gamx2*uden*uden
  uden = rho43*uu*uden
  sx  = -rho43*(a1cx+uu)/two13
  v1x = o43*(sx+two53*uden)/rho
  v2x = -two53*uden/gr
  !
  RETURN
  !
END SUBROUTINE optx
!
!
!---------------------------------------------------------------
SUBROUTINE wcx( rho, grho, sx, v1x, v2x )
  !---------------------------------------------------------------
  !!  Wu-Cohen exchange (without Slater exchange):
  !!  Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)
  !
  USE kind_l,   ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: kf, agrho, s1, s2, es2, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1, x1, x2, x3, &
              dxds1, dxds2, dxds3, sx_s
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0
  REAL(DP), PARAMETER :: third=1.d0 / 3.d0, c1=0.75d0/pi ,      &
                         c2=3.093667726280136d0, c5=4.d0*third, &
                         teneightyone = 0.123456790123d0
  ! parameters of the functional
  REAL(DP), PARAMETER :: k=0.804d0, mu=0.2195149727645171d0, &
                         cwc=0.00793746933516d0
  !
  agrho = SQRT(grho)
  kf  = c2 * rho**third
  dsg = 0.5d0 / kf
  s1  = agrho * dsg / rho
  s2  = s1 * s1
  es2 = EXP(-s2)
  ds  = - c5 * s1
  !
  !   Energy
  ! x = 10/81 s^2 + (mu - 10/81) s^2 e^-s^2 + ln (1 + c s^4)
  x1 = teneightyone * s2
  x2 = (mu - teneightyone) * s2 * es2
  x3 = LOG(1.d0 + cwc * s2 * s2)
  f1 = (x1 + x2 + x3) / k
  f2 = 1.d0 + f1
  f3 = k / f2
  fx = k - f3
  exunif = - c1 * kf
  sx_s = exunif * fx
  !
  !   Potential
  dxunif = exunif * third
  dfx1 = f2 * f2
  dxds1 = teneightyone
  dxds2 = (mu - teneightyone) * es2 * (1.d0 - s2)
  dxds3 = 2.d0 * cwc * s2 / (1.d0 + cwc * s2 *s2)
  dfx = 2.d0 * s1 * (dxds1 + dxds2 + dxds3) / dfx1
  !
  v1x = sx_s + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho
  sx  = sx_s * rho
  !
  RETURN
  !
END SUBROUTINE wcx
!
!
!-----------------------------------------------------------------------
SUBROUTINE pbexsr( rho, grho, sxsr, v1xsr, v2xsr, omega, in_err )
  !---------------------------------------------------------------------
  ! INCLUDE 'cnst.inc'
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: omega
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sxsr, v1xsr, v2xsr
  INTEGER :: in_err
  !
  ! ... local variables
  !
  REAL(DP) :: rs, vx, aa, rr, ex, s2, s, d1x, d2x, fx, dsdn, dsdg
  REAL(DP), PARAMETER :: small=1.D-20, smal2=1.D-08
  REAL(DP), PARAMETER :: us=0.161620459673995492D0, ax=-0.738558766382022406D0, &
                         um=0.2195149727645171D0, uk=0.8040D0, ul=um/uk
  REAL(DP), PARAMETER :: f1=-1.10783814957303361_DP, alpha=2.0_DP/3.0_DP
  !
  ! CALL XC(RHO,EX,EC,VX,VC)
  !
  rs = rho**(1.0_DP/3.0_DP)
  vx = (4.0_DP/3.0_DP)*f1*alpha*rs
  !
  ! aa = dmax1(grho,smal2)
  aa = grho
  ! rr = rho**(-4.0_DP/3.0_DP)
  rr = 1.0_DP/(rho*rs)
  ex = ax/rr
  s2 = aa*rr*rr*us*us
  !
  s = SQRT(s2)
  IF (s > 8.3D0) THEN
     s = 8.572844D0 - 18.796223D0/s2
  ENDIF
  !
  CALL wpbe_analy_erfc_approx_grad( rho, s, omega, fx, d1x, d2x, in_err )
  !
  sxsr  = ex*fx                        ! - ex
  dsdn  = -4.D0/3.D0*s/rho
  v1xsr = vx*fx + (dsdn*d2x+d1x)*ex    ! - VX
  dsdg  = us*rr
  v2xsr = ex*1.D0/SQRT(aa)*dsdg*d2x
  !
  ! NOTE, here sx is the total energy density,
  ! not just the gradient correction energy density as e.g. in pbex()
  ! And the same goes for the potentials V1X, V2X
  !
  RETURN
  !
END SUBROUTINE pbexsr
!
!
!
!-----------------------------------------------------------------------     
      SUBROUTINE axsr( IXC, RHO, GRHO, sx, V1X, V2X, OMEGA, IN_ERR )
!-----------------------------------------------------------------------     
!*** [Per Hyldgaard, No warranties. adapted from the pbesrx version above]
!-----------------------------------------------------------------------
!
!      INCLUDE 'cnst.inc'
      !
      use kind_l, ONLY : DP
      !
      IMPLICIT NONE
      !
      !$acc routine seq
      !
      INTEGER :: IXC, IN_ERR
      REAL(DP):: RHO, GRHO, V1X, V2X, OMEGA
      REAL(DP), PARAMETER :: SMALL=1.D-20, SMAL2=1.D-08
      REAL(DP), PARAMETER :: US=0.161620459673995492D0, &
              AX=-0.738558766382022406D0, &
              UM=0.2195149727645171D0,UK=0.8040D0,UL=UM/UK
      REAL(DP), PARAMETER :: f1 = -1.10783814957303361_DP, alpha = 2.0_DP/3.0_DP
      REAL(DP):: RS, VX, FX, AA, RR, EX, S2, S, D1X, D2X, SX, DSDN, DSDG
!     ==--------------------------------------------------------------==
      
!      CALL XC(RHO,EX,EC,VX,VC)
      RS = RHO**(1.0_DP/3.0_DP)
      VX = (4.0_DP/3.0_DP)*f1*alpha*RS

!      AA    = DMAX1(GRHO,SMAL2)
      AA    = GRHO
!      RR    = RHO**(-4.0_DP/3.0_DP)
      RR    = 1.0_DP/(RHO*RS)
      EX    = AX/RR
      S2    = AA*RR*RR*US*US

      S = SQRT(S2)
      IF(S.GT.8.3D0) THEN
        S = 8.572844D0 - 18.796223D0/S2
      ENDIF
      CALL wggax_analy_erfc(RHO,S,IXC,OMEGA,FX,D1X,D2X,IN_ERR)
      sx = EX*FX        ! - EX
      DSDN = -4.D0/3.D0*S/RHO
      V1X = VX*FX + (DSDN*D2X+D1X)*EX   ! - VX
      DSDG = US*RR
      V2X = EX*1.D0/SQRT(AA)*DSDG*D2X

!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE axsr
!
!
!-----------------------------------------------------------------------     
      SUBROUTINE wggax_analy_erfc(rho,s,nggatyp,omega,Fx_wgga, &
                                  dfxdn,dfxds,in_err)
!--------------------------------------------------------------------
!
!     Short-ranged wGGA Enhancement Factor (from erfc, analytical with
!     gradients)
!     
!     This codes the HJS analytical-xc-hole idea, JCP 128, 194 105 (2008).
!
!     Copyright Per Hyldgaard, GPL, No warranty, 2016-
!
!     Inputs:
!     rho     - electron density
!     s       - scaled graident of electron density
!     omega = default, short-range screening parameter
!     nggatyp = 1 is refitted to PBEx -- via analytical hole
!     nggatyp = 2 is fitted to PBEsolx - via analytical hole
!     nggatyp = 3 is fitted to cx13 - via analytical hole
!     nggatyp = 4 is fitted to rPW86 - via analytical hole
!     nggatyp = 5 is fitted to PHexplore - via analytical hole
!     nggatyp = 6 is fitted to test-reserve - via analytical hole
!     ....
!     nggatyp = 7 fitted to PBEx -- HJS fit to analytical hole
!     nggatyp = 8 fitted to PBEsolx -- HJS fit to analytical hole
!
!     Returns:
!     Fx_wgga - Analytic version of short-range gga enhancement factor
!     dfxdn   - Derivative of Fx_wgga with respect to density rho
!     dfxds   - Derivative of Fx_wgga with respect to scaled gradient s.
!
!--------------------------------------------------------------------

      use kind_l,               only: DP
      
      Implicit None
      
      !$acc routine seq
      
      REAL(DP), PARAMETER :: pi=3.14159265358979323846d0

      Real(dp) :: rho,s,omega,Fx_wgga,dfxdn,dfxds
      integer :: in_err
      integer :: nggatyp

      Real(dp) :: Abar,B,C,D,E
      parameter (Abar = 0.757211D0, B = -0.106364D0, C = -0.118649D0)
      parameter (D =  0.609650D0, E = -0.0477963D0)

      Real(dp) :: One, Two, Three, Four, Five, Six, Seven, Eight, Nine
      parameter(One=1.0D0,Two=2.0D0,Three=3.0D0,Four=4.0D0,Five=5.0D0)
      parameter(Six=6.0D0,Seven=7.0D0,Eight=8.0D0,Nine=9.0D0)

      Real(dp) :: f12, f32, f52, f72, f13, f23, f25, f45, f65, f54, f18, f38, f58
      parameter(f12=0.5D0,f32=Three*f12,f52=Five*f12,f72=Seven*f12)
      parameter (f13=One/Three,f23=Two*f13)
      parameter (f25=Two/Five, f45=Two*f25, f65=Three*f25,f54=Five/Four)
      parameter(f18=One/Eight,f38=Three*f18,f58=Five*f18)
      Real(dp) :: pi2, pisqrt
      parameter( pi2=pi**Two, pisqrt=pi**f12)

      Real(dp) :: s0val,s0sq
      parameter (s0val=2.0D0,s0sq=s0val*s0val)

      Real(dp) :: coef1,coef2,coef3
      parameter (coef1=(Four/Nine)*B,coef2=(Four/Nine),coef3=Eight/Nine)

!      Real(dp) ::, dimension(6) :: a2,a3,a4,a5,a6,a7
!      Real(dp) ::, dimension(6) :: b1,b2,b3,b4
!      Real(dp) ::, dimension(6) :: b5,b6,b7,b8,b9
      Real(dp), dimension(8) :: a2,a3,a4,a5,a6,a7
      Real(dp), dimension(8) :: b1,b2,b3,b4
      Real(dp), dimension(8) :: b5,b6,b7,b8,b9

!     HJS-type/ pbe-x pbesol-x cx13 rPW86 explore TestReserve / 
!  Last two columns are parameters from the original HJS fits for PBE/PBesol
      data a2 /  0.0154999D0,   4.58809D-3,   2.43873D-3,   6.00962D-7, &
                 0.0139602D0,   4.56198D-3, 0.0159941D0,   0.0047333D0 /
      data a3 / -0.0361006D0,  -8.57842D-3,  -4.15263D-3,  0.0402647D0, &
                -0.0363194D0,  -8.70003D-3, 0.0852995D0,   0.0403304D0 /
      data a4 /  0.0379567D0,   7.29562D-3,   2.58261D-3, -0.0353219D0, &
                 0.0469970D0,   7.36958D-3, -0.1603680D0,  -0.0574615D0 /
      data a5 / -0.0186715D0,  -3.20195D-3,   1.23940D-6,  0.0116112D0, &
                -0.0317508D0,  -3.02436D-3, 0.1526450D0,   0.0435395D0 /
      data a6 /  1.74264D-3,   6.04936D-4,  -7.58225D-4,  -1.55532D-4, &
                 8.26494D-3,   3.86773D-4, -0.0971263D0,  -0.0216251D0 /
      data a7 /  1.90765D-3,   2.16112D-5,   2.76378D-4,   5.03561D-5, &
                 1.45383D-3,   9.43741D-5, 0.0422061D0,   0.0063721D0 /
      data b1 / -2.7062566D0, -2.1449453D0, -2.2030319D0, -1.8779594D0, &
                -3.0623921D0, -2.2089330D0, 5.3331900D0,   8.5205600D0 /
      data b2 /  3.3316842D0,  2.0901104D0,  2.1759315D0,  1.5198811D0, &
                 4.3601225D0,  2.1968353D0, -12.4780000D0, -13.9885000D0 /
      data b3 / -2.3871819D0, -1.1935421D0, -1.2997841D0, -0.5383109D0, &
                -3.7025379D0, -1.2662249D0, 11.0988000D0,   9.2858300D0 /
      data b4 /  1.1197810D0,  0.4476392D0,  0.5347267D0,  0.1352399D0, &
                 2.0707006D0,  0.4689964D0, -5.1101300D0,  -3.2728700D0 /
      data b5 / -0.3606638D0, -0.1172367D0, -0.1588798D0, -0.0428465D0, &
                -0.7578009D0, -0.1165714D0, 1.7146800D0,   0.8434990D0 /
      data b6 /  0.0841990D0,  0.0231625D0,  0.0367329D0,  0.0117903D0, &
                 0.1666493D0,  0.0207188D0, -0.6103800D0,  -0.2355430D0 /
      data b7 / -0.0114719D0,  -3.52782D-3,  -7.73178D-3,   3.37908D-3, &
                -0.0178278D0,  -2.97718D-3, 0.3075550D0,   0.0847074D0 /
      data b8 /  1.69283D-3,   5.39942D-4,   1.26670D-3,  -4.93453D-5, &
                 7.58236D-3,   5.98226D-4, -0.0770547D0,  -0.0171561D0 /
      data b9 /  1.50540D-3,   1.58065D-5,   8.04175D-7,   7.09955D-6, &
                 7.01937D-5,   4.67972D-6, 0.0334840D0,   0.0050552D0 /
!     End HJS-type parameters: JPCM 34, 025902 (2022)

      integer :: i
      Real(dp) :: s2,s3,s4,s5,s6,s7,s8,s9
      Real(dp) :: hnom,hdenom,dhnomds,dhdenomds
      Real(dp) :: hs,dhds

      Real(dp) :: lam,eta,zeta,dzetads
      Real(dp) :: xi,phi,psi
      Real(dp) :: alpha,dalphadn,dalphads
      Real(dp) :: beta,dbetadn,dbetads

      Real(dp) :: chi,dchidn,dchids
      Real(dp) :: chiP1, chiP1p, dchiP1dn, dchiP1ds
      Real(dp) :: chiP2, chiP2p, dchiP2dn, dchiP2ds
      Real(dp) :: chiP3, chiP3p, dchiP3dn, dchiP3ds

      Real(dp) :: dampfac,cfbars,dcfbards
      Real(dp) :: egbars,degbards

      Real(dp) :: kf,ny,ny2,dnydn
      
      kf    = (Three*pi2*rho) ** f13
      ny= omega/kf
      ny2=ny*ny
      dnydn= -f13*ny/rho
      
!      if ((nggatyp.ge.1).or.(nggatyp.le.6)) then
      if ((nggatyp>=1).or.(nggatyp<=8)) then
         i = nggatyp
      else
         in_err = 3  ! wgga_analy_erfc: yet to be coded Wcx part
         return
      endif

      s2=s*s
      s3=s2*s
      s4=s2*s2
      s5=s2*s3
      s6=s3*s3
      s7=s4*s3
      s8=s4*s4
      s9=s4*s5

      hnom = a2(i)*s2+a3(i)*s3+a4(i)*s4
      hnom = hnom +a5(i)*s5+a6(i)*s6+a7(i)*s7
      dhnomds = Two*a2(i)*s+Three*a3(i)*s2+Four*a4(i)*s3
      dhnomds = dhnomds +Five*a5(i)*s4+Six*a6(i)*s5+Seven*a7(i)*s6
      hdenom = One+b1(i)*s+b2(i)*s2+b3(i)*s3+b4(i)*s4+b5(i)*s5
      hdenom = hdenom +b6(i)*s6+b7(i)*s7+b8(i)*s8+b9(i)*s9
      dhdenomds = b1(i)+Two*b2(i)*s+Three*b3(i)*s2+Four*b4(i)*s3
      dhdenomds = dhdenomds +Five*b5(i)*s4+Six*b6(i)*s5
      dhdenomds = dhdenomds +Seven*b7(i)*s6+Eight*b8(i)*s7
      dhdenomds = dhdenomds +Nine*b9(i)*s8
      hs=hnom/hdenom
      dhds=dhnomds/hdenom-hnom*dhdenomds/hdenom/hdenom

      zeta = s2*hs
      dzetads = Two*s*hs+s2*dhds
      lam=zeta+D
      eta=zeta+Abar

      dampfac = (One+s2/s0sq)
      cfbars=C-s2/dampfac/27.0D0 - zeta/Two
      dcfbards=-Two*s/dampfac/dampfac/27.0D0 - dzetads/Two

      egbars=-f25*cfbars*lam-(f45/Three)*B*lam**Two
      egbars=egbars-f65*Abar*lam**Three
      egbars=egbars-f45*pisqrt*lam**f72
      egbars=egbars-Three*f45*(zeta**f12-eta**f12)*lam**f72

      degbards=-f25*(dcfbards*lam+cfbars*dzetads)
      degbards=degbards-f23*f45*B*dzetads*lam
      degbards=degbards-Nine*f25*Abar*dzetads*lam**Two
      degbards=degbards-Seven*f25*pisqrt*dzetads*lam**f52
      degbards=degbards-Seven*f65*dzetads*(zeta**f12-eta**f12)*lam**f52
      degbards=degbards-f65*dzetads*(zeta**(-f12)-eta**(-f12))*lam**f72

      phi=(lam+ny2)**(f12)
      psi=(eta+ny2)**(f12)
      xi=(zeta+ny2)**(f12)

      alpha=Two*ny*(xi-psi)
      dalphadn=dnydn*Two*(xi-psi+ny2/xi-ny2/psi)
      dalphads=dzetads*(ny/xi-ny/psi)

      beta=Two*zeta*log((ny+xi)/(ny+phi))-Two*eta*log((ny+psi)/(ny+phi))
      dbetadn=Abar/phi
      dbetadn=dbetadn+zeta/xi-eta/psi
      dbetadn=Two*dbetadn*dnydn
      dbetads=Abar/(ny+phi)/phi+Two*log((ny+xi)/(ny+psi))
      dbetads=dbetads+zeta/(ny+xi)/xi-eta/(ny+psi)/psi
      dbetads=dbetads*dzetads

      chi=ny/phi
      dchidn=dnydn*lam/phi**Three
      dchids=-f12*chi*dzetads/phi/phi

      chiP1=One-chi
      chiP1p=-One
      dchiP1dn=chiP1p*dchidn
      dchiP1ds=chiP1p*dchids

      chiP2=One-f32*chi+f12*chi**Three
      chiP2p=-f32*(One-chi**Two)
      dchiP2dn=chiP2p*dchidn
      dchiP2ds=chiP2p*dchids

      chiP3=One-Five*f38*chi+f54*chi**Three-f38*chi**Five
      chiP3p=-Five*f38+Three*f54*chi**Two-Five*f38*chi**Four
      dchiP3dn=chiP3p*dchidn
      dchiP3ds=chiP3p*dchids

      Fx_wgga=Abar
      Fx_wgga=Fx_wgga-coef1*chiP1/lam
      Fx_wgga=Fx_wgga-coef2*cfbars*chiP2/lam**Two
      Fx_wgga=Fx_wgga-coef3*egbars*chiP3/lam**Three
      Fx_wgga=Fx_wgga+alpha+beta

      dfxdn=-coef1*dchiP1dn/lam
      dfxdn=dfxdn-coef2*cfbars*dchiP2dn/lam**Two
      dfxdn=dfxdn-coef3*egbars*dchiP3dn/lam**Three
      dfxdn=dfxdn+dalphadn+dbetadn

      dfxds=-coef1*(dchiP1ds/lam-chiP1*dzetads/lam**Two)
      dfxds=dfxds-coef2*(dcfbards*chiP2+cfbars*dchiP2ds)/lam**Two
      dfxds=dfxds+coef2*cfbars*chiP2*(Two*dzetads/lam**Three)
      dfxds=dfxds-coef3*(degbards*chiP3+egbars*dchiP3ds)/lam**Three
      dfxds=dfxds+coef3*egbars*chiP3*(Three*dzetads/lam**Four)
      dfxds=dfxds+dalphads+dbetads

!     ==--------------------------------------------------------------==
      RETURN
      END SUBROUTINE wggax_analy_erfc
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
SUBROUTINE rPW86( rho, grho, sx, v1x, v2x )
  !---------------------------------------------------------------------
  !! PRB 33, 8800 (1986) and J. Chem. Theory comp. 5, 2754 (2009).
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: s, s_2, s_3, s_4, s_5, s_6, fs, grad_rho, df_ds
  REAL(DP), PARAMETER :: a=1.851_DP, b=17.33_DP, c=0.163_DP, &
                         s_prefactor=6.18733545256027_DP,    &
                         Ax=-0.738558766382022_DP, four_thirds=4._DP/3._DP
  !
  grad_rho = SQRT(grho)
  !
  s = grad_rho/(s_prefactor*rho**(four_thirds))
  !
  s_2 = s**2
  s_3 = s_2 * s
  s_4 = s_2**2
  s_5 = s_3 * s_2
  s_6 = s_2 * s_4
  !
  ! Calculation of energy
  fs = (1 + a*s_2 + b*s_4 + c*s_6)**(1._DP/15._DP)
  sx = Ax * rho**(four_thirds) * (fs -1.0_DP)
  !
  ! Calculation of the potential
  df_ds = (1._DP/(15._DP*fs**(14.0_DP)))*(2*a*s + 4*b*s_3 + 6*c*s_5)
  !
  v1x = Ax*(four_thirds)*(rho**(1._DP/3._DP)*(fs -1.0_DP) &
        -grad_rho/(s_prefactor * rho)*df_ds)
  !
  v2x = Ax * df_ds/(s_prefactor*grad_rho)
  !
END SUBROUTINE rPW86
!
!
!-----------------------------------------------------------------
SUBROUTINE c09x( rho, grho, sx, v1x, v2x )
  !---------------------------------------------------------------
  !! Cooper '09 exchange for vdW-DF (without Slater exchange):
  !! V. R. Cooper, Phys. Rev. B 81, 161104(R) (2010)
  !
  !! Developed thanks to the contribution of
  !! Ikutaro Hamada - ikutaro@wpi-aimr.tohoku.ac.jp
  !! WPI-Advanced Institute of Materials Research, Tohoku University
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: kf, agrho, s1, s2, sx_s, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1, dfx2
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0
  REAL(DP), PARAMETER :: third=1._DP/3._DP, c1=0.75_DP/pi, &
                         c2=3.093667726280136_DP, c5=4._DP*third
  ! parameters of the functional
  REAL(DP) :: kappa, mu, alpha
  DATA kappa / 1.245_DP  /, &
       mu    / 0.0617_DP /, &
       alpha / 0.0483_DP /
  !
  agrho = SQRT(grho)
  kf = c2 * rho**third
  dsg = 0.5_DP / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  !
  ! ... Energy
  !
  f1 = EXP( - alpha * s2 )
  f2 = EXP( - alpha * s2 / 2.0_DP )
  f3 = mu * s2 * f1
  fx = f3 + kappa * ( 1.0_DP - f2 )
  exunif = - c1 * kf
  sx_s = exunif * fx
  !
  ! ... Potential
  !
  dxunif = exunif * third
  dfx1 = 2.0_DP * mu * s1 * ( 1.0_DP - alpha * s2 ) * f1
  dfx2 = kappa * alpha * s1 * f2
  dfx = dfx1 + dfx2
  v1x = sx_s + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho
  !
  sx  = sx_s * rho
  !
  RETURN
  !
END SUBROUTINE c09x
!
!
!---------------------------------------------------------------
SUBROUTINE sogga( rho, grho, sx, v1x, v2x )
  !-------------------------------------------------------------
  !! SOGGA exchange
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  ! input: charge and abs gradient
  ! output: energy and potential
  !
  ! ... local variables
  !
  REAL(DP) :: rho43, xs, xs2, dxs2_drho, dxs2_dgrho2
  REAL(DP) :: CX, denom, C1, C2, Fso, Fpbe, ex, Fx, dFx_dxs2, dex_drho
  !
  REAL(DP), PARAMETER :: one  = 1.0_DP, two   = 2.0_DP, three = 3.0_DP,        &
  &                      four = 4.0_DP, eight = 8.0_DP,                        &
  &                      f13 = one/three, f23 = two/three,   f43 = four/three, &
  &                      f34 = three/four,f83 = eight/three, f12 = one/two
  !
  REAL(DP), PARAMETER :: mu=0.12346_DP, kapa=0.552_DP
  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0
  !
  ! Cx LDA
  CX    =  f34 * (three/pi)**f13
  denom =  four * (three*pi**two)**f23
  C1    =  mu / denom
  C2    =  mu / (kapa * denom)
  !
  rho43 = rho**f43
  xs  = grho / rho43
  xs2 = xs * xs
  !
  dxs2_drho   = -f83 * xs2 / rho
  dxs2_dgrho2 = one /rho**f83
  !
  ex       = - CX * rho43
  dex_drho = - f43 * CX * rho**f13
  !
  Fso  = kapa * (one - EXP(-C2*xs2))
  Fpbe = C1 * xs2 / (one + C2*xs2)
  !
  Fx       = f12 * (Fpbe + Fso)
  dFx_dxs2 = f12 * (C1 / ((one + C2*xs2)**2) + C1*EXP(-C2*xs2))
  !
  !   Energy
  sx = Fx * ex
  !
  !   Potential
  v1x = dex_drho * Fx  +  ex * dFx_dxs2 * dxs2_drho
  v2x = two * ex * dFx_dxs2 * dxs2_dgrho2
  !
END SUBROUTINE sogga
!
!
!-------------------------------------------------------------------------
SUBROUTINE pbexgau( rho, grho, sxsr, v1xsr, v2xsr, alpha_gau )
  !-----------------------------------------------------------------------
  !! PBEX gaussian.
  !
  USE kind_l,  ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: alpha_gau
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sxsr, v1xsr, v2xsr
  !
  ! ... local variables
  !
  REAL(DP) :: rs, vx, aa, rr, ex, s2, s, d1x, d2x, fx, dsdn, dsdg
  !
  REAL(DP), PARAMETER :: small=1.D-20, smal2=1.D-08
  REAL(DP), PARAMETER :: us=0.161620459673995492D0, ax=-0.738558766382022406D0, &
                         um=0.2195149727645171D0, uk=0.8040D0, ul=um/uk
  REAL(DP), PARAMETER :: f1 = -1.10783814957303361_DP, alpha = 2.0_DP/3.0_DP
  !
  rs = rho**(1.0_DP/3.0_DP)
  vx = (4.0_DP/3.0_DP)*f1*alpha*rs
  aa = grho
  rr = 1.0_DP/(rho*rs)
  ex = ax/rr
  ! AX is 3/4/PI*(3*PI*PI)**(1/3). This is the same as -c1*c2 in pbex().
  s2 = aa*rr*rr*us*us
  s = SQRT(s2)
  IF (s > 10.D0) THEN
     s = 10.D0
  ENDIF
  CALL pbe_gauscheme( rho, s, alpha_gau, fx, d1x, d2x )
  sxsr = ex*fx                        ! - EX
  dsdn = -4.D0/3.D0*s/rho
  v1xsr = vx*fx + (dsdn*d2x+d1x)*ex   ! - VX
  dsdg = us*rr
  v2xsr = ex*1.D0/SQRT(aa)*dsdg*d2x
  !
  ! NOTE, here sx is the total energy density,
  ! not just the gradient correction energy density as e.g. in pbex()
  ! And the same goes for the potentials V1X, V2X
  !
  RETURN
  !
END SUBROUTINE pbexgau
    !
    !-----------------------------------------------------------------------
SUBROUTINE pbe_gauscheme( rho, s, alpha_gau, Fx, dFxdr, dFxds )
       !--------------------------------------------------------------------
       !
       USE kind_l, ONLY: DP
       !
       IMPLICIT NONE
       !
       !$acc routine seq
       !
       REAL(dp) :: rho,s,alpha_gau,Fx,dFxdr,dFxds
       ! input: charge and squared gradient and alpha_gau
       ! output: GGA enhancement factor of gau-PBE
       ! output: d(Fx)/d(s), d(Fx)/d(rho)
       !
       REAL(dp) :: Kx, Nx
       ! PBE96 GGA enhancement factor
       ! GGA enhancement factor of Gaussian Function
       !
       REAL(dp) :: bx, cx, PI, sqrtpial, Prefac, term_PBE, Third, KsF
       REAL(dp) :: d1sdr, d1Kxds, d1Kxdr, d1bxdr, d1bxds, d1bxdKx, &
              d1Nxdbx,d1Nxdr, d1Nxds
       !
       REAL(dp) :: Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
       !
       SAVE Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
       DATA Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten &
         / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /
       !
       REAL(dp) :: k , mu
       DATA k / 0.804d0 / , mu / 0.21951d0 /
       ! parameters of PBE functional
       !
       Third = One/Three
       PI = ACOS(-One)
       KsF = (Three*PI*PI*rho)**Third
       sqrtpial = SQRT(PI/alpha_gau)
       Prefac = Two * SQRT(PI/alpha_gau) / Three
       !
       ! PBE96 GGA enhancement factor part
       term_PBE = One / (One + s*s*mu/k)
       Kx =  One + k - k * term_PBE
       !
       ! GGA enhancement factor of Gaussian Function part
       bx = SQRT(Kx*alpha_gau) / KsF
       !
       ! cx = exp(-One/Four/bx/bx) - One
       IF (ABS(One/bx/bx) < 1.0D-4) THEN
          cx = TayExp(-One/bx/bx)
       ELSE
          cx = EXP(-One/bx/bx) - One
       ENDIF
       !
       Nx = bx * Prefac * ( SQRT(PI) * ERF(One/bx) + & 
        (bx - Two*bx*bx*bx)*cx - Two*bx )
       !
       ! for convergence
       IF (ABS(Nx) < 1.0D-15) THEN
         Nx = Zero
       ELSEIF ((One - ABS(Nx)) < 1.0D-15) THEN
         Nx = One
       ELSE
         Nx = Nx
       ENDIF
       ! for convergence end
       !
       Fx =  Kx * Nx
       !
       ! 1st derivatives
       d1sdr = - Four / Three * s / rho
       !
       d1Kxds = Two * s * mu * term_PBE * term_PBE
       d1Kxdr = d1Kxds * d1sdr
       d1bxdKx = bx / (Two* Kx)
       !
       d1bxdr = - bx /(Three*rho) + d1Kxdr * d1bxdKx
       !
       d1bxds =  d1bxdKx * d1Kxds
       !
       d1Nxdbx =  Nx/bx - Prefac * bx * Three * &
                   ( cx*(One + Two*bx*bx) + Two )
       !
       d1Nxdr = d1Nxdbx * d1bxdr
       d1Nxds = d1Nxdbx * d1bxds
       !
       dFxdr = d1Kxdr * Nx + Kx * d1Nxdr
       dFxds = d1Kxds * Nx + Kx * d1Nxds
       !
       RETURN
       !
END SUBROUTINE pbe_gauscheme
!
!
!-------------------------------------------------
FUNCTION TayExp(X)
  !-------------------------------------------
  USE kind_l,   ONLY: DP
  IMPLICIT NONE
  !$acc routine seq
  REAL(DP), INTENT(IN) :: X
  REAL(DP) :: TAYEXP
  INTEGER :: NTERM,I
  REAL(DP) :: SUMVAL,IVAL,COEF
  PARAMETER (NTERM=16)
  !
  SUMVAL = X
  IVAL = X
  COEF = 1.0D0
  DO 10 I = 2, NTERM
     COEF = COEF * I
     IVAL = IVAL * (X / COEF)
     SUMVAL = SUMVAL + IVAL
10     CONTINUE
  TAYEXP = SUMVAL
  !
  RETURN
  !
END FUNCTION TayExp
!
!
!
!-------------------------------------------------------------------------
SUBROUTINE PW86( rho, grho, sx, v1x, v2x )
  !-----------------------------------------------------------------------
  !! Perdew-Wang 1986 exchange gradient correction: PRB 33, 8800 (1986)
  !
  USE kind_l,  ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: s, s_2, s_3, s_4, s_5, s_6, fs, grad_rho, df_ds
  REAL(DP), PARAMETER :: a=1.296_DP, b=14._DP, c=0.2_DP,   &
                         s_prefactor=6.18733545256027_DP, &
                         Ax=-0.738558766382022_DP, four_thirds=4._DP/3._DP
  !
  grad_rho = SQRT(grho)
  !
  s = grad_rho / ( s_prefactor*rho**(four_thirds) )
  !
  s_2 = s**2
  s_3 = s_2 * s
  s_4 = s_2**2
  s_5 = s_3 * s_2
  s_6 = s_2 * s_4
  !
  ! Calculation of energy
  fs = (1 + a*s_2 + b*s_4 + c*s_6)**(1._DP/15._DP)
  sx = Ax * rho**(four_thirds) * (fs-1._DP)
  !
  ! Calculation of the potential
  df_ds = (1._DP/(15._DP*fs**(14._DP)))*(2*a*s + 4*b*s_3 + 6*c*s_5)
  !
  v1x = Ax*(four_thirds)*( rho**(1._DP/3._DP)*(fs-1._DP) &
            -grad_rho/(s_prefactor * rho)*df_ds )
  !
  v2x = Ax * df_ds/(s_prefactor*grad_rho)
  !
END SUBROUTINE PW86
!
!
!-----------------------------------------------------------------------
SUBROUTINE becke86b( rho, grho, sx, v1x, v2x )
  !-----------------------------------------------------------------------
  !! Becke 1986 gradient correction to exchange
  !! A.D. Becke, J. Chem. Phys. 85 (1986) 7184
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: arho, agrho
  REAL(DP) :: sgp1, sgp1_45, sgp1_95
  REAL(DP) :: rdg2_43, rdg2_73, rdg2_83, rdg2_4, rdg4_5
  REAL(DP), PARAMETER :: beta=0.00375_DP, gamma=0.007_DP
  !
  arho  = 0.5_DP  * rho
  agrho = 0.25_DP * grho
  !
  rdg2_43 = agrho / arho**(4d0/3d0)
  rdg2_73 = rdg2_43 / arho
  rdg2_83 = rdg2_43 * rdg2_43 / agrho
  rdg2_4 = rdg2_43 * rdg2_83 / agrho
  rdg4_5 = rdg2_73 * rdg2_83
  !
  sgp1 = 1d0 + gamma * rdg2_83
  sgp1_45 = sgp1**(-4d0/5d0)
  sgp1_95 = sgp1_45 / sgp1
  !
  sx  = -2d0 * beta * agrho / arho**(4d0/3d0) * sgp1_45
  v1x = -beta * (-4d0/3d0*rdg2_73*sgp1_45 + 32d0/15d0*gamma*rdg4_5*sgp1_95)
  v2x = -beta * (sgp1_45*rdg2_43/agrho - 4d0/5d0 *gamma*rdg2_4*sgp1_95)
  !
END SUBROUTINE becke86b
!
!
!---------------------------------------------------------------
SUBROUTINE b86b( rho, grho, iflag, sx, v1x, v2x )
  !-------------------------------------------------------------
  !! Becke exchange (without Slater exchange):
  !! iflag=1: A. D. Becke, J. Chem. Phys. 85, 7184 (1986) (B86b)
  !! iflag=2: J. Klimes, Phys. Rev. B 83, 195131 (2011). (OptB86b)
  !! iflag=3: I. Hamada, Phys. Rev. B 89, 121103(R) (B86R)
  !! iflag=4: D. Chakraborty, K. Berland, and T. Thonhauser, JCTC 16, 5893 (2020)
  !
  !! Ikutaro Hamada - HAMADA.Ikutaro@nims.go.jp
  !! National Institute for Materials Science
  !
  USE kind_l,     ONLY : DP
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  INTEGER, INTENT(IN) :: iflag
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: kf, agrho, s1, s2, sx_s, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  REAL(DP) :: dxunif, dfx, f1, f2, f3, dfx1
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  REAL(DP), PARAMETER :: pi=3.14159265358979323846d0
  REAL(DP), PARAMETER :: third=1._DP/3._DP, c1=0.75_DP/pi, &
                         c2=3.093667726280136_DP, c5=4._DP*third
  ! parameters of the functional
  REAL(DP) :: k(4), mu(4)
  DATA k / 0.5757_DP, 1.0000_DP, 0.711357_DP, 0.58_DP /, &
       mu/ 0.2449_DP, 0.1234_DP, 0.1234_DP, 0.12345679012345679_DP /
  !
  agrho = SQRT(grho)
  kf = c2 * rho**third
  dsg = 0.5_DP / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  !
  ! ... Energy
  !
  f1 = mu(iflag)*s2
  f2 = 1._DP + mu(iflag)*s2/k(iflag)
  f3 = f2**(4._DP/5._DP)
  fx = f1/f3
  exunif = - c1 * kf
  sx_s = exunif * fx
  !
  ! ... Potential
  !
  dxunif = exunif * third
  dfx1 = 1._DP + (1._DP/5._DP)*mu(iflag)*s2 / k(iflag)
  dfx  = 2._DP * mu(iflag) * s1 * dfx1 / (f2 * f3)
  v1x = sx_s + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho
  sx = sx_s * rho
  !
  RETURN
  !
END SUBROUTINE b86b
!
!
!-----------------------------------------------------------------------
SUBROUTINE cx13( rho, grho, sx, v1x, v2x )
  !-----------------------------------------------------------------------
  !! The new exchange partner for a vdW-DF1-cx suggested
  !! by K. Berland and P. Hyldgaard, see PRB 89, 035412 (2014),
  !! to test the plasmon nature of the vdW-DF1 inner functional.
  !
  USE kind_l, ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT) :: sx, v1x, v2x
  !
  ! ... local variables
  !
  REAL(DP) :: s, s_2, s_3, s_4, s_5, s_6, fs, fs_rPW86, df_rPW86_ds, grad_rho, df_ds
  REAL(DP), PARAMETER :: alp=0.021789_DP, beta=1.15_DP, a=1.851_DP, b=17.33_DP, &
                         c=0.163_DP, mu_LM=0.09434_DP,    &
                         s_prefactor=6.18733545256027_DP, &
                         Ax = -0.738558766382022_DP, four_thirds = 4._DP/3._DP
  !
  grad_rho = SQRT(grho)
  !
  s = grad_rho/(s_prefactor*rho**(four_thirds))
  !
  s_2 = s*s
  s_3 = s_2 * s
  s_4 = s_2 * s_2
  s_5 = s_3 * s_2
  s_6 = s_2 * s_2 *s_2
  !
  ! ... Energy
  fs_rPW86 = (1._DP + a*s_2 + b*s_4 + c*s_6)**(1._DP/15._DP)
  fs = 1._DP/(1._DP + alp*s_6) * (1._DP + mu_LM *s_2) &
       + alp*s_6/(beta+alp*s_6)*fs_rPW86
  !
  sx = Ax * rho**(four_thirds) * (fs-1._DP)
  !
  ! ... Potential
  df_rPW86_ds = (1._DP/(15._DP*fs_rPW86**(14._DP)))*(2*a*s + 4*b*s_3 + 6*c*s_5)
  !
  df_ds = 1._DP/(1._DP+alp*s_6)**2*( 2._DP*mu_LM*s*(1._DP+alp*s_6) &
            - 6._DP*alp*s_5*( 1._DP+mu_LM*s_2) )                   &
          + alp*s_6/(beta+alp*s_6)*df_rPW86_ds                     &
          + 6._DP*alp*s_5*fs_rPW86/(beta+alp*s_6)*(1._DP-alp*s_6/(beta + alp*s_6))
  !
  v1x = Ax*(four_thirds)*(rho**(1._DP/3._DP)*(fs-1._DP) &
        -grad_rho/(s_prefactor * rho)*df_ds)
  v2x = Ax * df_ds/(s_prefactor*grad_rho)
  !
END SUBROUTINE cx13
!
!
!
! ===========> SPIN <===========
!
!-----------------------------------------------------------------------
SUBROUTINE becke88_spin( rho_up, rho_dw, grho_up, grho_dw, sx_up, sx_dw, v1x_up, v1x_dw, v2x_up, v2x_dw )
  !-----------------------------------------------------------------------
  !! Becke exchange: A.D. Becke, PRA 38, 3098 (1988) - Spin polarized case
  !
  USE kind_l,    ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho_up, rho_dw
  !! charge
  REAL(DP), INTENT(IN) :: grho_up, grho_dw
  !! gradient
  REAL(DP), INTENT(OUT) :: sx_up, sx_dw
  !! the up and down energies
  REAL(DP), INTENT(OUT) :: v1x_up, v1x_dw
  !! first part of the potential
  REAL(DP), INTENT(OUT) :: v2x_up, v2x_dw
  !! second part of the potential
  !
  ! ... local variables
  !
  !INTEGER :: is
  REAL(DP), PARAMETER :: beta = 0.0042_DP, third = 1._DP/3._DP
  REAL(DP) :: rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee
  !
  !
  !DO is = 1, 2
     rho13 = rho_up**third
     rho43 = rho13**4
     xs  = SQRT(grho_up) / rho43
     xs2 = xs * xs
     sa2b8 = SQRT(1.0d0 + xs2)
     shm1  = LOG(xs + sa2b8)
     dd  = 1.0d0 + 6.0d0 * beta * xs * shm1
     dd2 = dd * dd
     ee = 6.0d0 * beta * xs2 / sa2b8 - 1.d0
     sx_up  = grho_up / rho43 * (-beta/dd)
     v1x_up = -(4.d0/3.d0) * xs2 * beta * rho13 * ee / dd2
     v2x_up = beta * (ee-dd) / (rho43*dd2)

     rho13 = rho_dw**third
     rho43 = rho13**4
     xs  = SQRT(grho_dw) / rho43
     xs2 = xs * xs
     sa2b8 = SQRT(1.0d0 + xs2)
     shm1  = LOG(xs + sa2b8)
     dd  = 1.0d0 + 6.0d0 * beta * xs * shm1
     dd2 = dd * dd
     ee = 6.0d0 * beta * xs2 / sa2b8 - 1.d0
     sx_dw  = grho_dw / rho43 * (-beta/dd)
     v1x_dw = -(4.d0/3.d0) * xs2 * beta * rho13 * ee / dd2
     v2x_dw = beta * (ee-dd) / (rho43*dd2)
  !ENDDO
  !
  RETURN
  !
END SUBROUTINE becke88_spin
!
!
!-----------------------------------------------------------------------------
SUBROUTINE wpbe_analy_erfc_approx_grad( rho, s, omega, Fx_wpbe, d1rfx, d1sfx, in_err )
      !-----------------------------------------------------------------------
      !! wPBE Enhancement Factor (erfc approx.,analytical, gradients).
      !
      USE kind_l,    ONLY: DP
      !
      IMPLICIT NONE
      !
      !$acc routine seq
      !
      REAL(DP) rho,s,omega,Fx_wpbe,d1sfx,d1rfx
      INTEGER in_err
      !
      REAL(DP) f12,f13,f14,f18,f23,f43,f32,f72,f34,f94,f1516,f98
      REAL(DP) pi,pi2,pi_23,srpi
      REAL(DP) Three_13
      !
      REAL(DP) ea1,ea2,ea3,ea4,ea5,ea6,ea7,ea8
      REAL(DP) eb1
      REAL(DP) A,B,C,D,E
      REAL(DP) Ha1,Ha2,Ha3,Ha4,Ha5
      REAL(DP) Fc1,Fc2
      REAL(DP) EGa1,EGa2,EGa3
      REAL(DP) EGscut,wcutoff,expfcutoff
      !
      REAL(DP) xkf, xkfrho
      REAL(DP) w,w2,w3,w4,w5,w6,w7,w8
      REAL(DP) d1rw
      REAL(DP) A2,A3,A4,A12,A32,A52,A72
      REAL(DP) X
      REAL(DP) s2,s3,s4,s5,s6
      !
      REAL(DP) H,F
      REAL(DP) Hnum,Hden,d1sHnum,d1sHden
      REAL(DP) d1sH,d1sF
      REAL(DP) G_a,G_b,EG
      REAL(DP) d1sG_a,d1sG_b,d1sEG
      !
      REAL(DP) Hsbw,Hsbw2,Hsbw3,Hsbw4,Hsbw12,Hsbw32,Hsbw52,Hsbw72
      REAL(DP) DHsbw,DHsbw2,DHsbw3,DHsbw4,DHsbw5
      REAL(DP) DHsbw12,DHsbw32,DHsbw52,DHsbw72,DHsbw92
      REAL(DP) d1sHsbw,d1rHsbw
      REAL(DP) HsbwA94,HsbwA9412
      REAL(DP) HsbwA942,HsbwA943,HsbwA945
      REAL(DP) piexperf,expei
      REAL(DP) piexperfd1,expeid1
      REAL(DP) d1spiexperf,d1sexpei
      REAL(DP) d1rpiexperf,d1rexpei
      REAL(DP) expei1,expei2,expei3,expei4
      REAL(DP) exint
      !
      REAL(DP) DHs,DHs2,DHs3,DHs4,DHs72,DHs92,DHsw,DHsw2,DHsw52,DHsw72
      REAL(DP) d1sDHs,d1rDHsw
      !
      REAL(DP) np1,np2
      REAL(DP) d1rnp1,d1rnp2
      REAL(DP) t1,t2t9,t10,t10d1
      REAL(DP) f2,f3,f4,f5,f6,f7,f8,f9
      REAL(DP) f2d1,f3d1,f4d1,f5d1,f6d1,f8d1,f9d1
      REAL(DP) d1sf2,d1sf3,d1sf4,d1sf5,d1sf6,d1sf7,d1sf8,d1sf9
      REAL(DP) d1rf2,d1rf3,d1rf4,d1rf5,d1rf6,d1rf7,d1rf8,d1rf9
      REAL(DP) d1st1,d1rt1
      REAL(DP) d1st2t9,d1rt2t9
      REAL(DP) d1st10,d1rt10
      REAL(DP) d1sterm1,d1rterm1,term1d1
      REAL(DP) d1sterm2
      REAL(DP) d1sterm3,d1rterm3
      REAL(DP) d1sterm4,d1rterm4
      REAL(DP) d1sterm5,d1rterm5
      !
      REAL(DP) term1,term2,term3,term4,term5
      !
      ! REAL(DP) ei
      !
      REAL(DP) Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
      REAL(DP) Fifteen,Sixteen
      REAL(DP) r12,r64,r36,r81,r256,r384,r864,r1944,r4374
      REAL(DP) r20,r25,r27,r48,r120,r128,r144,r288,r324,r512,r729
      REAL(DP) r30,r32,r75,r243,r2187,r6561,r40,r105,r54,r135
      REAL(DP) r1215,r15309
      !
      SAVE Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
      DATA Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten &
        / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /
      SAVE Fifteen,Sixteen
      DATA Fifteen,Sixteen / 1.5D1, 1.6D1 /
      SAVE r36,r64,r81,r256,r384,r864,r1944,r4374
      DATA r36,r64,r81,r256,r384,r864,r1944,r4374 &
        / 3.6D1,6.4D1,8.1D1,2.56D2,3.84D2,8.64D2,1.944D3,4.374D3 /
      SAVE r27,r48,r120,r128,r144,r288,r324,r512,r729
      DATA r27,r48,r120,r128,r144,r288,r324,r512,r729 &
        / 2.7D1,4.8D1,1.2D2,1.28D2,1.44D2,2.88D2,3.24D2,5.12D2,7.29D2 /
      SAVE r20,r32,r243,r2187,r6561,r40
      DATA r20,r32,r243,r2187,r6561,r40 &
        / 2.0d1,3.2D1,2.43D2,2.187D3,6.561D3,4.0d1 /
      SAVE r12,r25,r30,r54,r75,r105,r135,r1215,r15309
      DATA r12,r25,r30,r54,r75,r105,r135,r1215,r15309 &
        / 1.2D1,2.5d1,3.0d1,5.4D1,7.5d1,1.05D2,1.35D2,1.215D3,1.5309D4 /
      !
      ! ... General constants
      !
      f12    = 0.5d0
      f13    = One/Three
      f14    = 0.25d0
      f18    = 0.125d0
      !
      f23    = Two * f13
      f43    = Two * f23
      !
      f32    = 1.5d0
      f72    = 3.5d0
      f34    = 0.75d0
      f94    = 2.25d0
      f98    = 1.125d0
      f1516  = Fifteen / Sixteen
      !
      pi     = ACOS(-One)
      pi2    = pi*pi
      pi_23  = pi2**f13
      srpi   = SQRT(pi)
      !
      Three_13 = Three**f13
      !
      ! Constants from fit
      !
      ea1 = -1.128223946706117d0
      ea2 = 1.452736265762971d0
      ea3 = -1.243162299390327d0
      ea4 = 0.971824836115601d0
      ea5 = -0.568861079687373d0
      ea6 = 0.246880514820192d0
      ea7 = -0.065032363850763d0
      ea8 = 0.008401793031216d0
      !
      eb1 = 1.455915450052607d0
      !
      ! Constants for PBE hole
      !
      A      =  1.0161144d0
      B      = -3.7170836d-1
      C      = -7.7215461d-2
      D      =  5.7786348d-1
      E      = -5.1955731d-2
      X      = - Eight/Nine
      !
      ! Constants for fit of H(s) (PBE)
      !
      Ha1    = 9.79681d-3
      Ha2    = 4.10834d-2
      Ha3    = 1.87440d-1
      Ha4    = 1.20824d-3
      Ha5    = 3.47188d-2
      !
      ! Constants for F(H) (PBE)
      !
      Fc1    = 6.4753871d0
      Fc2    = 4.7965830d-1
      !
      ! Constants for polynomial expansion for EG for small s
      !
      EGa1   = -2.628417880d-2
      EGa2   = -7.117647788d-2
      EGa3   =  8.534541323d-2
      !
      ! Constants for large x expansion of exp(x)*ei(-x)
      !
      expei1 = 4.03640D0
      expei2 = 1.15198D0
      expei3 = 5.03627D0
      expei4 = 4.19160D0
      !
      ! Cutoff criterion below which to use polynomial expansion
      !
      EGscut     = 8.0d-2
      wcutoff    = 1.4D1
      expfcutoff = 7.0D2
      !
      ! Calculate prelim variables
      !
      xkf    = (Three*pi2*rho) ** f13
      xkfrho = xkf * rho
      !
      A2 = A*A
      A3 = A2*A
      A4 = A3*A
      A12 = SQRT(A)
      A32 = A12*A
      A52 = A32*A
      A72 = A52*A
      !
      w      = omega / xkf
      w2    = w * w
      w3    = w2 * w
      w4    = w2 * w2
      w5    = w3 * w2
      w6    = w5 * w
      w7    = w6 * w
      w8    = w7 * w
      !
      d1rw  = -(One/(Three*rho))*w
      !
      X      = - Eight/Nine
      !
      s2     = s*s
      s3     = s2*s
      s4     = s2*s2
      s5     = s4*s
      s6     = s5*s
      !
      ! Calculate wPBE enhancement factor
      !
      Hnum    = Ha1*s2 + Ha2*s4
      Hden    = One + Ha3*s4 + Ha4*s5 + Ha5*s6
      !
      H       = Hnum/Hden
      !
      d1sHnum = Two*Ha1*s + Four*Ha2*s3
      d1sHden = Four*Ha3*s3 + Five*Ha4*s4 + Six*Ha5*s5
      !
      d1sH    = (Hden*d1sHnum - Hnum*d1sHden) / (Hden*Hden)
      !
      F      = Fc1*H + Fc2
      d1sF   = Fc1*d1sH
      !
      ! Change exponent of Gaussian if we're using the simple approx.
      !
      IF (w > wcutoff) eb1 = 2.0d0
      !
      ! Calculate helper variables (should be moved later on...)
      !
      Hsbw = s2*H + eb1*w2
      Hsbw2 = Hsbw*Hsbw
      Hsbw3 = Hsbw2*Hsbw
      Hsbw4 = Hsbw3*Hsbw
      Hsbw12 = SQRT(Hsbw)
      Hsbw32 = Hsbw12*Hsbw
      Hsbw52 = Hsbw32*Hsbw
      Hsbw72 = Hsbw52*Hsbw
      !
      d1sHsbw  = d1sH*s2 + Two*s*H
      d1rHsbw  = Two*eb1*d1rw*w
      !
      DHsbw = D + s2*H + eb1*w2
      DHsbw2 = DHsbw*DHsbw
      DHsbw3 = DHsbw2*DHsbw
      DHsbw4 = DHsbw3*DHsbw
      DHsbw5 = DHsbw4*DHsbw
      DHsbw12 = SQRT(DHsbw)
      DHsbw32 = DHsbw12*DHsbw
      DHsbw52 = DHsbw32*DHsbw
      DHsbw72 = DHsbw52*DHsbw
      DHsbw92 = DHsbw72*DHsbw
      !
      HsbwA94   = f94 * Hsbw / A
      HsbwA942  = HsbwA94*HsbwA94
      HsbwA943  = HsbwA942*HsbwA94
      HsbwA945  = HsbwA943*HsbwA942
      HsbwA9412 = SQRT(HsbwA94)
      !
      DHs    = D + s2*H
      DHs2   = DHs*DHs
      DHs3   = DHs2*DHs
      DHs4   = DHs3*DHs
      DHs72  = DHs3*SQRT(DHs)
      DHs92  = DHs72*DHs
      !
      d1sDHs = Two*s*H + s2*d1sH
      !
      DHsw   = DHs + w2
      DHsw2  = DHsw*DHsw
      DHsw52 = SQRT(DHsw)*DHsw2
      DHsw72 = DHsw52*DHsw
      !
      d1rDHsw = Two*d1rw*w
      !
      IF (s > EGscut) THEN
        !
        G_a    = srpi * (Fifteen*E + Six*C*(One+F*s2)*DHs + &
                         Four*B*(DHs2) + Eight*A*(DHs3))    &
                      * (One / (Sixteen * DHs72))           &
                       - f34*pi*SQRT(A) * EXP(f94*H*s2/A) * &
                         (One - ERF(f32*s*SQRT(H/A)))
        !
        d1sG_a = (One/r32)*srpi *                           &
                 ((r36*(Two*H + d1sH*s) / (A12*SQRT(H/A)))  &
                  + (One/DHs92) *                           &
                     (-Eight*A*d1sDHs*DHs3 - r105*d1sDHs*E  &
                      -r30*C*d1sDHs*DHs*(One+s2*F)          &
                      +r12*DHs2*(-B*d1sDHs + C*s*(d1sF*s + Two*F)))  &
                  - ((r54*EXP(f94*H*s2/A)*srpi*s*(Two*H+d1sH*s)*     &
                     ERFC(f32*SQRT(H/A)*s))                       &
                     / A12))
        !
        G_b    = (f1516 * srpi * s2) / DHs72
        !
        d1sG_b = (Fifteen*srpi*s*(Four*DHs - Seven*d1sDHs*s)) &
                 / (r32*DHs92)
        !
        EG     = - (f34*pi + G_a) / G_b
        !
        d1sEG  = (-Four*d1sG_a*G_b + d1sG_b*(Four*G_a + Three*pi)) &
                 / (Four*G_b*G_b)
        !
      ELSE
        !
        EG    = EGa1 + EGa2*s2 + EGa3*s4
        d1sEG = Two*EGa2*s + Four*EGa3*s3
        !
      ENDIF
      !
      ! Calculate the terms needed in any case
      !
      term2 =       (DHs2*B + DHs*C + Two*E + DHs*s2*C*F + Two*s2*EG) / &
                    (Two*DHs3)
      !
      d1sterm2 = (-Six*d1sDHs*(EG*s2 + E)                     &
                  + DHs2 * (-d1sDHs*B + s*C*(d1sF*s + Two*F)) &
                  + Two*DHs * (Two*EG*s - d1sDHs*C            &
                  + s2 * (d1sEG - d1sDHs*C*F)))               &
                 / (Two*DHs4)

      term3 = - w  * (Four*DHsw2*B + Six*DHsw*C + Fifteen*E &
                      + Six*DHsw*s2*C*F + Fifteen*s2*EG) /  &
                     (Eight*DHs*DHsw52)
      !
      d1sterm3 = w * (Two*d1sDHs*DHsw * (Four*DHsw2*B         &
                         + Six*DHsw*C + Fifteen*E             &
                         + Three*s2*(Five*EG + Two*DHsw*C*F)) &
                      + DHs * (r75*d1sDHs*(EG*s2 + E)         &
                         + Four*DHsw2*(d1sDHs*B               &
                              - Three*s*C*(d1sF*s + Two*F))   &
                         - Six*DHsw*(-Three*d1sDHs*C          &
                              + s*(Ten*EG + Five*d1sEG*s      &
                                  - Three*d1sDHs*s*C*F))))    &
                 / (Sixteen*DHs2*DHsw72)
      !
      d1rterm3 = (-Two*d1rw*DHsw * (Four*DHsw2*B              &
                         + Six*DHsw*C + Fifteen*E             &
                         + Three*s2*(Five*EG + Two*DHsw*C*F)) &
                      + w * d1rDHsw * (r75*(EG*s2 + E)        &
                         + Two*DHsw*(Two*DHsw*B + Nine*C      &
                                     + Nine*s2*C*F)))         &
                 / (Sixteen*DHs*DHsw72)

      term4 = - w3 * (DHsw*C + Five*E + DHsw*s2*C*F + Five*s2*EG) /  &
                     (Two*DHs2*DHsw52)
      !
      d1sterm4 = (w3 * (Four*d1sDHs*DHsw * (DHsw*C + Five*E   &
                             + s2 * (Five*EG + DHsw*C*F))     &
                        + DHs * (r25*d1sDHs*(EG*s2 + E)       &
                             - Two*DHsw2*s*C*(d1sF*s + Two*F) &
                             + DHsw * (Three*d1sDHs*C + s*(-r20*EG  &
                                   - Ten*d1sEG*s              &
                                   + Three*d1sDHs*s*C*F)))))  &
                 / (Four*DHs3*DHsw72)
      !
      d1rterm4 = (w2 * (-Six*d1rw*DHsw * (DHsw*C + Five*E   &
                             + s2 * (Five*EG + DHsw*C*F))   &
                        + w * d1rDHsw * (r25*(EG*s2 + E) +  &
                             Three*DHsw*C*(One + s2*F))))  &
                 / (Four*DHs2*DHsw72)
      !
      term5 = - w5 * (E + s2*EG) / &
                     (DHs3*DHsw52)
      !
      d1sterm5 = (w5 * (Six*d1sDHs*DHsw*(EG*s2 + E)               &
                        + DHs * (-Two*DHsw*s * (Two*EG + d1sEG*s) &
                             + Five*d1sDHs * (EG*s2 + E))))       &
                 / (Two*DHs4*DHsw72)
      !
      d1rterm5 = (w4 * Five*(EG*s2 + E) * (-Two*d1rw*DHsw   &
                                           + d1rDHsw * w))  &
                 / (Two*DHs3*DHsw72)
      !
      !
      IF ((s > 0.0d0).OR.(w > 0.0d0)) THEN
        !
        t10    = (f12)*A*LOG(Hsbw / DHsbw)
        t10d1  = f12*A*(One/Hsbw - One/DHsbw)
        d1st10 = d1sHsbw*t10d1
        d1rt10 = d1rHsbw*t10d1
        !
      ENDIF
      !
      ! Calculate exp(x)*f(x) depending on size of x
      !
      IF (HsbwA94 < expfcutoff) THEN
        !
        piexperf = pi*EXP(HsbwA94)*ERFC(HsbwA9412)
        ! expei    = Exp(HsbwA94)*Ei(-HsbwA94)
        CALL expint(1,HsbwA94,exint,in_err)
        expei    = EXP(HsbwA94)*(-exint)
      ELSE
        !
        ! print *,rho,s," LARGE HsbwA94"
        !
        piexperf = pi*(One/(srpi*HsbwA9412)          &
                   - One/(Two*SQRT(pi*HsbwA943))     &
                   + Three/(Four*SQRT(pi*HsbwA945)))
        !
        expei  = - (One/HsbwA94) *                         &
                   (HsbwA942 + expei1*HsbwA94 + expei2) /  &
                   (HsbwA942 + expei3*HsbwA94 + expei4)

      ENDIF
      !
      ! Calculate the derivatives (based on the orig. expression)
      ! --> Is this ok? ==> seems to be ok...
      !
      piexperfd1  = - (Three*srpi*SQRT(Hsbw/A))/(Two*Hsbw)  &
                    + (Nine*piexperf)/(Four*A)
      d1spiexperf = d1sHsbw*piexperfd1
      d1rpiexperf = d1rHsbw*piexperfd1

      expeid1  = f14*(Four/Hsbw + (Nine*expei)/A)
      d1sexpei = d1sHsbw*expeid1
      d1rexpei = d1rHsbw*expeid1
      !
      IF (w == Zero) THEN
        !
        ! Fall back to original expression for the PBE hole
        !
        t1 = -f12*A*expei
        d1st1 = -f12*A*d1sexpei
        d1rt1 = -f12*A*d1rexpei
        !
        ! write(*,*) s, t1, t10, d1st1,d1rt1,d1rt10
        !
        IF (s > 0.0D0) THEN
          !
          term1    = t1 + t10
          d1sterm1 = d1st1 + d1st10
          d1rterm1 = d1rt1 + d1rt10
          !
          Fx_wpbe = X * (term1 + term2)
          !
          d1sfx = X * (d1sterm1 + d1sterm2)
          d1rfx = X * d1rterm1
          !
        ELSE
          !
          Fx_wpbe = 1.0d0
          !
          ! TODO    This is checked to be true for term1
          !         How about the other terms???
          !
          d1sfx   = 0.0d0
          d1rfx   = 0.0d0
          !
        ENDIF
        !
        !
      ELSEIF (w > wcutoff) THEN
        !
        ! Use simple Gaussian approximation for large w
        !
        ! print *,rho,s," LARGE w"
        !
        term1 = -f12*A*(expei+LOG(DHsbw)-LOG(Hsbw))

        term1d1  = - A/(Two*DHsbw) - f98*expei
        d1sterm1 = d1sHsbw*term1d1
        d1rterm1 = d1rHsbw*term1d1

        Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5)

        d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3  &
                              + d1sterm4 + d1sterm5)

        d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)
        !
      ELSE
         !
         ! For everything else, use the full blown expression
         !
         ! First, we calculate the polynomials for the first term
         !
         np1    = -f32*ea1*A12*w + r27*ea3*w3/(Eight*A12)     &
                  - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52)
        !
        d1rnp1 = - f32*ea1*d1rw*A12 + (r81*ea3*d1rw*w2)/(Eight*A12) &
                 - (r1215*ea5*d1rw*w4)/(r32*A32)                    &
                 + (r15309*ea7*d1rw*w6)/(r128*A52)
        !
        np2 = -A + f94*ea2*w2 - r81*ea4*w4/(Sixteen*A)        &
              + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3)
        !
        !
        d1rnp2 =   f12*(Nine*ea2*d1rw*w)         &
                 - (r81*ea4*d1rw*w3)/(Four*A)    &
                 + (r2187*ea6*d1rw*w5)/(r32*A2)  &
                 - (r6561*ea8*d1rw*w7)/(r32*A3)
        !
        ! The first term is
        !
        t1    = f12*(np1*piexperf + np2*expei)
        d1st1 = f12*(d1spiexperf*np1 + d1sexpei*np2)
        d1rt1 = f12*(d1rnp2*expei + d1rpiexperf*np1 +  &
                     d1rexpei*np2 + d1rnp1*piexperf)
        !
        ! The factors for the main polynomoal in w and their derivatives
        !
        f2    = (f12)*ea1*srpi*A / DHsbw12
        f2d1  = - ea1*srpi*A / (Four*DHsbw32)
        d1sf2 = d1sHsbw*f2d1
        d1rf2 = d1rHsbw*f2d1
        !
        f3    = (f12)*ea2*A / DHsbw
        f3d1  = - ea2*A / (Two*DHsbw2)
        d1sf3 = d1sHsbw*f3d1
        d1rf3 = d1rHsbw*f3d1
        !
        f4    =  ea3*srpi*(-f98 / Hsbw12     &
                 + f14*A / DHsbw32)
        f4d1  = ea3*srpi*((Nine/(Sixteen*Hsbw32))-   &
                          (Three*A/(Eight*DHsbw52)))
        d1sf4 = d1sHsbw*f4d1
        d1rf4 = d1rHsbw*f4d1
        !
        f5    = ea4*(One/r128) * (-r144*(One/Hsbw)   &
                 + r64*(One/DHsbw2)*A)
        f5d1  = ea4*((f98/Hsbw2)-(A/DHsbw3))
        d1sf5 = d1sHsbw*f5d1
        d1rf5 = d1rHsbw*f5d1
        !
        f6    = ea5*(Three*srpi*(Three*DHsbw52*(Nine*Hsbw-Two*A) &
                 + Four*Hsbw32*A2))                              &
                 / (r32*DHsbw52*Hsbw32*A)
        f6d1  = ea5*srpi*((r27/(r32*Hsbw52))-        &
                    (r81/(r64*Hsbw32*A))-            &
                    ((Fifteen*A)/(Sixteen*DHsbw72)))
        d1sf6 = d1sHsbw*f6d1
        d1rf6 = d1rHsbw*f6d1
        !
        f7    = ea6*(((r32*A)/DHsbw3                 &
                 + (-r36 + (r81*s2*H)/A)/Hsbw2)) / r32
        d1sf7 = ea6*(Three*(r27*d1sH*DHsbw4*Hsbw*s2 +           &
                Eight*d1sHsbw*A*(Three*DHsbw4 - Four*Hsbw3*A) + &
                r54*DHsbw4*s*(Hsbw - d1sHsbw*s)*H))/            &
                (r32*DHsbw4*Hsbw3*A)
        d1rf7 = ea6*d1rHsbw*((f94/Hsbw3)-((Three*A)/DHsbw4)     &
                           -((r81*s2*H)/(Sixteen*Hsbw3*A)))
        !
        f8    = ea7*(-Three*srpi*(-r40*Hsbw52*A3                &
                 +Nine*DHsbw72*(r27*Hsbw2-Six*Hsbw*A+Four*A2))) &
                 / (r128 * DHsbw72*Hsbw52*A2)
        f8d1  = ea7*srpi*((r135/(r64*Hsbw72)) + (r729/(r256*Hsbw32*A2))  &
                         -(r243/(r128*Hsbw52*A))                         &
                         -((r105*A)/(r32*DHsbw92)))
        d1sf8 = d1sHsbw*f8d1
        d1rf8 = d1rHsbw*f8d1
        !
        f9    = (r324*ea6*eb1*DHsbw4*Hsbw*A                      &
                + ea8*(r384*Hsbw3*A3 + DHsbw4*(-r729*Hsbw2       &
                + r324*Hsbw*A - r288*A2))) / (r128*DHsbw4*Hsbw3*A2)
        f9d1  = -((r81*ea6*eb1)/(Sixteen*Hsbw3*A))               &
                + ea8*((r27/(Four*Hsbw4))+(r729/(r128*Hsbw2*A2)) &
                      -(r81/(Sixteen*Hsbw3*A))                   &
                      -((r12*A/DHsbw5)))
        d1sf9 = d1sHsbw*f9d1
        d1rf9 = d1rHsbw*f9d1
        !
        t2t9    = f2*w  + f3*w2 + f4*w3 + f5*w4 + f6*w5          &
                        + f7*w6 + f8*w7 + f9*w8
        d1st2t9 = d1sf2*w + d1sf3*w2 + d1sf4*w3 + d1sf5*w4       &
                          + d1sf6*w5 + d1sf7*w6 + d1sf8*w7       &
                          + d1sf9*w8
        d1rt2t9 = d1rw*f2 + d1rf2*w + Two*d1rw*f3*w   &
                  + d1rf3*w2 + Three*d1rw*f4*w2       &
                  + d1rf4*w3 + Four*d1rw*f5*w3        &
                  + d1rf5*w4 + Five*d1rw*f6*w4        &
                  + d1rf6*w5 + Six*d1rw*f7*w5         &
                  + d1rf7*w6 + Seven*d1rw*f8*w6       &
                  + d1rf8*w7 + Eight*d1rw*f9*w7 + d1rf9*w8
        !
        ! The final value of term1 for 0 < omega < wcutoff is:
        !
        term1 = t1 + t2t9 + t10
        !
        d1sterm1 = d1st1 + d1st2t9 + d1st10
        d1rterm1 = d1rt1 + d1rt2t9 + d1rt10
        !
        ! The final value for the enhancement factor and its
        ! derivatives is:
        !
        Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5)
        !
        d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3    &
                              + d1sterm4 + d1sterm5)
        !
        d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)
        !
      ENDIF

END SUBROUTINE wpbe_analy_erfc_approx_grad
!
!------------------------------------------------------------------
SUBROUTINE EXPINT(n, x, exin, in_err)
!-----------------------------------------------------------------------
!! Evaluates the exponential integral \(E_n(x)\). 
!! Inspired by Numerical Recipes.
! Parameters: maxit is the maximum allowed number of iterations,
! eps is the desired relative error, not smaller than the machine precision,
! big is a number near the largest representable floating-point number,
!
      USE kind_l,               ONLY: DP
      IMPLICIT NONE
      !$acc routine seq
      INTEGER, INTENT(IN) :: n
      REAL(DP), INTENT(IN) :: x
      REAL(DP), INTENT(OUT) :: exin
      INTEGER :: in_err
      INTEGER, parameter :: maxit=200
      REAL(DP), parameter :: eps=1E-12, big=huge(x)*eps
      REAL(DP), parameter :: euler = 0.577215664901532860606512d0
!     EPS=1E-9, FPMIN=1E-30

      INTEGER :: i, nm1, k
      REAL(DP) :: a,b,c,d,del,fact,h,iarsum

      IF (.NOT. ((n >= 0).AND.(x >= 0.0).AND.((x > 0.0).OR.(n > 1)))) THEN
         in_err = 1   ! expint: bad arguments
         RETURN
      END IF
      
      IF (n == 0) THEN
         exin= exp(-x)/x
         RETURN
      END IF
      nm1 = n-1
      IF (x == 0.0d0) THEN
         exin = 1.0d0/nm1
      ELSE IF (x > 1.0d0) THEN
         b = x+n
         c = big
         d = 1.0d0/b
         h = d
         DO i=1,maxit
            a = -i*(nm1+i)
            b = b+2.0d0
            d = 1.0d0/(a*d+b)
            c = b+a/c
            del = c*d
            h = h*del
            IF (ABS(del-1.0d0) <= EPS) EXIT
         END DO
         IF (i > maxit) THEN
           in_err = 2   ! expint: continued fraction failed
           RETURN
         ENDIF
         exin = h*EXP(-x)
      ELSE
         IF (nm1 /= 0) THEN
            exin = 1.0d0/nm1
         ELSE
            exin = -LOG(x)-euler
         END IF
         fact = 1.0d0
         do i=1,maxit
            fact = -fact*x/i
            IF (i /= nm1) THEN
               del = -fact/(i-nm1)
            ELSE

               iarsum = 0.0d0
               do k=1,nm1
                  iarsum = iarsum + 1.0d0/k
               end do

               del = fact*(-LOG(x)-euler+iarsum)
!               del = fact*(-LOG(x)-euler+sum(1.0d0/arth(1,1,nm1)))
            END IF
            exin = exin + del
            IF (ABS(del) < ABS(exin)*eps) EXIT
         END DO
         IF (i > maxit) THEN
           in_err = 2   ! expint: series failed
           RETURN
         ENDIF
      END IF
END SUBROUTINE EXPINT
!
END MODULE

