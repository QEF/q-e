!-----------------------------------------------------------------------
SUBROUTINE slater_ext( rs, ex, vx )                    !<GPU:DEVICE>
  !---------------------------------------------------------------------
  !! Slater exchange with alpha=2/3
  !
  USE kind_l,  ONLY: DP
  !
  IMPLICIT NONE
  !!
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(OUT) :: ex
  !! Exchange energy (per unit volume)
  REAL(DP), INTENT(OUT) :: vx
  !! Exchange potential
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER   :: f = -0.687247939924714_DP, alpha = 2.0_DP/3.0_DP
  !                        f = -9/8*(3/2pi)^(2/3)
  ex = f * alpha / rs
  vx = 4._DP / 3._DP * f * alpha / rs
  !
  RETURN
  !
END SUBROUTINE slater_ext
!
!-----------------------------------------------------------------------
SUBROUTINE slater_spin_ext( rho, zeta, ex, vx_up, vx_dw )                 !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! Slater exchange with alpha=2/3, spin-polarized case.
  !
  USE kind_l, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho
  !! total charge density
  REAL(DP), INTENT(IN) :: zeta
  !! zeta = (rho_up - rho_dw) / rho_tot
  REAL(DP), INTENT(OUT) :: ex
  !! exchange energy
  REAL(DP), INTENT(OUT) :: vx_up, vx_dw
  !! exchange potential (up, down)
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: f = -1.10783814957303361d0, alpha = 2.0d0/3.0d0
  !                      f = -9/8*(3/pi)^(1/3)
  REAL(DP), PARAMETER :: third = 1.d0/3.d0, p43 = 4.d0/3.d0
  REAL(DP) :: exup, exdw, rho13
  !
  !
  rho13 = ( (1.d0 + zeta)*rho )**third
  exup = f * alpha * rho13
  vx_up = p43 * f * alpha * rho13
  !
  rho13 = ( (1.d0 - zeta)*rho )**third
  exdw = f * alpha * rho13
  vx_dw = p43 * f * alpha * rho13
  !
  ex = 0.5d0 * ( (1.d0 + zeta)*exup + (1.d0 - zeta)*exdw)
  !
  RETURN
  !
END SUBROUTINE slater_spin_ext
!
!
!-------------------------------------------------------------------------
SUBROUTINE pz_ext( rs, iflag, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! LDA parametrization from Monte Carlo DATA:
  !
  !! * iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981);
  !! * iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994).
  !
  USE kind_l,  ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iflag                        !<GPU:VALUE>
  !! see routine comments
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: vc
  !! correlation potential
  !
  ! ... local variables
  !
  REAL(DP) :: a(2), b(2), c(2), d(2), gc(2), b1(2), b2(2)
  REAL(DP) :: lnrs, rs12, ox, dox
  !
  DATA a / 0.0311d0, 0.031091d0 /, b / -0.048d0,  -0.046644d0 /, &
       c / 0.0020d0, 0.00419d0  /, d / -0.0116d0, -0.00983d0  /
  !
  DATA gc / -0.1423d0, -0.103756d0 /, b1 / 1.0529d0, 0.56371d0 /, &
       b2 /  0.3334d0,  0.27358d0  /
  !
  IF ( rs < 1.0d0 ) THEN
     !
     ! high density formula
     lnrs = LOG(rs)
     !
     ec = a(iflag)*lnrs + b(iflag) + c(iflag)*rs*lnrs + d(iflag)*rs
     !
     vc = a(iflag)*lnrs + ( b(iflag) - a(iflag)/3.d0 ) + 2.d0/3.d0 * &
          c(iflag)*rs*lnrs + ( 2.d0*d(iflag) - c(iflag) )/3.d0*rs
  ELSE
     !
     ! interpolation formula
     rs12 = SQRT(rs)
     !
     ox  = 1.d0 + b1(iflag)*rs12 + b2(iflag)*rs
     dox = 1.d0 + 7.d0/6.d0*b1(iflag)*rs12 + 4.d0/3.d0*b2(iflag)*rs
     !
     ec = gc(iflag)/ox
     vc = ec*dox/ox
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE pz_ext
!

!-----------------------------------------------------------------------
SUBROUTINE pz_polarized_ext( rs, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! J.P. Perdew and A. Zunger, PRB 23, 5048 (1981).
  !! spin-polarized energy and potential.
  !
  USE kind_l, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: vc
  !! correlation potential
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: a=0.01555d0, b=-0.0269d0, c=0.0007d0, &
                         d=-0.0048d0, gc=-0.0843d0, b1=1.3981d0, b2=0.2611d0
  REAL(DP) :: lnrs, rs12, ox, dox
  REAL(DP), PARAMETER :: xcprefact=0.022575584d0, pi34=0.6203504908994d0
  ! REAL(DP) :: betha, etha, csi, prefact
  !
  !
  IF ( rs < 1.0d0 ) THEN
     ! high density formula
     lnrs = LOG(rs)
     ec = a * lnrs + b + c * rs * lnrs + d * rs
     vc = a * lnrs + (b - a / 3.d0) + 2.d0 / 3.d0 * c * rs * lnrs + &
          (2.d0 * d-c) / 3.d0 * rs
  ELSE
     ! interpolation formula
     rs12 = SQRT(rs)
     ox = 1.d0 + b1 * rs12 + b2 * rs
     dox = 1.d0 + 7.d0 / 6.d0 * b1 * rs12 + 4.d0 / 3.d0 * b2 * rs
     ec = gc / ox
     vc = ec * dox / ox
  ENDIF
  !
  !  IF ( lxc_rel ) THEN
  !     betha = prefact * pi34 / rs
  !     etha = DSQRT( 1 + betha**2 )
  !     csi = betha + etha
  !     prefact = 1.0D0 - (3.0D0/2.0D0) * ( (betha*etha - LOG(csi))/betha**2 )**2
  !     ec = ec * prefact
  !     vc = vc * prefact
  !  ENDIF
  !
  RETURN
  !
END SUBROUTINE pz_polarized_ext
!
!
!-----------------------------------------------------------------------
SUBROUTINE pw_ext( rs, iflag, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! * iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !! * iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  INTEGER, INTENT(IN)  :: iflag                 !<GPU:VALUE>
  !! see routine comments
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: vc
  !! correlation potential
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: a=0.031091d0, b1=7.5957d0, b2=3.5876d0, c0=a, &
                         c1=0.046644d0, c2=0.00664d0, c3=0.01043d0, d0=0.4335d0, &
                         d1=1.4408d0
  REAL(DP) :: lnrs, rs12, rs32, rs2, om, dom, olog
  REAL(DP) :: a1(2), b3(2), b4(2)
  DATA a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /
  !
  ! high- and low-density formulae implemented but not used in PW case
  ! (reason: inconsistencies in PBE/PW91 functionals).
  !
  IF ( rs < 1d0 .AND. iflag == 2 ) THEN
     !
     ! high density formula
     lnrs = LOG(rs)
     ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs
     vc = c0 * lnrs - (c1 + c0 / 3.d0) + 2.d0 / 3.d0 * c2 * rs * &
               lnrs - (2.d0 * c3 + c2) / 3.d0 * rs
     !
  ELSEIF ( rs > 100.d0 .AND. iflag == 2 ) THEN
     !
     ! low density formula
     ec = - d0 / rs + d1 / rs**1.5d0
     vc = - 4.d0 / 3.d0 * d0 / rs + 1.5d0 * d1 / rs**1.5d0
     !
  ELSE
     !
     ! interpolation formula
     rs12 = SQRT(rs)
     rs32 = rs * rs12
     rs2  = rs**2
     om   = 2.d0*a*( b1*rs12 + b2*rs + b3(iflag) * rs32 + b4(iflag)*rs2 )
     dom  = 2.d0*a*( 0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3(iflag) * &
            rs32 + 2.d0 * b4(iflag) * rs2 )
     olog = LOG( 1.d0 + 1.0d0 / om )
     !
     ec = - 2.d0 * a * (1.d0 + a1(iflag) * rs) * olog
     vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1(iflag) * rs) &
              * olog - 2.d0 / 3.d0 * a * (1.d0 + a1(iflag) * rs) * dom / &
              (om * (om + 1.d0) )
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE pw_ext

!
!-----------------------------------------------------------------------
SUBROUTINE pw_spin_ext( rs, zeta, ec, vc_up, vc_dw )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! J.P. Perdew and Y. Wang, PRB 45, 13244 (1992).
  !
  USE kind_l, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(IN) :: zeta
  !! zeta = (rho_up - rho_dw)/rho_tot
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: vc_up, vc_dw
  !! correlation potential (up, down)
  !
  ! ... local variables
  !
  REAL(DP) :: rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz
  REAL(DP) :: om, dom, olog, epwc, vpwc
  REAL(DP) :: omp, domp, ologp, epwcp, vpwcp
  REAL(DP) :: oma, doma, ologa, alpha, vpwca
  !
  ! xc parameters, unpolarised
  REAL(DP), PARAMETER :: a = 0.031091d0, a1 = 0.21370d0, b1 = 7.5957d0, b2 = &
           3.5876d0, b3 = 1.6382d0, b4 = 0.49294d0, c0 = a, c1 = 0.046644d0, &
           c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, d1 = 1.4408d0
  ! xc parameters, polarised
  REAL(DP), PARAMETER :: ap = 0.015545d0, a1p = 0.20548d0, b1p = 14.1189d0, b2p &
               = 6.1977d0, b3p = 3.3662d0, b4p = 0.62517d0, c0p = ap, c1p =     &
              0.025599d0, c2p = 0.00319d0, c3p = 0.00384d0, d0p = 0.3287d0, d1p &
              = 1.7697d0
  ! xc PARAMETERs, antiferro
  REAL(DP), PARAMETER :: aa = 0.016887d0, a1a = 0.11125d0, b1a = 10.357d0, b2a = &
               3.6231d0, b3a = 0.88026d0, b4a = 0.49671d0, c0a = aa, c1a =       &
               0.035475d0, c2a = 0.00188d0, c3a = 0.00521d0, d0a = 0.2240d0, d1a &
               = 0.3969d0
  REAL(DP), PARAMETER :: fz0 = 1.709921d0
  !
  !     if (rs < 0.5d0) then
  ! high density formula (not implemented)
  !
  !     elseif (rs > 100.d0) then
  ! low density formula (not implemented)
  !
  !     else
  ! interpolation formula
  !
  zeta2 = zeta * zeta
  zeta3 = zeta2 * zeta
  zeta4 = zeta3 * zeta
  rs12 = SQRT(rs)
  rs32 = rs * rs12
  rs2 = rs**2
  !
  ! unpolarised
  om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
  dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 * rs32 &
       + 2.d0 * b4 * rs2)
  olog = LOG(1.d0 + 1.0d0 / om)
  epwc = - 2.d0 * a * (1.d0 + a1 * rs) * olog
  vpwc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 * rs) * olog - 2.d0 / &
         3.d0 * a * (1.d0 + a1 * rs) * dom / (om * (om + 1.d0) )
  !
  ! polarized
  omp  = 2.d0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
  domp = 2.d0 * ap * (0.5d0 * b1p * rs12 + b2p * rs + 1.5d0 * b3p * &
         rs32 + 2.d0 * b4p * rs2)
  ologp = LOG(1.d0 + 1.0d0 / omp)
  epwcp = - 2.d0 * ap * (1.d0 + a1p * rs) * ologp
  vpwcp = - 2.d0 * ap * (1.d0 + 2.d0 / 3.d0 * a1p * rs) * ologp - &
          2.d0 / 3.d0 * ap * (1.d0 + a1p * rs) * domp / (omp * (omp + 1.d0))
  !
  ! antiferro
  oma = 2.d0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
  doma = 2.d0 * aa * ( 0.5d0 * b1a * rs12 + b2a * rs + 1.5d0 * b3a * &
         rs32 + 2.d0 * b4a * rs2 )
  ologa = LOG( 1.d0 + 1.0d0/oma )
  alpha = 2.d0 * aa * (1.d0 + a1a*rs) * ologa
  vpwca = + 2.d0 * aa * (1.d0 + 2.d0/3.d0 * a1a * rs) * ologa + &
          2.d0 / 3.d0 * aa * (1.d0 + a1a*rs) * doma / (oma * (oma + 1.d0))
  !
  !
  fz = ( (1.d0 + zeta)**(4.d0 / 3.d0) + (1.d0 - zeta)**(4.d0 / &
          3.d0) - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2.d0)
  dfz = ( (1.d0 + zeta)**(1.d0 / 3.d0) - (1.d0 - zeta)**(1.d0 / &
          3.d0) ) * 4.d0 / (3.d0 * (2.d0** (4.d0 / 3.d0) - 2.d0) )
  !
  !
  ec = epwc + alpha * fz * (1.d0 - zeta4) / fz0 + (epwcp - epwc) &
              * fz * zeta4
  !
  vc_up = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
                 * fz * zeta4 + (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
                 zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
                 * (1.d0 - zeta)
  !
  vc_dw = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
                 * fz * zeta4 - (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
                 zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
                 * (1.d0 + zeta)
  !
  RETURN
  !
END SUBROUTINE pw_spin_ext

!-----------------------------------------------------------------------
SUBROUTINE lyp_ext( rs, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! C. Lee, W. Yang, and R.G. Parr, PRB 37, 785 (1988).
  !! LDA part only.
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: vc
  !! correlation potential
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: a=0.04918d0, b=0.132d0*2.87123400018819108d0
  REAL(DP), PARAMETER :: pi43=1.61199195401647d0
  !                      pi43=(4pi/3)^(1/3)
  REAL(DP), PARAMETER :: c=0.2533d0*pi43, d=0.349d0*pi43
  REAL(DP) :: ecrs, ox
  !
  !
  ecrs = b*EXP( -c*rs )
  ox = 1.d0 / (1.d0 + d*rs)
  !
  ec = - a*ox*(1.d0 + ecrs)
  vc = ec - rs/3.d0*a*ox*( d*ox + ecrs*(d*ox + c) )
  !
  RETURN
  !
END SUBROUTINE lyp_ext


!-----------------------------------------------------------------------------
SUBROUTINE lsd_lyp_ext( rho, zeta, elyp, vlyp_up, vlyp_dw )                    !<GPU:DEVICE>
  !==--------------------------------------------------------------==
  !==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
  !==  THIS IS ONLY THE LDA PART                                   ==
  !==--------------------------------------------------------------==
  !
  USE kind_l,       ONLY: DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho
  !! total charge density
  REAL(DP), INTENT(IN) :: zeta
  !! zeta = (rho_up - rho_dw)/rho_tot
  REAL(DP), INTENT(OUT) :: elyp
  !! correlation energy
  REAL(DP), INTENT(OUT) :: vlyp_up, vlyp_dw
  !! correlation potential (up, down)
  !
  ! ... local variables
  !
  REAL(DP) :: ra,rb,rm3,dr,e1,or,dor,e2,de1a,de1b,de2a,de2b
  REAL(DP), PARAMETER :: small=1.D-24, a=0.04918D0, b=0.132D0, &
                         c=0.2533D0,   d=0.349D0,  cf=2.87123400018819108D0
  !==--------------------------------------------------------------==
  ra = rho*0.5D0*(1.D0+zeta)
  ra = MAX(ra,small)
  rb = rho*0.5D0*(1.D0-zeta)
  rb = MAX(rb,small)
  !
  rm3= rho**(-1.D0/3.D0)
  dr = (1.D0+d*rm3)
  !
  e1 = 4.D0*a*ra*rb/rho/dr
  or = EXP(-c*rm3)/dr*rm3**11.D0
  dor= -1.D0/3.D0*rm3**4*or*(11.D0/rm3-c-d/dr)
  e2 = 2.D0**(11.D0/3.D0)*cf*a*b*or*ra*rb*(ra**(8.d0/3.d0)+ rb**(8.d0/3.d0))
  !
  elyp = (-e1-e2)/rho
  !
  de1a = -e1*(1.D0/3.D0*d*rm3**4/dr+1./ra-1./rho)
  de1b = -e1*(1.D0/3.D0*d*rm3**4/dr+1./rb-1./rho)
  de2a = -2.D0**(11.D0/3.D0)*cf*a*b*(dor*ra*rb*(ra**(8.d0/3.d0)+ &
          rb**(8.d0/3.d0))+or*rb*(11.d0/3.d0*ra**(8.d0/3.d0)+ &
          rb**(8.d0/3.d0)))
  de2b = -2.D0**(11.D0/3.D0)*cf*a*b*(dor*ra*rb*(ra**(8.d0/3.d0)+ &
          rb**(8.d0/3.d0))+or*ra*(11.d0/3.d0*rb**(8.d0/3.d0)+ &
          ra**(8.d0/3.d0)))
  !
  vlyp_up = de1a + de2a
  vlyp_dw = de1b + de2b
  !==--------------------------------------------------------------==
  !
  RETURN
  !
END SUBROUTINE lsd_lyp_ext

!---------------------------------------------------------------
SUBROUTINE pbec_ext( rho, grho, iflag, sc, v1c, v2c )                    !<GPU:DEVICE>
  !---------------------------------------------------------------
  !! PBE correlation (without LDA part)
  !
  !! * iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  !! * iflag=2: J.P.Perdew et al., PRL 100, 136406 (2008).
  !! * iflag=3: L. Chiodo et al, PRL 108, 126402 (2012)  (PBEQ2D)
  !
  USE kind_l,    ONLY: DP
  USE corr_lda_l, ONLY: pw_l
  USE corr_gga_l, ONLY: cpbe2d_l
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: iflag              !<GPU:VALUE>
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
  CALL pw_l( rs, 1, ec, vc )                      !<GPU:pw_l=>pw_l_d>
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
     CALL cpbe2d_l( rho, grho, sc2D, v1c2D, v2c2D )       !<GPU:cpbe2d=>cpbe2d_d>
     sc  = sc  + sc2D
     v1c = v1c + v1c2D
     v2c = v2c + v2c2D
  ENDIF
  !
  RETURN
  !
END SUBROUTINE pbec_ext
!
!
!-------------------------------------------------------------------
SUBROUTINE pbec_spin_ext( rho, zeta, grho, iflag, sc, v1c_up, v1c_dw, v2c )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------
  !! PBE correlation (without LDA part) - spin-polarized.
  !
  !! * iflag = 1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996);
  !! * iflag = 2: J.P.Perdew et al., PRL 100, 136406 (2008).
  !
  USE corr_lda_l, ONLY: pw_spin_l
  USE kind_l,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iflag        !<GPU:VALUE>
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
  CALL pw_spin_l( rs, zeta, ec, vc_up, vc_dn )                                 !<GPU:pw_spin=>pw_spin_d>
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
END SUBROUTINE pbec_spin_ext
!
!
!------------------------------------------------------------------------
SUBROUTINE lsd_glyp_ext( rho_in_up, rho_in_dw, grho_up, grho_dw, grho_ud, sc, v1c_up, v1c_dw, v2c_up, v2c_dw, v2c_ud )                     !<GPU:DEVICE>
  !----------------------------------------------------------------------
  !! Lee, Yang, Parr: gradient correction part.
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
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
END SUBROUTINE lsd_glyp_ext
