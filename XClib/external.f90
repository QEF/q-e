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




