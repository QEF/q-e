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

!-------------------------------------------------------------------------
SUBROUTINE tpsscxc_ext( rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! TPSS metaGGA corrections for exchange and correlation - Hartree a.u.
  !
  !! Definition:  E_x = \int E_x(rho,grho) dr
  !
  USE kind_l,      ONLY : DP
  USE metagga,     ONLY : metax, metac
  !
  IMPLICIT NONE
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
  CALL metax( rho, grho, tau, sx, v1x, v2x, v3x )                   !<GPU:metax=>metax_d>
  ! correlation
  CALL metac( rho, grho, tau, sc, v1c, v2c, v3c )                   !<GPU:metac=>metac_d>
  !
  RETURN
  !
END SUBROUTINE tpsscxc_ext
