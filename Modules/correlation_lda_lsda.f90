!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE corr_lda !<GPU:corr_lda=>corr_lda_gpu>

USE ldaxc_interfaces, ONLY: pz

CONTAINS

!
!
! ... LSDA
!
!-----------------------------------------------------------------------
SUBROUTINE pz_polarized( rs, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! J.P. Perdew and A. Zunger, PRB 23, 5048 (1981).
  !! spin-polarized energy and potential.
  !
  USE kinds, ONLY : DP
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
END SUBROUTINE pz_polarized
!
!
!-----------------------------------------------------------------------
SUBROUTINE pz_spin( rs, zeta, ec, vc_up, vc_dw )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! Perdew and Zunger, PRB 23, 5048 (1981). Spin polarized case.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(IN) :: zeta
  !! zeta = (rho_up - rho_dw) / rho_tot
  REAL(DP), INTENT(OUT) :: ec
  !! correlation energy
  REAL(DP), INTENT(OUT) :: vc_up, vc_dw
  !! correlation potential (up, down)
  !
  ! ... local variables
  !
  REAL(DP) :: ecu, vcu, ecp, vcp
  REAL(DP) :: fz, dfz
  REAL(DP), PARAMETER :: p43=4.0d0/3.d0, third=1.d0/3.d0
  !
  ! unpolarized part (Perdew-Zunger formula)
  CALL pz( rs, 1, ecu, vcu )                                        !<GPU:pz=>pz_d>
  !
  ! polarization contribution
  CALL pz_polarized( rs, ecp, vcp )                                 !<GPU:pz_polarized=>pz_polarized_d>
  !
  fz = ( (1.0d0 + zeta)**p43 + (1.d0 - zeta)**p43 - 2.d0) / &
         (2.d0**p43 - 2.d0)
  dfz = p43 * ( (1.0d0 + zeta)**third- (1.d0 - zeta)**third) &
          / (2.d0**p43 - 2.d0)
  !
  ec = ecu + fz * (ecp - ecu)
  vc_up = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * ( 1.d0 - zeta)
  vc_dw = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (-1.d0 - zeta)
  !
  RETURN
  !
END SUBROUTINE pz_spin
!
!-------------------------------------------------------------------------------
SUBROUTINE vwn_spin( rs, zeta, ec, vc_up, vc_dw )                    !<GPU:DEVICE>
   !------------------------------------------------------------------------------
   !! S.H. Vosko, L. Wilk, and M. Nusair.  Spin polarized case.
   !
   USE kinds,     ONLY: DP
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN) :: rs
   !! Wigner-Seitz radius
   REAL(DP), INTENT(IN) :: zeta
   !! zeta = (rho_up - rho_dw) / rho_tot
   REAL(DP), INTENT(OUT) :: ec
   !! correlation energy
   REAL(DP), INTENT(OUT) :: vc_up, vc_dw
   !! correlation potential (up, down)
   !
   ! ... local variables
   !
   REAL(DP) :: zeta3, zeta4, trup, trdw, trup13, trdw13, fz, dfz, fzz4
   REAL(DP) :: SQRTrs, ecP, ecF, ac, De, vcP, vcF, dac, dec1, dec2
   REAL(DP) :: cfz, cfz1, cfz2, iddfz0
   !
   ! parameters:   e_c/para,    e_c/ferro,     alpha_c
   REAL(DP), PARAMETER :: &
      A(3)  = (/ 0.0310907_DP, 0.01554535_DP, -0.01688686394039_DP /), &
      x0(3) = (/ -0.10498_DP, -0.32500_DP, -0.0047584_DP /), &
      b(3)  = (/3.72744_DP, 7.06042_DP, 1.13107_DP /), &
      c(3)  = (/ 12.9352_DP, 18.0578_DP, 13.0045_DP /),&
      Q(3)  = (/ 6.15199081975908_DP, 4.73092690956011_DP, 7.12310891781812_DP /), &
      tbQ(3) = (/ 1.21178334272806_DP, 2.98479352354082_DP, 0.31757762321188_DP /), &
      fx0(3) = (/ 12.5549141492_DP, 15.8687885_DP, 12.99914055888256_DP /), &
      bx0fx0(3) = (/ -0.03116760867894_DP, -0.14460061018521_DP, -0.00041403379428_DP /)
   !
   ! N.B.: A is expressed in Hartree
   ! Q = SQRT(4*c - b^2)
   ! tbQ = 2*b/Q
   ! fx0 = X(x_0) = x_0^2 + b*x_0 + c
   ! bx0fx0 = b*x_0/X(x_0)
   !
   !
   ! coefficients for f(z), df/dz, ddf/ddz(0)
   cfz = 2.0_DP**(4.0_DP/3.0_DP) - 2.0_DP
   cfz1 = 1.0_DP / cfz
   cfz2 = 4.0_DP/3.0_DP * cfz1
   iddfz0 = 9.0_DP / 8.0_DP *cfz
   !
   SQRTrs = SQRT(rs)
   zeta3 = zeta**3
   zeta4 = zeta3*zeta
   trup = 1.0_DP + zeta
   trdw = 1.0_DP - zeta
   trup13 = trup**(1.0_DP/3.0_DP)
   trdw13 = trdw**(1.0_DP/3.0_DP)
   fz = cfz1 * (trup13*trup + trdw13*trdw - 2.0_DP) ! f(zeta)
   dfz = cfz2 * (trup13 - trdw13)                   ! df / dzeta
   !
   CALL padefit_ParSet1( sqrtrs, 1, ecP, vcP )              ! ecF = e_c Paramagnetic     !<GPU:padefit_ParSet1=>padefit_ParSet1_d>
   CALL padefit_ParSet1( sqrtrs, 2, ecF, vcF )              ! ecP = e_c Ferromagnetic    !<GPU:padefit_ParSet1=>padefit_ParSet1_d>
   CALL padefit_ParSet1( sqrtrs, 3, ac,  dac )              ! ac = "spin stiffness"      !<GPU:padefit_ParSet1=>padefit_ParSet1_d>
   !
   ac = ac * iddfz0
   dac = dac * iddfz0
   De = ecF - ecP - ac                              ! e_c[F]-e_c[P]-alpha_c/(ddf/ddz(z=0))
   fzz4 = fz * zeta4
   ec = ecP + ac * fz  + De * fzz4
   !
   dec1 = vcP + dac*fz + (vcF - vcP - dac) * fzz4   ! e_c - (r_s/3)*(de_c/dr_s)
   dec2 = ac*dfz + De*(4.0_DP*zeta3*fz + zeta4*dfz) ! de_c/dzeta
   !
   vc_up = dec1 + (1.0_DP - zeta)*dec2              ! v_c[s] = e_c - (r_s/3)*(de_c/dr_s)
   vc_dw = dec1 - (1.0_DP + zeta)*dec2              !          + [sign(s)-zeta]*(de_c/dzeta)
   !
END SUBROUTINE
!
!
!----
SUBROUTINE padefit_ParSet1( x, i, fit, dfit )                           !<GPU:DEVICE>
   !----
   !! It implements formula [4.4] in:
   !! S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
   !
   USE kinds, ONLY: DP
   !
   IMPLICIT NONE
   !
   REAL(DP) :: x
   !! input: x is SQRT(r_s)
   INTEGER :: i
   !! input: i is the index of the fit
   REAL(DP) :: fit
   !! output: Pade fit calculated in x [eq. 4.4]
   REAL(DP) :: dfit
   !! output: dfit/drho = fit - (rs/3)*dfit/drs = ec - (x/6)*dfit/dx
   !
   ! ... local variables
   !
   REAL(DP), PARAMETER :: &
   A(3)  = (/ 0.0310907d0, 0.01554535d0, -0.01688686394039d0 /), &
   x0(3) = (/ -0.10498d0, -0.32500d0, -0.0047584d0 /), &
   b(3)  = (/3.72744d0, 7.06042d0, 1.13107d0 /), &
   c(3)  = (/ 12.9352d0, 18.0578d0, 13.0045d0 /),&
   Q(3)  = (/ 6.15199081975908d0, 4.73092690956011d0, 7.12310891781812d0 /), &
   tbQ(3) = (/ 1.21178334272806d0, 2.98479352354082d0, 0.31757762321188d0 /), &
   fx0(3) = (/ 12.5549141492d0, 15.8687885d0, 12.99914055888256d0 /), &
   bx0fx0(3) = (/ -0.03116760867894d0, -0.14460061018521d0, -0.00041403379428d0 /)

   REAL(DP) :: sqx, xx0, Qtxb, atg, fx
   REAL(DP) :: txb, txbfx, itxbQ
   !
   sqx = x * x                          ! x^2 = r_s
   xx0 = x - x0(i)                      ! x - x_0
   Qtxb = Q(i) / (2.0_DP*x + b(i))      ! Q / (2x+b)
   atg = ATAN(Qtxb)                     ! tan^-1(Q/(2x+b))
   fx = sqx + b(i)*x + c(i)             ! X(x) = x^2 + b*x + c
   !
   fit = A(i) * (  LOG(sqx/fx) + tbQ(i)*atg - &
       bx0fx0(i) * ( LOG(xx0*xx0/fx) + (tbQ(i) + 4.0_DP*x0(i)/Q(i)) * atg )  )
   !
   txb = 2.0_DP*x + b(i)
   txbfx = txb / fx
   itxbQ = 1.0_DP / (txb*txb + Q(i)*Q(i))
   !
   dfit = fit - A(i) / 3.0_DP + A(i)*x/6.0_DP * (  txbfx + 4.0_DP*b(i)*itxbQ + &
           bx0fx0(i) * ( 2.0_DP/xx0 - txbfx - 4.0_DP*(b(i)+2.0_DP*x0(i))*itxbQ )  )
   !
END SUBROUTINE
SUBROUTINE padefit_ParSet2( x, i, fit, dfit )                           !<GPU:DEVICE>
   !----
   !! It implements formula [4.4] in:
   !! S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
   !
   USE kinds, ONLY: DP
   !
   IMPLICIT NONE
   !
   REAL(DP) :: x
   !! input: x is SQRT(r_s)
   INTEGER :: i
   !! input: i is the index of the fit
   REAL(DP) :: fit
   !! output: Pade fit calculated in x [eq. 4.4]
   REAL(DP) :: dfit
   !! output: dfit/drho = fit - (rs/3)*dfit/drs = ec - (x/6)*dfit/dx
   !
   ! ... local variables
   !
   REAL(DP), PARAMETER :: &
      A(3)  = (/  0.0310907_DP, 0.01554535_DP, -0.01688686394039_DP /), &
      x0(3) = (/ -0.409286_DP, -0.743294_DP,   -0.228344_DP /), &
      b(3)  = (/ 13.0720_DP,   20.1231_DP,      1.06835_DP  /), &
      c(3)  = (/ 42.7198_DP,  101.578_DP,      11.4813_DP /),&
      Q(3)  = (/  0.044899888641577_DP,      1.171685277708971_DP,  6.692072046645942_DP /), &
      tbQ(3) = (/ 582.273159042780890_DP,   34.348984975465861_DP,  0.319288254087299_DP /), &
      fx0(3) = (/ 37.537128437796000_DP,    87.173106479036008_DP, 11.289489669936000_DP /), &
      bx0fx0(3) = (/ -0.142530524167984_DP, -0.171582499414508_DP, -0.021608710360898_DP /)

   REAL(DP) :: sqx, xx0, Qtxb, atg, fx
   REAL(DP) :: txb, txbfx, itxbQ
   !
   sqx = x * x                          ! x^2 = r_s
   xx0 = x - x0(i)                      ! x - x_0
   Qtxb = Q(i) / (2.0_DP*x + b(i))      ! Q / (2x+b)
   atg = ATAN(Qtxb)                     ! tan^-1(Q/(2x+b))
   fx = sqx + b(i)*x + c(i)             ! X(x) = x^2 + b*x + c
   !
   fit = A(i) * (  LOG(sqx/fx) + tbQ(i)*atg - &
       bx0fx0(i) * ( LOG(xx0*xx0/fx) + (tbQ(i) + 4.0_DP*x0(i)/Q(i)) * atg )  )
   !
   txb = 2.0_DP*x + b(i)
   txbfx = txb / fx
   itxbQ = 1.0_DP / (txb*txb + Q(i)*Q(i))
   !
   dfit = fit - A(i) / 3.0_DP + A(i)*x/6.0_DP * (  txbfx + 4.0_DP*b(i)*itxbQ + &
           bx0fx0(i) * ( 2.0_DP/xx0 - txbfx - 4.0_DP*(b(i)+2.0_DP*x0(i))*itxbQ )  )
   !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------------------------
SUBROUTINE vwn1_rpa_spin( rs, zeta, ec, vc_up, vc_dw )                    !<GPU:DEVICE>
   !---------------------------------------------------------------------------------------
   !
   USE kinds, ONLY: DP
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
   REAL(DP) :: zeta3, zeta4, trup, trdw, trup13, trdw13, fz, dfz, fzz4
   REAL(DP) :: SQRTrs, ecP, ecF, ac, De, vcP, vcF, dac, dec1, dec2
   REAL(DP) :: cfz, cfz1, cfz2, iddfz0
   !
   ! PARAMETERs:   e_c/para,    e_c/ferro,     alpha_c
   REAL(DP), PARAMETER :: &
      A(3)  = (/  0.0310907_DP, 0.01554535_DP, -0.01688686394039_DP /), &
      x0(3) = (/ -0.409286_DP, -0.743294_DP,   -0.228344_DP /), &
      b(3)  = (/ 13.0720_DP,   20.1231_DP,      1.06835_DP  /), &
      c(3)  = (/ 42.7198_DP,  101.578_DP,      11.4813_DP /),&
      Q(3)  = (/  0.044899888641577_DP,      1.171685277708971_DP,  6.692072046645942_DP /), &
      tbQ(3) = (/ 582.273159042780890_DP,   34.348984975465861_DP,  0.319288254087299_DP /), &
      fx0(3) = (/ 37.537128437796000_DP,    87.173106479036008_DP, 11.289489669936000_DP /), &
      bx0fx0(3) = (/ -0.142530524167984_DP, -0.171582499414508_DP, -0.021608710360898_DP /)
   ! N.B.: A is expressed in Hartree
   ! Q = SQRT(4*c - b^2)
   ! tbQ = 2*b/Q
   ! fx0 = X(x_0) = x_0^2 + b*x_0 + c
   ! bx0fx0 = b*x_0/X(x_0)
   !
   !
   ! coefficients for f(z), df/dz, ddf/ddz(0)
   cfz = 2.0_DP**(4.0_DP/3.0_DP) - 2.0_DP
   cfz1 = 1.0_DP / cfz
   cfz2 = 4.0_DP/3.0_DP * cfz1
   iddfz0 = 9.0_DP / 8.0_DP *cfz
   !
   SQRTrs = SQRT(rs)
   zeta3 = zeta**3
   zeta4 = zeta3*zeta
   trup = 1.0_DP + zeta
   trdw = 1.0_DP - zeta
   trup13 = trup**(1.0_DP/3.0_DP)
   trdw13 = trdw**(1.0_DP/3.0_DP)
   fz = cfz1 * (trup13*trup + trdw13*trdw - 2.0_DP)     ! f(zeta)
   dfz = cfz2 * (trup13 - trdw13)                       ! df / dzeta
   !
   CALL padefit_ParSet2( sqrtrs, 1, ecP, vcP ) ! ecF = e_c Paramagnetic    !<GPU:padefit_ParSet2=>padefit_ParSet2_d>
   CALL padefit_ParSet2( sqrtrs, 2, ecF, vcF ) ! ecP = e_c Ferromagnetic   !<GPU:padefit_ParSet2=>padefit_ParSet2_d>
   CALL padefit_ParSet2( sqrtrs, 3, ac,  dac ) ! ac = "spin stiffness"     !<GPU:padefit_ParSet2=>padefit_ParSet2_d>
   !
   ac = ac * iddfz0
   dac = dac * iddfz0
   De = ecF - ecP - ac                                  ! e_c[F]-e_c[P]-alpha_c/(ddf/ddz(z=0))
   fzz4 = fz * zeta4
   ec = ecP + ac * fz  + De * fzz4
   !
   dec1 = vcP + dac*fz + (vcF - vcP - dac) * fzz4       ! e_c - (r_s/3)*(de_c/dr_s)
   dec2 = ac*dfz + De*(4.0_DP*zeta3*fz + zeta4*dfz)     ! de_c/dzeta
   !
   vc_up = dec1 + (1.0_DP - zeta)*dec2                  ! v_c[s] = e_c - (r_s/3)*(de_c/dr_s)
   vc_dw = dec1 - (1.0_DP + zeta)*dec2                  !          +[sign(s)-zeta]*(de_c/dzeta)
   !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE pw_spin( rs, zeta, ec, vc_up, vc_dw )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! J.P. Perdew and Y. Wang, PRB 45, 13244 (1992).
  !
  USE kinds, ONLY : DP
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
END SUBROUTINE pw_spin
!
!
!-----------------------------------------------------------------------------
SUBROUTINE lsd_lyp( rho, zeta, elyp, vlyp_up, vlyp_dw )                    !<GPU:DEVICE>
  !==--------------------------------------------------------------==
  !==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
  !==  THIS IS ONLY THE LDA PART                                   ==
  !==--------------------------------------------------------------==
  !
  USE kinds,       ONLY: DP
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
END SUBROUTINE lsd_lyp

END MODULE
