!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE corr_lda_l !<GPU:corr_lda=>corr_lda_gpu>

 CONTAINS

!-------------------------------------------------------------------------
SUBROUTINE pz_l( rs, iflag, ec, vc )                    !<GPU:DEVICE>
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
END SUBROUTINE pz_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE pzKZK_l( rs, ec, vc, vol )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! LDA parametrization from Monte Carlo DATA:
  !
  !! * iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !! * iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
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
  REAL(DP) :: vol                                        !<GPU:VALUE>
  !! volume element
  !
  ! ... local variables
  !
  INTEGER  :: iflag, kr
  REAL(DP) :: ec0(2), vc0(2), ec0p
  REAL(DP) :: a(2), b(2), c(2), d(2), gc(2), b1(2), b2(2)
  REAL(DP) :: lnrs, rs12, ox, dox, lnrsk, rsk
  REAL(DP) :: a1, grs, g1, g2, g3, g4, dL, gh, gl, grsp
  REAL(DP) :: f3, f2, f1, f0, pi
  REAL(DP) :: D1, D2, D3, P1, P2, ry2h
  !
  DATA a / 0.0311_DP, 0.031091_DP /, b / -0.048_DP,  -0.046644_DP /, &
       c / 0.0020_DP, 0.00419_DP  /, d / -0.0116_DP, -0.00983_DP /
  !
  DATA gc / -0.1423_DP, -0.103756_DP /, b1 / 1.0529_DP, 0.56371_DP /, &
       b2 / 0.3334_DP, 0.27358_DP /
  !
  DATA a1 / -2.2037_DP /, g1 / 0.1182_DP/, g2 / 1.1656_DP/, g3 / -5.2884_DP/, &
       g4 / -1.1233_DP /
  !
  DATA ry2h / 0.5_DP /
  !
  iflag = 1
  pi = 4.d0 * ATAN(1.d0)
  dL = vol**(1.d0/3.d0)
  gh = 0.5d0 * dL / (2.d0 * pi)**(1.d0/3.d0)
  gl = dL * (3.d0 / 2.d0 / pi)**(1.d0/3.d0)
  !
  rsk = gh
  !
  DO kr = 1, 2
     !
     lnrsk = LOG(rsk)
     IF (rsk < 1.0d0) THEN
        ! high density formula
        ec0(kr) = a(iflag) *lnrsk + b(iflag) + c(iflag) * rsk * lnrsk + d(iflag) * rsk
        vc0(kr) = a(iflag) * lnrsk + (b(iflag) - a(iflag) / 3.d0) + 2.d0 / &
                  3.d0 * c (iflag) * rsk * lnrsk + (2.d0 * d (iflag) - c (iflag) ) &
                  / 3.d0 * rsk
        !
     ELSE
        ! interpolation formula
        rs12 = SQRT(rsk)
        ox  = 1.d0 + b1 (iflag) * rs12 + b2 (iflag) * rsk
        dox = 1.d0 + 7.d0 / 6.d0 * b1 (iflag) * rs12 + 4.d0 / 3.d0 * &
              b2 (iflag) * rsk
        ec0(kr) = gc (iflag) / ox
        vc0(kr) = ec0(kr) * dox / ox
        !
     ENDIF
     !
     grs  = g1 * rsk * lnrsk + g2 * rsk + g3 * rsk**1.5d0 + g4 * rsk**2.d0
     grsp = g1 * lnrsk + g1 + g2 + 1.5d0 * g3 * rsk**0.5d0 + 2.d0 * g4 * rsk
     !
     ec0(kr)  = ec0(kr) + (-a1 * rsk / dL**2.d0 + grs / dL**3.d0) * ry2h
     vc0(kr)  = vc0(kr) + (-2.d0 * a1 * rsk / dL**2.d0 / 3.d0 + &
           grs / dL**3.d0 -  grsp * rsk / 3.d0 / dL**3.d0) * ry2h
     !
     rsk = rs
     !
  ENDDO
  !
  lnrs = LOG(rs)
  !
  IF (rs <= gh) THEN
     ec = ec0(2)
     vc = vc0(2)
  ELSE
     IF ( rs <= gl ) THEN
        ec0p = 3.d0 * (ec0(1) - vc0(1)) / gh
        P1 = 3.d0 *  ec0(1) - gh * ec0p
        P2 = ec0p
        D1 = gl - gh
        D2 = gl**2.d0 - gh**2.d0
        D3 = gl**3.d0 - gh**3.d0
        f2 = 2.d0 * gl**2.d0 * P2 * D1 + D2 * P1
        f2 = f2 / (-(2.d0*gl*D1)**2.d0 + 4.d0*gl*D1*D2 - D2**2.d0 )
        f3 = - (P2 + 2.d0*D1*f2) / (3.d0 * D2)
        f1 = - (P1 + D2 * f2) / (2.d0 * D1)
        f0 = - gl * (gl * f2 + 2.d0 * f1) / 3.d0
        !
        ec = f3 * rs**3.d0 + f2 * rs**2.d0 + f1 * rs + f0
        vc = f2 * rs**2.d0 / 3.d0 + f1 * 2.d0 * rs / 3.d0 + f0
     ELSE
        ec = 0.d0
        vc = 0.d0
     ENDIF
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE pzKZK_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE vwn_l( rs, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980).
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
  REAL(DP), PARAMETER :: a = 0.0310907d0, b = 3.72744d0, &
                         c = 12.9352d0, x0 = -0.10498d0
  REAL(DP) :: q, f1, f2, f3, rs12, fx, qx, tx, tt
  !
  !
  q  = SQRT( 4.d0*c - b*b )
  f1 = 2.d0*b/q
  f2 = b*x0 / ( x0*x0 + b*x0 + c )
  f3 = 2.d0*( 2.d0 * x0 + b ) / q
  !
  rs12 = SQRT(rs)
  fx = rs + b*rs12 + c
  qx = ATAN( q / (2.d0*rs12 + b) )
  !
  ec = a * ( LOG( rs/fx ) + f1*qx - f2*( LOG( (rs12 - x0)**2 / fx) &
             + f3*qx) )
  !
  tx = 2.d0*rs12 + b
  tt = tx*tx + q*q
  !
  vc = ec - rs12*a / 6.d0*( 2.d0 / rs12 - tx/fx - 4.d0*b/tt - f2 * &
            (2.d0 / (rs12 - x0) - tx/fx - 4.d0*(2.d0 * x0 + b)/tt) )
  !
  RETURN
  !
END SUBROUTINE vwn_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE vwn1_rpa_l( rs, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980).
  !
  USE kind_l,       ONLY: DP
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
  REAL(DP), PARAMETER :: a = 0.0310907_DP, b = 13.0720_DP, &
                         c = 42.7198_DP, x0 = -0.409286_DP
  REAL(DP) :: q, f1, f2, f3, rs12, fx, qx, tx, tt
  !
  q  = SQRT(4.d0*c - b*b)
  f1 = 2.d0*b/q
  f2 = b*x0 / (x0*x0 + b*x0 + c)
  f3 = 2.d0 * (2.d0 * x0 + b) / q
  !
  rs12 = SQRT(rs)
  fx = rs + b * rs12 + c
  qx = ATAN(q / (2.d0 * rs12 + b) )
  !
  ec = a*( LOG(rs/fx) + f1*qx - f2*(LOG((rs12 - x0)**2/fx) + f3*qx) )
  !
  tx = 2.d0*rs12 + b
  tt = tx*tx + q*q
  !
  vc = ec - rs12*a/6.d0*( 2.d0/rs12 - tx/fx - 4.d0*b/tt - f2 * &
            (2.d0/(rs12 - x0) - tx/fx - 4.d0*(2.d0*x0 + b)/tt) )
  !
  !
  RETURN
  !
END SUBROUTINE vwn1_rpa_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE lyp_l( rs, ec, vc )                    !<GPU:DEVICE>
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
END SUBROUTINE lyp_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE pw_l( rs, iflag, ec, vc )                    !<GPU:DEVICE>
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
END SUBROUTINE pw_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE wignerc_l( rs, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! Wigner correlation.
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
  REAL(DP) :: rho13 !rho13=rho^(1/3)
  REAL(DP), PARAMETER :: pi34 = 0.6203504908994d0
  !                      pi34 = (3/4pi)^(1/3)
  !
  !
  rho13 = pi34 / rs
  vc = - rho13 * ( (0.943656d0 + 8.8963d0 * rho13) / (1.d0 + &
         12.57d0 * rho13)**2 )
  ec = - 0.738d0 * rho13 * ( 0.959d0 / (1.d0 + 12.57d0 * rho13) )
  !
  RETURN
  !
END SUBROUTINE wignerc_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE hl_l( rs, ec, vc )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! L. Hedin and  B.I. Lundqvist,  J. Phys. C 4, 2064 (1971).
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
  REAL(DP) :: a, x
  !
  !
  a = LOG(1.0d0 + 21.d0/rs)
  x = rs / 21.0d0
  ec = a + (x**3 * a - x*x) + x/2.d0 - 1.0d0/3.0d0
  ec = - 0.0225d0 * ec
  vc = - 0.0225d0 * a
  !
  RETURN
  !
END SUBROUTINE hl_l
!
!
!-----------------------------------------------------------------------
SUBROUTINE gl_l( rs, ec, vc )                    !<GPU:DEVICE>
  !---------------------------------------------------------------------
  !! O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976).
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
  REAL(DP) :: x
  REAL(DP), PARAMETER :: c=0.0333d0, r=11.4d0
  !                      c=0.0203,   r=15.9   for the paramagnetic case
  !
  x = rs / r
  vc = - c * LOG(1.d0 + 1.d0 / x)
  ec = - c * ( (1.d0 + x**3) * LOG(1.d0 + 1.d0 / x) - 1.0d0 / &
                                     3.0d0 + x * (0.5d0 - x) )
  !
  RETURN
  !
END SUBROUTINE gl_l
!
!
END MODULE
