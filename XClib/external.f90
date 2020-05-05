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
