!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE pz( length, rs, iflag, ec, vc )
  !-----------------------------------------------------------------------
  !     LDA PARAMETERization from Monte Carlo DATA
  !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  USE kinds,  ONLY: DP
! #if defined(__LIBXC)
!   USE xc_f90_types_m
!   USE xc_f90_lib_m
! #endif
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length, iflag
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
! #if defined(__LIBXC)
!   REAL(DP), DIMENSION(length) :: rho 
!   REAL(DP), PARAMETER   :: pi34 = 0.6203504908994d0 ! pi34=(3/4pi)^(1/3)
!   INTEGER  :: func_id = 9   ! Perdew & Zunger
!   INTEGER  :: size
!   TYPE(xc_f90_pointer_t) :: xc_func
!   TYPE(xc_f90_pointer_t) :: xc_info
!   !
!   size=length
!   !
!   IF (iflag == 1)  func_id = 9   ! Perdew & Zunger
!   IF (iflag == 2)  func_id = 11  ! Ortiz & Ballone (PZ)
!   !
!   rho = (pi34/rs)**3
!   CALL xc_f90_func_init( xc_func, xc_info, func_id, XC_UNPOLARIZED )
!   CALL xc_f90_lda_exc_vxc( xc_func, size, rho(1), ec(1), vc(1) )
!   CALL xc_f90_func_end( xc_func )
! #else
  !
  REAL(DP) :: a(2), b(2), c(2), d(2), gc(2), b1(2), b2(2)
  REAL(DP) :: lnrs, rs12, ox, dox
  INTEGER  :: i
  !
  DATA a / 0.0311d0, 0.031091d0 /, b / -0.048d0,  -0.046644d0 /, &
       c / 0.0020d0, 0.00419d0  /, d / -0.0116d0, -0.00983d0  /
  !
  DATA gc / -0.1423d0, -0.103756d0 /, b1 / 1.0529d0, 0.56371d0 /, &
       b2 /  0.3334d0,  0.27358d0  /
  !
  DO i = 1, length
     IF ( rs(i) < 1.0d0 ) THEN
       !
       ! high density formula
       lnrs = LOG(rs(i))
       !
       ec(i) = a(iflag)*lnrs + b(iflag) + c(iflag)*rs(i)*lnrs + d(iflag)*rs(i)
       !
       vc(i) = a(iflag)*lnrs + ( b(iflag) - a(iflag)/3.d0 ) + 2.d0/3.d0 * &
               c(iflag)*rs(i)*lnrs + ( 2.d0*d(iflag) - c(iflag) )/3.d0*rs(i)
     ELSE
       !
       ! interpolation formula
       rs12 = SQRT(rs(i))
       !
       ox  = 1.d0 + b1(iflag)*rs12 + b2(iflag)*rs(i)
       dox = 1.d0 + 7.d0/6.d0*b1(iflag)*rs12 + 4.d0/3.d0*b2(iflag)*rs(i)
       !
       ec(i) = gc(iflag)/ox
       vc(i) = ec(i)*dox/ox
       !
     ENDIF
  ENDDO
! #endif
  !
  RETURN
  !
END SUBROUTINE pz
!
!
!-----------------------------------------------------------------------
SUBROUTINE pzKZK( length, rs, ec, vc, vol )
  !-----------------------------------------------------------------------
  !     LDA PARAMETERization from Monte Carlo DATA
  !     iflag=1: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  USE kinds,      ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  INTEGER  :: iflag, i
  REAL(DP) :: ec0(2), vc0(2), ec0p 
  REAL(DP) :: a(2), b(2), c(2), d(2), gc(2), b1(2), b2(2)
  REAL(DP) :: lnrs, rs12, ox, dox, lnrsk, rsk
  REAL(DP) :: a1, grs, g1, g2, g3, g4, dL, vol, gh, gl, grsp
  REAL(DP) :: f3, f2, f1, f0, pi
  REAL(DP) :: D1, D2, D3, P1, P2, ry2h
  !
  DATA a / 0.0311_dp, 0.031091_dp /, b / -0.048_dp,  -0.046644_dp /, &
       c / 0.0020_dp, 0.00419_dp  /, d / -0.0116_dp, -0.00983_dp /
  !
  DATA gc / -0.1423_dp, -0.103756_dp /, b1 / 1.0529_dp, 0.56371_dp /, &
       b2 / 0.3334_dp, 0.27358_dp /
  !
  DATA a1 / -2.2037_dp /, g1 / 0.1182_dp/, g2 / 1.1656_dp/, g3 / -5.2884_dp/, &
       g4 / -1.1233_dp /
  !
  DATA ry2h / 0.5_dp /
  !
  iflag = 1
  pi = 4.d0*ATAN(1.d0)
  dL = vol**(1.d0/3.d0)
  gh = 0.5d0 * dL / (2.d0 *pi)**(1.d0/3.d0)
  gl = dL * (3.d0 / 2.d0 / pi)**(1.d0/3.d0)
  !
  !-------------------------
  lnrsk = LOG(rsk)
  IF ( rsk < 1.0d0 ) THEN
     ! high density formula
     ec0(1) = a(iflag)*lnrsk + b(iflag) + c(iflag)*rsk*lnrsk + d(iflag)*rsk
     vc0(1) = a(iflag)*lnrsk + ( b(iflag) - a(iflag)/3.d0) + 2.d0/3.d0* &
               c(iflag)*rsk*lnrsk + ( 2.d0*d(iflag) - c(iflag) )/3.d0*rsk
  ELSE
     ! interpolation formula
     rs12 = SQRT(rsk)
     ox  = 1.d0 + b1(iflag)*rs12 + b2(iflag)*rsk
     dox = 1.d0 + 7.d0/6.d0*b1(iflag)*rs12 + 4.d0/3.d0*b2(iflag)*rsk
     ec0(1) = gc(iflag)/ox
     vc0(1) = ec0(1)*dox/ox
  ENDIF
  !
  grs  = g1*rsk*lnrsk + g2*rsk + g3*rsk**1.5d0 + g4*rsk**2.d0
  grsp = g1*lnrsk + g1 + g2 + 1.5d0*g3*rsk**0.5d0 + 2.d0*g4*rsk
  !
  ec0(1) = ec0(1) + ( -a1*rsk/dL**2.d0 + grs/dL**3.d0 )*ry2h
  vc0(1) = vc0(1) + ( -2.d0*a1*rsk/dL**2.d0/3.d0 + &
            grs/dL**3.d0 - grsp*rsk/3.d0/dL**3.d0 )*ry2h
  !--------------------------
  !
  !
  DO i = 1, length
     !
     rsk = rs(i)
     lnrsk = LOG(rsk)
     IF ( rsk < 1.0d0 ) THEN
        !
        ! high density formula
        ec0(2) = a(iflag)*lnrsk + b(iflag) + c(iflag)*rsk*lnrsk + d(iflag)*rsk
        !
        vc0(2) = a(iflag)*lnrsk + ( b(iflag) - a(iflag)/3.d0) + 2.d0/3.d0* &
                  c(iflag)*rsk*lnrsk + ( 2.d0*d(iflag) - c(iflag) )/3.d0*rsk
        !
     ELSE
        !
        ! interpolation formula
        rs12 = SQRT(rsk)
        !
        ox  = 1.d0 + b1(iflag)*rs12 + b2(iflag)*rsk
        dox = 1.d0 + 7.d0/6.d0*b1(iflag)*rs12 + 4.d0/3.d0*b2(iflag)*rsk
        !
        ec0(2) = gc(iflag)/ox
        vc0(2) = ec0(2)*dox/ox
        !
     ENDIF
     !
     grs  = g1*rsk*lnrsk + g2*rsk + g3*rsk**1.5d0 + g4*rsk**2.d0
     grsp = g1*lnrsk + g1 + g2 + 1.5d0*g3*rsk**0.5d0 + 2.d0*g4*rsk
     !
     ec0(2) = ec0(2) + ( -a1*rsk/dL**2.d0 + grs/dL**3.d0 )*ry2h
     vc0(2) = vc0(2) + ( -2.d0*a1*rsk/dL**2.d0/3.d0 + &
               grs/dL**3.d0 - grsp*rsk/3.d0/dL**3.d0 )*ry2h
     !
     !
     lnrs = LOG(rs(i))
     IF ( rs(i) < gh ) THEN
        !
        ec(i) = ec0(2)
        vc(i) = vc0(2)
        !
     ELSE
        !
        IF ( rs(i) < gl ) THEN
           ec0p = 3.d0*( ec0(1) - vc0(1) )/gh
           P1 = 3.d0*ec0(1) - gh*ec0p
           P2 = ec0p
           D1 = gl - gh
           D2 = gl**2.d0 - gh**2.d0
           D3 = gl**3.d0 - gh**3.d0
           f2 = 2.d0*gl**2.d0*P2*D1 + D2*P1
           f2 = f2 / ( -( 2.d0*gl*D1 )**2.d0 + 4.d0*gl*D1*D2 - D2**2.d0 )
           f3 = -( P2 + 2.d0*D1*f2 )/( 3.d0*D2 )
           f1 = -( P1 + D2*f2 ) / ( 2.d0*D1 )
           f0 = -gl*( gl*f2 + 2.d0*f1 )/3.d0
           !
           ec(i) = f3*rs(i)**3.d0 + f2*rs(i)**2.d0 + f1*rs(i) + f0
           vc(i) = f2*rs(i)**2.d0/3.d0 + f1*2.d0*rs(i)/3.d0 + f0
        ELSE
           ec(i) = 0.d0
           vc(i) = 0.d0
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  !
  RETURN
  !
END SUBROUTINE pzKZK
!
!
!-----------------------------------------------------------------------
SUBROUTINE vwn( length, rs, ec, vc )
  !-----------------------------------------------------------------------
  !     S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
  !
  USE kinds,      ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  REAL(DP) :: a, b, c, x0
  PARAMETER (a = 0.0310907d0, b = 3.72744d0, c = 12.9352d0, x0 = -0.10498d0)
  REAL(DP) :: q, f1, f2, f3, rs12, fx, qx, tx, tt
  INTEGER  :: i
  !
  q  = SQRT( 4.d0*c - b*b )
  f1 = 2.d0*b/q
  f2 = b*x0 / ( x0*x0 + b*x0 + c )
  f3 = 2.d0*( 2.d0 * x0 + b ) / q
  !
  DO i = 1, length
     !
     rs12 = SQRT(rs(i))
     fx = rs(i) + b*rs12 + c
     qx = ATAN( q / (2.d0*rs12 + b) )
     !
     ec(i) = a * (LOG( rs(i)/fx ) + f1*qx - f2*( LOG( (rs12 - x0) **2 / fx) &
                 + f3 * qx) )
     !
     tx = 2.d0*rs12 + b
     tt = tx*tx + q*q
     !
     vc(i) = ec(i) - rs12*a / 6.d0*( 2.d0 / rs12 - tx/fx - 4.d0*b/tt - f2 * &
                  (2.d0 / (rs12 - x0) - tx/fx - 4.d0*(2.d0 * x0 + b)/tt) )
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE vwn
!
!
!-----------------------------------------------------------------------
SUBROUTINE vwn1_rpa( length, rs, ec, vc )
  !-----------------------------------------------------------------------
  !     S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
  !
  USE kinds,      ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  REAL(DP) :: a, b, c, x0
  PARAMETER (a = 0.0310907_dp, b = 13.0720_dp, c = 42.7198_dp, x0 = -0.409286_dp)
  REAL(DP) :: q, f1, f2, f3, rs12, fx, qx, tx, tt
  INTEGER  :: i
  !
  q  = SQRT(4.d0*c - b*b)
  f1 = 2.d0*b/q
  f2 = b*x0 / (x0*x0 + b*x0 + c)
  f3 = 2.d0 * (2.d0 * x0 + b) / q
  !
  DO i = 1, length
     !
     rs12 = SQRT(rs(i))
     fx = rs(i) + b * rs12 + c
     qx = ATAN(q / (2.d0 * rs12 + b) )
     !
     ec(i) = a*( LOG(rs(i)/fx) + f1*qx - f2*(LOG((rs12 - x0)**2/fx) + f3*qx) )
     !
     tx = 2.d0*rs12 + b
     tt = tx*tx + q*q
     !
     vc(i) = ec(i) - rs12*a/6.d0*( 2.d0/rs12 - tx/fx - 4.d0*b/tt - f2 * &
                  (2.d0/(rs12 - x0) - tx/fx - 4.d0*(2.d0*x0 + b)/tt) )
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE vwn1_rpa
!
!
!-----------------------------------------------------------------------
SUBROUTINE lyp( length, rs, ec, vc )
  !-----------------------------------------------------------------------
  !     C. Lee, W. Yang, and R.G. Parr, PRB 37, 785 (1988)
  !     LDA part only
  !
  USE kinds,      ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  REAL(DP) :: a, b, c, d, pi43
  PARAMETER (a = 0.04918d0, b = 0.132d0 * 2.87123400018819108d0)
           ! pi43 = (4pi/3)^(1/3)
  PARAMETER (pi43 = 1.61199195401647d0, c = 0.2533d0 * pi43, &
                                        d = 0.349d0  * pi43)
  REAL(DP) :: ecrs, ox
  INTEGER  :: i
  !
  DO i = 1, length
     ecrs = b*EXP( -c*rs(i) )
     ox = 1.d0 / (1.d0 + d*rs(i))
     !
     ec(i) = - a*ox*(1.d0 + ecrs)
     vc(i) = ec(i) - rs(i)/3.d0*a*ox*( d*ox + ecrs*(d*ox + c) )
  ENDDO
  !
  RETURN
  !
END SUBROUTINE lyp
!
!
!-----------------------------------------------------------------------
SUBROUTINE pw( length, rs, iflag, ec, vc )
  !-----------------------------------------------------------------------
  !     iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  USE kinds,      ONLY: DP
! #if defined(__LIBXC)
!   USE xc_f90_types_m
!   USE xc_f90_lib_m
! #endif
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: length, iflag
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
! #if defined(__LIBXC)
!   REAL(DP), DIMENSION(length) :: rho
!   REAL(DP), PARAMETER    :: pi34 = 0.6203504908994d0 ! pi34=(3/4pi)^(1/3)
!   INTEGER :: func_id = 12   ! Perdew & Zunger
!   INTEGER :: size
!   TYPE(xc_f90_pointer_t) :: xc_func
!   TYPE(xc_f90_pointer_t) :: xc_info
!   !
!   size = length
!   !
!   IF (iflag == 1)  func_id = 12  ! Perdew & Wang
!   IF (iflag == 2)  func_id = 14  ! Ortiz & Ballone (PW)
!   !
!   rho = (pi34/rs)**3
!   CALL xc_f90_func_init(xc_func, xc_info, func_id, XC_UNPOLARIZED)    
!   CALL xc_f90_lda_exc_vxc(xc_func,size , rho(1), ec(1), vc(1))
!   CALL xc_f90_func_end(xc_func)  
! #else
  !
  REAL(DP) :: a, b1, b2, c0, c1, c2, c3, d0, d1
  PARAMETER (a = 0.031091d0, b1 = 7.5957d0, b2 = 3.5876d0, c0 = a, &
            c1 = 0.046644d0, c2 = 0.00664d0,c3 = 0.01043d0,d0 = 0.4335d0,&
            d1 = 1.4408d0)
  REAL(DP) :: lnrs, rs12, rs32, rs2, om, dom, oLOG
  REAL(DP) :: a1(2), b3(2), b4(2)
  DATA a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /
  INTEGER  :: i
  !
  ! high- and low-density formulae implemented but not used in PW case
  ! (reason: inconsistencies in PBE/PW91 functionals)
  !
  DO i = 1, length
     !
     IF ( rs(i) < 1d0 .and. iflag == 2 ) THEN
        !
        ! high density formula
        lnrs = LOG(rs(i))
        ec(i) = c0 * lnrs - c1 + c2 * rs(i) * lnrs - c3 * rs(i)
        vc(i) = c0 * lnrs - (c1 + c0 / 3.d0) + 2.d0 / 3.d0 * c2 * rs(i) * &
                lnrs - (2.d0 * c3 + c2) / 3.d0 * rs(i)
        !
     ELSEIF ( rs(i) > 100.d0 .and. iflag == 2 ) THEN
        !
        ! low density formula
        ec(i) = - d0 / rs(i) + d1 / rs(i)**1.5d0
        vc(i) = - 4.d0 / 3.d0 * d0 / rs(i) + 1.5d0 * d1 / rs(i)**1.5d0
        !
     ELSE
        !
        ! interpolation formula
        rs12 = SQRT(rs(i))
        rs32 = rs(i) * rs12
        rs2  = rs(i)**2
        om   = 2.d0*a*( b1*rs12 + b2*rs(i) + b3(iflag) * rs32 + b4(iflag)*rs2 )
        dom  = 2.d0*a*( 0.5d0 * b1 * rs12 + b2 * rs(i) + 1.5d0 * b3(iflag) * &
               rs32 + 2.d0 * b4(iflag) * rs2 )
        olog = LOG( 1.d0 + 1.0d0 / om )
        !
        ec(i)   = - 2.d0 * a * (1.d0 + a1(iflag) * rs(i)) * olog
        vc(i)   = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1(iflag) * rs(i)) &
                 * olog - 2.d0 / 3.d0 * a * (1.d0 + a1(iflag) * rs(i)) * dom / &
                 (om * (om + 1.d0) )
        !
     ENDIF
     !
  ENDDO
  !
! #endif
  !
  RETURN
  !
END SUBROUTINE pw
!
!
!-----------------------------------------------------------------------
SUBROUTINE wignerc( length, rs, ec, vc )
  !-----------------------------------------------------------------------
  !        Wigner correlation
  !
  USE kinds,      ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  REAL(DP) :: pi34, rho13
  PARAMETER (pi34 = 0.6203504908994d0)
  !          pi34 = (3/4pi)^(1/3), rho13=rho^(1/3)
  INTEGER  :: i
  !
  DO i = 1, length
     rho13 = pi34 / rs(i)
     vc(i) = - rho13 * ( (0.943656d0 + 8.8963d0 * rho13) / (1.d0 + &
              12.57d0 * rho13)**2)
     ec(i) = - 0.738d0 * rho13 * (0.959d0 / (1.d0 + 12.57d0 * rho13) )
  ENDDO
  !
  RETURN
  !
END SUBROUTINE wignerc
!
!
!-----------------------------------------------------------------------
SUBROUTINE hl( length, rs, ec, vc )
  !-----------------------------------------------------------------------
  !     L. Hedin and  B.I. Lundqvist,  J. Phys. C 4, 2064 (1971)
  !
  USE kinds,      ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  REAL(DP) :: a, x
  INTEGER  :: i
  !
  DO i = 1, length
     a = LOG(1.0d0 + 21.d0 / rs(i))
     x = rs(i) / 21.0d0
     ec(i) = a + (x**3 * a - x*x) + x / 2.d0 - 1.0d0 / 3.0d0
     ec(i) = - 0.0225d0 * ec(i)
     vc(i) = - 0.0225d0 * a
  ENDDO
  !
  RETURN
  !
END SUBROUTINE hl
!
!-----------------------------------------------------------------------
SUBROUTINE gl( length, rs, ec, vc )
  !---------------------------------------------------------------------
  !  O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976)
  !
  USE kinds,      ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  REAL(DP) :: c, r, x
  PARAMETER (c = 0.0333d0, r = 11.4d0)
  ! c=0.0203, r=15.9 for the paramagnetic case
  INTEGER :: i
  !
  DO i = 1, length 
     x = rs(i) / r
     vc(i) = - c * LOG(1.d0 + 1.d0 / x)
     ec(i) = - c * ( (1.d0 + x**3) * LOG(1.d0 + 1.d0 / x) - 1.0d0 / &
                                           3.0d0 + x * (0.5d0 - x) )
  ENDDO
  !
  RETURN
  !
END SUBROUTINE gl
!
!
! ... LSDA
!
!-----------------------------------------------------------------------
SUBROUTINE pz_polarized( length, rs, ec, vc )
  !-----------------------------------------------------------------------
  !     J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
  !     spin-polarized energy and potential
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rs
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec, vc
  !
  REAL(DP) :: a, b, c, d, gc, b1, b2
  PARAMETER (a = 0.01555d0, b = - 0.0269d0, c = 0.0007d0, d = &
             - 0.0048d0, gc = - 0.0843d0, b1 = 1.3981d0, b2 = 0.2611d0)
  REAL(DP) :: lnrs, rs12, ox, dox
  REAL(DP), PARAMETER :: xcprefact = 0.022575584d0, pi34 = 0.6203504908994d0 
  ! REAL(DP) :: betha, etha, csi, prefact
  INTEGER :: i
  !
  DO i = 1, length
    IF ( rs(i) < 1.0d0 ) THEN
       ! high density formula
       lnrs = LOG(rs(i))
       ec(i) = a * lnrs + b + c * rs(i) * lnrs + d * rs(i)
       vc(i) = a * lnrs + (b - a / 3.d0) + 2.d0 / 3.d0 * c * rs(i) * lnrs + &
            (2.d0 * d-c) / 3.d0 * rs(i)
    ELSE
       ! interpolation formula
       rs12 = SQRT(rs(i))
       ox = 1.d0 + b1 * rs12 + b2 * rs(i)
       dox = 1.d0 + 7.d0 / 6.d0 * b1 * rs12 + 4.d0 / 3.d0 * b2 * rs(i)
       ec(i) = gc / ox
       vc(i) = ec(i) * dox / ox
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
  ENDDO
  !
  RETURN
  !
END SUBROUTINE pz_polarized
!
!
!-----------------------------------------------------------------------
SUBROUTINE pz_spin( length, rs, zeta, ec, vc )
  !-----------------------------------------------------------------------
  !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length)   :: rs, zeta
  REAL(DP), INTENT(OUT), DIMENSION(length)   :: ec
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vc
  !
  REAL(DP), DIMENSION(length) :: ecu, vcu, ecp, vcp
  REAL(DP) :: fz, dfz
  REAL(DP) :: p43, third
  PARAMETER (p43 = 4.0d0 / 3.d0, third = 1.d0 / 3.d0)
  INTEGER :: i
  !
  ! unpolarized part (Perdew-Zunger formula)
  CALL pz( length, rs, 1, ecu, vcu )
  !
  ! polarization contribution
  call pz_polarized ( length, rs, ecp, vcp ) 
  !
  do i = 1, length
    fz = ( (1.0d0 + zeta(i)) **p43 + (1.d0 - zeta(i)) **p43 - 2.d0) / &
         (2.d0**p43 - 2.d0)
    dfz = p43 * ( (1.0d0 + zeta(i)) **third- (1.d0 - zeta(i)) **third) &
         / (2.d0**p43 - 2.d0)
    !
    ec(i) = ecu(i) + fz * (ecp(i) - ecu(i))
    vc(i,1) = vcu(i) + fz * (vcp(i) - vcu(i)) + (ecp(i) - ecu(i)) * dfz * (1.d0 - zeta(i))
    vc(i,2) = vcu(i) + fz * (vcp(i) - vcu(i)) + (ecp(i) - ecu(i)) * dfz * (-1.d0- zeta(i))
  enddo
    !
  return
END SUBROUTINE pz_spin
!
!---------
SUBROUTINE vwn_spin( length, rs_v, zeta_v, ec, vc )
   !
   USE kinds, ONLY: DP
   !
   IMPLICIT NONE
   !
   ! PARAMETERs:   e_c/para,    e_c/ferro,     alpha_c
   REAL(DP), PARAMETER :: &
      A(3)  = (/ 0.0310907_dp, 0.01554535_dp, -0.01688686394039_dp /), &
      x0(3) = (/ -0.10498_dp, -0.32500_dp, -0.0047584_dp /), &
      b(3)  = (/3.72744_dp, 7.06042_dp, 1.13107_dp /), &
      c(3)  = (/ 12.9352_dp, 18.0578_dp, 13.0045_dp /),&
      Q(3)  = (/ 6.15199081975908_dp, 4.73092690956011_dp, 7.12310891781812_dp /), &
      tbQ(3) = (/ 1.21178334272806_dp, 2.98479352354082_dp, 0.31757762321188_dp /), &
      fx0(3) = (/ 12.5549141492_dp, 15.8687885_dp, 12.99914055888256_dp /), &
      bx0fx0(3) = (/ -0.03116760867894_dp, -0.14460061018521_dp, -0.00041403379428_dp /)
   ! N.B.: A is expressed in Hartree
   ! Q = SQRT(4*c - b^2)
   ! tbQ = 2*b/Q
   ! fx0 = X(x_0) = x_0^2 + b*x_0 + c
   ! bx0fx0 = b*x_0/X(x_0)
   !
   INTEGER, INTENT(IN) :: length
   REAL(DP), INTENT(IN),  DIMENSION(length)   :: rs_v, zeta_v
   REAL(DP), INTENT(OUT), DIMENSION(length)   :: ec
   REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vc
   !
   REAL(DP) :: rs,zeta,zeta3, zeta4, trup, trdw, trup13, trdw13, fz, dfz, fzz4
   REAL(DP) :: SQRTrs, ecP, ecF, ac, De, vcP, vcF, dac, dec1, dec2
   REAL(DP) :: cfz, cfz1, cfz2, iddfz0
   INTEGER  :: i
   !
   ! coefficients for f(z), df/dz, ddf/ddz(0)
   cfz = 2.0_dp**(4.0_dp/3.0_dp) - 2.0_dp
   cfz1 = 1.0_dp / cfz
   cfz2 = 4.0_dp/3.0_dp * cfz1
   iddfz0 = 9.0_dp / 8.0_dp *cfz
   !
   DO i=1,length
      !
      zeta=zeta_v(i) ;  rs=rs_v(i)
      !
      SQRTrs = SQRT(rs)
      zeta3 = zeta**3
      zeta4 = zeta3*zeta
      trup = 1.0_dp + zeta
      trdw = 1.0_dp - zeta
      trup13 = trup**(1.0_dp/3.0_dp)
      trdw13 = trdw**(1.0_dp/3.0_dp)
      fz = cfz1 * (trup13*trup + trdw13*trdw - 2.0_dp)         ! f(zeta)
      dfz = cfz2 * (trup13 - trdw13)     ! d f / d zeta
      !
      CALL padefit( sqrtrs, 1, ecP, vcP )            ! ecF = e_c Paramagnetic
      CALL padefit( sqrtrs, 2, ecF, vcF )            ! ecP = e_c Ferromagnetic
      CALL padefit( sqrtrs, 3, ac, dac )             ! ac = "spin stiffness"
      !
      ac = ac * iddfz0
      dac = dac * iddfz0
      De = ecF - ecP - ac ! e_c[F] - e_c[P] - alpha_c/(ddf/ddz(z=0))
      fzz4 = fz * zeta4
      ec(i) = ecP + ac * fz  + De * fzz4
      !
      dec1 = vcP + dac*fz + (vcF - vcP - dac) * fzz4       ! e_c - (r_s/3)*(de_c/dr_s)
      dec2 = ac*dfz + De*(4.0_dp*zeta3*fz + zeta4*dfz)       ! de_c/dzeta
      !
      ! v_c[s] = e_c - (r_s/3)*(de_c/dr_s) + [sign(s)-zeta]*(de_c/dzeta)
      vc(i,1) = dec1 + (1.0_dp - zeta)*dec2
      vc(i,2) = dec1 - (1.0_dp + zeta)*dec2
   ENDDO
   !
   !
  CONTAINS
   !---
   SUBROUTINE padefit( x, i, fit, dfit )
      !----
      ! implements formula [4.4] in:
      ! S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
      !
      USE kinds, ONLY: DP
      !
      IMPLICIT NONE
      !
      ! input
      REAL(DP) :: x     ! x is SQRT(r_s) 
      INTEGER :: i      ! i is the index of the fit

      ! output
      REAL(DP) :: fit, dfit
      ! Pade fit calculated in x and its derivative w.r.t. rho
      ! rs = inv((rho*)^(1/3)) = x^2
      ! fit  [eq. 4.4]
      ! dfit/drho = fit - (rs/3)*dfit/drs = ec - (x/6)*dfit/dx
      !
      ! local
      REAL(DP) :: sqx, xx0, Qtxb, atg, fx
      REAL(DP) :: txb, txbfx, itxbQ
      !
      sqx = x * x                          ! x^2 = r_s
      xx0 = x - x0(i)                      ! x - x_0
      Qtxb = Q(i) / (2.0_dp*x + b(i))      ! Q / (2x+b)
      atg = ATAN(Qtxb)                     ! tan^-1(Q/(2x+b))
      fx = sqx + b(i)*x + c(i)             ! X(x) = x^2 + b*x + c
      !
      fit = A(i) * (  LOG(sqx/fx) + tbQ(i)*atg - &
            bx0fx0(i) * ( LOG(xx0*xx0/fx) + (tbQ(i) + 4.0_dp*x0(i)/Q(i)) * atg )  )
      !
      txb = 2.0_dp*x + b(i)
      txbfx = txb / fx
      itxbQ = 1.0_dp / (txb*txb + Q(i)*Q(i))
      !
      dfit = fit - A(i) / 3.0_dp + A(i)*x/6.0_dp * (  txbfx + 4.0_dp*b(i)*itxbQ + &
              bx0fx0(i) * ( 2.0_dp/xx0 - txbfx - 4.0_dp*(b(i)+2.0_dp*x0(i))*itxbQ )  )
      !
   END SUBROUTINE
   !
END SUBROUTINE

!---------
SUBROUTINE vwn1_rpa_spin(length, rs_v, zeta_v, ec, vc)
   !
   USE kinds, ONLY: DP
   !
   IMPLICIT NONE
   !
   ! PARAMETERs:   e_c/para,    e_c/ferro,     alpha_c
   REAL(DP), PARAMETER :: &
      A(3)  = (/  0.0310907_dp, 0.01554535_dp, -0.01688686394039_dp /), &
      x0(3) = (/ -0.409286_dp, -0.743294_dp,   -0.228344_dp /), &
      b(3)  = (/ 13.0720_dp,   20.1231_dp,      1.06835_dp  /), &
      c(3)  = (/ 42.7198_dp,  101.578_dp,      11.4813_dp /),&
      Q(3)  = (/  0.044899888641577_dp,      1.171685277708971_dp,  6.692072046645942_dp /), &
      tbQ(3) = (/ 582.273159042780890_dp,   34.348984975465861_dp,  0.319288254087299_dp /), &
      fx0(3) = (/ 37.537128437796000_dp,    87.173106479036008_dp, 11.289489669936000_dp /), &
      bx0fx0(3) = (/ -0.142530524167984_dp, -0.171582499414508_dp, -0.021608710360898_dp /)
   ! N.B.: A is expressed in Hartree
   ! Q = SQRT(4*c - b^2)
   ! tbQ = 2*b/Q
   ! fx0 = X(x_0) = x_0^2 + b*x_0 + c
   ! bx0fx0 = b*x_0/X(x_0)
   !
   INTEGER,  INTENT(IN) :: length
   REAL(DP), INTENT(IN),  DIMENSION(length)   :: rs_v, zeta_v
   REAL(DP), INTENT(OUT), DIMENSION(length)   :: ec
   REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vc
   !
   ! local
   REAL(DP) :: rs,zeta,zeta3, zeta4, trup, trdw, trup13, trdw13, fz, dfz, fzz4 
   REAL(DP) :: SQRTrs, ecP, ecF, ac, De, vcP, vcF, dac, dec1, dec2
   REAL(DP) :: cfz, cfz1, cfz2, iddfz0
   INTEGER  :: i
   !
   ! coefficients for f(z), df/dz, ddf/ddz(0)
   cfz = 2.0_dp**(4.0_dp/3.0_dp) - 2.0_dp
   cfz1 = 1.0_dp / cfz
   cfz2 = 4.0_dp/3.0_dp * cfz1
   iddfz0 = 9.0_dp / 8.0_dp *cfz
   !
   DO i=1,length
     rs=rs_v(i)
     zeta=zeta_v(i)
     !
     SQRTrs = SQRT(rs)
     zeta3 = zeta**3
     zeta4 = zeta3*zeta
     trup = 1.0_dp + zeta
     trdw = 1.0_dp - zeta
     trup13 = trup**(1.0_dp/3.0_dp)
     trdw13 = trdw**(1.0_dp/3.0_dp)
     fz = cfz1 * (trup13*trup + trdw13*trdw - 2.0_dp)         ! f(zeta)
     dfz = cfz2 * (trup13 - trdw13)     ! d f / d zeta
     !
     CALL padefit( sqrtrs, 1, ecP, vcP )            ! ecF = e_c Paramagnetic
     CALL padefit( sqrtrs, 2, ecF, vcF )            ! ecP = e_c Ferromagnetic
     CALL padefit( sqrtrs, 3, ac, dac )             ! ac = "spin stiffness"
     !
     ac = ac * iddfz0
     dac = dac * iddfz0
     De = ecF - ecP - ac ! e_c[F] - e_c[P] - alpha_c/(ddf/ddz(z=0))
     fzz4 = fz * zeta4
     ec(i) = ecP + ac * fz  + De * fzz4
     !
     dec1 = vcP + dac*fz + (vcF - vcP - dac) * fzz4       ! e_c - (r_s/3)*(de_c/dr_s)
     dec2 = ac*dfz + De*(4.0_dp*zeta3*fz + zeta4*dfz)       ! de_c/dzeta
     !
     ! v_c[s] = e_c - (r_s/3)*(de_c/dr_s) + [sign(s)-zeta]*(de_c/dzeta)
     vc(i,1) = dec1 + (1.0_dp - zeta)*dec2
     vc(i,2) = dec1 - (1.0_dp + zeta)*dec2
   ENDDO
   !
   !
  CONTAINS
   !---
   SUBROUTINE padefit(x, i, fit, dfit)
      !----
      ! implements formula [4.4] in:
      ! S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
      !
      USE kinds, ONLY: DP
      !
      IMPLICIT NONE
      !
      ! input
      REAL(DP) :: x     ! x is SQRT(r_s) 
      INTEGER :: i      ! i is the index of the fit

      ! output
      REAL(DP) :: fit, dfit
      ! Pade fit calculated in x and its derivative w.r.t. rho
      ! rs = inv((rho*)^(1/3)) = x^2
      ! fit  [eq. 4.4]
      ! dfit/drho = fit - (rs/3)*dfit/drs = ec - (x/6)*dfit/dx

      ! local
      REAL(DP) :: sqx, xx0, Qtxb, atg, fx
      REAL(DP) :: txb, txbfx, itxbQ

      sqx = x * x                          ! x^2 = r_s
      xx0 = x - x0(i)                      ! x - x_0
      Qtxb = Q(i) / (2.0_dp*x + b(i))      ! Q / (2x+b)
      atg = ATAN(Qtxb)                     ! tan^-1(Q/(2x+b))
      fx = sqx + b(i)*x + c(i)             ! X(x) = x^2 + b*x + c

      fit = A(i) * (  LOG(sqx/fx) + tbQ(i)*atg - &
            bx0fx0(i) * ( LOG(xx0*xx0/fx) + (tbQ(i) + 4.0_dp*x0(i)/Q(i)) * atg )  )

      txb = 2.0_dp*x + b(i)
      txbfx = txb / fx
      itxbQ = 1.0_dp / (txb*txb + Q(i)*Q(i))

      dfit = fit - A(i) / 3.0_dp + A(i)*x/6.0_dp * (  txbfx + 4.0_dp*b(i)*itxbQ + &
              bx0fx0(i) * ( 2.0_dp/xx0 - txbfx - 4.0_dp*(b(i)+2.0_dp*x0(i))*itxbQ )  )

   END SUBROUTINE
   !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE pw_spin( length, rs_v, zeta_v, ec, vc )
  !-----------------------------------------------------------------------
  !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length)   :: rs_v, zeta_v
  REAL(DP), INTENT(OUT), DIMENSION(length)   :: ec
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vc
  !
  ! xc PARAMETERs, unpolarised
  REAL(DP) :: rs, zeta
  REAL(DP) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
  PARAMETER (a = 0.031091d0, a1 = 0.21370d0, b1 = 7.5957d0, b2 = &
       3.5876d0, b3 = 1.6382d0, b4 = 0.49294d0, c0 = a, c1 = 0.046644d0, &
       c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, d1 = 1.4408d0)
  ! xc PARAMETERs, polarised
  REAL(DP) :: ap, a1p, b1p, b2p, b3p, b4p, c0p, c1p, c2p, c3p, d0p, &
       d1p
  PARAMETER (ap = 0.015545d0, a1p = 0.20548d0, b1p = 14.1189d0, b2p &
       = 6.1977d0, b3p = 3.3662d0, b4p = 0.62517d0, c0p = ap, c1p = &
       0.025599d0, c2p = 0.00319d0, c3p = 0.00384d0, d0p = 0.3287d0, d1p &
       = 1.7697d0)
  ! xc PARAMETERs, antiferro
  REAL(DP) :: aa, a1a, b1a, b2a, b3a, b4a, c0a, c1a, c2a, c3a, d0a, &
       d1a
  PARAMETER (aa = 0.016887d0, a1a = 0.11125d0, b1a = 10.357d0, b2a = &
       3.6231d0, b3a = 0.88026d0, b4a = 0.49671d0, c0a = aa, c1a = &
       0.035475d0, c2a = 0.00188d0, c3a = 0.00521d0, d0a = 0.2240d0, d1a &
       = 0.3969d0)
  REAL(DP) :: fz0
  PARAMETER (fz0 = 1.709921d0)
  REAL(DP) :: rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz
  REAL(DP) :: om, dom, olog, epwc, vpwc
  REAL(DP) :: omp, domp, ologp, epwcp, vpwcp
  REAL(DP) :: oma, doma, ologa, alpha, vpwca
  INTEGER :: i
  !
  !     if(rs.lt.0.5d0) then
  ! high density formula (not implemented)
  !
  !     else if(rs.gt.100.d0) then
  ! low density formula  (not implemented)
  !
  !     else
  ! interpolation formula
  
  DO i = 1,length
    zeta=zeta_v(i)
    rs=rs_v(i)
    !
    zeta2 = zeta * zeta
    zeta3 = zeta2 * zeta
    zeta4 = zeta3 * zeta
    rs12 = SQRT (rs)
    rs32 = rs * rs12
    rs2 = rs**2
    ! unpolarised
    om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
    dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 * rs32 &
         + 2.d0 * b4 * rs2)
    olog = LOG (1.d0 + 1.0d0 / om)
    epwc = - 2.d0 * a * (1.d0 + a1 * rs) * olog
    vpwc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 * rs) * olog - 2.d0 / &
          3.d0 * a * (1.d0 + a1 * rs) * dom / (om * (om + 1.d0) )
    ! polarized
    omp  = 2.d0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
    domp = 2.d0 * ap * (0.5d0 * b1p * rs12 + b2p * rs + 1.5d0 * b3p * &
          rs32 + 2.d0 * b4p * rs2)
    ologp = LOG (1.d0 + 1.0d0 / omp)
    epwcp = - 2.d0 * ap * (1.d0 + a1p * rs) * ologp
    vpwcp = - 2.d0 * ap * (1.d0 + 2.d0 / 3.d0 * a1p * rs) * ologp - &
       2.d0 / 3.d0 * ap * (1.d0 + a1p * rs) * domp / (omp * (omp + 1.d0))
    ! antiferro
    oma = 2.d0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
    doma = 2.d0 * aa * (0.5d0 * b1a * rs12 + b2a * rs + 1.5d0 * b3a * &
            rs32 + 2.d0 * b4a * rs2)
    ologa = LOG (1.d0 + 1.0d0 / oma)
    alpha = 2.d0 * aa * (1.d0 + a1a * rs) * ologa
    vpwca = + 2.d0 * aa * (1.d0 + 2.d0 / 3.d0 * a1a * rs) * ologa + &
            2.d0 / 3.d0 * aa * (1.d0 + a1a * rs) * doma / (oma * (oma + 1.d0))
    !
    fz = ( (1.d0 + zeta) ** (4.d0 / 3.d0) + (1.d0 - zeta) ** (4.d0 / &
            3.d0) - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2.d0)
    dfz = ( (1.d0 + zeta) ** (1.d0 / 3.d0) - (1.d0 - zeta) ** (1.d0 / &
            3.d0) ) * 4.d0 / (3.d0 * (2.d0** (4.d0 / 3.d0) - 2.d0) )
    !
    ec(i) = epwc + alpha * fz * (1.d0 - zeta4) / fz0 + (epwcp - epwc) &
              * fz * zeta4
    !
    vc(i,1) = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
               * fz * zeta4 + (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
               zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
               * (1.d0 - zeta)

    vc(i,2) = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
               * fz * zeta4 - (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
               zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
               * (1.d0 + zeta)
    !
  ENDDO  
  !        endif
  !
  RETURN
  !
END SUBROUTINE pw_spin
!
!
!-----------------------------------------------------------------------------
SUBROUTINE lsd_lyp( length, rho, zeta, elyp, vlyp )
  !==--------------------------------------------------------------==
  !==  C. LEE, W. YANG, AND R.G. PARR, PRB 37, 785 (1988)          ==
  !==  THIS IS ONLY THE LDA PART                                   ==
  !==--------------------------------------------------------------==
  !
  USE kinds,       ONLY: DP
  !
  IMPLICIT NONE
  !
  ! ... arguments
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN),  DIMENSION(length)   :: rho, zeta
  REAL(DP), INTENT(OUT), DIMENSION(length)   :: elyp
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vlyp
  !
  ! ... locals
  !
  REAL(DP) :: ra,rb,rm3,dr,e1,or,dor,e2,de1a,de1b,de2a,de2b
  REAL(DP), PARAMETER :: small=1.D-24, a=0.04918D0, b=0.132D0, &
                         c=0.2533D0,   d=0.349D0,  cf=2.87123400018819108D0
  INTEGER :: i
  !==--------------------------------------------------------------==
  DO i = 1, length
     !
     ra = rho(i)*0.5D0*(1.D0+zeta(i))
     ra = MAX(ra,small)
     rb = rho(i)*0.5D0*(1.D0-zeta(i))
     rb = MAX(rb,small)
     !
     rm3= rho(i)**(-1.D0/3.D0)
     dr = (1.D0+d*rm3)
     !
     e1 = 4.D0*a*ra*rb/rho(i)/dr
     or = EXP(-c*rm3)/dr*rm3**11.D0
     dor= -1.D0/3.D0*rm3**4*or*(11.D0/rm3-c-d/dr)
     e2 = 2.D0**(11.D0/3.D0)*cf*a*b*or*ra*rb*(ra**(8.d0/3.d0)+ rb**(8.d0/3.d0))
     !
     elyp(i) = (-e1-e2)/rho(i)
     !
     de1a = -e1*(1.D0/3.D0*d*rm3**4/dr+1./ra-1./rho(i))
     de1b = -e1*(1.D0/3.D0*d*rm3**4/dr+1./rb-1./rho(i))
     de2a = -2.D0**(11.D0/3.D0)*cf*a*b*(dor*ra*rb*(ra**(8.d0/3.d0)+ &
             rb**(8.d0/3.d0))+or*rb*(11.d0/3.d0*ra**(8.d0/3.d0)+ &
             rb**(8.d0/3.d0)))
     de2b = -2.D0**(11.D0/3.D0)*cf*a*b*(dor*ra*rb*(ra**(8.d0/3.d0)+ &
             rb**(8.d0/3.d0))+or*ra*(11.d0/3.d0*rb**(8.d0/3.d0)+ &
             ra**(8.d0/3.d0)))
     !
     vlyp(i,1) = de1a + de2a
     vlyp(i,2) = de1b + de2b
     !
  ENDDO
  !==--------------------------------------------------------------==
  !
  RETURN
  !
END SUBROUTINE lsd_lyp
