!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
MODULE exch_lda
!----------------------------------------------------------------------
!! LDA exchange functionals.
!
CONTAINS
!
!-----------------------------------------------------------------------
SUBROUTINE slater( rs, ex, vx )
  !---------------------------------------------------------------------
  !! Slater exchange with alpha=2/3
  !
  USE kind_l,  ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
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
END SUBROUTINE slater
!
!
!-----------------------------------------------------------------------
SUBROUTINE slater1( rs, ex, vx )
  !---------------------------------------------------------------------
  !! Slater exchange with alpha=1, corresponding to -1.374/r_s Ry.
  !! Used to recover old results.
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(OUT) :: ex
  !! Exchange energy (per unit volume)
  REAL(DP), INTENT(OUT) :: vx
  !! Exchange potential
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER   :: f = -0.687247939924714d0, alpha = 1.0d0
  !
  ex = f * alpha / rs
  vx = 4.d0 / 3.d0 * f * alpha / rs
  !
  RETURN
  !
END SUBROUTINE slater1
!
!
!-----------------------------------------------------------------------
SUBROUTINE slater_rxc( rs, ex, vx )
  !---------------------------------------------------------------------
  !! Slater exchange with alpha=2/3 and Relativistic exchange.
  !
  USE kind_l,      ONLY: DP
  USE constants_l,  ONLY: pi, c_au
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(OUT) :: ex
  !! Exchange energy (per unit volume)
  REAL(DP), INTENT(OUT) :: vx
  !! Exchange potential
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER   :: zero=0.d0, one=1.d0, pfive=0.5d0, &
                           opf=1.5d0 !, C014=0.014D0
  REAL(DP) :: trd, ftrd, tftm, a0, alp, z, fz, fzp, vxp, xp, &
              beta, sb, alb, c014
  !
  trd  = one/3.d0
  ftrd = 4.d0*trd
  tftm = 2**ftrd-2.d0
  A0   = (4.d0/(9.d0*PI))**trd
  C014 = 1.0_DP/a0/c_au
  ! --- X-alpha PARAMETER:
  alp = 2.d0 * trd
  !
  z  = zero
  fz = zero
  fzp= zero
  !
  vxp = -3.d0*alp/( 2.d0*PI*A0*rs )
  xp  = 3.d0*vxp/4.d0
  beta= C014 / rs
  sb  = SQRT(1.d0+beta*beta)
  alb = LOG(beta+sb)
  vxp = vxp * ( -pfive + opf * alb / (beta*sb) )
  xp  = xp * ( one-opf*((beta*sb-alb)/beta**2)**2 )
  !  vxf = 2**trd*vxp
  !  exf = 2**trd*xp
  vx = vxp
  ex = xp
  !
  RETURN
  !
END SUBROUTINE slater_rxc
!
!
!-----------------------------------------------------------------------
SUBROUTINE slaterKZK( rs, ex, vx, vol )
  !---------------------------------------------------------------------
  !! Slater exchange with alpha=2/3, Kwee, Zhang and Krakauer KE
  !! correction.
  !
  USE kind_l,      ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rs
  !! Wigner-Seitz radius
  REAL(DP), INTENT(OUT) :: ex
  !! Exchange energy (per unit volume)
  REAL(DP), INTENT(OUT) :: vx
  !! Exchange potential
  REAL(DP) :: vol
  !! Finite size volume element
  !
  ! ... local variables
  !
  REAL(DP) :: dL, ga, pi, a0
  REAL(DP), PARAMETER :: a1 = -2.2037d0, &
                         a2 = 0.4710d0, a3 = -0.015d0, ry2h = 0.5d0
  REAL(DP), PARAMETER :: f = -0.687247939924714d0, alpha = 2.0d0/3.0d0
  !                      f = -9/8*(3/2pi)^(2/3)
  !
  pi = 4.d0 * ATAN(1.d0)
  a0 = f * alpha * 2.d0
  !
  dL = vol**(1.d0/3.d0)
  ga = 0.5d0 * dL *(3.d0 /pi)**(1.d0/3.d0)
  !
  IF ( rs < ga ) THEN
     ex = a0 / rs + a1 * rs / dL**2.d0 + a2 * rs**2.d0 / dL**3.d0
     vx = (4.d0 * a0 / rs + 2.d0 * a1 * rs / dL**2.d0 + &
              a2 * rs**2.d0 / dL**3.d0 ) / 3.d0
  ELSE
     ex = a0 / ga + a1 * ga / dL**2.d0 + a2 * ga**2.d0 / dL**3.d0 ! solids
     vx = ex
     ! ex = a3 * dL**5.d0 / rs**6.d0                           ! molecules
     ! vx = 3.d0 * ex
  ENDIF
  !
  ex = ry2h * ex    ! Ry to Hartree
  vx = ry2h * vx
  !
  RETURN
  !
END SUBROUTINE slaterKZK
!
!
!  ... LSDA
!
!-----------------------------------------------------------------------
SUBROUTINE slater_spin( rho, zeta, ex, vx_up, vx_dw )
  !-----------------------------------------------------------------------
  !! Slater exchange with alpha=2/3, spin-polarized case.
  !
  USE kind_l, ONLY : DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
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
END SUBROUTINE slater_spin
!
!
!-----------------------------------------------------------------------
SUBROUTINE slater_rxc_spin( rho, z, ex, vx_up, vx_dw )
  !-----------------------------------------------------------------------
  !! Slater exchange with alpha=2/3, relativistic exchange case.
  !! Spin-polarized case.
  !
  USE kind_l,       ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
  !
  REAL(DP), INTENT(IN) :: rho
  !! total charge density
  REAL(DP), INTENT(IN) :: z
  !! z = (rho_up - rho_dw) / rho_tot
  REAL(DP), INTENT(OUT) :: ex
  !! exchange energy
  REAL(DP), INTENT(OUT) :: vx_up, vx_dw
  !! exchange potential (up, down)
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: zero=0.D0, one=1.D0, pfive=.5D0, &
                         opf=1.5D0, C014=0.014D0, &
                         pi = 3.14159265358979323846D0
  REAL(DP) :: rs, trd, ftrd, tftm, a0, alp, fz, fzp, vxp, xp, &
              beta, sb, alb, vxf, exf
  !
  trd = one/3.d0
  ftrd = 4.d0*trd
  tftm = 2**ftrd-2.d0
  A0 = (4.d0/(9.d0*PI))**trd
  !
  !      X-alpha PARAMETER:
  alp = 2.d0 * trd
  !
  !
  IF ( rho <=  zero ) THEN
     ex = zero
     vx_up = zero ; vx_dw=zero
     RETURN
  ELSE
     fz = ((1.d0+z)**ftrd+(1.d0-Z)**ftrd-2.d0)/tftm
     fzp = ftrd*((1.d0+Z)**trd-(1.d0-Z)**trd)/tftm
  ENDIF
  !
  RS = (3.d0 / (4.d0*PI*rho) )**trd
  vxp = -3.d0*alp/(2.d0*PI*A0*RS)
  XP = 3.d0*vxp/4.d0
  !
  beta = C014/RS
  SB = SQRT(1.d0+beta*beta)
  alb = LOG(beta+SB)
  vxp = vxp * (-pfive + opf * alb / (beta*SB))
  xp = xp * (one-opf*((beta*SB-alb)/beta**2)**2)

  vxf = 2.d0**trd*vxp
  exf = 2.d0**trd*xp
  vx_up = vxp + fz*(vxf-vxp) + (1.d0-z)*fzp*(exf-xp)
  vx_dw = vxp + fz*(vxf-vxp) - (1.d0+z)*fzp*(exf-xp)
  ex    = xp  + fz*(exf-xp)
  !
  RETURN
  !
END SUBROUTINE slater_rxc_spin
!
!
!-----------------------------------------------------------------------
SUBROUTINE slater1_spin( rho, zeta, ex, vx_up, vx_dw )
  !-----------------------------------------------------------------------
  !! Slater exchange with alpha=2/3, spin-polarized case
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  !$acc routine seq
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
  REAL(DP), PARAMETER :: f = - 1.10783814957303361d0, alpha = 1.0d0, &
                         third = 1.d0 / 3.d0, p43 = 4.d0 / 3.d0
  !                      f = -9/8*(3/pi)^(1/3)
  REAL(DP) :: exup, exdw, rho13
  !
  !
  rho13 = ( (1.d0 + zeta) * rho)**third
  exup = f * alpha * rho13
  vx_up = p43 * f * alpha * rho13
  !
  rho13 = ( (1.d0 - zeta) * rho)**third
  exdw = f * alpha * rho13
  vx_dw = p43 * f * alpha * rho13
  !
  ex = 0.5d0 * ( (1.d0 + zeta) * exup + (1.d0 - zeta) * exdw)
  !
  !
  RETURN
  !
END SUBROUTINE slater1_spin
!
!
END MODULE exch_lda
