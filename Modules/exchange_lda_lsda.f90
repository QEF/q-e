!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


MODULE exch_lda  !<GPU:exch_lda=>exch_lda_gpu>
!
!! Module containing LDA functionals.
!
CONTAINS
!-----------------------------------------------------------------------
!
!  ... LSDA
!
!-----------------------------------------------------------------------
SUBROUTINE slater_spin( rho, zeta, ex, vx_up, vx_dw )                 !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! Slater exchange with alpha=2/3, spin-polarized case.
  !
  USE kinds, ONLY : DP
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
END SUBROUTINE slater_spin
!
!
!-----------------------------------------------------------------------
SUBROUTINE slater_rxc_spin( rho, z, ex, vx_up, vx_dw )                 !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! Slater exchange with alpha=2/3, relativistic exchange case.
  !! Spin-polarized case.
  !
  USE kinds,       ONLY: DP
  !
  IMPLICIT NONE
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
SUBROUTINE slater1_spin( rho, zeta, ex, vx_up, vx_dw )                 !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !     Slater exchange with alpha=2/3, spin-polarized case
  !
  USE kinds, ONLY: DP
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

END MODULE
