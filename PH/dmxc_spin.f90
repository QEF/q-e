!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dmxc_spin (rhoup, rhodw, dmuxc_uu, dmuxc_ud, dmuxc_du, &
     dmuxc_dd)
  !-----------------------------------------------------------------------
  !  derivative of the xc potential with respect to the local density
  !
  !
  USE kinds, only : DP
  use funct
  implicit none
  ! I/O variables


  real(kind=DP) :: rhoup, rhodw, dmuxc_uu, dmuxc_ud, dmuxc_du, dmuxc_dd
  ! input: the up charge density
  ! input: the dw charge density
  ! output: the up-up derivative of the xc en
  ! output: the up-dw derivative of the xc en
  ! output: the dw-up derivative of the xc en
  ! output: the dw-dw derivative of the xc en
  ! DFT functional

  !-auxiliary variables

  real(kind=DP) :: rhotot, rs, zeta, fz, fz1, fz2, ex, vx, ecu, ecp, vcu, &
       vcp, dmcu, dmcp, dpz, dpz_polarized, aa, bb, cc
  real(kind=DP) :: dr, dz, ec, vxupm, vxdwm, vcupm, vcdwm, rho, vxupp, &
       vxdwp, vcupp, vcdwp


  integer :: iflg
  !-parameters
  real(kind=DP) :: small, e2, pi34, third, p43, p49, m23
  parameter (small = 1.d-30, e2 = 2.d0)

  parameter (pi34 = 0.75d0 / 3.141592653589793d+00, third = 1.d0 / &
       3.d0, p43 = 4.d0 / 3.d0, p49 = 4.d0 / 9.d0, m23 = - 2.d0 / 3.d0)
  dmuxc_uu = 0.d0
  dmuxc_du = 0.d0
  dmuxc_ud = 0.d0

  dmuxc_dd = 0.d0
  rhotot = rhoup + rhodw
  if (rhotot.le.small) return
  zeta = (rhoup - rhodw) / rhotot

  if (abs (zeta) .gt.1.d0) return
  if (iexch.eq.1.and.icorr.eq.1) then
     !
     !    first case: analytical derivative available
     !
     !..exchange
     rs = (pi34 / (2.d0 * rhoup) ) **third
     call slater (rs, ex, vx)
     dmuxc_uu = vx / (3.d0 * rhoup)
     rs = (pi34 / (2.d0 * rhodw) ) **third
     call slater (rs, ex, vx)
     dmuxc_dd = vx / (3.d0 * rhodw)
     !..correlation
     rs = (pi34 / rhotot) **third
     iflg = 2
     if (rs.lt.1.0d0) iflg = 1
     dmcu = dpz (rs, iflg)
     dmcp = dpz_polarized (rs, iflg)
     call pz (rs, 1, ecu, vcu)

     call pz_polarized (rs, ecp, vcp)
     fz = ( (1.d0 + zeta) **p43 + (1.d0 - zeta) **p43 - 2.d0) &
          / (2.d0**p43 - 2.d0)
     fz1 = p43 * ( (1.d0 + zeta) **third- (1.d0 - zeta) **third) &
          / (2.d0**p43 - 2.d0)

     fz2 = p49 * ( (1.d0 + zeta) **m23 + (1.d0 - zeta) **m23) &
          / (2.d0**p43 - 2.d0)
     aa = dmcu + fz * (dmcp - dmcu)
     bb = 2.d0 * fz1 * (vcp - vcu - (ecp - ecu) ) / rhotot

     cc = fz2 * (ecp - ecu) / rhotot
     dmuxc_uu = dmuxc_uu + aa + (1.d0 - zeta) * bb + (1.d0 - zeta) &
          **2 * cc
     dmuxc_du = dmuxc_du + aa + ( - zeta) * bb + (zeta**2 - 1.d0) &
          * cc
     dmuxc_ud = dmuxc_du

     dmuxc_dd = dmuxc_dd+aa - (1.d0 + zeta) * bb + (1.d0 + zeta) ** &
          2 * cc

  else
     rho = rhoup + rhodw
     dr = min (1.d-6, 1.d-4 * rho)
     call xc_spin (rho - dr, zeta, ex, ec, vxupm, vxdwm, vcupm, &
          vcdwm)

     call xc_spin (rho + dr, zeta, ex, ec, vxupp, vxdwp, vcupp, &
          vcdwp)
     dmuxc_uu = (vxupp + vcupp - vxupm - vcupm) / (2.d0 * dr)
     dmuxc_ud = dmuxc_uu
     dmuxc_dd = (vxdwp + vcdwp - vxdwm - vcdwm) / (2.d0 * dr)

     dmuxc_du = dmuxc_dd
     dz = min (1.d-6, 1.d-4 * abs (zeta) )
     call xc_spin (rho, zeta - dz, ex, ec, vxupm, vxdwm, vcupm, &
          vcdwm)

     call xc_spin (rho, zeta + dz, ex, ec, vxupp, vxdwp, vcupp, &
          vcdwp)
     dmuxc_uu = dmuxc_uu + (vxupp + vcupp - vxupm - vcupm) * &
          (1.d0 - zeta) / rho / (2.d0 * dz)
     dmuxc_ud = dmuxc_ud- (vxupp + vcupp - vxupm - vcupm) * (1.d0 + &
          zeta) / rho / (2.d0 * dz)
     dmuxc_du = dmuxc_du + (vxdwp + vcdwp - vxdwm - vcdwm) * &
          (1.d0 - zeta) / rho / (2.d0 * dz)
     dmuxc_dd = dmuxc_dd- (vxdwp + vcdwp - vxdwm - vcdwm) * (1.d0 + &
          zeta) / rho / (2.d0 * dz)
  endif
  !
  ! scales to rydberg units
  !
  dmuxc_uu = e2 * dmuxc_uu
  dmuxc_du = e2 * dmuxc_du
  dmuxc_ud = e2 * dmuxc_ud
  dmuxc_dd = e2 * dmuxc_dd
  !
  return

end subroutine dmxc_spin
!-----------------------------------------------------------------------
function dpz_polarized (rs, iflg)
!-----------------------------------------------------------------------
!  derivative of the correlation potential with respect to the local den
!  Perdew and Zunger parameterization of the C.A. functional
!
USE kinds, only : DP
implicit none
real(kind=DP) :: rs, dpz_polarized
                             ! input : the value of rs
                             ! output: the derivative of the corr. poten

integer :: iflg
                             ! input : flag to choose the functional for

real(kind=DP) :: b1, b2, a1, a2, gc, a, b, c, d, pi, fpi
                             !\
                             ! \
                             !  \
                             !   \
                             !    parameter which define the functional
                             !
                             !
                             !
                             !  /
                             ! /
                             !/
parameter (a = 0.01555d0, b = - 0.0269d0, c = 0.0007d0, d = &
 - 0.0048d0, gc = - 0.0843d0, b1 = 1.3981d0, b2 = 0.2611d0, a1 = &
 7.0d0 * b1 / 6.d0, a2 = 4.d0 * b2 / 3.d0, pi = 3.14159265358979d0, &
 fpi = 4.d0 * pi)

real(kind=DP) :: x, den, dmx, dmrs
                              ! auxiliary variable
                              ! auxiliary variable
                              ! auxiliary variable
                              ! auxiliary variable
if (iflg.eq.1) then
   dmrs = a / rs + 2.d0 / 3.d0 * c * (log (rs) + 1.d0) + (2.d0 * &
    d-c) / 3.d0
else
   x = sqrt (rs)
   den = 1.d0 + x * (b1 + x * b2)
   dmx = gc * ( (a1 + 2.d0 * a2 * x) * den - 2.d0 * (b1 + 2.d0 * &
    b2 * x) * (1.d0 + x * (a1 + x * a2) ) ) / den**3
   dmrs = 0.5d0 * dmx / x
endif
!

dpz_polarized = - fpi * rs**4.d0 / 9.d0 * dmrs
return
end function dpz_polarized
