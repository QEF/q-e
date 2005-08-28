!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
function dmxc (rho)
  !-----------------------------------------------------------------------
  !
  !  derivative of the xc potential with respect to the local density
  !
  USE kinds, only : DP
  use funct
  implicit none
  ! I/O variables

  real(DP) :: rho, dmxc
  ! input: the charge density ( positive )
  ! output: the derivative of the xc potential


  real(DP) :: dr, vxp, vcp, vxm, vcm, ex, ec
  ! delta rrho for numerical derivatives
  ! the potentials for + charge
  ! the potentials for - charge
  ! the energy
  ! DFT functional
  ! auxiliary variables
  real(DP) :: vx, rs, dpz


  integer :: iflg
  ! parameters
  real(DP) :: small, e2, pi34, third
  parameter (small = 1.d-30, e2 = 2.d0)

  parameter (pi34 = 0.75d0 / 3.141592653589793d+00, third = 1.d0 / &
       3.d0)
  dmxc = 0.d0
  if (rho.le.small) then
     return

  endif
  !
  !    first case: analytical derivatives available
  !
  if (iexch.eq.1.and.icorr.eq.1) then
     rs = (pi34 / rho) **third
     !..exchange
     call slater (rs, ex, vx)
     dmxc = vx / (3.d0 * rho)
     !..correlation
     iflg = 2
     if (rs.lt.1.0d0) iflg = 1
     dmxc = dmxc + dpz (rs, iflg)
  else
     !
     !     second case: numerical derivatives
     !
     dr = min (1.d-6, 1.d-4 * rho)
     call xc (rho + dr, ex, ec, vxp, vcp)
     call xc (rho - dr, ex, ec, vxm, vcm)
     dmxc = (vxp + vcp - vxm - vcm) / (2.d0 * dr)
  endif
  !
  ! scales to rydberg units
  !

  dmxc = e2 * dmxc
  return

end function dmxc
!-----------------------------------------------------------------------
function dpz (rs, iflg)
!-----------------------------------------------------------------------
!  derivative of the correlation potential with respect to the local den
!  Perdew and Zunger parameterization of the C.A. functional
!
USE kinds, only : DP
implicit none
real(DP) :: rs, dpz
                             ! input : the value of rs
                             ! output: the derivative of the corr. poten

integer :: iflg
                             ! input : flag to choose the functional for

real(DP) :: b1, b2, a1, a2, gc, a, b, c, d, pi, fpi
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
parameter (a = 0.0311d0, b = - 0.048d0, c = 0.0020d0, d = - &
 0.0116d0, gc = - 0.1423d0, b1 = 1.0529d0, b2 = 0.3334d0, a1 = &
 7.0d0 * b1 / 6.d0, a2 = 4.d0 * b2 / 3.d0, pi = 3.14159265358979d0, &
 fpi = 4.d0 * pi)

real(DP) :: x, den, dmx, dmrs
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

dpz = - fpi * rs**4.d0 / 9.d0 * dmrs
return
end function dpz
