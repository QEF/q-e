!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
function d2mxc (rho)
  !-----------------------------------------------------------------------
  !
  !  second derivative of the xc potential with respect to the local density
  !  Perdew and Zunger parameterization of the Ceperley-Alder functional
  !
  USE kinds, ONLY: DP
  implicit none
  !
  real (DP) :: rho, d2mxc
  ! input: the charge density ( positive )
  ! output: the second derivative of the xc potent

  real (DP) :: b1, b2, gc, a, b, c, d, pi, thofpi_3, fpioth_3, &
       thopi_3, tm1, tm2, tm3, tm4, tm5, tm6
  ! parameters defining the functionals
  !
  !    pi
  !    (3/4/pi)^0.333
  !    (4*pi/3)^0.333
  !      (3/pi)^0.333
  !    35.d0*b1,
  !    76.d0*b1*b1 + 64.d0*b2
  !    35.d0*b1*b1*b1 + 234.d0*b1*b2
  !    140.d0*b2*b1*b1 + 176.d0*b2*b2
  !    175.d0*b1*b2*b2
  !    64.d0*b2*b2*b2

  parameter (b1 = 1.0529d0, b2 = 0.3334d0, gc = - 0.1423d0, a = &
       0.0311d0, b = - 0.0480d0, c = 0.0020d0, d = - 0.0116d0, pi = &
       3.14159265358979d0, fpioth_3 = 1.61199195401647d0, thofpi_3 = &
       0.620350490899400d0, thopi_3 = 0.98474502184270d0, tm1 = &
       36.85150d0, tm2 = 105.59107916d0, tm3 = 122.996139546115d0, tm4 = &
       71.30831794516d0, tm5 = 20.4812455967d0, tm6 = 2.371792877056d0)

  real (DP) :: rs, x, den

  rs = thofpi_3 * (1.d0 / rho) **0.3333333333333333d0
  if (rs.ge.1.d0) then
     x = sqrt (rs)
     den = 1.d0 + x * b1 + b2 * x**2
     d2mxc = - gc * (tm1 * x + tm2 * x**2 + tm3 * x**3 + tm4 * x**4 &
          + tm5 * x**5 + tm6 * x**6) / ( (rho**2) * (den**4) * 216.d0)
  else
     d2mxc = (9.d0 * a + (6.d0 * c + 8.d0 * d) * rs + 8.d0 * c * rs &
          * log (rs) ) / (rho**2) / 27.d0

  endif

  rs = rs * fpioth_3
  d2mxc = d2mxc + (2.d0 / 9.d0 * thopi_3 * rs**5)

  d2mxc = 2.d0 * d2mxc
  return
end function d2mxc
