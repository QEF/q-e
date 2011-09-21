!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
function w_1gauss (x, n)
!-----------------------------------------------------------------------
!
!     the derivative of w0gauss:
!
! --> (n=-99): second derivative of Fermi-Dirac function
!
USE kinds, ONLY : DP
USE constants, ONLY : sqrtpm1
!
implicit none
real (DP) :: w_1gauss, x
                       ! output: the value of the function
                      ! input: the point where to compute the function

integer :: n
                      ! input: the order of the smearing function
!
!    here the local variables
!
real (DP) :: a, arg, hp, hd, aux1, aux2
                      ! the coefficients a_n
                      ! the argument of the exponential
                      ! the hermite function
                      ! the hermite function
                      ! auxiliary variable
                      ! auxiliary variable
integer :: i, ni
                      ! counter on n values
                      ! counter on 2n values
! Fermi-Dirac smearing
if (n.eq. - 99) then
   aux1 = exp (x)
   aux2 = exp ( - x)
   w_1gauss = (aux2 - aux1) / (2.d0 + aux1 + aux2) **2
   return
endif
!
arg = min (200.d0, x**2)
w_1gauss = - 2.d0 * x * exp ( - arg) * sqrtpm1
if (n.eq.0) return
hd = exp ( - arg)
hp = 2.d0 * x * exp ( - arg)
ni = 1
a = sqrtpm1
do i = 1, n
hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
ni = ni + 1
a = - a / (DBLE (i) * 4.0d0)
hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
ni = ni + 1
w_1gauss = w_1gauss - a * hp
enddo
return



end function w_1gauss
