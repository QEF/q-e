!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
function w1gauss (x, n)
  !-----------------------------------------------------------------------
  !
  !    w1gauss(x,n) = \int_{-\infty}^x   y delta(y) dy
  !    where delta(x) is the current approximation for the delta function,
  !    as obtained from w0gauss(x,n)
  !
  ! --> (n>=0) : Methfessel-Paxton case
  !
  ! --> (n=-1): Cold smearing (Marzari-Vanderbilt)
  !     w1gauss = 1/sqrt(2*pi)*(x-1/sqrt(2))*exp(-(x-1/sqrt(2))**2)
  !
  ! --> (n=-99): Fermi-Dirac case. In this case w1gauss corresponds
  !     to the negative of the electronic entropy.
  !
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  implicit none
  real(DP) :: w1gauss, x
  ! output: the value of the function
  ! input: the point where to compute the function

  integer :: n
  ! input: the order of the smearing function
  !
  !    here the local variables
  !

  real(DP) :: a, hp, arg, hpm1, hd, f, onemf, xp
  ! the coefficients a_n
  ! the hermite function
  ! the argument of the exponential
  ! the hermite function
  ! the hermite function
  ! Fermi-Dirac occupation number
  ! 1 - f
  ! auxiliary variable (cold smearing)

  integer :: i, ni
  ! counter on n values
  ! counter on 2n values

  ! Fermi-Dirac smearing
  if (n.eq. - 99) then
     if (abs (x) .le.36.0) then
        f = 1.0d0 / (1.0d0 + exp ( - x) )
        onemf = 1.0d0 - f
        w1gauss = f * log (f) + onemf * log (onemf)
        ! in order to avoid problems for large values of x
     else
        ! neglect w1gauss when abs(w1gauss) < 1.0d-14
        w1gauss = 0.0d0
     endif
     return

  endif
  ! Cold smearing
  if (n.eq. - 1) then
     xp = x - 1.0d0 / sqrt (2.0d0)
     arg = min (200.d0, xp**2)
     w1gauss = 1.0d0 / sqrt (2.0d0 * pi) * xp * exp ( - arg)
     return

  endif
  ! Methfessel-Paxton
  arg = min (200.d0, x**2)
  w1gauss = - 0.5d0 * exp ( - arg) / sqrt (pi)
  if (n.eq.0) return
  hd = 0.d0
  hp = exp ( - arg)
  ni = 0
  a = 1.d0 / sqrt (pi)
  do i = 1, n
     hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
     ni = ni + 1
     hpm1 = hp
     hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
     ni = ni + 1
     a = - a / (DBLE (i) * 4.0d0)
     w1gauss = w1gauss - a * (0.5d0 * hp + DBLE (ni) * hpm1)
  enddo
  return


end function w1gauss
