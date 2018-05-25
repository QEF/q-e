!
! Copyright (C) 2004-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE splinelib_gpum
  !---------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: splint_eq_gpu
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    ! This subroutine is equivalent to calling splint for a list of x values with
    ! an equispaced xdata interpolation grid with spacing dq, xdata = [0, dq, 2dq....]
    !------------------------------------------------------------------------
    SUBROUTINE splint_eq_gpu( dq, ydata, d2y, xlist, s)
      !------------------------------------------------------------------------
      !
#if defined(__CUDA)
      USE cudafor
#endif
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: ydata(:), d2y(:), xlist(:)
      REAL(DP), INTENT(OUT) :: s(:)
      REAL(DP), INTENT(IN) :: dq
      !
      INTEGER  :: xdim, n, ig, khi, klo
      REAL(DP) :: a, b, h, xlo, xhi, x
      !
#if defined(__CUDA)
      attributes(device) :: ydata, d2y, xlist, s
#endif
      !
      xdim = size(ydata)
      n  = size(xlist)

      !$cuf kernel do(1) <<<*,*>>>
      do ig = 1, n
        !
        x = xlist(ig)
        klo = MAX( MIN( int(x/dq) + 1, ( xdim - 1 ) ), 1 )
        !
        khi = klo + 1
        !
        xlo = (klo - 1) * dq
        xhi = (khi - 1) * dq
        h = xhi - xlo
        !
        a = ( xhi - x ) / h
        b = ( x - xlo ) / h
        !
        s(ig) = a * ydata(klo) + b * ydata(khi) + &
                 ( ( a*a*a - a ) * d2y(klo) + ( b*b*b - b ) * d2y(khi) ) * &
                 ( h*h ) / 6.0_DP
      end do

      END SUBROUTINE splint_eq_gpu

END MODULE splinelib_gpum
