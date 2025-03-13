!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE suscept_g0(n, g, x, dx, ddx)
  !---------------------------------------------------------------------------
  !
  ! ... calculate 1st and 2nd derivative of susceptibility at G = 0.
  !
  USE constants, ONLY : eps12
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: n
  REAL(DP), INTENT(IN)  :: g(1:*)
  REAL(DP), INTENT(IN)  :: x(1:*)
  REAL(DP), INTENT(OUT) :: dx
  REAL(DP), INTENT(OUT) :: ddx
  !
  REAL(DP) :: denom
  REAL(DP) :: g2, g3
  REAL(DP) :: x2, x3
  !
  IF (n > 2) THEN
    g2    = g(2) - g(1)
    g3    = g(3) - g(1)
    x2    = x(2) - x(1)
    x3    = x(3) - x(1)
    denom = g2 * g3 * (g3 - g2)
    IF (ABS(denom) > eps12) THEN
      dx  = -(g2*g2*x3 - g3*g3*x2) / denom
      ddx = 2.0_DP * (g2*x3 - g3*x2) / denom
    ELSE
      dx  = 0.0_DP
      ddx = 0.0_DP
    END IF
    !
  ELSE IF (n > 1) THEN
    denom = g(2) - g(1)
    IF (ABS(denom) > eps12) THEN
      dx  = (x(2) - x(1)) / denom
      ddx = 0.0_DP
    ELSE
      dx  = 0.0_DP
      ddx = 0.0_DP
    END IF
    !
  ELSE
    dx  = 0.0_DP
    ddx = 0.0_DP
  END IF
  !
END SUBROUTINE suscept_g0
