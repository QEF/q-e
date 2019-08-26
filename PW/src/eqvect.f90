!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
LOGICAL FUNCTION eqvect( x, y, f, accep )
  !-----------------------------------------------------------------------
  !! This function tests if the difference x-y-f is an integer.
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: x(3)
  !! first 3d vector in crystal axis
  REAL(DP), INTENT(IN) :: y(3)
  !! second 3d vector in crystal axis
  REAL(DP), INTENT(IN) :: f(3)
  !! fractionary translation
  REAL(DP), INTENT(IN) :: accep
  !! threshold of acceptability
  !
  eqvect = ABS( x(1)-y(1)-f(1) - NINT(x(1)-y(1)-f(1)) ) < accep .AND. &
           ABS( x(2)-y(2)-f(2) - NINT(x(2)-y(2)-f(2)) ) < accep .AND. &
           ABS( x(3)-y(3)-f(3) - NINT(x(3)-y(3)-f(3)) ) < accep
  !
  RETURN
  !
END FUNCTION eqvect
