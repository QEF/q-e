!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE volume (alat, a1, a2, a3, omega)
  !---------------------------------------------------------------------
  !
  !     Compute the volume of the unit cell defined by 3 vectors
  !     a1, a2, a3, given in units of "alat" (alat may be 1):
  !        omega = alat^3 * [ a1 . (a2 x a3) ]
  !     ( . = scalar product, x = vector product )
  !
  USE kinds, ONLY: dp
  IMPLICIT NONE
  !
  REAL(dp), INTENT(IN) :: alat, a1(3), a2(3), a3(3)
  REAL(dp), INTENT(OUT) :: omega
  !
  omega = a1(1) * ( a2(2)*a3(3)-a2(3)*a3(2) ) - &
          a1(2) * ( a2(1)*a3(3)-a2(3)*a3(1) ) + &
          a1(3) * ( a2(1)*a3(2)-a2(2)*a3(1) )
  !
  IF ( omega < 0.0_dp) THEN
     call infomsg('volume','axis vectors are left-handed')
     omega = ABS (omega)
  END IF
  !
  IF ( alat < 1.0_dp) call infomsg('volume','strange lattice parameter')
  omega = omega * alat**3
  !
  RETURN
  !
END SUBROUTINE volume
