!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------
SUBROUTINE setdqf( nqf, qfcoef, mesh, r, l, drho )
  !-----------------------------------------------------------------------
  !! Computes the derivative of the Q function, \(dQ/dr\), from its 
  !! polynomial expansion (valid for r < rinner).
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: nqf
  !! number of polynomial coefficients
  INTEGER, INTENT(in) :: l
  !! angular momentum
  INTEGER, INTENT(in) :: mesh
  !! number of mesh point
  REAL(dp), INTENT(in) :: r(mesh)
  !! the radial mesh
  REAL(dp), INTENT(in) :: qfcoef(nqf)
  !! the coefficients defining Q
  REAL(dp), INTENT(out) :: drho(mesh)
  !! \(dQ(r)/dr\)
  !
  INTEGER  :: ir, i
  !
  DO ir = 1, mesh
     !
     drho(ir) = 0.0_dp
     DO i = max( 1, 2-l ), nqf
        drho(ir) = drho(ir) + qfcoef(i)*r(ir)**(2*i-3+l)*(2*i-2+l)
     ENDDO
     !
  END DO
  !
END SUBROUTINE setdqf
