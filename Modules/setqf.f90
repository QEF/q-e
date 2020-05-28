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
  !
  ! ... Computes the derivative of the Q function, dQ/dr,
  ! ... from its polynomial expansion (valid for r < rinner)
  ! ... On input: nqf = number of polynomial coefficients
  ! ...    qfcoef(nqf)= the coefficients defining Q
  ! ...          mesh = number of mesh point
  ! ...        r(mesh)= the radial mesh
  ! ...             l = angular momentum
  ! ... On output:
  ! ...      drho(mesh)= dQ(r)/dr
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(in):: nqf, l, mesh
  REAL(dp), INTENT(in) :: r(mesh), qfcoef(nqf)
  REAL(dp), INTENT(out) :: drho(mesh)
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
