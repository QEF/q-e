!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE setqfnew( nqf, qfcoef, mesh, r, l, n, rho )
  !-----------------------------------------------------------------------
  !
  ! ... Computes the Q function from its polynomial expansion (r < rinner)
  ! ... On input: nqf = number of polynomial coefficients
  ! ...    qfcoef(nqf)= the coefficients defining Q
  ! ...          mesh = number of mesh point
  ! ...        r(mesh)= the radial mesh
  ! ...             l = angular momentum
  ! ...             n = additional exponent, result is multiplied by r^n
  ! ... On output:
  ! ...      rho(mesh)= r^n * Q(r)
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(in):: nqf, l, mesh, n
  REAL(dp), INTENT(in) :: r(mesh), qfcoef(nqf)
  REAL(dp), INTENT(out) :: rho(mesh)
  !
  INTEGER  :: ir, i
  REAL(dp) :: rr
  !
  DO ir = 1, mesh
     rr = r(ir)**2
     rho(ir) = qfcoef(1)
     DO i = 2, nqf
        rho(ir) = rho(ir) + qfcoef(i)*rr**(i-1)
     ENDDO
     rho(ir) = rho(ir)*r(ir)**(l+n)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE setqfnew
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
