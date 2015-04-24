!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE setqfnew( nqf, qfcoef, mesh, r, ltot, n, rho )
  !-----------------------------------------------------------------------
  !
  ! ... Computes the Q function from its polynomial expansion (r < rinner)
  ! ... On input: nqf = number of polynomial coefficients
  ! ...    qfcoef(nqf)= the coefficients defining Q
  ! ...          mesh = number of mesh point
  ! ...        r(mesh)= the radial mesh
  ! ...          ltot = angular momentum
  ! ...             n = additional exponent, result is multiplied by r^n
  ! ... On output:
  ! ...      rho(mesh)= r^n * Q(r)
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(in):: nqf, ltot, mesh, n
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
     rho(ir) = rho(ir)*r(ir)**(ltot+n)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE setqfnew
!
!------------------------------------------------------------------------
SUBROUTINE setdqf( nqf, qfcoef, mesh, r, ltot, drho )
  !-----------------------------------------------------------------------
  !
  ! ... Computes the derivative of the Q function, dQ/dr,
  ! ... from its polynomial expansion (valid for r < rinner)
  ! ... On input: nqf = number of polynomial coefficients
  ! ...    qfcoef(nqf)= the coefficients defining Q
  ! ...          mesh = number of mesh point
  ! ...        r(mesh)= the radial mesh
  ! ...          ltot = angular momentum
  ! ... On output:
  ! ...      drho(mesh)= r^n * dQ(r)/dr
  !
  USE kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(in):: nqf, ltot, mesh
  REAL(dp), INTENT(in) :: r(mesh), qfcoef(nqf)
  REAL(dp), INTENT(out) :: drho(mesh)
  !
  INTEGER  :: ir, i
  REAL(dp) :: rr
  !
  DO ir = 1, mesh
     !
     rr = r(ir)*r(ir)
     drho(ir) = 0.D0
     DO i = max( 1, 2-ltot ), nqf
        drho(ir) = drho(ir) + qfcoef(i)*rr**(i-2+ltot)*(i-1+ltot)
     ENDDO
     !
  END DO
  !
END SUBROUTINE setdqf
