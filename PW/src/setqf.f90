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
  ! ... This routine computes the first part of the Q function up to rinner.
  ! ... On input: nqf = number of coefficients
  ! ...    qfcoef(nqf)= the coefficients defining Q
  ! ...          mesh = number of mesh point
  ! ...        r(mesh)= the radial mesh
  ! ...          ltot = angular momentum
  ! ...             n = exponent
  ! ... On output:
  ! ...      rho(mesh)= Q(r)*r^n
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
