!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE fcp_hessian(hess)
  !----------------------------------------------------------------------------
  !
  ! ... calculate inverse of the Hessian:
  ! ...
  ! ...     d^2E/dN^2 = d(ef)/dN = 1/DOS(ef)
  !
  USE ener,  ONLY : ef
  USE kinds, ONLY : DP
  USE klist, ONLY : nkstot, wk, degauss, ngauss
  USE wvfct, ONLY : nbnd, et
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: hess
  !
  INTEGER :: ik
  INTEGER :: ibnd
  !
  REAL(DP), EXTERNAL :: w0gauss
  !
  hess = 0.0_DP
  !
  DO ik = 1, nkstot
     !
     DO ibnd = 1, nbnd
        !
        hess = hess + wk (ik) * &
             & w0gauss((ef - et(ibnd, ik)) / degauss, ngauss) / degauss
        !
     END DO
     !
  END DO
  !
END SUBROUTINE fcp_hessian
