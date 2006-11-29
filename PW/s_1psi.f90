!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE s_1psi( npwx, n, psi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... spsi = S*psi for one wavefunction
  ! ... Wrapper routine - calls ccalbec and s_psi
  !
  USE kinds,  ONLY : DP
  USE uspp,   ONLY : vkb, nkb
  USE wvfct,  ONLY : gamma_only
  USE becmod, ONLY : becp, rbecp, becp_nc
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER          :: npwx, n
  COMPLEX(DP) :: psi(n), spsi(n)
  !
  !
  CALL start_clock( 's_1psi' )
  !
  IF ( gamma_only ) THEN
     !
     CALL ccalbec( nkb, npwx, n, 1, rbecp, vkb, psi )
     !
  ELSE
     !
     IF ( noncolin ) THEN
        !
        CALL ccalbec_nc( nkb, npwx, n, npol, 1, becp_nc, vkb, psi )
        !
     ELSE
        !
        CALL ccalbec( nkb, npwx, n, 1, becp, vkb, psi )
        !
     END IF
     !
  END IF 
  !
  IF ( noncolin ) THEN
     !
     CALL s_psi_nc( npwx, n, 1, psi, spsi )
     !
  ELSE
     !
     CALL s_psi( npwx, n, 1, psi, spsi )
     !
  END IF
  !
  CALL stop_clock( 's_1psi' )
  !
  RETURN
  !
END SUBROUTINE s_1psi
