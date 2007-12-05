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
  ! ... Wrapper routine - calls calbec and s_psi
  !
  USE kinds,  ONLY : DP
  USE uspp,   ONLY : vkb, nkb
  USE becmod, ONLY : becp, rbecp, becp_nc, calbec
  USE control_flags,    ONLY : gamma_only 
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER          :: npwx, n
  COMPLEX(DP) :: psi(npwx*npol,1), spsi(npwx*npol)
  !
  !
  CALL start_clock( 's_1psi' )
  !
  IF ( gamma_only ) THEN
     !
     CALL calbec( n, vkb, psi, rbecp )
     CALL s_psi( npwx, n, 1, psi, spsi )
     !
  ELSE IF ( noncolin ) THEN
     !
     CALL calbec( n, vkb, psi, becp_nc )
     CALL s_psi_nc( npwx, n, 1, psi, spsi )
     !
  ELSE
     !
     CALL calbec( n, vkb, psi, becp )
     CALL s_psi( npwx, n, 1, psi, spsi )
     !
  END IF
  !
  CALL stop_clock( 's_1psi' )
  !
  RETURN
  !
END SUBROUTINE s_1psi
