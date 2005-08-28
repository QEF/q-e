!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE s_1psi( lda, n, psi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... spsi = S*psi for one wavefunction
  ! ... Wrapper routine - calls ccalbec and s_psi
  !
  USE kinds,  ONLY : DP
  USE wvfct,  ONLY : npwx
  USE uspp,   ONLY : vkb, nkb
  USE wvfct,  ONLY : gamma_only
  USE becmod, ONLY : becp, rbecp
  !
  IMPLICIT NONE
  !
  INTEGER          :: lda, n
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
     CALL ccalbec( nkb, npwx, n, 1, becp, vkb, psi )
     !
  END IF 
  !
  CALL s_psi( lda, n, 1, psi, spsi )
  !
  CALL stop_clock( 's_1psi' )
  !
  RETURN
  !
END SUBROUTINE s_1psi
