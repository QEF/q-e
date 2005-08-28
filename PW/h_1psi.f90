!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE h_1psi( lda, n, psi, hpsi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the Hamiltonian and the S matrix
  ! ... to a vector psi and puts the result in hpsi and spsi
  ! ... Wrapper routine - calls h_psi and s_psi
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER           :: lda, n
  COMPLEX (DP) :: psi(n), hpsi(n), spsi(n)
  !
  !
  CALL start_clock( 'h_1psi' )
  !
  CALL h_psi( lda, n, 1, psi, hpsi )
  CALL s_psi( lda, n, 1, psi, spsi )
  !
  CALL stop_clock( 'h_1psi' )
  !
  RETURN
  !
END SUBROUTINE h_1psi
