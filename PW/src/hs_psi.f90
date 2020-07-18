!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE hs_psi( lda, n, m, psi, hpsi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the Hamiltonian and the S matrix
  ! ... to a vector psi and puts the result in hpsi and spsi
  ! ... Wrapper routine - calls h_psi_ and s_psi_
  !
  ! ... No bgrp parallelization here !
  !
  USE kinds,  ONLY: DP
  USE noncollin_module, ONLY: npol 
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX (DP) :: psi(lda*npol, m), hpsi(lda*npol,m), spsi(lda*npol,m)
  !
  CALL start_clock( 'hs_psi' )
  ! 
  CALL h_psi_ ( lda, n, m, psi, hpsi ) ! apply H to m wfcs (no bgrp parallelization here)
  CALL s_psi_ ( lda, n, m, psi, spsi ) ! apply S to m wfcs (no bgrp parallelization here)
  !
  CALL stop_clock( 'hs_psi' )
  !
  RETURN
  !
END SUBROUTINE hs_psi
