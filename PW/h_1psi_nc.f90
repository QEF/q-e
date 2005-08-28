!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine h_1psi_nc (lda, n, npol, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !    This routine applies the Hamiltonian and the S matrix
  !    to a vector psi and puts the result in hpsi and spsi
  !    Wrapper routine - calls h_psi
  USE kinds, only: DP
  implicit none
  !
  integer :: lda, n, npol
  complex(DP) :: psi (lda,npol), hpsi (lda,npol), spsi(lda,npol)
  !
  call start_clock ('h_1psi')
  call h_psi_nc (lda, n, 1, psi, hpsi)
  call s_psi_nc (lda, n, 1, psi, spsi)
  call stop_clock ('h_1psi')
  return
end subroutine h_1psi_nc

