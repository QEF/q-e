!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine h_1psi (lda, n, psi, hpsi, spsi)
  !-----------------------------------------------------------------------
  !
  !    This routine applies the Hamiltonian and the S matrix
  !    to a vector psi and puts the result in hpsi and spsi
  !    Wrapper routine - calls h_psi
  use parameters, only: DP
  implicit none
  !
  integer :: lda, n
  complex(kind=DP) :: psi (n), hpsi (n), spsi(n)
  !
  call start_clock ('h_1psi')
  call h_psi (lda, n, 1, psi, hpsi, spsi)
  call stop_clock ('h_1psi')
  return
end subroutine h_1psi

