!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine s_1psi (lda, n, psi, spsi)
  !-----------------------------------------------------------------------
  !
  !    spsi = S*psi for one wavefunction
  !    Wrapper routine - calls ccalbec and s_psi
  !
  use parameters, only: DP
  use becmod, only: becp
  use wvfct, only: npwx
  use us, only: vkb, nkb
  implicit none
  !
  integer :: lda, n
  complex(kind=DP) :: psi (n), spsi (n)
  !
  call start_clock ('s_1psi')
  call ccalbec (nkb, npwx, n, 1, becp, vkb, psi)
  call s_psi (lda, n, 1, psi, spsi)
  call stop_clock ('s_1psi')
  return
end subroutine s_1psi
