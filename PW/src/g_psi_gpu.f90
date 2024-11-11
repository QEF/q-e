!
! Copyright (C) 2001-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define TEST_NEW_PRECONDITIONING
!
!-----------------------------------------------------------------------
subroutine g_1psi (lda, n, psi, e)
  !-----------------------------------------------------------------------
  !
  !    As g_psi, for a single wavefunction
  !
  USE kinds,     ONLY: dp
  USE g_psi_mod, ONLY : h_diag, s_diag
  USE noncollin_module, ONLY : npol

  implicit none

  integer, intent(in) :: lda
  !! leading dimension of psi
  integer, intent(in) :: n
  !! real dimension of psi
  complex(DP), intent(inout) :: psi (lda, npol)
  !! a single wavefunction psi
  real(DP), intent(in) :: e
  !! Kohn -Sham eigenvalue
  !
  !    Local variables
  !
  real(DP), parameter :: eps = 1.0d-4
  !! a small number
  real(DP) :: scala = 1.0_dp
  !! scaling factor
  real(DP) :: x, denm
  !! various factors
  integer :: i, ipol
  ! do loop counters
  !
  call start_clock_gpu ('g_psi')
  !
  !$acc parallel loop collapse(2) present(h_diag, s_diag, psi)
  do ipol=1,npol
     do i = 1, n
#ifdef TEST_NEW_PRECONDITIONING
        x = (h_diag(i,ipol) - e * s_diag(i,ipol))*scala
        denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))/scala
#else
        !
        ! denm = g2+v(g=0)-e(k)*s  with  |denm| >= eps
        !
        denm = h_diag (i,ipol) - e * s_diag (i,ipol)
        if (abs (denm) < eps) denm = sign (eps, denm)
#endif
        psi (i, ipol) = psi (i, ipol) / denm
     enddo
  enddo
  !
  call stop_clock ('g_psi')
  !
end subroutine g_1psi
