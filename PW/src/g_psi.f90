!
! Copyright (C) 2001-2024 Quantum ESPRESO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define TEST_NEW_PRECONDITIONING
!
!-----------------------------------------------------------------------
subroutine g_psi (lda, n, m, npol, psi, e)
  !---------------------------------------------------------------------
  USE kinds
  USE g_psi_mod, ONLY : h_diag, s_diag
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: lda
  !! leading dimension of psi
  INTEGER,intent(in) :: n
  !! real dimension of psi
  INTEGER,intent(in) :: m
  !! number of bands
  INTEGER, intent(in) :: npol
  !! number of polarizations
  COMPLEX(DP), intent(inout) :: psi(lda, npol, m)
  !! wavefunctions (psi)
  REAL(DP), intent(in) :: e(m)
  !! Kohn-Sham eigenvalues
  !
  !    Local variables
  !
  real(DP), parameter :: eps = 1.0d-4
  !! a small number
  real(DP) :: scala = 1.0_dp
  !! scaling factor
  real(DP) :: x, denm
  !! various factors
  integer :: k, i, ipol
  ! do loop counters
  !
  call start_clock_gpu ('g_psi')
  !
  !$acc parallel loop collapse(3) present(h_diag, s_diag, psi, e)
  do ipol=1,npol
     do k = 1, m
        do i = 1, n
#ifdef TEST_NEW_PRECONDITIONING
           x = (h_diag(i,ipol) - e(k)*s_diag(i,ipol))*scala
           denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))/scala
#else
           !
           ! denm = g2+v(g=0)-e(k)*s  with  |denm| >= eps
           !
           denm = h_diag (i,ipol) - e (k) * s_diag (i,ipol)
           if (abs (denm) < eps) denm = sign (eps, denm)
#endif
           psi (i, ipol, k) = psi (i, ipol, k) / denm
        enddo
     enddo
  enddo
  call stop_clock ('g_psi')
  !
end subroutine g_psi
