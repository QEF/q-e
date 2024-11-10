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
  implicit none
  integer :: lda, n, m, npol, ipol
  ! input: the leading dimension of psi
  ! input: the real dimension of psi
  ! input: the number of bands
  ! input: the number of coordinates of psi
  ! local variable: counter of coordinates of psi
  real(DP) :: e (m)
  ! input: the eigenvectors
  complex(DP) :: psi (lda, npol, m)
  ! inp/out: the psi vector
  !
  !    Local variables
  !
  real(DP), parameter :: eps = 1.0d-4
  ! a small number
  real(DP) :: x, scala, denm
  integer :: k, i
  ! counter on psi functions
  ! counter on G vectors
  !
  call start_clock_gpu ('g_psi')
  !
  !$acc data present(h_diag, s_diag, psi, e)
  scala = 1.d0
  !$acc parallel loop collapse(3)
  do ipol=1,npol
     do k = 1, m
        do i = 1, n
#ifdef TEST_NEW_PRECONDITIONING
           x = (h_diag(i,ipol) - e(k)*s_diag(i,ipol))*scala
           denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))/scala
           psi (i, ipol, k) = psi (i, ipol, k) / denm
#else
           !
           ! denm = g2+v(g=0)-e(k)*s  with  |denm| >= eps
           !
           denm = h_diag (i,ipol) - e (k) * s_diag (i,ipol)
           if (abs (denm) < eps) denm = sign (eps, denm)
           psi (i, ipol, k) = psi (i, ipol, k) / denm
#endif
        enddo
     enddo
  enddo
  !$acc end data
  call stop_clock ('g_psi')
  return
end subroutine g_psi
!-----------------------------------------------------------------------
subroutine g_1psi (lda, n, psi, e)
  !-----------------------------------------------------------------------
  !
  !    As g_psi, for a single wavefunction
  !
  USE kinds
  USE noncollin_module,     ONLY : npol

  implicit none

  integer, intent(in) :: lda, & ! input: the leading dimension of psi
                         n      ! input: the real dimension of psi
  complex(DP) :: psi (lda, npol) ! inp/out: the psi vector
  real(DP), intent(in) :: e     ! input: the eigenvector
  !
  ! convert scalar e to size-1 vector [e] to exactly match g_psi argument type
  CALL g_psi (lda, n, 1, npol, psi, [e])

  return
  end
