!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define TEST_NEW_PRECONDITIONING
!
!-----------------------------------------------------------------------
subroutine g_psi (lda, n, m, npol, psi, e)
  !-----------------------------------------------------------------------
  !
  !    This routine computes an estimate of the inverse Hamiltonian
  !    and applies it to m wavefunctions
  !
  USE kinds
  USE g_psi_mod
  USE g_psi_mod_gpum, ONLY : using_h_diag, using_s_diag
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
  !
  integer :: k, i
  ! counter on psi functions
  ! counter on G vectors
  integer, parameter :: blocksize = 256
  integer :: iblock, numblock
  ! chunking parameters
  !
  call using_h_diag(0); call using_s_diag(0)
  call start_clock ('g_psi')
  !
  ! compute the number of chuncks
  numblock  = (n+blocksize-1)/blocksize
  !
#ifdef TEST_NEW_PRECONDITIONING
  scala = 1.d0
  !$omp parallel do collapse(3) private(x, denm)
  do k = 1, m
     do ipol=1, npol
        do iblock = 1, numblock
           do i = (iblock-1)*blocksize+1, MIN(iblock*blocksize, n)
              x = (h_diag(i,ipol) - e(k)*s_diag(i,ipol))*scala
              denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))/scala
              psi (i, ipol, k) = psi (i, ipol, k) / denm
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do
#else
  !$omp parallel do collapse(3) private(denm)
  do ipol=1,npol
     do k = 1, m
        do iblock = 1, numblock
           do i = (iblock-1)*blocksize+1, MIN(iblock*blocksize, n)
              denm = h_diag (i,ipol) - e (k) * s_diag (i,ipol)
              !
              ! denm = g2+v(g=0) - e(k)
              !
                 if (abs (denm) < eps) denm = sign (eps, denm)
              !
              ! denm = sign( max( abs(denm),eps ), denm )
              !
              psi (i, ipol, k) = psi (i, ipol, k) / denm
           enddo
        enddo
     enddo
  enddo
  !$omp end parallel do
#endif

  call stop_clock ('g_psi')
  return
end subroutine g_psi
