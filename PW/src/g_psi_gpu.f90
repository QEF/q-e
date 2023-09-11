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
subroutine g_psi_gpu (lda, n, m, npol, psi_d, e_d)
  USE kinds
#if defined(__CUDA)
  USE cudafor
#endif
  USE g_psi_mod_gpum, ONLY : h_diag_d, s_diag_d, using_h_diag_d, using_s_diag_d
  !
  implicit none
  integer :: lda, n, m, npol, ipol
  ! input: the leading dimension of psi
  ! input: the real dimension of psi
  ! input: the number of bands
  ! input: the number of coordinates of psi
  ! local variable: counter of coordinates of psi
  real(DP) :: e_d (m)
  ! input: the eigenvectors
  complex(DP) :: psi_d (lda, npol, m)
  ! inp/out: the psi vector
#if defined(__CUDA)
  attributes( device ) :: e_d, psi_d
#endif
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
  call using_h_diag_d(0)
  call using_s_diag_d(0)
  !
  call start_clock_gpu ('g_psi')
  !
#ifdef TEST_NEW_PRECONDITIONING
  scala = 1.d0

!$cuf kernel do(3) <<<*,*>>>
  do ipol=1,npol
     do k = 1, m
        do i = 1, n
           x = (h_diag_d(i,ipol) - e_d(k)*s_diag_d(i,ipol))*scala
           denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))/scala
           psi_d (i, ipol, k) = psi_d (i, ipol, k) / denm
        enddo
     enddo
  enddo
#else
!$cuf kernel do(3) <<<*,*>>>
  do ipol=1,npol
     do k = 1, m
        do i = 1, n
           denm = h_diag_d (i,ipol) - e_d (k) * s_diag_d (i,ipol)
        !
        ! denm = g2+v(g=0) - e(k)
        !
           if (abs (denm) < eps) denm = sign (eps, denm)
        !
        ! denm = sign( max( abs(denm),eps ), denm )
        !
           psi_d (i, ipol, k) = psi_d (i, ipol, k) / denm
        enddo
     enddo
  enddo
#endif

  call stop_clock_gpu ('g_psi')
  return
end subroutine g_psi_gpu
!-----------------------------------------------------------------------
subroutine g_1psi_gpu (lda, n, psi_d, e_d)
  !-----------------------------------------------------------------------
  !
  !    This routine computes an estimate of the inverse Hamiltonian
  !    and applies it to one wavefunction
  !
  USE kinds
  USE noncollin_module,     ONLY : npol
  USE iso_c_binding
#if defined(__CUDA)
  USE cudafor
#endif

  implicit none

  integer :: lda, & ! input: the leading dimension of psi
             n      ! input: the real dimension of psi
  complex(DP) :: psi_d (lda, npol) ! inp/out: the psi vector
  real(DP), target :: e_d     ! input: the eigenvectors
  real(DP), dimension(:), pointer :: e_d_vec
#if defined(__CUDA)
  attributes(device) :: psi_d, e_d, e_d_vec
#endif
  !
  call start_clock ('g_1psi')

  ! cast scalar to size 1 vector to exactly match g_psi_gpu argument type
#if defined(__CUDA)
  call C_F_POINTER(C_DEVLOC(e_d), e_d_vec, [1])
#else
  call C_F_POINTER(C_LOC(e_d), e_d_vec, [1])
#endif
  CALL g_psi_gpu (lda, n, 1, npol, psi_d, e_d_vec)

  call stop_clock ('g_1psi')

  return

end subroutine g_1psi_gpu
