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
#if defined(__CUDA)
subroutine g_psi_gpu (lda, n, m, npol, psi_d, e_d)
  !-----------------------------------------------------------------------
  !
  !    This routine computes an estimate of the inverse Hamiltonian
  !    and applies it to m wavefunctions
  !
  USE kinds
  USE cudafor
  USE g_psi_mod, ONLY : h_diag, s_diag
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
  ! PUT BUFFER HERE!!!! buffer buFeFERRR
  
  !
  attributes( device ) :: e_d, psi_d
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

  COMPLEX(DP), ALLOCATABLE :: psi(:,:,:)
  REAL(DP), ALLOCATABLE :: e(:)

  
  ALLOCATE(psi(lda,npol,m), e(m) )
  
  psi = psi_d; e = e_d
  CALL g_psi(lda, n, m, npol, psi, e)
  psi_d = psi; e_d = e
  
  DEALLOCATE(psi, e )

  return
end subroutine g_psi_gpu
#endif 
