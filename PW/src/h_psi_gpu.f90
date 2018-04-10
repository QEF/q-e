
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
#if defined(__CUDA)
SUBROUTINE h_psi_gpu( lda, n, m, psi_d, hpsi_d )
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m
  COMPLEX(DP), DEVICE, INTENT(IN)  :: psi_d(lda*npol,m) 
  COMPLEX(DP), DEVICE, INTENT(OUT) :: hpsi_d(lda*npol,m)


  COMPLEX(DP), ALLOCATABLE :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:)

  
  ALLOCATE(psi(lda*npol,m), hpsi(lda*npol,m) )
  
  psi = psi_d
  CALL h_psi( lda, n, m, psi, hpsi )
  hpsi_d = hpsi
  
  DEALLOCATE(psi, hpsi )
  
END SUBROUTINE h_psi_gpu
#endif
