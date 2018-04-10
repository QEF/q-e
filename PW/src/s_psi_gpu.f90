!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
#if defined(__CUDA)
SUBROUTINE s_psi_gpu( lda, n, m, psi_d, spsi_d )
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : npol
  USE cudafor
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m
  COMPLEX(DP), DEVICE, INTENT(IN)  :: psi_d(lda*npol,m) 
  COMPLEX(DP), DEVICE, INTENT(OUT) :: spsi_d(lda*npol,m)


  COMPLEX(DP), ALLOCATABLE  :: psi(:,:)
  COMPLEX(DP), ALLOCATABLE  :: spsi(:,:)

  
  ALLOCATE(psi(lda*npol,m), spsi(lda*npol,m) )
  
  psi = psi_d
  CALL s_psi( lda, n, m, psi, spsi )
  spsi_d = spsi
  
  DEALLOCATE(psi, spsi )
  
END SUBROUTINE s_psi_gpu
#endif
