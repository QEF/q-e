!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------
SUBROUTINE pbcg_psi (lda, n, m, psi, h_diag, flag)
  !-----------------------------------------------------------------
  !
  !    This routine gives a preconditioning to the linear system solver.
  !    The preconditioning is diagonal in reciprocal space
  !
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE

  INTEGER :: lda, n, m, flag
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors
  ! input: flag=1 use h_diag, flag=-1 use CONJG(h_diag)
  COMPLEX(kind=DP) :: psi (lda, m)
  ! inp/out: the vector to be preconditioned

  COMPLEX(kind=DP) :: h_diag (lda, m)
  ! input: the preconditioning vector

  INTEGER :: k, i
  ! counter on bands
  ! counter on the elements of the vector
  !
  DO k = 1, m
     DO i = 1, n
       IF (flag == 1) THEN
         psi (i, k) = psi (i, k) * h_diag (i, k)
       ELSEIF (flag == -1) THEN
         psi (i, k) = psi (i, k) * conjg(h_diag (i, k))
       ELSE
         PRINT*, 'flag is neither 1 nor -1. Stop'
       ENDIF
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE pbcg_psi
