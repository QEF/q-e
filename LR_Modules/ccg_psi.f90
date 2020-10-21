!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------
subroutine ccg_psi (lda, n, m, psi, h_diag, flag)
  !-----------------------------------------------------------------
  !
  !  This routine gives a preconditioning to the linear system solver.
  !  The preconditioning is diagonal in reciprocal space
  !
  USE kinds, only : DP
  USE noncollin_module, only : noncolin, npol
  implicit none

  integer :: lda, n, m, flag
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors
  ! input: flag=1 use h_diag, flag=-1 use conjg(h_diag)

  complex(DP) :: psi (lda*npol, m)
  ! inp/out: the vector to be preconditioned

  complex(DP) :: h_diag (lda*npol, m)
  ! input: the preconditioning vector

  integer :: k, i
  ! counter on bands
  ! counter on the elements of the vector
  !
  IF (flag == 1) THEN
     do k = 1, m
        do i = 1, n
           psi (i, k) = psi (i, k) * h_diag (i, k)
        enddo
        IF (noncolin) THEN
           do i = 1, n
              psi (i+lda, k) = psi (i+lda, k) * h_diag (i+lda, k)
           enddo
        END IF
     enddo
  ELSEIF (flag == -1) THEN
     do k = 1, m
        do i = 1, n
           psi (i, k) = psi (i, k) * CONJG(h_diag (i, k))
        enddo
        IF (noncolin) THEN
           do i = 1, n
              psi (i+lda, k) = psi (i+lda, k) * CONJG(h_diag (i+lda, k))
           enddo
        END IF
     enddo
  END IF
  !
  return
  !
end subroutine ccg_psi
