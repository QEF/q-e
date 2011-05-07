!
! Copyright (C) 2001 PWSCF group
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
  !    This routine gives a preconditioning to the linear system solver.
  !    The preconditioning is diagonal in reciprocal space
  !
  !
  USE kinds, only : DP
  USE noncollin_module, ONLY : noncolin, npol
  implicit none

  integer :: lda, n, m, flag
  ! input: the leading dimension of the psi vector
  ! input: the real dimension of the vector
  ! input: the number of vectors
  ! input: flag=1 use h_diag, flag=-1 use conjg(h_diag)
  complex(kind=DP) :: psi (lda*npol, m)
  ! inp/out: the vector to be preconditioned

  complex(kind=DP) :: h_diag (lda*npol, m)
  ! input: the preconditioning vector

  integer :: k, i
  ! counter on bands
  ! counter on the elements of the vector
  !
  do k = 1, m
     do i = 1, n
       if (flag .eq. 1) then
         psi (i, k) = psi (i, k) * h_diag (i, k)
       else if (flag .eq. -1) then
         psi (i, k) = psi (i, k) * CONJG(h_diag (i, k))
       else
         print*, 'flag is neither 1 nor -1. Stop'
       endif  
     enddo
     IF (noncolin) THEN
        do i = 1, n
           if (flag .eq. 1) then
              psi (i+lda, k) = psi (i+lda, k) * h_diag (i+lda, k)
           else if (flag .eq. -1) then
              psi (i+lda, k) = psi (i+lda, k) * CONJG(h_diag (i+lda, k))
           else
              print*, 'flag is neither 1 nor -1. Stop'
           endif  
        end do
     END IF
  enddo
  return
end subroutine ccg_psi
