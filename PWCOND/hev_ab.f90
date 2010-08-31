!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine hev_ab (n, amt, lda, eigen, veigen, el, eh, m)
!
! It solves EP: A x = lambda X with hermitean matrix A.
!
  USE kinds, only : DP
  implicit none
  character :: lrange
  integer ::                                  &
       n,      &  ! dimension of the matrix to be diagonalized
       m,      &  ! input:  flag
                  ! output: number of eigenvalues in the [el,eh]
       il, ih, &  ! number of lowest and highest eigenvalues
       lda,    &  ! leading dimension of A
       lwork,  &  ! aux. var.
       info       ! -1 --> all eigenvalues are computed

  integer, allocatable :: iwork(:), ifail(:)
  real(DP) ::  &
       eigen(n),    & ! eigenvalues
       el, eh,      & ! interval for eigenvalue searching
       abstol         ! accuracy for eigenvalues
  real(DP), allocatable  :: rwork(:)
  complex(DP), allocatable :: work(:)
  complex(DP) ::  &
       amt(lda, n),    &   ! A
       veigen(lda, n)      ! X

  lwork = 16*n

  allocate( work( lwork ) )
  allocate( rwork( 7*n ) )
  allocate( iwork( 5*n ) )
  allocate( ifail( n ) )

  abstol=0.d0

!
! If m=-1 find all the solutions, otherwise only el<lambda<eh
!
  if (m.eq.-1) then
     lrange='A'
  else
     lrange='V'
  endif

  call ZHEEVX('V', lrange, 'U', n, amt, lda, el, eh, il, ih, abstol,  &
               m, eigen, veigen, n, work, lwork, rwork, iwork, ifail, &
               info)
  call errore ('hev_ab','info =/= 0',abs(info))

  deallocate(rwork)
  deallocate(work)
  deallocate(iwork)
  deallocate(ifail)

  return
end subroutine hev_ab
