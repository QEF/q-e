!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine gk_sort (k, ngm, g, ecut, ngk, igk, gk)
  !-----------------------------------------------------------------------
  !
  ! sorts k+g in order of increasing magnitude, up to ecut
  ! NB: this version will yield the same ordering for different ecut
  !     and the same ordering in all machines
  !
  use parameters
  implicit none
  !
  !      Here the dummy variables
  !
  integer :: ngm, ngk, igk(ngk)
  ! input: the number of g vectors
  ! output: the number of k+G vectors inside
  ! output: the correspondence k+G <-> G
  real(kind=DP) :: k (3), g (3, ngm), ecut, gk (ngk)
  ! input: the k point
  ! input: the coordinates of G vectors
  ! input: the cut-off energy
  ! output: the moduli of k+G
  !
  !     one parameter
  !
  real(kind=DP) :: eps
  ! a small number
  parameter (eps = 1.d-8)
  !
  !      here the local variables
  !
  integer :: ng, nk
  ! counter on   G vectors
  ! counter on k+G vectors

  real(kind=DP) :: q, q2x, d (3), dnorm
  ! |k+G|^2
  ! upper bound for |G|
  ! d vector for ordering
  ! the norm of d
  !
  !    set the d vector for unique ordering
  !
  d (1) = 0.5465246754d0
  d (2) = 0.8765365676d0
  d (3) = 0.1524654376d0
  dnorm = sqrt (d (1) * d (1) + d (2) * d (2) + d (3) * d (3) )
  d (1) = d (1) / dnorm
  d (2) = d (2) / dnorm
  d (3) = d (3) / dnorm
  !
  !    first we count the number of k+G vectors inside the cut-off sphere
  !
  q2x = (sqrt (k (1) **2 + k (2) **2 + k (3) **2) + sqrt (ecut) ) ** &
       2
  ngk = 0
  do ng = 1, ngm
     q = (k (1) + g (1, ng) ) **2 + (k (2) + g (2, ng) ) **2 + (k (3) &
          + g (3, ng) ) **2
     !
     ! here if |k+G|^2 <= Ecut
     !
     if (q.le.ecut) then
        ngk = ngk + 1
        !
        ! gk is a fake quantity giving the same ordering on all machines
        !
        if (q.gt.eps) then
           gk (ngk) = 1.d4 * q + ( (k (1) + g (1, ng) ) * d (1) &
                + (k (2) + g (2, ng) ) * d (2) + (k (3) + g (3, ng) ) &
                * d (3) ) / sqrt (q)
        else
           gk (ngk) = 0.d0
        endif
        !   set the initial value of index array
        igk (ngk) = ng
     else
        !  if |G| > |k| + sqrt(Ecut)  stop search and order vectors
        if (g (1, ng) **2 + g (2, ng) **2 + g (3, ng) **2.gt.q2x) goto &
             10
     endif
  enddo
  call errore ('gk_sort', 'unexpected exit from do-loop', - 1)
  !
  ! order vector gk keeping initial position in index
  !
10 call hpsort (ngk, gk, igk)
  !
  !    now order true |k+G|
  !
  do nk = 1, ngk
     gk (nk) = (k (1) + g (1, igk (nk) ) ) **2 + (k (2) + g (2, igk ( &
          nk) ) ) **2 + (k (3) + g (3, igk (nk) ) ) **2
  enddo
  return
end subroutine gk_sort

!-----------------------------------------------------------------------
subroutine gk_l2gmap (ngm, ig_l2g, ngk, igk, igk_l2g)
  !-----------------------------------------------------------------------
  !
  ! This subroutine maps local G+k index to the global G vector index
  ! the mapping is used to collect wavefunctions subsets distributed
  ! across processors.
  ! Written by Carlo Cavazzoni
  !
  use parameters
  implicit none
  !
  !      Here the dummy variables
  !
  integer :: ngm, ngk, igk(ngk), ig_l2g(ngm)   ! input
  integer :: igk_l2g(ngk)                      ! output
  integer :: nk

  ! input: mapping between local and global G vector index
  do nk = 1, ngk
     igk_l2g(nk) = ig_l2g( igk(nk) )
  enddo
  return
end subroutine gk_l2gmap
