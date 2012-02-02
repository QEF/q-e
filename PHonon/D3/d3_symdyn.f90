!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine d3_symdyn (d3dyn, u, ug0, xq, s, invs, rtau, irt, irgq, &
     at, bg, nsymq, nat, irotmq, minus_q, npert_i, npert_f)
  !-----------------------------------------------------------------------
  !
  !    This routine symmetrize the dynamical matrix written in the basis
  !    of the modes
  !
  !
  USE kinds, only : DP
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum

  implicit none
  integer :: nat, s (3, 3, 48), irt (48, nat), irgq (48), invs (48), &
       nsymq, npert_i, npert_f, irotmq
  ! input: the number of atoms
  ! input: the symmetry matrices
  ! input: the rotated of each atom
  ! input: the small group of q
  ! input: the inverse of each matrix
  ! input: the order of the small gro
  ! input: the symmetry q -> -q+G

  real (DP) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3)
  ! input: the coordinates of q
  ! input: the R associated at each r
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors

  logical :: minus_q
  ! input: if true symmetry sends q->

  complex (DP) :: d3dyn (3 * nat, 3 * nat, 3 * nat), &
       ug0 (3 * nat, 3 * nat), u (3 * nat, 3 * nat)
  ! inp/out: matrix to symmetr
  ! input: the q=0 patterns
  ! input: the patterns

  integer :: i, j, i1, icart, jcart, kcart, na, nb, nc, mu, nu, om
  ! counters

  complex (DP) :: work, wrk (3, 3)
  ! auxiliary variables
  complex (DP), allocatable :: phi (:,:,:,:,:,:)
  ! the dynamical matrix

  allocate  (phi( 3, 3, 3, nat, nat, nat))
  !
  ! First we transform in the cartesian coordinates
  !
  phi = (0.d0, 0.d0)
  do i1 = npert_i, npert_f
     nc = (i1 - 1) / 3 + 1
     kcart = i1 - 3 * (nc - 1)
     do i = 1, 3 * nat
        na = (i - 1) / 3 + 1
        icart = i - 3 * (na - 1)
        do j = 1, 3 * nat
           nb = (j - 1) / 3 + 1
           jcart = j - 3 * (nb - 1)
           work = (0.d0, 0.d0)
           do om = 1, 3 * nat
              do mu = 1, 3 * nat
                 do nu = 1, 3 * nat
                    work = work + CONJG(ug0 (i1, om) ) * u (i, mu) * &
                         d3dyn (om, mu, nu) * CONJG(u (j, nu) )
                 enddo
              enddo
           enddo
           phi (kcart, icart, jcart, nc, na, nb) = work
        enddo
     enddo
  enddo
#ifdef __MPI
  call mp_sum( phi, inter_pool_comm )
#endif
  !
  ! Then we transform to the crystal axis
  !
  do nc = 1, nat
     do na = 1, nat
        do nb = 1, nat
           call trntnsc_3 (phi (1, 1, 1, nc, na, nb), at, bg, - 1)
        enddo
     enddo
  enddo
  !
  !   And we symmetrize in this basis
  !
  call d3_symdynph (xq, phi, s, invs, rtau, irt, irgq, nsymq, nat, &
       irotmq, minus_q)
  !
  !  Back to cartesian coordinates
  !
  do nc = 1, nat
     do na = 1, nat
        do nb = 1, nat
           call trntnsc_3 (phi (1, 1, 1, nc, na, nb), at, bg, + 1)
        enddo
     enddo
  enddo
  !
  !  rewrite the dynamical matrix on the array dyn with dimension 3nat x 3
  !
  do i1 = 1, 3 * nat
     nc = (i1 - 1) / 3 + 1
     kcart = i1 - 3 * (nc - 1)
     do i = 1, 3 * nat
        na = (i - 1) / 3 + 1
        icart = i - 3 * (na - 1)
        do j = 1, 3 * nat
           nb = (j - 1) / 3 + 1
           jcart = j - 3 * (nb - 1)
           d3dyn (i1, i, j) = phi (kcart, icart, jcart, nc, na, nb)
        enddo
     enddo

  enddo
  deallocate (phi)

  return

end subroutine d3_symdyn
