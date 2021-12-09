!
! Copyright (C) 2001-2012 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine symdyn_munu_new( dyn, u, xq, s, invs, rtau, irt, at, &
                            bg, nsymq, nat, irotmq, minus_q )
  !-----------------------------------------------------------------------
  !! This routine symmetrize the dynamical matrix written in the basis
  !! of the modes.
  
  USE kinds, only : DP
  
  implicit none
  
  integer :: nat
  !1 input: the number of atoms
  integer :: s(3,3,48)
  !! input: the symmetry matrices
  integer :: irt(48,nat)
  !! input: the rotated of each atom
  integer :: invs(48)
  !! input: the inverse of each matrix
  integer :: nsymq
  !! input: the order of the small group
  integer :: irotmq
  !! input: the small group of q.  
  !! input: the symmetry q -> -q+G
  real(DP) :: xq(3)
  !! input: the coordinates of q
  real(DP) :: rtau(3,48,nat)
  !! input: the R associated at each r
  real(DP) :: at(3,3)
  !! input: direct lattice vectors
  real(DP) :: bg(3,3)
  !! input: reciprocal lattice vectors  
  logical :: minus_q
  !! input: if true symmetry sends q-> -q+G
  complex(DP) :: dyn(3*nat,3*nat)
  !! inp/out: matrix to symmetrize
  complex(DP) :: u(3*nat,3*nat)
  !! input: the patterns
  !
  ! ... local variables
  !
  integer :: i, j, icart, jcart, na, nb, mu, nu
  ! counter on modes
  ! counter on modes
  ! counter on cartesian coordinates
  ! counter on cartesian coordinates
  ! counter on atoms
  ! counter on atoms
  ! counter on modes
  ! counter on modes

  complex(DP) :: work, phi (3, 3, nat, nat)
  ! auxiliary variable
  ! the dynamical matrix
  !
  ! First we transform in the cartesian coordinates
  !
  CALL dyn_pattern_to_cart(nat, u, dyn, phi)
  !
  ! Then we transform to the crystal axis
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, - 1)
     enddo
  enddo
  !
  !   And we symmetrize in this basis
  !
  call symdynph_gq_new (xq, phi, s, invs, rtau, irt, nsymq, nat, &
       irotmq, minus_q)
  !
  !  Back to cartesian coordinates
  !
  do na = 1, nat
     do nb = 1, nat
        call trntnsc (phi (1, 1, na, nb), at, bg, + 1)
     enddo
  enddo
  !
  !  rewrite the dynamical matrix on the array dyn with dimension 3nat x 3nat
  !
  CALL compact_dyn(nat, dyn, phi)

  return
end subroutine symdyn_munu_new
