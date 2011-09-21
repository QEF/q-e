!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, &
     irt, irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, &
     npert, nirr, gi, gimq, iverbosity)
  !---------------------------------------------------------------------
  !
  !     This routine substitute set_irr when there are no symmetries.
  !     The irreducible representations are all one dimensional and
  !     we set them to the displacement of a single atom in one direction
  !
  USE kinds, only : DP
  implicit none
  !
  !   first the dummy variables
  !

  integer :: nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       iverbosity, npert (3 * nat), irgq (48), nsymq, irotmq, nirr, npertx
  ! input: the number of atoms
  ! input: the number of symmetries
  ! input: the symmetry matrices
  ! input: the inverse of each matrix
  ! input: the rotated of each atom
  ! input: write control
  ! output: the dimension of each represe
  ! output: the small group of q
  ! output: the order of the small group
  ! output: the symmetry sending q -> -q+
  ! output: the number of irr. representa

  real(DP) :: xq (3), rtau (3, 48, nat), at (3, 3), bg (3, 3), &
       gi (3, 48), gimq (3)
  ! input: the q point
  ! input: the R associated to each tau
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  ! output: [S(irotq)*q - q]
  ! output: [S(irotmq)*q + q]

  complex(DP) :: u(3*nat, 3*nat), t(npertx, npertx, 48, 3*nat),&
       tmq (npertx, npertx, 3 * nat)
  ! output: the pattern vectors
  ! output: the symmetry matrices
  ! output: the matrice sending q -> -q+G

  logical :: minus_q
  ! output: if true one symmetry send q -> -q+G
  integer :: imode
  ! counter on modes
  !
  !    set the information on the symmetry group
  !
  call smallgq (xq,at,bg,s,nsym,irgq,nsymq,irotmq,minus_q,gi,gimq)
  !
  !     set the modes
  !
  u (:,:) = (0.d0, 0.d0)
  do imode = 1, 3 * nat
     u (imode, imode) = (1.d0, 0.d0)
  enddo
  nirr = 3 * nat
  do imode = 1, 3 * nat
     npert (imode) = 1
  enddo
  !
  !   And we compute the matrices which represent the symmetry transformat
  !   in the basis of the displacements
  !
  t(:, :, :, :) = (0.d0, 0.d0)
  do imode = 1, 3 * nat
     t (1, 1, 1, imode) = (1.d0, 0.d0)
  enddo

  tmq (:, :, :) = (0.d0, 0.d0)
  if (minus_q) then
     tmq (1, 1, :) = (1.d0, 0.d0)
  end if

  return
end subroutine set_irr_nosym
