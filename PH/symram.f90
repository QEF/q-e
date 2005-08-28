!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symram (phi, nsym, s, nat, irt)
  !-----------------------------------------------------------------------
  !
  ! Symmetrizes the Raman tensor.
  ! The tensor in input is a real tensor in crystal coordinates.
  ! The first two indexes correspond to the electric fields; the third
  !   to atomic displacements
  !
#include "f_defs.h"
  use kinds, only : DP
  implicit none

  integer :: nsym, s (3, 3, 48), nat, irt (48, nat)
  ! input: the number of symmetries
  ! input: the rotation matrix
  ! input: the number of atoms
  ! input: correspondence between rotated atoms

  real(DP) :: phi (3, 3, 3, nat)
  ! matrix to symmetrize

  integer :: isym, i, j, k, l, m, n, na, sna
  ! counter on symmetries
  ! counter on axis
  ! counter on atoms
  ! the rotated atom

  real(DP) :: work (3, 3, 3, nat)
  ! auxiliary work space
  !
  if (nsym == 1) return
  !
  work(:,:,:,:) = 0.d0
  !
  do na = 1, nat
     do isym = 1, nsym
        sna = irt (isym, na)
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
           do l = 1, 3
           do m = 1, 3
           do n = 1, 3

              work (i, j, k, na) = work (i, j, k, na) +    &
                                   s (i, l, isym) *        &
                                   s (j, m, isym) *        &
                                   s (k, n, isym) * phi (l, m, n, sna)
           enddo
           enddo
           enddo
        enddo
        enddo
        enddo
     enddo
  enddo
  !
  phi(:,:,:,:) =   work(:,:,:,:) / DBLE (nsym)
  !
  return
end subroutine symram
