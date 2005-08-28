!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine sym_elop (phi, nsym, s, nat, irt)
  !-----------------------------------------------------------------------
  !
  ! Symmetrizes the Electr-optic tensor
  ! The tensor in input is a real tensor in crystal coordinates.
  !
#include "f_defs.h"
  use kinds, only : DP
  implicit none

  integer :: nsym, s (3, 3, 48), nat, irt (48, nat)
  ! input: the number of symmetries
  ! input: the rotation matrix
  ! input: the number of atoms
  ! input: correspondence between rotated atoms

  real(DP) :: phi (3, 3, 3)
  ! matrix to symmetrize

  integer :: isym, i, j, k, l, m, n, na, sna
  ! counter on symmetries
  ! counter on axis
  ! counter on atoms
  ! the rotated atom

  real(DP) :: work (3, 3, 3)
  ! auxiliary space

  if (nsym.eq.1) return
  work (:,:,:) = 0.d0
  do isym = 1, nsym
     do i = 1, 3
     do j = 1, 3
     do k = 1, 3
        do l = 1, 3
        do m = 1, 3
        do n = 1, 3
           work (i, j, k) = work (i, j, k) +        &
                            s (i, l, isym) *        &
                            s (j, m, isym) *        &
                            s (k, n, isym) * phi (l, m, n)
        enddo
        enddo
        enddo
     enddo
     enddo
     enddo
  enddo

  call DSCAL (27, 1.d0 / DBLE (nsym), work, 1)
  call DCOPY (27, work, 1, phi, 1)

  return
end subroutine sym_elop

