!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine symz (phi, nsym, s, nat, irt)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  USE kinds, only : DP
  implicit none

  integer :: nsym, s (3, 3, 48), nat, irt (48, nat)
  ! input: the number of symmetries
  ! input: the rotation matrix
  ! input: the number of atoms
  ! input: correspondence between rotated atoms

  real(kind=DP) :: phi (3, 3, nat)
  ! matrix to symmetrize

  integer :: isym, i, j, k, l, na, sna
  ! counter on symmetries
  ! counter on points
  ! counter on atoms
  ! the rotated atom

  real(kind=DP) :: work (3, 3, nat)
  ! auxiliary space

  !

  if (nsym.eq.1) return
  call setv (9 * nat, 0.d0, work, 1)
  !
  do na = 1, nat
     do isym = 1, nsym
        sna = irt (isym, na)
        do i = 1, 3
           do j = 1, 3
              do k = 1, 3
                 do l = 1, 3
                    work (i, j, na) = work (i, j, na) + s (i, k, isym) * s (j, l, &
                         isym) * phi (k, l, sna)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  call DSCAL (9 * nat, 1.d0 / float (nsym), work, 1)
  call DCOPY (9 * nat, work, 1, phi, 1)
  !
  return
end subroutine symz
