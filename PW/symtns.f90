!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symtns (phi, nsym, s)
!-----------------------------------------------------------------------
!
!     symmetrize a tensor in the basis of crystallographic axis
!
#include "machine.h"
use parameters
implicit none
integer :: nsym, s (3, 3, 48), isym, i, j, k, l
real(kind=DP) :: phi (3, 3), work (3, 3)
external DSCAL, DCOPY
!
if (nsym.eq.1) return
call setv (9, 0.d0, work, 1)
!
do isym = 1, nsym
do i = 1, 3
do j = 1, 3
do k = 1, 3
do l = 1, 3
work (i, j) = work (i, j) + s (i, k, isym) * s (j, l, isym) &
 * phi (k, l)
enddo
enddo
enddo
enddo
enddo
!
call DSCAL (9, 1.d0 / nsym, work, 1)
call DCOPY (9, work, 1, phi, 1)
!

return
end subroutine symtns

