!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifdef T3D_BENCHLIB
!================================================================
subroutine scopy_t3e (n, a, ia, b, ib)
!----------------------------------------------------------------
! optimized scopy for t3e system (SC 2/99)
!
use parameters
implicit none
integer :: i, ni
integer,  intent (in) ::n
integer,  intent (in) ::ia
integer,  intent (in) ::ib
real(kind=DP) :: b (n)
real(kind=DP) :: a (n)

integer :: lputp
include "mpp/shmem.fh"
intrinsic my_pe

external lputp

if ( (ia.eq.1) .and. (ib.eq.1) ) then
! benchlib ...
   do i = 1, n, 480
   ni = min (480, n - i + 1)
   call lgetv (a (i), 1, ni)
   call lputv (b (i), 1, ni)
   enddo

  123    if (lputp () .ne.0) goto 123

else
   call scopy (n, a, ia, b, ib)

endif
return

end subroutine scopy_t3e
#else
subroutine dummy
return

end subroutine dummy

#endif
