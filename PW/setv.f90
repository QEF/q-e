!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine setv (n, sa, sx, incx)  
!----------------------------------------------------------------------
!
!     This routine sets the array sx to the value sa.
!
use parameters
implicit none  
real(kind=DP) :: sa, sx ( * )  
                           ! The value to be set
                           ! The array to be set
integer :: n, incx, i  
                           ! The number of elements which are changed
                           ! The distance between changed elements
                           ! The positions counter
#ifdef T3D_BENCHLIB
integer :: lputp  
external lputp  
!
! t3e optimization (SC 2/2/99)
!
if (incx.eq.1) then  
   call lsetv (sx, 1, n, sa)  
  123    if (lputp () .ne.0) goto 123  
else  

!DIR$ CACHE_BYPASS  SX
   do i = 1, (n - 1) * incx + 1, incx  
   sx (i) = sa  
   enddo  
endif  
#else
do i = 1, (n - 1) * incx + 1, incx  
   sx (i) = sa  
enddo  
#endif
return  
end subroutine setv

