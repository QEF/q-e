!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine simpson (mesh, func, rab, asum)  
!-----------------------------------------------------------------------
!
!     simpson's rule integrator for function stored on the
!     radial logarithmic mesh
!
use parameters
implicit none  
integer :: i, mesh  


real(kind=DP) :: rab (mesh), func (mesh), f1, f2, f3, r12, asum  
!     routine assumes that mesh is an odd number so run check
!     if ( mesh+1 - ( (mesh+1) / 2 ) * 2 .ne. 1 ) then
!       write(*,*) '***error in subroutine radlg'
!       write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1
!       stop
!     endif
asum = 0.0d0  
r12 = 1.0d0 / 12.0d0  

f3 = func (1) * rab (1) * r12  
do i = 2, mesh - 1, 2  
f1 = f3  
f2 = func (i) * rab (i) * r12  
f3 = func (i + 1) * rab (i + 1) * r12  
asum = asum + 4.0d0 * f1 + 16.0d0 * f2 + 4.0d0 * f3  

enddo  
return  
end subroutine simpson
