!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine updathes (nax3, nat3, oldforce, force, hessm1, dtau)
   !-----------------------------------------------------------------------
   !
   ! updates the inverse hessian
   !
#include "f_defs.h"
   USE kinds
   implicit none
   integer :: nat3, nax3

   real(kind=DP) :: hessm1(nax3,nat3), dtau(nat3), force(nat3), oldforce(nat3)
   integer :: i, j
                                                    ! work arrays
   real(kind=DP) :: fac1, fac2, DDOT
   real(kind=DP), allocatable ::  hdg (:), u (:)
   external DAXPY, DDOT, DGEMV

   allocate ( hdg(nat3),u(nat3) )

   call DAXPY (nat3, - 1.d0, force, 1, oldforce, 1)
   hdg(:) = 0.d0
   call DGEMV ('n', nat3, nat3, 1.d0, hessm1, nax3, oldforce, 1, 0.d0, hdg, 1)
   !      call matv (nax3,nat3,hessm1,oldforce,hdg)

   fac1 = DDOT (nat3, oldforce, 1, dtau, 1)
   fac2 = DDOT (nat3, oldforce, 1, hdg, 1)
   do i = 1, nat3
      u(i) = dtau (i) / fac1 - hdg (i) / fac2
   enddo
   do i = 1, nat3
      do j = 1, nat3
         hessm1(i,j) = hessm1(i,j) + dtau(i) * dtau(j) / fac1 - &
                       hdg(i) * hdg(j) / fac2 + u(i) * u(j) * fac2
      enddo
   enddo
   deallocate (hdg,u)
   return
end subroutine updathes

