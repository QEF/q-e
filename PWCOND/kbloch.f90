!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine kbloch(ntot, val)
!
! It obtaines complex k in 1st B.Z. for Bloch states from
! lambda=\exp{ikd}, d is the unit cell length.
! The result is in the units (2\pi/d)
!
#include "machine.h"
  USE kinds, only : DP
  use constants, only : tpi
  implicit none

  integer ::    &
         ntot,  &  ! number of Bloch states
         in
  real(kind=DP) :: d, tpiba, rho, g1, g2, k1, k2
  complex(kind=DP) ::  &
         val(ntot) ! complex k values

  do in=1, ntot
     g1=DREAL(val(in))
     g2=DIMAG(val(in))
     rho=DSQRT(g1**2+g2**2)
     k1=DACOS(g1/rho)
     k2=-DLOG(rho)       
     if (g2.le.0.d0) k1=tpi-k1 
     k1=k1/tpi
     k2=k2/tpi
     k1=k1-1.d0*INT(k1) 
     if (k1.gt.0.5d0) k1=k1-1.d0        
     val(in)=CMPLX(k1,k2) 
  !   WRITE( stdout,'(i5, 2f12.7)') in, DREAL(val(in)), DIMAG(val(in))
  enddo          

  return
end subroutine kbloch
