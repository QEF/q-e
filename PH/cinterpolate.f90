!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine cinterpolate (v, vs, iflag)
!
!     This subroutine interpolates :
!     from the smooth mesh (vs) to a  thicker mesh (v)  (iflag>0)
!        vs is unchanged on output
!     from the  thick mesh (v ) to a smoother mesh (vs) (iflag<=0)
!        v  is unchanged on output
!     V and Vs are complex and in real space . V and Vs may coincide
!
#include"machine.h"
use pwcom
complex(kind=DP) :: v (nrxx), vs (nrxxs)
                                 ! function on thick mesh
                                 ! function on smooth mesh

integer :: iflag
                                 ! gives the direction of the interpolat

complex(kind=DP), allocatable :: aux (:), auxs (:)
                                 ! work array on thick mesh
                                 ! work array on smooth mesh

integer :: ig

call start_clock ('cinterpolate')
if (iflag.le.0) then
!
!    from thick to smooth
!
   if (doublegrid) then
      allocate (aux ( nrxx))    
      call ZCOPY (nrxx, v, 1, aux, 1)
      call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
      call setv (2 * nrxxs, 0.d0, vs, 1)
      do ig = 1, ngms
      vs (nls (ig) ) = aux (nl (ig) )
      enddo
      call cft3s (vs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1)
      deallocate (aux)
   else
      call ZCOPY (nrxx, v, 1, vs, 1)
   endif
else
!
!   from smooth to thick
!
   if (doublegrid) then
      allocate (auxs ( nrxxs))    
      call ZCOPY (nrxxs, vs, 1, auxs, 1)
      call cft3s (auxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, &
       - 1)
      call setv (2 * nrxx, 0.d0, v, 1)
      do ig = 1, ngms
      v (nl (ig) ) = auxs (nls (ig) )
      enddo
      call cft3 (v, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
      deallocate (auxs)
   else
      call ZCOPY (nrxx, vs, 1, v, 1)
   endif
endif
call stop_clock ('cinterpolate')
return
end subroutine cinterpolate
