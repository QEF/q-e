!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine interpolate (v, vs, iflag)  
  !
  !     This subroutine interpolates :
  !     from the smooth mesh (vs) to a  thicker mesh (v)  (iflag>0)
  !        vs is unchanged on output
  !     from the  thick mesh (v ) to a smoother mesh (vs) (iflag<=0)
  !        v  is unchanged on output
  !     V and Vs are real and in real space . V and Vs may coincide
  !
#include"machine.h"
  use pwcom  
  use allocate
  implicit none
  real(kind=DP) :: v (nrxx), vs (nrxxs)  
  ! function on thick mesh
  ! function on smooth mesh

  complex(kind=DP),pointer :: aux (:), auxs (:)  
  ! work array on thick mesh
  ! work array on smooth mesh

  integer :: iflag  
  ! gives the direction of the interpolat

  integer :: ig, ir  

  call start_clock ('interpolate')  
  if (iflag.le.0) then  
     !
     !    from thick to smooth
     !
     if (doublegrid) then  
        call mallocate(aux, nrxx)  
        call mallocate(auxs,nrxxs)  
        do ir = 1, nrxx  
           aux (ir) = v (ir)  
        enddo
        call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)  
        call setv (2 * nrxxs, 0.d0, auxs, 1)  
        do ig = 1, ngms  
           auxs (nls (ig) ) = aux (nl (ig) )  
        enddo
        call cft3s (auxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1)  
        do ir = 1, nrxxs  
           vs (ir) = auxs (ir)  
        enddo
        call mfree (auxs)  
        call mfree (aux)  
     else  
        call DCOPY (nrxx, v, 1, vs, 1)  
     endif
  else  
     !
     !   from smooth to thick
     !
     if (doublegrid) then  
        call mallocate(aux, nrxx)  
        call mallocate(auxs,nrxxs)    
        do ir = 1, nrxxs  
           auxs (ir) = vs (ir)  
        enddo
        call cft3s (auxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 1)
        call setv (2 * nrxx, 0.d0, aux, 1)  
        do ig = 1, ngms  
           aux (nl (ig) ) = auxs (nls (ig) )  
        enddo
        call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)  
        do ir = 1, nrxx  
           v (ir) = aux (ir)  
        enddo
        call mfree (auxs)  
        call mfree (aux)  
     else  
        call DCOPY (nrxx, vs, 1, v, 1)  
     endif
  endif
  call stop_clock ('interpolate')  
  return  
end subroutine interpolate

