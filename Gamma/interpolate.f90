!
! Copyright (C) 2003 PWSCF group
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
  use gamma
  implicit none
  real(kind=DP) :: v (nrxx), vs (nrxxs)
  ! function on thick mesh
  ! function on smooth mesh

  complex(kind=DP), allocatable :: aux (:), auxs (:)
  ! work array on thick mesh
  ! work array on smooth mesh

  integer :: iflag
  ! gives the direction of the interpolation

  integer :: ig, ir

  call start_clock ('interpolate')
  if (iflag.le.0) then
     !
     !    from thick to smooth
     !
     if (doublegrid) then
        allocate (aux( nrxx))    
        allocate (auxs(nrxxs))    
        aux (:) = v (:)
        call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
        auxs (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           auxs (nls (ig) ) = aux (nl (ig) )
           auxs (nlsm(ig) ) = aux (nlm(ig) )
        enddo
        call cft3s (auxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1)
        vs (:) = auxs (:)
        deallocate (auxs)
        deallocate (aux)
     else
        do ir = 1, nrxx
           vs (ir) = v (ir)
        enddo
     endif
  else
     !
     !   from smooth to thick
     !
     if (doublegrid) then
        allocate (aux( nrxx))    
        allocate (auxs(nrxxs))    
        auxs (:) = vs (:)
        call cft3s (auxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 1)
        aux (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           aux (nl (ig) ) = auxs (nls (ig) )
           aux (nlm(ig) ) = auxs (nlsm(ig) )
        enddo
        call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        v (:) = aux (:)
        deallocate (auxs)
        deallocate (aux)
     else
        do ir = 1, nrxx
           v (ir) = vs (ir)
        enddo
     endif
  endif
  call stop_clock ('interpolate')
  return
end subroutine interpolate

