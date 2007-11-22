!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE kinds, ONLY: DP
  USE gvect,  ONLY: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, nlm
  USE gsmooth,ONLY: nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,nrxxs,ngms, &
       nls, nlsm, doublegrid
  USE control_flags, ONLY: gamma_only
  !
  implicit none
  real(DP) :: v (nrxx), vs (nrxxs)
  ! function on thick mesh
  ! function on smooth mesh

  complex(DP), allocatable :: aux (:), auxs (:)
  ! work array on thick mesh
  ! work array on smooth mesh

  integer :: iflag
  ! gives the direction of the interpolation

  integer :: ig, ir

  call start_clock ('interpolate')
  if (iflag <= 0) then
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
        enddo
        if (gamma_only) then
           do ig = 1, ngms
              auxs (nlsm(ig) ) = aux (nlm(ig) )
           enddo
        end if
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
        enddo
        if (gamma_only) then
           do ig = 1, ngms
              aux (nlm(ig) ) = auxs (nlsm(ig) )
           enddo
        end if
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
  USE kinds, ONLY: DP
  USE gvect,  ONLY: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nl, nlm
  USE gsmooth,ONLY: nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,nrxxs,ngms, &
       nls, nlsm, doublegrid
  USE control_flags, ONLY: gamma_only
  !
  IMPLICIT NONE
  complex(DP) :: v (nrxx), vs (nrxxs)
  ! function on thick mesh
  ! function on smooth mesh

  integer :: iflag
  ! gives the direction of the interpolation

  complex(DP), allocatable :: aux (:), auxs (:)
  ! work array on thick mesh
  ! work array on smooth mesh

  integer :: ig

  if (gamma_only) call errore ('cinterpolate','not allowed', 1)
  call start_clock ('cinterpolate')
  if (iflag <= 0) then
     !
     !    from thick to smooth
     !
     if (doublegrid) then
        allocate (aux ( nrxx))    
        aux (:) = v(:) 
        call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
        vs (:) = (0.d0, 0.d0)
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
        auxs (:) = vs(:)
        call cft3s (auxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1)
        v (:) = (0.d0, 0.d0)
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
