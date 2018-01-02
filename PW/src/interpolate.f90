!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine fft_interpolate_real (v_in, dfft_in, v_out, dfft_out )
  !
  !   This subroutine interpolates an array  v_in   defined on fft grit  dfft_in
  !                           to   an array  v_out  defined on fft grid  dfft_out
  !   v_in and v_out are assumed to be real arrays and may concide
  !
  USE kinds,          ONLY: DP
  USE control_flags,  ONLY : gamma_only
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : fwfft, invfft
  ! I/O variables
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
  REAL(DP),INTENT(IN)  :: v_in (dfft_in%nnr)
  REAL(DP),INTENT(OUT) :: v_out (dfft_out%nnr)
  ! local variables
  INTEGER :: ngm
  COMPLEX(DP), ALLOCATABLE :: aux_in (:), aux_out (:)

  call start_clock ('interpolate')

  IF (dfft_out%grid_id == dfft_in%grid_id) THEN

     v_out (:) = v_in (:)

  ELSE

     ALLOCATE (aux_in( dfft_in%nnr), aux_out(dfft_out%nnr))

     aux_in (:) = v_in(:)

     CALL fwfft ('Rho', aux_in, dfft_in)

     aux_out(:) = (0.d0, 0.d0)

     ngm = min(dfft_in%ngm, dfft_out%ngm)

     aux_out (dfft_out%nl (1:ngm) ) = aux_in (dfft_in%nl (1:ngm) )
     IF (gamma_only) aux_out (dfft_out%nlm (1:ngm) ) = aux_in (dfft_in%nlm (1:ngm) )

     CALL invfft ('Rho', aux_out, dfft_out)

     v_out (:) = aux_out (:)

     DEALLOCATE (aux_in, aux_out)

  END IF

  call stop_clock ('interpolate')

  return

end subroutine fft_interpolate_real
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
  USE gvecs,ONLY: ngms, doublegrid
  USE control_flags, ONLY: gamma_only
  USE fft_base,      ONLY : dfftp, dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  complex(DP) :: v (dfftp%nnr), vs (dffts%nnr)
  ! function on thick mesh
  ! function on smooth mesh

  integer :: iflag
  ! gives the direction of the interpolation

  complex(DP), allocatable :: aux (:), auxs (:)
  ! work array on thick mesh
  ! work array on smooth mesh

  integer :: ig

  if (gamma_only) call errore ('cinterpolate','not allowed', 1)
  call start_clock ('interpolate')
  if (iflag <= 0) then
     !
     !    from thick to smooth
     !
     if (doublegrid) then
        allocate (aux ( dfftp%nnr))    
        aux (:) = v(:) 
        CALL fwfft ('Rho', aux, dfftp)
        vs (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           vs (dffts%nl (ig) ) = aux (dfftp%nl (ig) )
        enddo
        CALL invfft ('Rho', vs, dffts)
        deallocate (aux)
     else
        call zcopy (dfftp%nnr, v, 1, vs, 1)
     endif
  else
     !
     !   from smooth to thick
     !
     if (doublegrid) then
        allocate (auxs (dffts%nnr))    
        auxs (:) = vs(:)
        CALL fwfft ('Rho', auxs, dffts)
        v (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           v (dfftp%nl (ig) ) = auxs (dffts%nl (ig) )
        enddo
        CALL invfft ('Rho', v, dfftp)
        deallocate (auxs)
     else
        call zcopy (dfftp%nnr, vs, 1, v, 1)
     endif
  endif
  call stop_clock ('interpolate')
  return
end subroutine cinterpolate
