!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine fft_interpolate_real (dfft_in, v_in, dfft_out, v_out )
  !
  !   This subroutine interpolates an array  v_in   defined on fft grid  dfft_in
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
subroutine fft_interpolate_complex (dfft_in, v_in, dfft_out, v_out )
  !
  !   This subroutine interpolates an array  v_in   defined on fft grid  dfft_in
  !                           to   an array  v_out  defined on fft grid  dfft_out
  !   v_in and v_out are assumed to be complex arrays and may concide
  !
  USE kinds,          ONLY: DP
  USE control_flags,  ONLY : gamma_only
  USE fft_types,      ONLY : fft_type_descriptor
  USE fft_interfaces, ONLY : fwfft, invfft
  ! I/O variables
  TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
  COMPLEX(DP),INTENT(IN)  :: v_in (dfft_in%nnr)
  COMPLEX(DP),INTENT(OUT) :: v_out (dfft_out%nnr)
  ! local variables
  INTEGER :: ngm
  COMPLEX(DP), ALLOCATABLE :: aux_in (:)

  if (gamma_only) call errore ('cinterpolate','not allowed', 1)
  call start_clock ('interpolate')

  IF (dfft_out%grid_id == dfft_in%grid_id) THEN

     v_out (:) = v_in (:)

  ELSE

     ALLOCATE (aux_in( dfft_in%nnr))

     aux_in (:) = v_in(:)

     CALL fwfft ('Rho', aux_in, dfft_in)

     v_out(:) = (0.d0, 0.d0)

     ngm = min(dfft_in%ngm, dfft_out%ngm)

     v_out (dfft_out%nl (1:ngm) ) = aux_in (dfft_in%nl (1:ngm) )

     CALL invfft ('Rho', v_out, dfft_out)

     DEALLOCATE (aux_in)

  END IF

  call stop_clock ('interpolate')

  return

end subroutine fft_interpolate_complex
!
