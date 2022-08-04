!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
!
MODULE fft_wave
  !
  !! This module contains wrapper to FFT and inverse FFTs of w.f.
  !
  USE kinds,           ONLY: DP
  USE fft_interfaces,  ONLY: invfft
  USE fft_types,       ONLY: fft_type_descriptor
  USE control_flags,   ONLY: gamma_only
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wave_g2r
  !
CONTAINS
  !
  !
  !----------------------------------------------------------------------
  SUBROUTINE wave_g2r( f_in, f_out, dfft, dim2, igk )
    !--------------------------------------------------------------------
    !
    USE fft_helper_subroutines, ONLY: c2psi_gamma, c2psi_k
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    INTEGER, INTENT(IN) :: dim2
    COMPLEX(DP) :: f_in(:,:)
    COMPLEX(DP) :: f_out(:)
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !
    INTEGER :: i2, npw, numblock
    INTEGER :: j, idx, ioff, ntgrp, right_nnr
    INTEGER, PARAMETER :: blocksize = 256
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !
    npw = SIZE(f_in(:,1))
    !
    IF (gamma_only) THEN
      IF ( dim2/=2 ) CALL c2psi_gamma( dfft, f_out, f_in(:,1) )
      IF ( dim2==2 ) CALL c2psi_gamma( dfft, f_out, f_in(:,1), f_in(:,2) )
    ELSE
      !$acc data present_or_copyin(igk)
#if defined(__CUDA)
      CALL c2psi_k( dfft, f_out, f_in(:,1), igk, npw )
#else
      !$omp parallel
      CALL threaded_barrier_memset( f_out, 0.D0, dfft%nnr*2 )
      CALL c2psi_k( dfft, f_out, f_in(:,1), igk, npw )
      !$omp end parallel
#endif
      !$acc end data
    ENDIF
    !
    !$acc host_data use_device( f_out )
    CALL invfft( 'Wave', f_out, dfft )
    !$acc end host_data
    !
    !$acc end data
    !
    RETURN
    !
  END SUBROUTINE wave_g2r
  !
  !
END MODULE fft_wave
