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
  !! This module contains wrapper routines that enclose the calls to FFTXlib needed
  !! by QE. The 
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
#if defined(__CUDA)
    USE fft_helper_subroutines, ONLY: c2psi_gamma_gpu, c2psi_k_gpu
#else
    USE fft_helper_subroutines, ONLY: c2psi_gamma_cpu, c2psi_k_cpu
#endif
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
#if defined(__CUDA) && defined(_OPENACC)
    INTEGER, POINTER, DEVICE :: nl_d(:), nlm_d(:)
    !
    nl_d  => dfft%nl_d
    nlm_d => dfft%nlm_d
#endif
    !
    !$acc data present_or_copyin(f_in,igk) present_or_copyout(f_out)
    !
    npw = SIZE(f_in(:,1))
    !
    IF (gamma_only) THEN
#if defined(__CUDA)
      IF ( dim2/=2 ) CALL c2psi_gamma_gpu( dfft, f_out, f_in(:,1) )
      IF ( dim2==2 ) CALL c2psi_gamma_gpu( dfft, f_out, f_in(:,1), f_in(:,2) )
#else
      IF ( dim2/=2 ) CALL c2psi_gamma_cpu( dfft, f_out, f_in(:,1) )
      IF ( dim2==2 ) CALL c2psi_gamma_cpu( dfft, f_out, f_in(:,1), f_in(:,2) )
#endif
    ELSE
#if defined(__CUDA)
      CALL c2psi_k_gpu( dfft, f_out, f_in(:,1), igk, npw )
#else
      !$omp parallel
      CALL threaded_barrier_memset( f_out, 0.D0, dfft%nnr*2 )
      CALL c2psi_k_cpu( dfft, f_out, f_in(:,1), igk, npw )
      !$omp end parallel
#endif
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
