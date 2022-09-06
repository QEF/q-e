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
  USE fft_interfaces,  ONLY: fwfft, invfft
  USE fft_types,       ONLY: fft_type_descriptor
  USE control_flags,   ONLY: gamma_only
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wave_r2g, wave_g2r, tgwave_g2r, tgwave_r2g
  !
CONTAINS
  !
  !
  !----------------------------------------------------------------------
  SUBROUTINE wave_r2g( f_in, f_out, dfft, igk )
    !--------------------------------------------------------------------
    !! Wave function FFT from R to G-space.
    !
    USE fft_helper_subroutines,  ONLY: fftx_psi2c_gamma, fftx_psi2c_k
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP), INTENT(IN)  :: f_in(:)
    COMPLEX(DP), INTENT(OUT) :: f_out(:,:)
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !
    COMPLEX(DP), ALLOCATABLE :: psic(:)
    INTEGER :: dim2, nrxxs
    !
    nrxxs = SIZE(f_in)
    dim2 = SIZE(f_out(1,:))
    !
    ALLOCATE( psic(nrxxs) )
    psic = f_in
    !
    CALL fwfft( 'Wave', psic, dfft )
    !
    IF (gamma_only) THEN
      IF (dim2==1) CALL fftx_psi2c_gamma( dfft, psic, f_out(:,1) )
      IF (dim2==2) CALL fftx_psi2c_gamma( dfft, psic, f_out(:,1), f_out(:,2) )
    ELSE
      CALL fftx_psi2c_k( dfft, psic, f_out(:,1), igk )
    ENDIF
    !
    DEALLOCATE( psic )
    !
    RETURN
    !
  END SUBROUTINE wave_r2g
  !
  !
  !----------------------------------------------------------------------
  SUBROUTINE wave_g2r( f_in, f_out, dfft, igk, howmany_set )
    !--------------------------------------------------------------------
    !! Wave function FFT from G to R-space.
    !
    USE fft_helper_subroutines, ONLY: c2psi_gamma, c2psi_k
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP) :: f_in(:,:)
    COMPLEX(DP) :: f_out(:)
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(3)
    !
    INTEGER :: npw, dim2
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !
    npw  = SIZE(f_in(:,1))
    dim2 = SIZE(f_in(1,:))
    !
    IF (gamma_only) THEN
      IF (dim2/=2) CALL c2psi_gamma( dfft, f_out, f_in(:,1) )
      IF (dim2==2) CALL c2psi_gamma( dfft, f_out, f_in(:,1), f_in(:,2) )
    ELSE
      !$acc data present_or_copyin(igk)
      IF (PRESENT(howmany_set)) THEN     !only when ACC is active
        npw = howmany_set(3)
        CALL c2psi_k( dfft, f_out, f_in, igk, npw, howmany_set )
      ELSE
        CALL c2psi_k( dfft, f_out, f_in, igk, npw )
      ENDIF
      !$acc end data
    ENDIF
    !
    !$acc host_data use_device( f_out )
    IF (PRESENT(howmany_set)) THEN
      CALL invfft( 'Wave', f_out, dfft, howmany=howmany_set(1) )
    ELSE
      CALL invfft( 'Wave', f_out, dfft )
    ENDIF
    !$acc end host_data
    !
    !$acc end data
    !
    RETURN
    !
  END SUBROUTINE wave_g2r
  !
  !
  !----------------------------------------------------------------------
  SUBROUTINE tgwave_g2r( f_in, f_out, dfft, ibnd, ibnd_end, igk )
    !--------------------------------------------------------------------
    !! Wave function FFT from G to R-space. Task-group version.
    !
    USE fft_helper_subroutines,  ONLY: c2psi_gamma_tg, c2psi_k_tg
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP) :: f_in(:,:)
    COMPLEX(DP) :: f_out(:)
    INTEGER, INTENT(IN) :: ibnd, ibnd_end
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !
    INTEGER :: npw
    !
    !$acc data present_or_copyin(f_in,igk) present_or_copyout(f_out)
    !
    npw = SIZE(f_in(:,1))
    !
    !$acc kernels
    f_out(:) = (0.D0,0.D0)
    !$acc end kernels
    !
    IF (gamma_only) THEN
      CALL c2psi_gamma_tg( dfft, f_out, f_in, ibnd, ibnd_end )
    ELSE
      CALL c2psi_k_tg( dfft, f_out, f_in, igk, npw, ibnd, ibnd_end )
    ENDIF
    !
    !$acc host_data use_device(f_out)
    CALL invfft( 'tgWave', f_out, dfft )
    !$acc end host_data
    !
    !$acc end data
    !
    RETURN
    !
  END SUBROUTINE tgwave_g2r
  !
  !
  !----------------------------------------------------------------------
  SUBROUTINE tgwave_r2g( f_in, f_out, dfft, n, ibnd, ibnd_end, igk )
    !--------------------------------------------------------------------
    !! Wave function FFT from R to G-space. Task-group version.
    !
    USE fft_helper_subroutines,  ONLY: psi2c_gamma_tg, psi2c_k_tg
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP), INTENT(IN)  :: f_in(:)
    COMPLEX(DP), INTENT(OUT) :: f_out(:,:)
    INTEGER, INTENT(IN) :: ibnd, ibnd_end, n
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !
    COMPLEX(DP), ALLOCATABLE :: psic(:)
    INTEGER :: nrxxs
    !
    nrxxs = SIZE(f_in)
    !
    ALLOCATE( psic(nrxxs) )
    psic = f_in
    !
    CALL fwfft( 'tgWave', psic, dfft )
    !
    IF (gamma_only) THEN
      CALL psi2c_gamma_tg( dfft, psic, f_out, n, ibnd, ibnd_end )
    ELSE
      CALL psi2c_k_tg( dfft, psic, f_out, igk, n, ibnd, ibnd_end )
    ENDIF
    !
    DEALLOCATE( psic )
    !
    RETURN
    !
  END SUBROUTINE tgwave_r2g
  !
  !
END MODULE fft_wave
