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
  USE control_flags,   ONLY: gamma_only, many_fft
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wave_fft_init, wave_fft_finalize
  PUBLIC :: wave_r2g, wave_g2r, tgwave_r2g, tgwave_g2r
  !
  ! ... workspace arrays
  !
  COMPLEX(DP), ALLOCATABLE :: wpsic(:), wpsic_m(:)
  !! wave-FFT workspace.
  COMPLEX(DP), ALLOCATABLE :: tg_wpsic(:)
  !! task-group wave-FFT workspace.
  !
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE wave_fft_init( dfft, do_many_fft )
    !----------------------------------------------------------------------
    !! Allocation of wave-FFT workspace array.
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    LOGICAL, INTENT(IN), OPTIONAL :: do_many_fft
    !
    INTEGER :: incr
    !
    IF (dfft%has_task_groups .AND. .NOT.ALLOCATED(tg_wpsic)) THEN
      ALLOCATE( tg_wpsic(dfft%nnr_tg) )
      !$acc enter data create(tg_wpsic)
    ENDIF
    !
    IF (.NOT. ALLOCATED(wpsic)) THEN
      ALLOCATE( wpsic(dfft%nnr) )
      !$acc enter data create(wpsic)
    ENDIF
    !
#if defined(__CUDA)
    IF ((.NOT.ALLOCATED(wpsic_m)).AND.PRESENT(do_many_fft)) THEN
      IF (do_many_fft) THEN
        incr = many_fft
        IF (gamma_only) incr = 2*many_fft
        ALLOCATE( wpsic_m(dfft%nnr*incr) )
        !$acc enter data create(wpsic_m)
      ENDIF
    ENDIF
#endif
    !
    RETURN
    !
  END SUBROUTINE wave_fft_init
  !
  !---------------------------------------------------------------------------
  SUBROUTINE wave_fft_finalize( dfft )
    !------------------------------------------------------------------------
    !! Deallocation of wave-FFT workspace array.
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    !
    IF (dfft%has_task_groups .AND. ALLOCATED(tg_wpsic)) THEN
      !$acc exit data delete(tg_wpsic)
      DEALLOCATE( tg_wpsic )
    ENDIF
    !
    IF (ALLOCATED(wpsic)) THEN
      !$acc exit data delete(wpsic)
      DEALLOCATE( wpsic )
    ENDIF
    !
    IF (ALLOCATED(wpsic_m)) THEN
      !$acc exit data delete(wpsic_m)
      DEALLOCATE( wpsic_m )
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE wave_fft_finalize
  !
  !----------------------------------------------------------------------
  SUBROUTINE wave_r2g( f_in, f_out, dfft, igk, howmany_set )
    !--------------------------------------------------------------------
    !! Wave function FFT from R to G-space.
    !
    USE fft_helper_subroutines,  ONLY: fftx_psi2c_gamma, fftx_psi2c_k
    USE control_flags,           ONLY: many_fft
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP), INTENT(IN)  :: f_in(:)
    COMPLEX(DP), INTENT(OUT) :: f_out(:,:)
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(2)
    !
    INTEGER :: dim1, dim2
    !
    dim1 = SIZE(f_in(:))
    dim2 = SIZE(f_out(1,:))
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    IF (PRESENT(howmany_set)) THEN
      !$acc kernels
      wpsic_m(1:dim1) = f_in
      !$acc end kernels
    ELSE
      !$acc kernels
      wpsic(1:dim1) = f_in
      !$acc end kernels
    ENDIF
    !
    IF (PRESENT(howmany_set)) THEN
      !$acc host_data use_device(wpsic_m)
      CALL fwfft( 'Wave', wpsic_m, dfft, howmany=howmany_set(1) )
      !$acc end host_data
    ELSE
      !$acc host_data use_device(wpsic)
      CALL fwfft( 'Wave', wpsic, dfft )
      !$acc end host_data
    ENDIF
    !
    IF (gamma_only) THEN
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_psi2c_gamma( dfft, wpsic_m, f_out, howmany_set=howmany_set )
      ELSE
        IF (dim2==1) CALL fftx_psi2c_gamma( dfft, wpsic, f_out(:,1:1) )
        IF (dim2==2) CALL fftx_psi2c_gamma( dfft, wpsic, f_out(:,1:1), &
                                                         vout2=f_out(:,2) )             
      ENDIF
    ELSE
      !$acc data present_or_copyin(igk)
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_psi2c_k( dfft, wpsic_m, f_out, igk, howmany_set )
      ELSE
        CALL fftx_psi2c_k( dfft, wpsic, f_out(:,1:1), igk )
      ENDIF
      !$acc end data
    ENDIF
    !
    !$acc end data
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
    USE fft_helper_subroutines, ONLY: fftx_c2psi_gamma, fftx_c2psi_k
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP), INTENT(IN) :: f_in(:,:)
    COMPLEX(DP) :: f_out(:)
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(2)
    !
    INTEGER :: npw, dim2
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !
    npw  = SIZE(f_in(:,1))
    dim2 = SIZE(f_in(1,:))
    !
    IF (gamma_only) THEN
      !
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_c2psi_gamma( dfft, f_out, f_in, howmany_set=howmany_set )
      ELSE
        IF (dim2/=2) CALL fftx_c2psi_gamma( dfft, f_out, f_in(:,1:1) )
        IF (dim2==2) CALL fftx_c2psi_gamma( dfft, f_out, f_in(:,1:1), ca=f_in(:,2) )
      ENDIF
      !
    ELSE
      !
      !$acc data present_or_copyin(igk)
      IF (PRESENT(howmany_set)) THEN     !only when ACC is active
        npw = howmany_set(2)
        CALL fftx_c2psi_k( dfft, f_out, f_in, igk, npw, howmany_set(1) )
      ELSE
        CALL fftx_c2psi_k( dfft, f_out, f_in, igk, npw )
      ENDIF
      !$acc end data
      !
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
  SUBROUTINE tgwave_g2r( f_in, f_out, dfft, n, igk )
    !--------------------------------------------------------------------
    !! Wave function FFT from G to R-space. Task-group version.
    !
    USE fft_helper_subroutines,  ONLY: fftx_c2psi_gamma_tg, fftx_c2psi_k_tg
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP) :: f_in(:,:)
    COMPLEX(DP) :: f_out(:)
    INTEGER, INTENT(IN) :: n
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !
    INTEGER :: npw, dbnd
    !
    !$acc data present_or_copyin(f_in,igk) present_or_copyout(f_out)
    !
    npw = SIZE(f_in(:,1))
    dbnd = SIZE(f_in(1,:))
    IF (n/=npw) npw = n
    !
    !$acc kernels
    f_out(:) = (0.D0,0.D0)
    !$acc end kernels
    !
    IF (gamma_only) THEN
      CALL fftx_c2psi_gamma_tg( dfft, f_out, f_in, npw, dbnd )
    ELSE
      CALL fftx_c2psi_k_tg( dfft, f_out, f_in, igk, npw, dbnd )
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
  SUBROUTINE tgwave_r2g( f_in, f_out, dfft, n, igk )
    !--------------------------------------------------------------------
    !! Wave function FFT from R to G-space. Task-group version.
    !
    USE fft_helper_subroutines,  ONLY: fftx_psi2c_gamma_tg, fftx_psi2c_k_tg
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    COMPLEX(DP), INTENT(IN)  :: f_in(:)
    COMPLEX(DP), INTENT(OUT) :: f_out(:,:)
    INTEGER, INTENT(IN) :: n
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !
    INTEGER :: dbnd
    !
    dbnd = SIZE(f_out(1,:))
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !$acc kernels
    tg_wpsic = f_in
    !$acc end kernels
    !
    !$acc host_data use_device(tg_wpsic)
    CALL fwfft( 'tgWave', tg_wpsic, dfft )
    !$acc end host_data
    !
    IF (gamma_only) THEN
      CALL fftx_psi2c_gamma_tg( dfft, tg_wpsic, f_out, n, dbnd )
    ELSE
      CALL fftx_psi2c_k_tg( dfft, tg_wpsic, f_out, igk, n, dbnd )
    ENDIF
    !
    !$acc end data
    !
    RETURN
    !
  END SUBROUTINE tgwave_r2g
  !
  !
END MODULE fft_wave
