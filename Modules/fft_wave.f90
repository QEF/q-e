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
  !! This module contains wrappers to FFT and inverse FFTs of the wave function,
  !! which it enclose the calls to g-vect/FFT-grid transposition routines too.
  !
  USE kinds,                  ONLY: DP
  USE fft_interfaces,         ONLY: fwfft, invfft
#if defined(__OPENMP_GPU)
  USE fft_interfaces,         ONLY: fwfft_y_omp, invfft_y_omp
  USE fft_helper_subroutines, ONLY: fftx_psi2c_k_omp, fftx_c2psi_k_omp
#endif
  USE fft_types,              ONLY: fft_type_descriptor
  USE control_flags,          ONLY: gamma_only, many_fft
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wave_r2g, wave_g2r, tgwave_r2g, tgwave_g2r
  !
  !
CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE wave_r2g( f_in, f_out, dfft, igk, howmany_set, omp_mod )
    !--------------------------------------------------------------------
    !! Wave function FFT from R to G-space.
    !
    USE fft_helper_subroutines,  ONLY: fftx_psi2c_gamma, fftx_psi2c_k, &
                                       fftx_psi2c_gamma_omp, fftx_psi2c_k_omp
    
    USE control_flags,           ONLY: many_fft
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    !! FFT descriptor
    COMPLEX(DP) :: f_in(:)
    !! input: r-space wave-function
    COMPLEX(DP), INTENT(OUT) :: f_out(:,:)
    !! output: g-space wave-function
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !! index of G corresponding to a given index of k+G
    INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(3)
    !! gpu-enabled only (many_fft>1 case):  
    !! (1) group_size;  
    !! (2) true dimension of psi;  
    !! (3) howmany (if gamma) or group_size (if k).
    INTEGER, OPTIONAL, INTENT(IN) :: omp_mod
    !! whether to execute the FFT on GPU with OMP5 offload
    !
    ! ... local variables
    INTEGER :: dim1, dim2
    LOGICAL :: omp_offload, omp_map
    !
    dim1 = SIZE(f_in(:))
    dim2 = SIZE(f_out(1,:))
    !
    omp_offload = .FALSE.
    omp_map     = .FALSE.
    IF (PRESENT(omp_mod)) THEN
      omp_offload = omp_mod>=0 ! run FFT on device (data already mapped)
      omp_map     = omp_mod>=1 ! map data and run FFT on device
    ENDIF 
    IF(omp_offload.AND.PRESENT(howmany_set)) CALL errore( 'wave_r2g','omp_offload &
                                                          &and many FFT NYI', 1 )
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !
    !$acc host_data use_device(f_in)
    IF (PRESENT(howmany_set)) THEN
      CALL fwfft( 'Wave', f_in, dfft, howmany=howmany_set(3) )
    ELSE
      IF(omp_offload) THEN
#if defined (__OPENMP_GPU)
        IF(omp_map) THEN
          !$omp target data map(tofrom:f_in)
          CALL fwfft_y_omp( 'Wave', f_in, dfft )
          !$omp end target data 
        ELSE
          CALL fwfft_y_omp( 'Wave', f_in, dfft )
        ENDIF
#endif
      ELSE
        CALL fwfft( 'Wave', f_in, dfft )
      ENDIF 
    ENDIF
    !$acc end host_data
    !
    IF (gamma_only) THEN
      !
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_psi2c_gamma( dfft, f_in, f_out, howmany_set=howmany_set(1:2) )
      ELSE
        IF (omp_offload) THEN
#if defined (__OPENMP_GPU)
          IF (omp_map) THEN
            !$omp target data map(to:f_in) map(from:f_out)
            IF (dim2==1) CALL fftx_psi2c_gamma_omp( dfft, f_in, f_out(:,1:1) )
            IF (dim2==2) CALL fftx_psi2c_gamma_omp( dfft, f_in, f_out(:,1:1), &
                                                          vout2=f_out(:,2) )
            !$omp end target data
          ELSE
            IF (dim2==1) CALL fftx_psi2c_gamma_omp( dfft, f_in, f_out(:,1:1) )
            IF (dim2==2) CALL fftx_psi2c_gamma_omp( dfft, f_in, f_out(:,1:1), &
                                                            vout2=f_out(:,2) )
          ENDIF
#endif
        ELSE
          IF (dim2==1) CALL fftx_psi2c_gamma( dfft, f_in, f_out(:,1:1) )
          IF (dim2==2) CALL fftx_psi2c_gamma( dfft, f_in, f_out(:,1:1), &
                                                      vout2=f_out(:,2) )
        ENDIF
      ENDIF
      !
    ELSE
      !$acc data present_or_copyin(igk)
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_psi2c_k( dfft, f_in, f_out, igk, howmany_set(1:2) )
      ELSE
        IF(omp_offload) THEN
#if defined (__OPENMP_GPU)
          IF(omp_map) THEN
            !$omp target data map(to:f_in,igk) map(tofrom:f_out)
            CALL fftx_psi2c_k_omp( dfft, f_in, f_out(:,1:1), igk )
            !$omp end target data
          ELSE
            CALL fftx_psi2c_k_omp( dfft, f_in, f_out(:,1:1), igk )
          ENDIF
#endif
        ELSE
          CALL fftx_psi2c_k( dfft, f_in, f_out(:,1:1), igk )
        ENDIF 
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
  SUBROUTINE wave_g2r( f_in, f_out, dfft, igk, howmany_set, omp_mod )
    !--------------------------------------------------------------------
    !! Wave function FFT from G to R-space.
    !
    USE fft_helper_subroutines, ONLY: fftx_c2psi_gamma, fftx_c2psi_k, &
                                      fftx_c2psi_gamma_omp, fftx_c2psi_k_omp
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: dfft
    !! FFT wave descriptor
    COMPLEX(DP), INTENT(IN) :: f_in(:,:)
    !! input: g-space wave-function
    COMPLEX(DP) :: f_out(:)
    !! output: r-space wave-function
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !! index of G corresponding to a given index of k+G
    INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(3)
    !! gpu-enabled only (many_fft>1 case):  
    !! (1) group_size;  
    !! (2) true dimension of psi;  
    !! (3) howmany (if gamma) or group_size (if k).
    INTEGER, OPTIONAL, INTENT(IN) :: omp_mod
    !! whether to execute the FFT on GPU with OMP5 offload
    !
    ! ... local variables
    INTEGER :: npw, dim2
    LOGICAL :: omp_offload, omp_map
    !
    omp_offload = .FALSE.
    omp_map     = .FALSE.
    IF (PRESENT(omp_mod)) THEN
      omp_offload = omp_mod>=0 ! run FFT on device (data already mapped)
      omp_map     = omp_mod>=1 ! map data and run FFT on device
    ENDIF
    IF (omp_offload.AND.PRESENT(howmany_set)) CALL errore('wave_r2g','omp_offload &
                                                           &and many FFT NYI', 1 )
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !
    npw  = SIZE(f_in(:,1))
    dim2 = SIZE(f_in(1,:))
    !
    IF (gamma_only) THEN
      !
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_c2psi_gamma( dfft, f_out, f_in, howmany_set=howmany_set(1:2) )
      ELSE
        IF(omp_offload) THEN
#if defined (__OPENMP_GPU)
          IF(omp_map) THEN
            !$omp target data map(to:f_in) map(from:f_out)
            IF (dim2/=2) CALL fftx_c2psi_gamma_omp( dfft, f_out, f_in(:,1:1) )
            IF (dim2==2) CALL fftx_c2psi_gamma_omp( dfft, f_out, f_in(:,1:1), &
                                                              ca=f_in(:,2) )
            !$omp end target data
          ELSE 
            IF (dim2/=2) CALL fftx_c2psi_gamma_omp( dfft, f_out, f_in(:,1:1) )
            IF (dim2==2) CALL fftx_c2psi_gamma_omp( dfft, f_out, f_in(:,1:1), &
                                                              ca=f_in(:,2) )
          ENDIF
#endif
        ELSE
          IF (dim2/=2) CALL fftx_c2psi_gamma( dfft, f_out, f_in(:,1:1) )
          IF (dim2==2) CALL fftx_c2psi_gamma( dfft, f_out, f_in(:,1:1), ca=f_in(:,2) )
        ENDIF
      ENDIF
      !
    ELSE
      !
      !$acc data present_or_copyin(igk)
      IF (PRESENT(howmany_set)) THEN     !only when ACC is active
        npw = howmany_set(2)
        CALL fftx_c2psi_k( dfft, f_out, f_in, igk, npw, howmany_set(1) )
      ELSE
        IF(omp_offload) THEN 
#if defined (__OPENMP_GPU)
          IF(omp_map) THEN
            !$omp target data map(to:f_in,igk) map(tofrom:f_out) 
            CALL fftx_c2psi_k_omp( dfft, f_out, f_in, igk, npw )
            !$omp end target data
          ELSE
            CALL fftx_c2psi_k_omp( dfft, f_out, f_in, igk, npw )
          ENDIF
#endif
        ELSE
          CALL fftx_c2psi_k( dfft, f_out, f_in, igk, npw )
        ENDIF 
      ENDIF
      !$acc end data
      !
    ENDIF
    !
    !$acc host_data use_device( f_out )
    IF (PRESENT(howmany_set)) THEN
      CALL invfft( 'Wave', f_out, dfft, howmany=howmany_set(3) )
    ELSE
      IF(omp_offload) THEN
#if defined (__OPENMP_GPU)
        IF(omp_map) THEN 
          !$omp target data map(tofrom:f_out)
          CALL invfft_y_omp( 'Wave', f_out, dfft, howmany=1 )
          !$omp end target data 
        ELSE
          CALL invfft_y_omp( 'Wave', f_out, dfft )
        ENDIF
#endif
      ELSE
        CALL invfft( 'Wave', f_out, dfft )
      ENDIF
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
    !! FFT descriptor
    COMPLEX(DP) :: f_in(:,:)
    !! input: wave in g-space - task group chunk
    COMPLEX(DP) :: f_out(:)
    !! output: wave in r-space - task group chunk
    INTEGER, INTENT(IN) :: n
    !! true dimension of f_in
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !! index of G corresponding to a given index of k+G
    !
    ! ... local variables
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
    !! FFT descriptor
    COMPLEX(DP) :: f_in(:)
    !! input: wave in g-space - task group chunk
    COMPLEX(DP), INTENT(OUT) :: f_out(:,:)
    !! output: wave in r-space - task group chunk
    INTEGER, INTENT(IN) :: n
    !! true dimension of f_out
    INTEGER, OPTIONAL, INTENT(IN) :: igk(:)
    !! index of G corresponding to a given index of k+G
    !
    ! ... local variables
    INTEGER :: dbnd
    !
    dbnd = SIZE(f_out(1,:))
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !
    !$acc host_data use_device(f_in)
    CALL fwfft( 'tgWave', f_in, dfft )
    !$acc end host_data
    !
    IF (gamma_only) THEN
      CALL fftx_psi2c_gamma_tg( dfft, f_in, f_out, n, dbnd )
    ELSE
      CALL fftx_psi2c_k_tg( dfft, f_in, f_out, igk, n, dbnd )
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
