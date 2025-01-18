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
  USE kinds,           ONLY: DP
  USE fft_interfaces,  ONLY: fwfft, invfft
  USE fft_types,       ONLY: fft_type_descriptor
  USE control_flags,   ONLY: gamma_only
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wave_r2g, wave_g2r, tgwave_r2g, tgwave_g2r, fwfft_wave, invfft_wave
  !
  !
CONTAINS
  !
  !----------------------------------------------------------------------
  SUBROUTINE wave_r2g( f_in, f_out, dfft, igk, howmany_set )
    !--------------------------------------------------------------------
    !! Wave function FFT from R to G-space.
    !
    USE fft_helper_subroutines,  ONLY: fftx_psi2c_gamma, fftx_psi2c_k
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
    !
    ! ... local variables
    INTEGER :: dim1, dim2
    !
    dim1 = SIZE(f_in(:))
    dim2 = SIZE(f_out(1,:))
    !
    !$acc data present_or_copyin(f_in) present_or_copyout(f_out)
    !
    !$acc host_data use_device(f_in)
    IF (PRESENT(howmany_set)) THEN
      CALL fwfft( 'Wave', f_in, dfft, howmany=howmany_set(3) )
    ELSE
      CALL fwfft( 'Wave', f_in, dfft )
    ENDIF
    !$acc end host_data
    !
    IF (gamma_only) THEN
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_psi2c_gamma( dfft, f_in, f_out, howmany_set=howmany_set(1:2) )
      ELSE
        IF (dim2==1) CALL fftx_psi2c_gamma( dfft, f_in, f_out(:,1:1) )
        IF (dim2==2) CALL fftx_psi2c_gamma( dfft, f_in, f_out(:,1:1), &
                                                         vout2=f_out(:,2) )             
      ENDIF
    ELSE
      !$acc data present_or_copyin(igk)
      IF (PRESENT(howmany_set)) THEN
        CALL fftx_psi2c_k( dfft, f_in, f_out, igk, howmany_set(1:2) )
      ELSE
        CALL fftx_psi2c_k( dfft, f_in, f_out(:,1:1), igk )
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
    !
    ! ... local variables
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
        CALL fftx_c2psi_gamma( dfft, f_out, f_in, howmany_set=howmany_set(1:2) )
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
      CALL invfft( 'Wave', f_out, dfft, howmany=howmany_set(3) )
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
SUBROUTINE fwfft_wave (npwx, npw, igk, evc_g, evc_r)
  !
  !! Wave function FFT from R to G-space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dffts
  USE fft_interfaces,   ONLY : fwfft
  USE noncollin_module, ONLY : noncolin, npol

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npwx, npw, igk(npw)
  COMPLEX(DP), INTENT(INOUT) :: evc_g (npwx*npol), evc_r (dffts%nnr,npol)
  !
  INTEGER :: ig, ik

#if defined(__CUDA) && defined(_OPENACC)
  INTEGER, POINTER, DEVICE :: nl(:)
  nl => dffts%nl_d
#else
  INTEGER, ALLOCATABLE :: nl(:)
  ALLOCATE( nl(dffts%ngm) )
  nl = dffts%nl
#endif

  !$acc host_data use_device(evc_r)
  CALL fwfft ('Wave', evc_r(:,1), dffts)
  !$acc end host_data
  !$acc parallel loop private(ik)
  DO ig = 1, npw
     ik = nl(igk(ig))
     evc_g (ig) = evc_g (ig) + evc_r (ik,1)
  ENDDO
  IF (noncolin) THEN
     !$acc host_data use_device(evc_r)
     CALL fwfft ('Wave', evc_r(:,2), dffts)
     !$acc end host_data
     !$acc parallel loop private(ik) 
     DO ig = 1, npw
        ik = nl(igk(ig))
        evc_g (ig+npwx) = evc_g (ig+npwx) + evc_r (ik,2)
     ENDDO
  ENDIF

#if !defined(__CUDA) || !defined(_OPENACC)
  DEALLOCATE(nl)
#endif
END SUBROUTINE fwfft_wave
  !
SUBROUTINE invfft_wave (npwx, npw, igk, evc_g, evc_r)
  !
  !! Wave function FFT from G to R-space.
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dffts
  USE fft_interfaces,   ONLY : invfft
  USE noncollin_module, ONLY : noncolin, npol

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npwx, npw, igk(npw)
  COMPLEX(DP), INTENT(IN) :: evc_g (npwx*npol)
  COMPLEX(DP), INTENT(OUT):: evc_r (dffts%nnr,npol)
  !
  INTEGER :: ig, ik

#if defined(__CUDA) && defined(_OPENACC)
  INTEGER, POINTER, DEVICE :: nl(:)
  nl => dffts%nl_d
#else
  INTEGER, ALLOCATABLE :: nl(:)
  ALLOCATE( nl(dffts%ngm) )
  nl = dffts%nl
#endif

  !$acc kernels
  evc_r(:,:) = (0.0_dp, 0.0_dp)
  !$acc end kernels
  !$acc parallel loop private(ik)
  DO ig = 1, npw
     ik = nl(igk(ig))
     evc_r (ik, 1) = evc_g (ig)
  ENDDO
  !$acc host_data use_device(evc_r)
  CALL invfft ('Wave', evc_r(:,1), dffts)
  !$acc end host_data
  IF (noncolin) THEN
     !$acc parallel loop private(ik)
     DO ig = 1, npw
        ik = nl(igk(ig))
        evc_r (ik, 2) = evc_g (ig+npwx)
     ENDDO
     !$acc host_data use_device(evc_r)
     CALL invfft ('Wave', evc_r(:,2), dffts)
     !$acc end host_data
  ENDIF

#if !defined(__CUDA) || !defined(_OPENACC)
  DEALLOCATE(nl)
#endif
END SUBROUTINE invfft_wave
  !
END MODULE fft_wave
