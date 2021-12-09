!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------!
! FFT scalar drivers Module - contains machine-dependent routines for      !
! internal FFTW, FFTW v.3, IBM ESSL, Intel DFTI
! (both 3d for serial execution and 1d+2d FFTs for parallel execution);    !
! legacy NEC ASL libraries (3d only, no parallel execution)                !
! CUDA FFT for NVidiia GPUs
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini, Jason Wood                !
! Last update Feb 2021
!--------------------------------------------------------------------------!

!=----------------------------------------------------------------------=!
   MODULE fft_scalar
!=----------------------------------------------------------------------=!

     USE fft_param
#if defined(__FFTW3)
     USE fft_scalar_fftw3
#elif defined(__DFTI)
     USE fft_scalar_dfti
#elif defined(__LINUX_ESSL)
     USE fft_scalar_essl
#elif defined(__SX6)
     USE fft_scalar_sx6
#elif defined(__FFTW)
     USE fft_scalar_fftw
#else
#error No fft_scalar backend selected!
#endif
#if defined(__CUDA)
     USE fft_scalar_cuFFT
#endif
     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds
#if defined(__CUDA)
     PUBLIC :: cft_1z_gpu, cft_2xy_gpu, cfft3d_gpu, cfft3ds_gpu
#endif

   END MODULE fft_scalar
