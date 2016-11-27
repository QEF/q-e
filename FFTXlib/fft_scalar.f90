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
! internal FFTW, FFTW v.3, IBM ESSL, Intel DFTI, ARMlib                    !
! (both 3d for serial execution and 1d+2d FFTs for parallel execution);    !
! legacy NEC ASL libraries (3d only, no parallel execution)                !
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini - Last update Jul 2015     !
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
#elif defined(__ARM_LIB)
     USE fft_scalar_arm
#else
     USE fft_scalar_fftw
#endif
       
     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds

   END MODULE fft_scalar
