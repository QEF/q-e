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
! FFTW, FFTW3, ESSL (both 3d for serial execution and 1d+2d FFTs for       !
! parallel execution; NEC ASL libraries (3d only, no parallel execution)   !
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini - Last update Jul 2015     !
!--------------------------------------------------------------------------!

#include "fft_defs.h"

#if defined(__FFTW3)

#include "fft_scalar.FFTW3.f90"

#elif defined(__DFTI)

#include "fft_scalar.DFTI.f90"

#elif defined(__LINUX_ESSL)

#include "fft_scalar.ESSL.f90"

#elif defined(__SX6)

#include "fft_scalar.SX6.f90"

#elif defined(__ARM_LIB)

#include "fft_scalar.ARM_LIB.f90"

#else

#include "fft_scalar.FFTW.f90"

#endif
