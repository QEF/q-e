!
! Copyright (C) 2006-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! FFT data Module.
! Written by Carlo Cavazzoni
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE fft_base
!=----------------------------------------------------------------------=!

        USE parallel_include

        USE fft_types, ONLY: fft_dlay_descriptor

        IMPLICIT NONE

        ! ... data structure containing all information
        ! ... about fft data distribution for a given
        ! ... potential grid, and its wave functions sub-grid.

        TYPE ( fft_dlay_descriptor ) :: dfftp ! descriptor for dense grid
             !  Dimensions of the 3D real and reciprocal space FFT grid
             !  relative to the charge density and potential ("dense" grid)
        TYPE ( fft_dlay_descriptor ) :: dffts ! descriptor for smooth grid
             !  Dimensions of the 3D real and reciprocal space
             !  FFT grid relative to the smooth part of the charge density
             !  (may differ from the full charge density grid for USPP )
        TYPE ( fft_dlay_descriptor ) :: dfftb ! descriptor for box grids
             !  Dimensions of the 3D real and reciprocal space
             !  FFT grid relative to the "small box" computation
             !  of the atomic augmentation part of the 
             !  charge density used in USPP (to speed up CPV iterations)
        TYPE ( fft_dlay_descriptor ) :: dfft3d

        SAVE

        PRIVATE

        PUBLIC :: dfftp, dffts, dfftb, dfft3d, fft_dlay_descriptor

!=----------------------------------------------------------------------=!
   END MODULE fft_base
!=----------------------------------------------------------------------=!
