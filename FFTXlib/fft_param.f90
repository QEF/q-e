!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  
  INTEGER, PARAMETER :: nfftx = 2049
  !! Number of different FFT tables that the module
  !!could keep into memory without reinitialization

  INTEGER, PARAMETER :: ndims = 10
  !!Max allowed fft dimension

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

  INTEGER, PARAMETER :: stdout = 6    ! unit connected to standard output
