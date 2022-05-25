!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE fft_param
  use iso_fortran_env, only : stderr=>ERROR_UNIT, stdout=>OUTPUT_UNIT
#if defined(__MPI)
#if defined(__MPI_MODULE)
  USE mpi
#else
  INCLUDE 'mpif.h'
#endif
#else
  INTEGER, PARAMETER :: MPI_COMM_WORLD =  0
  INTEGER, PARAMETER :: MPI_COMM_NULL  = -1
  INTEGER, PARAMETER :: MPI_COMM_SELF  = -2
#endif
  
  INTEGER, PARAMETER :: ndims = 20
  !! Number of different FFT tables that the module
  !!could keep into memory without reinitialization

  INTEGER, PARAMETER :: nfftx = 16385
  !!Max allowed fft dimension

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

  REAL(DP), PARAMETER :: eps8  = 1.0E-8_DP

END MODULE fft_param
