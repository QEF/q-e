!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE ppcg_param

#if defined(__MPI)
#if defined(__MPI_MODULE)
  USE mpi
#else
  INCLUDE 'mpif.h'
#endif
#endif
  
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: stdout = 6    ! unit connected to standard output

  LOGICAL :: gamma_only =.false.
  LOGICAL :: use_para_diag =.false.

END MODULE ppcg_param
