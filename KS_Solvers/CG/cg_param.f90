!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE cg_param

  USE parallel_include
  
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: stdout = 6    ! unit connected to standard output

  LOGICAL :: gamma_only =.false.
  LOGICAL :: use_para_diag =.false.

END MODULE cg_param
