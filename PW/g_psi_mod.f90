!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module g_psi_mod
  use parameters, only : DP
  implicit none
  ! These are the variables needed in g_psi
  real(kind=DP), allocatable :: h_diag (:),& ! diagonal part of the Hamiltonian
       s_diag (:)                        ! diagonal part of the overlap matrix
end module g_psi_mod
