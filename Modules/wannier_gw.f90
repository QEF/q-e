!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE wannier_gw
  !
  ! ... The variables needed for gww-gwl code (head.x)
  !
  USE kinds, ONLY: DP
  !
  SAVE
  !
  LOGICAL :: l_head!if true calculates the head of the symmetrized dielectric matrix -1
  INTEGER :: n_gauss!number of frequency steps for head calculation
  REAL(kind=DP) :: omega_gauss!period for frequency calculation
  INTEGER :: grid_type!0 GL -T,T 2 GL 0 T 3 Equally spaced 0 Omega
  INTEGER :: nsteps_lanczos!number of lanczos steps
    !options for grid_freq=5
  INTEGER :: second_grid_n!sub spacing for second grid
  INTEGER :: second_grid_i!max regular step using the second grid
  LOGICAL :: l_scissor!if true displaces occupied manifold of scissor
  REAL(kind=DP) :: scissor!see above

END MODULE wannier_gw
