!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE start_k
  !
  ! ... Basic variables for k-points generations, as read from input
  !
  USE kinds,      ONLY : DP
  !
  SAVE
  !
  ! ... uniform k-point grid parameters
  !
  INTEGER :: &
       nk1, nk2, nk3,   &! the special-point grid
       k1, k2, k3        ! the offset from the origin
  !
  !
  ! ... list of k-points and weights
  !
  INTEGER :: nks_start=0 ! number of initial k points
  REAL(DP), ALLOCATABLE :: wk_start(:)
    ! initial weight of k points
  REAL(DP), ALLOCATABLE :: xk_start(:,:)
    ! initial coordinates of k points

END MODULE start_k
