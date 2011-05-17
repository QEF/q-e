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
  
  CONTAINS
 
    LOGICAL FUNCTION init_start_grid ( nk1_, nk2_, nk3_, k1_, k2_, k3_ ) 
       INTEGER, INTENT (IN) :: nk1_, nk2_, nk3_, k1_, k2_, k3_
       !
       IF ( nk1_ > 0 ) nk1 = nk1_
       IF ( nk2_ > 0 ) nk2 = nk2_
       IF ( nk3_ > 0 ) nk3 = nk3_
       IF (  k1_ > 0 )  k1 = k1_
       IF (  k2_ > 0 )  k2 = k2_
       IF (  k3_ > 0 )  k3 = k3_
       init_start_grid = (nk1_*nk2_*nk3_ > 0)
       !
    END FUNCTION init_start_grid

END MODULE start_k
