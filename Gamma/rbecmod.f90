!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE rbecmod
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  REAL(KIND=DP), ALLOCATABLE :: &
       becp(:,:)    !  contains products of wavefunctions and beta
  !     
END MODULE rbecmod

