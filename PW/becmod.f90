!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE becmod
  USE kinds, ONLY :  DP
  !
  SAVE
  ! variables containing <beta|psi>
  REAL(KIND=DP), ALLOCATABLE :: &
       rbecp(:,:) !   <beta|psi> for real (at Gamma) wavefunctions 
  COMPLEX(KIND=DP), ALLOCATABLE ::  &
       becp (:,:) !  as above for complex wavefunctions
  !!!       becp_nc(:,:,:) !  as above for spinors
  !
END MODULE becmod


