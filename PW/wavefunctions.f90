!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
!
MODULE wavefunctions
  USE parameters, ONLY :  DP
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE, TARGET :: &
       evc(:,:)     ! wavefunctions in the PW basis
  !
  COMPLEX(KIND=DP) , ALLOCATABLE, TARGET :: &
       psic(:)      ! additional memory for FFT
  !
END MODULE wavefunctions
