!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
!
module wavefunctions
  use parameters, only: DP
  !
  complex(kind=DP) , allocatable, target :: &
       evc(:,:)     ! wavefunctions in the PW basis
  !
  complex(kind=DP) , allocatable, target :: &
       psic(:)     ! additional memory for FFT
end module wavefunctions
