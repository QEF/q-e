!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module becmod
  use parameters, only: DP
  complex(kind=DP) , allocatable ::  becp(:,:)
  ! contains products of wavefunctions and beta
end module becmod
!
