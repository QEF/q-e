!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine estimate (hessm1, nax3, nat, nat3)  
  !-----------------------------------------------------------------------
  !
  use parameters
  implicit none  
  integer :: nax3, nat3, nat, i  
  real(kind=DP) :: hessm1 (nax3, nat3)  

  external setv  
  call setv (nax3 * nat3, 0.d0, hessm1, 1)  
  do i = 1, nat3  
     hessm1 (i, i) = 1.d0  

  enddo
  return  
end subroutine estimate
