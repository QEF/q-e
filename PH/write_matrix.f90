!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine write_matrix (alpha, wdyn, nat)  
  !-----------------------------------------------------------------------
  use parameters, only : DP
  implicit none
  integer :: i, j, na, nb, nat  
  complex(kind=DP) :: wdyn (3, 3, nat, nat)  

  character (len=*) :: alpha  
  write (6, '(a)') alpha  
  do na = 1, nat  
     do nb = 1, nat  
        write (6, '(2i4)') na, nb  
        do i = 1, 3  
           write (6, '(6f10.5)') (wdyn (i, j, na, nb) , j = 1, 3)  
        enddo
     enddo

  enddo
  return  
end subroutine write_matrix
