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
  !! Write out \(\text{wdyn}\) matrix.
  !
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  !
  implicit none
  !
  character (len=*) :: alpha
  complex(DP) :: wdyn(3,3,nat,nat)
  integer :: nat
  !
  ! ... local variables
  !
  integer :: i, j, na, nb
  !
  WRITE( stdout, '(a)') alpha
  do na = 1, nat
     do nb = 1, nat
        WRITE( stdout, '(2i4)') na, nb
        do i = 1, 3
           WRITE( stdout, '(6f10.5)') (wdyn (i, j, na, nb) , j = 1, 3)
        enddo
     enddo

  enddo
  return
end subroutine write_matrix
