!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine tra_write_matrix (alpha, adyn, u, nat)
  !-----------------------------------------------------------------------
  !! This routine writes on output the symmetrized dynamical matrix in
  !! cartesian coordinates. The input matrix \(\text{adyn}\) is in the
  !! basis of the modes. On output \(\text{adyn}\) is unchanged.
  !
  USE io_global,    ONLY : stdout
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at, bg
  USE symm_base,    ONLY : s, irt, invs
  USE lr_symm_base, ONLY : rtau, nsymq, irotmq, minus_q
  USE qpoint,       ONLY : xq

  implicit none
  
  integer :: i, j, na, nb, nat
  complex(DP) :: adyn(3*nat,3*nat), u(3*nat,3*nat)
  complex(DP) :: auxdyn(3*nat,3*nat)
  character (len=*) :: alpha

  auxdyn=adyn
  CALL symdyn_munu_new (auxdyn, u, xq, s, invs, rtau, irt, at, bg, &
          nsymq, nat, irotmq, minus_q)

  WRITE( stdout, '(a)') alpha
  do na = 1, nat
     do nb = 1, nat
        WRITE( stdout, '(2i4)') na, nb
        do i = 1, 3
           WRITE( stdout, '(6f12.7)') (auxdyn(3*(na-1)+i, 3*(nb-1)+j),j=1,3)
        enddo
     enddo
  enddo
  return
end subroutine tra_write_matrix

!-----------------------------------------------------------------------
subroutine tra_write_matrix_no_sym (alpha, adyn, nat)
  !-----------------------------------------------------------------------
  !! This routine writes on output the symmetrized dynamical matrix in
  !! cartesian coordinates. The input matrix \(\text{adyn}\) is in the 
  !! basis of the modes. On output \(\text{adyn}\) is unchanged.
  
  USE io_global,    ONLY : stdout
  USE kinds,        ONLY : DP

  implicit none
  
  integer :: nat
  complex(DP) :: adyn(3*nat,3*nat)
  character (len=*) :: alpha
  
  ! ... local variables
  
  integer :: i, j, na, nb

  WRITE( stdout, '(a)') alpha
  do na = 1, nat
     do nb = 1, nat
        WRITE( stdout, '(2i4)') na, nb
        do i = 1, 3
           WRITE( stdout, '(6f12.7)') (adyn(3*(na-1)+i, 3*(nb-1)+j),j=1,3)
        enddo
     enddo
  enddo
  return
end subroutine tra_write_matrix_no_sym
