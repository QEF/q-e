!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!-----------------------------------------------------------------------
subroutine tra_write_matrix (alpha, adyn, u, nat)
  !-----------------------------------------------------------------------
#include "machine.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !
  !    This routine writes on output the dynamical matrix in the cartesian
  !    It transform it from the basis of the modes but on output adyn is
  !    unchanged
  !
  integer :: i, j, na, nb, nat
  integer :: icart, jcart, mu, nu
  complex(kind=DP) :: adyn (3 * nat, 3 * nat), u (3 * nat, 3 * nat)
  complex(kind=DP) :: wdyn (3, 3, nat, nat), work
  character (len=*) :: alpha
  WRITE( stdout, * ) nat
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icart = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcart = j - 3 * (nb - 1)
        work = (0.d0, 0.d0)
        do mu = 1, 3 * nat
           do nu = 1, 3 * nat
              work = work + u (i, mu) * adyn (mu, nu) * conjg (u (j, nu) )
           enddo
        enddo
        wdyn (icart, jcart, na, nb) = work
     enddo


  enddo
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
end subroutine tra_write_matrix
