!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine poolscatter (nsize, nkstot, f_in, nks, f_out)
  !-----------------------------------------------------------------------
#include "machine.h"
#ifdef __PARA
  use para
  use parameters, only : DP
  implicit none
  !
  !    This routine scatters a quantity (typically the eigenvalues)
  !    among the pools. On input, f_in is required only on the
  !    first node of the first pool. f_in and f_out may coincide.
  !    Not a smart implementation!
  !
  integer :: nsize, nkstot, nks
  ! first dimension of vectors f_in and f_out
  ! number of k-points per pool
  ! total number of k-points
  real (kind=DP) :: f_in (nsize, nkstot), f_out (nsize, nks)
  ! input (contains values for all k-poi
  ! output(only for k-points of mypool)


  integer :: rest, nbase
  ! the rest of the integer division nkstot/npo
  ! the position in the original list
  !
  ! copy from the first node of the first pool
  !        to the first node of all the other pools
  !
  if (me.eq.1) call poolbcast (nsize * nkstot, f_in)
  !
  ! distribute the vector on the first node of each pool
  !
  rest = nkstot/kunit - (nkstot / kunit / npool) * npool 
  nbase = nks * (mypool - 1)
  if (mypool.gt.rest) nbase = nbase + rest*kunit
  call DCOPY (nsize * nks, f_in (1, nbase+1), 1, f_out, 1)
  !
  ! copy from the first node  of every pool
  !        to the other nodes of every pool
  !
  call broadcast (nsize * nks, f_out)
#endif
  return


end subroutine poolscatter

