!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine divide_et_impera (xk, wk, isk, lsda, nkstot, nks)
  !-----------------------------------------------------------------------
#include "machine.h"
#ifdef __PARA
  USE kinds, only : DP
  use para
  implicit none
  !
  !    This routine divides the k points (with granularity kunit) among no
  !    and sets the variable nkstot equal to the total number of k-points
  !

  integer :: nkstot, nks, ik, isk (nks), rest, nbase
  ! total number of k-points
  ! number of k-points per pool
  ! counter on kpoints
  ! spin index of each kpoint (when lsda=.t.)
  ! the rest of the integer division nkstot/npo
  ! the position in the original list of the fi
  ! point that belongs to this pool - 1

  logical :: lsda
  ! logical for local spin density approx.
  real (kind=DP) :: xk (3,nks), wk (nks)
  ! k-points
  ! k-point weights
  !
  if (mod (nks, kunit) .ne.0) call errore ('d_&_i', &
       ' nks/kunit is not an integer', nks)
  !
  nkstot = nks
  nks = kunit * (nkstot / kunit / npool)
  if (nks.eq.0) call errore ('d_&_i', ' nks = 0 for some nodes', 1)
  rest = (nkstot - nks * npool) / kunit
  if (mypool.le.rest) nks = nks + kunit
  !
  ! calculates nbase = the position in the list of the first point that
  !                   belong to this npool - 1
  !
  nbase = nks * (mypool - 1)
  if (mypool.gt.rest) nbase = nbase+rest * kunit
  !
  ! displaces these points in the first positions of the list
  !
  if (nbase.gt.0) then
     do ik = 1, nks
        call DCOPY (3, xk (1, nbase+ik), 1, xk (1, ik), 1)
        wk (ik) = wk (nbase+ik)
     enddo
     if (lsda) then
        do ik = 1, nks
           isk (ik) = isk (nbase+ik)
        enddo
     endif

  endif
#endif
  return
end subroutine divide_et_impera

