!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------


subroutine poolextreme (ps, iflag)
  !-----------------------------------------------------------------------
  !
  !     Finds the maximum (iflag.gt.0) or the minimum (iflag.le.0) value o
  !     a real variable among the values distributed across the pools and
  !     returns this value to all pools.
  !
  !-----------------------------------------------------------------------
#include "machine.h"
#ifdef PARA
  use para
  use parameters, only : DP
  implicit none

  integer :: iflag, info
  real (kind=DP) :: ps, psr


  include 'mpif.h'
  if (npool.le.1) return
  !
  call mpi_barrier (MPI_COMM_WORLD, info)
  call error ('poolextreme', 'info<>0 at barrier', info)
  if (iflag.gt.0) then
     call mpi_allreduce (ps, psr, 1, MPI_REAL8, MPI_MAX, &
          MPI_COMM_ROW, info)
     call error ('poolextreme', 'info<>0 in allreduce1', info)
  else
     call mpi_allreduce (ps, psr, 1, MPI_REAL8, MPI_MIN, &
          MPI_COMM_ROW, info)
     call error ('poolextreme', 'info<>0 in allreduce2', info)
  endif

  ps = psr
#endif
  return

end subroutine poolextreme

