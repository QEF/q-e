!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------


subroutine extreme (ps, iflag)
  !-----------------------------------------------------------------------
  !
  !     Finds the maximum (iflag.gt.0) or the minimum (iflag.le.0) value o
  !     a real variable among the values distributed on a given pool
  !
  !-----------------------------------------------------------------------
#include "machine.h"
#ifdef PARA
  use para
  use parameters, only : DP
  implicit none
  include 'mpif.h'

  real (kind=DP) :: ps, psr

  integer :: iflag, info
  call mpi_barrier (MPI_COMM_WORLD, info)
  if (iflag.gt.0) then
     call mpi_allreduce (ps, psr, 1, MPI_REAL8, MPI_MAX, &
          MPI_COMM_WORLD, info)
  else
     call mpi_allreduce (ps, psr, 1, MPI_REAL8, MPI_MIN, &
          MPI_COMM_WORLD, info)

  endif

  ps = psr
#endif
  return

end subroutine extreme

