!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine poolreduce (size, ps)
  !-----------------------------------------------------------------------
  !
  !     Sums a distributed variable ps(size) over the pools.
  !     This MPI-only version uses a fixed-length buffer
  !
  !-----------------------------------------------------------------------
#include "machine.h"
#ifdef __PARA
  use para
#endif
  use parameters, only : DP
  implicit none
  integer :: size

  real (kind=DP) :: ps (size)
#ifdef __PARA

  include 'mpif.h'
  integer :: MAXB
  parameter (MAXB = 10000)
  real (8) :: buff (MAXB)

  integer :: info, nbuf, n
  if (size.le.0.or.npool.le.1) return
  call start_clock ('poolreduce')
  !
  !  MPI syncronize processes
  !

  call mpi_barrier (MPI_COMM_WORLD, info)

  call errore ('poolreduce', 'info<>0 at barrier', info)
  nbuf = size / MAXB
  do n = 1, nbuf
     call mpi_allreduce (ps (1 + (n - 1) * MAXB), buff, MAXB, &
          MPI_REAL8, MPI_SUM, MPI_COMM_ROW, info)
     call errore ('poolreduce', 'info<>0 at allreduce1', info)
     call DCOPY (MAXB, buff, 1, ps (1 + (n - 1) * MAXB), 1)
  enddo
  if (size-nbuf * MAXB.gt.0) then
     call mpi_allreduce (ps (1 + nbuf * MAXB), buff, size-nbuf * &
          MAXB, MPI_REAL8, MPI_SUM, MPI_COMM_ROW, info)
     call errore ('poolreduce', 'info<>0 at allreduce2', info)
     call DCOPY (size-nbuf * MAXB, buff, 1, ps (1 + nbuf * MAXB), &
          1)

  endif
  call stop_clock ('poolreduce')
#endif
  return

end subroutine poolreduce

