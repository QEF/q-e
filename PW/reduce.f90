!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine reduce (size, ps)
  !-----------------------------------------------------------------------
  !
  !     sums a distributed variable s(size) over the processors.
  !     This version uses a fixed-length buffer of appropriate (?) size
  !                  uses shmem for the t3d/t3e, MPI otherwhise
  !
#undef SHMEM
  !
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
  integer :: info, n, nbuf
#define MAXB 10000
  real (kind=DP) :: buff (MAXB)
#ifdef SHMEM
  include 'mpp/shmem.fh'
  integer :: pWrkSync, pWrkData, start
  common / SH_SYNC / pWrkSync (SHMEM_BARRIER_SYNC_SIZE)
  common / SH_DATA / pWrkData (1024 * 1024)
  data pWrkData / 1048576 * 0 /
  data pWrkSync / SHMEM_BARRIER_SYNC_SIZE * SHMEM_SYNC_VALUE /
  !DIR$ CACHE_ALIGN /SH_SYNC/
  !DIR$ CACHE_ALIGN /SH_DATA/
#endif
  if (nprocp.le.1) return
  if (size.le.0) return
  call start_clock ('reduce')
  !
  !  syncronize processes - maybe unneeded on t3d but necessary on t3e !!!
  !
  call mpi_barrier (MPI_COMM_POOL, info)

  call errore ('reduce', 'error in barrier', info)

  nbuf = size / MAXB
#ifdef SHMEM
  start = (mypool - 1) * nprocp
#endif
  do n = 1, nbuf
#ifdef SHMEM
     call SHMEM_REAL8_SUM_TO_ALL (buff, ps (1 + (n - 1) * MAXB), &
          MAXB, start, 0, nprocp, pWrkData, pWrkSync)
#else
     call mpi_allreduce (ps (1 + (n - 1) * MAXB), buff, MAXB, &
          MPI_REAL8, MPI_SUM, MPI_COMM_POOL, info)
     call errore ('reduce', 'error in allreduce1', info)
#endif
     call DCOPY (MAXB, buff, 1, ps (1 + (n - 1) * MAXB), 1)
  enddo
  !
  !    possible remaining elements < maxb
  !
  if (size-nbuf * MAXB.gt.0) then
#ifdef SHMEM
     call SHMEM_REAL8_SUM_TO_ALL (buff, ps (1 + nbuf * MAXB), &
          size-nbuf * MAXB, start, 0, nprocp, pWrkData, pWrkSync)
#else
     call mpi_allreduce (ps (1 + nbuf * MAXB), buff, size-nbuf * &
          MAXB, MPI_REAL8, MPI_SUM, MPI_COMM_POOL, info)
     call errore ('reduce', 'error in allreduce2', info)
#endif
     call DCOPY (size-nbuf * MAXB, buff, 1, ps (1 + nbuf * MAXB), &
          1)

  endif
  call stop_clock ('reduce')
#endif
  return
end subroutine reduce
!
!-----------------------------------------------------------------------
subroutine ireduce (size, is)
  !-----------------------------------------------------------------------
  !
  !     sums a distributed variable is(size) over the processors.
  !
#include "machine.h"
#ifdef __PARA
  use para
#endif
  implicit none

  integer :: size, is (size)
#ifdef __PARA
  include 'mpif.h'
  integer :: info, n, m, nbuf
#define MAXI 500

  integer :: buff (MAXI)
  if (nprocp.le.1) return
  if (size.le.0) return
  !
  !  syncronize processes
  !
  call mpi_barrier (MPI_COMM_POOL, info)

  call errore ('reduce', 'error in barrier', info)

  nbuf = size / MAXI
  do n = 1, nbuf
     call mpi_allreduce (is (1 + (n - 1) * MAXI), buff, MAXI, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_POOL, info)
     call errore ('ireduce', 'error in allreduce 1', info)
     do m = 1, MAXI
        is (m + (n - 1) * MAXI) = buff (m)
     enddo
  enddo
  !
  !    possible remaining elements < MAXI
  !
  if (size-nbuf * MAXI.gt.0) then
     call mpi_allreduce (is (1 + nbuf * MAXI), buff, size-nbuf * &
          MAXI, MPI_INTEGER, MPI_SUM, MPI_COMM_POOL, info)
     call errore ('reduce', 'error in allreduce 2', info)
     do m = 1, size-nbuf * MAXI
        is (m + nbuf * MAXI) = buff (m)
     enddo

  endif
#endif
  return
end subroutine ireduce

