!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine fft_scatter1 (f_in, nrx3, nxx_, f_aux, ncp_, npp_, sign)
  !-----------------------------------------------------------------------
  !
  ! transpose the fft grid across nodes
  ! a) From columns to planes (sign > 0)
  !
  !    "columns" (or "pencil") representation:
  !    processor "me" has ncp_(me) contiguous columns along z
  !    Each column has nrx3 elements for a fft of order nr3
  !    nrx3 can be =nr3+1 in order to reduce memory conflicts.
  !
  !    The transpose take places in two steps:
  !    1) on each processor the columns are divided into slices along z
  !       that are stored contiguously. On processor "me", slices for
  !       processor "proc" are npp_(proc)*ncp_(me) big
  !    2) all processors communicate to exchange slices
  !       (all columns with z in the slice belonging to "me"
  !        must be received, all the others must be sent to "proc")
  !    Finally one gets the "planes" representation:
  !    processor "me" has npp_(me) complete xy planes
  !
  !  b) From planes to columns (sign < 0)
  !
  !  Quite the same in the opposite direction
  !
  !  The output is overwritten on f_in ; f_aux is used as work space
  !
#include "machine.h"
#ifdef __PARA
  use para
  use parameters, only : DP
  implicit none

  include 'mpif.h'
  integer :: nrx3, nxx_, sign, ncp_ (maxproc), npp_ (maxproc)

  real (kind=DP) :: f_in (2*nxx_), f_aux (2*nxx_)
  integer :: dest, from, k, offset1 (maxproc), sendcount (maxproc), &
       sdispls (maxproc), recvcount (maxproc), rdispls (maxproc), &
       proc, ierr
  !
  !
  if (nprocp.eq.1) return
  !
  call start_clock ('fft_scatter')
  !
  ! sendcount(proc): amount of data processor "me" must send to processor
  ! recvcount(proc): amount of data processor "me" must receive from
  !
  do proc = 1, nprocp
     sendcount (proc) = 2 * npp_ (proc) * ncp_ (me)
     recvcount (proc) = 2 * npp_ (me) * ncp_ (proc)
  enddo
  !
  ! offset1(proc) is used to locate the slices to be sent to proc
  ! sdispls(proc)+1 is the beginning of data that must be sent to proc
  ! rdispls(proc)+1 is the beginning of data that must be received from pr
  !
  offset1 (1) = 1
  sdispls (1) = 0
  rdispls (1) = 0
  do proc = 2, nprocp
     offset1 (proc) = offset1 (proc - 1) + 2 * npp_ (proc - 1)
     sdispls (proc) = sdispls (proc - 1) + sendcount (proc - 1)
     rdispls (proc) = rdispls (proc - 1) + recvcount (proc - 1)
  enddo
  !

  ierr = 0
  if (sign.gt.0) then
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     do proc = 1, nprocp
        from = offset1 (proc)
        dest = 1 + sdispls (proc)
        do k = 1, ncp_ (me)
           call DCOPY (2*npp_(proc), f_in (from + 2*(k-1)*nrx3),       1, &
                                     f_aux(dest + 2*(k-1)*npp_(proc)), 1)
        enddo
     enddo
     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     f_in(:) = (0.d0,0.d0)
     !
     ! step two: communication
     !
     call mpi_barrier (MPI_COMM_POOL, ierr)
     call mpi_alltoallv (f_aux, sendcount, sdispls, MPI_REAL8, f_in, &
          recvcount, rdispls, MPI_REAL8, MPI_COMM_POOL, ierr)
     call errore ('fft_scatter', 'info<>0', ierr)
     !
  else
     !
     !  "backward" scatter from planes to columns
     !
     !  step two: communication
     !
     call mpi_barrier (MPI_COMM_POOL, ierr)
     call mpi_alltoallv (f_in, recvcount, rdispls, MPI_REAL8, f_aux, &
          sendcount, sdispls, MPI_REAL8, MPI_COMM_POOL, ierr)
     call errore ('fft_scatter', 'info<>0', ierr)
     !
     !  step one: store contiguously the columns
     !
     f_in(:) = (0.d0,0.d0)
     !
     do proc = 1, nprocp
        from = 1 + sdispls (proc)
        dest = offset1 (proc)
        do k = 1, ncp_ (me)
           call DCOPY(2*npp_(proc), f_aux(from + 2*(k-1)*npp_(proc)), 1, &
                                    f_in (dest + 2*(k-1)*nrx3),       1)
        enddo

     enddo

  endif
  call stop_clock ('fft_scatter')
#endif
  return

end subroutine fft_scatter1
