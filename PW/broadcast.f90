!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine broadcast (ndata, data)
  !-----------------------------------------------------------------------
  !
  ! this routine broadcasts the real*8 vector 'data' from the first proces
  ! of each pool to all the other processors of the pool
  !
#ifdef __PARA
  use para
#endif
  implicit none
  ! on INPUT
  integer :: ndata
  real (8) :: data (ndata)
#ifdef __PARA
  include 'mpif.h'

  integer :: root, ierr
  if (nprocp.le.1) return
  !
  ! root is the "broadcasting" processors.
  ! N.B. in MPI the processors are labeled from 0 to nproc-1
  !
  root = 0
  call mpi_barrier (MPI_COMM_POOL, ierr)
  call mpi_bcast (data, ndata, MPI_REAL8, root, MPI_COMM_POOL, ierr)
  call errore ('broadcast', 'info<>0', ierr)
#endif
  return
end subroutine broadcast

