!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine init_pool  
  !-----------------------------------------------------------------------
  !
  !    This routine initialize the pool
  !
#ifdef PARA
  !
  use para
  use mp_global , only: mp_global_group_start

  !
  implicit none  
  !
  ! MPI pool division in pools
  !
  include 'mpif.h'  
  integer :: ierr, rank  
  !
  ! set "mypool" and reset "me"
  !

  rank = me  
  mypool = (me-1) / nprocp + 1  

  me = me- (mypool - 1) * nprocp  
  call mpi_barrier (MPI_COMM_WORLD, ierr)  
  call error ('init_pool', 'at the input barrier', ierr)  
  call mpi_comm_split (MPI_COMM_WORLD, mypool, rank, MPI_COMM_POOL, ierr)

  call error ('init_pool', 'MPI_COMM_POOL is wrong', ierr)  

  if (npool.le.1) return  
  call mpi_barrier (MPI_COMM_WORLD, ierr)  
  call error ('init_pool', 'at the second barrier', ierr)  
  call mpi_comm_split (MPI_COMM_WORLD, me, rank, MPI_COMM_ROW, ierr)  

  call error ('init_pool', 'MPI_COMM_ROW is wrong', ierr)  
  !
  ! Initialize globally accessible pool variables
  !
  call mp_global_group_start( (me-1), (mypool-1), MPI_COMM_POOL, &
       MPI_COMM_ROW, nprocp)

#endif
  return  
end subroutine init_pool
