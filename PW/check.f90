!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------


subroutine check (size, ps)
  !-----------------------------------------------------------------------
#include "machine.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: size
  real (kind=DP) :: ps (size)
#ifdef __PARA
  integer :: info, i
  real (kind=DP), allocatable :: massimo (:), minimo (:)
  real(kind=DP) ::  chisq


  include 'mpif.h'

  allocate (massimo( size))    
  allocate (minimo ( size))    

  call mpi_barrier (MPI_COMM_WORLD, info)
  call errore ('check', 'at the initial barrier', info)
  call mpi_allreduce (ps, massimo, size, MPI_REAL8, MPI_MAX, &
       MPI_COMM_WORLD, info)
  call errore ('check', 'at the first allreduce', info)
  call mpi_allreduce (ps, minimo, size, MPI_REAL8, MPI_MIN, &
       MPI_COMM_WORLD, info)

  call errore ('check', 'at the second allreduce', info)
  chisq = 0.d0
  do i = 1, size
     chisq = chisq + massimo (i) - minimo (i)
  enddo
  if (chisq.ne.0.d0) then
     !        call errore('check','WARNING, using first proc. data',-1)

     WRITE( stdout, * ) '*** WARNING, using first proc. data ***'
     WRITE( stdout, '(5x,"chisq = ",1pe9.2)') chisq
     call mpi_bcast (ps, size, MPI_REAL8, 0, MPI_COMM_WORLD, info)
     call errore ('check', 'at the first broadcast', info)

  endif
  deallocate (minimo)
  deallocate (massimo)
#endif
  return

end subroutine check

