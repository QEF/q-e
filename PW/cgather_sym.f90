!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine cgather_sym (f_in, f_out)
  !-----------------------------------------------------------------------
  ! gather complex data for symmetrization (in phonon code)
  ! COMPLEX*16 f_in = distributed variable (nrxx)
  ! COMPLEX*16 f_out= gathered variable (nrx1*nrx2*nrx3)
  !
#ifdef PARA
#include "machine.h"
  use para
  implicit none

  complex(kind=8) :: f_in (nxx), f_out ( * )
  include 'mpif.h'


  integer :: root, proc, info, displs (nprocp), recvcount (nprocp)

  call start_clock ('cgather')
  root = 0
  do proc = 1, nprocp
     recvcount (proc) = 2 * ncplane * npp (proc)
     if (proc.eq.1) then
        displs (proc) = 0
     else
        displs (proc) = displs (proc - 1) + recvcount (proc - 1)
     endif

  enddo
  call mpi_barrier (MPI_COMM_POOL, info)
  call mpi_allgatherv (f_in, recvcount (me), MPI_REAL8, f_out, &
       recvcount, displs, MPI_REAL8, MPI_COMM_POOL, info)

  call error ('cgather', 'info<>0', info)

  call mpi_barrier (MPI_COMM_WORLD, info)
  call stop_clock ('cgather')
#endif
  return
end subroutine cgather_sym

