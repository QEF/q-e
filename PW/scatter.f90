!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine scatter (f_in, f_out)
  !-----------------------------------------------------------------------
  ! scatters data from the first processor of every pool
  !
  ! REAL*8 f_in = gathered variable (nrx1*nrx2*nrx3)
  ! REAL*8 f_out= distributed variable (nxx)
  !
#ifdef __PARA
#include "machine.h"
  use para
  implicit none

  real (8) :: f_in ( * ), f_out (nxx)
  include 'mpif.h'


  integer :: root, proc, info, displs (nprocp), sendcount (nprocp)
  root = 0
  call start_clock ('scatter')
  do proc = 1, nprocp
     sendcount (proc) = ncplane * npp (proc)
     if (proc.eq.1) then
        displs (proc) = 0
     else
        displs (proc) = displs (proc - 1) + sendcount (proc - 1)
     endif

  enddo
  call mpi_barrier (MPI_COMM_POOL, info)
  call mpi_scatterv (f_in, sendcount, displs, MPI_REAL8, f_out, &
       sendcount (me), MPI_REAL8, root, MPI_COMM_POOL, info)
  call errore ('gather_pot', 'info<>0', info)
  if (sendcount(me).ne.nxx) f_out(sencount(me)+1:nxx) = 0.d0

  call stop_clock ('scatter')
#endif
  return
end subroutine scatter

