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
  use parameters, only : DP
#ifdef PARA
  use para
  use allocate
#endif
  implicit none  
  integer :: size  
  real (kind=DP) :: ps (size)  
#ifdef PARA
  integer :: info, i  
  real (kind=DP), pointer :: massimo (:), minimo (:)
  real(kind=DP) ::  chisq  


  include 'mpif.h'  

  call mallocate(massimo, size)  
  call mallocate(minimo , size)  

  call mpi_barrier (MPI_COMM_WORLD, info)  
  call error ('check', 'at the initial barrier', info)  
  call mpi_allreduce (ps, massimo, size, MPI_REAL8, MPI_MAX, &
       MPI_COMM_WORLD, info)
  call error ('check', 'at the first allreduce', info)  
  call mpi_allreduce (ps, minimo, size, MPI_REAL8, MPI_MIN, &
       MPI_COMM_WORLD, info)

  call error ('check', 'at the second allreduce', info)  
  chisq = 0.d0  
  do i = 1, size  
     chisq = chisq + massimo (i) - minimo (i)  
  enddo
  if (chisq.ne.0.d0) then  
     !        call error('check','WARNING, using first proc. data',-1)

     write (6, * ) '*** WARNING, using first proc. data ***'  
     write (6, '(5x,"chisq = ",1pe9.2)') chisq  
     call mpi_bcast (ps, size, MPI_REAL8, 0, MPI_COMM_WORLD, info)  
     call error ('check', 'at the first broadcast', info)  

  endif
  call mfree (minimo)  
  call mfree (massimo)  
#endif
  return  

end subroutine check

