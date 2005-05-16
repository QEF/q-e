!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

subroutine errore(a,b,n)
!-----------------------------------------------------------------------
      use io_global, only: stdout
#ifdef __PARA
      use parallel_include
#endif
      character(len=*) a,b
      integer n, ierr
!
      WRITE( stdout,1) a,b,n
    1 format(//' program ',a,':',a,'.',8x,i8,8x,'stop')
#ifdef __MPI
      call mpi_abort( MPI_COMM_WORLD, ierr, ierr)
      call mpi_finalize(ierr)
#endif
!
      stop
end subroutine
