!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


       SUBROUTINE SYNCRONIZE 
#if defined __MPI
       IMPLICIT NONE
       include 'mpif.h'
       INTEGER IERR
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif
       RETURN 
       END 
