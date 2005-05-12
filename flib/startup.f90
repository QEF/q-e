!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


       SUBROUTINE PARALLEL_STARTUP(NPROC,MPIME) 
       IMPLICIT NONE

#if defined __MPI
       include 'mpif.h'
#endif

       INTEGER NPROC
       INTEGER MPIME, I, ERR

       MPIME = 0
       NPROC = 1

#if defined __MPI
!      ---------------------------
!      INITIALIZE MPI ENVIRONEMENT 
!      ---------------------------
       CALL MPI_INIT(ERR)  
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ERR)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPIME,ERR)
#endif

       RETURN
       END SUBROUTINE parallel_startup
