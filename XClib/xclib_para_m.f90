!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------
MODULE xclib_parallel_include
!-----------------------------------------------------
!! MPI stuff
!
#if defined (__MPI)
        !
        !     Include file for MPI
        !
#if defined (__MPI_MODULE)
        USE mpi
#else
        INCLUDE 'mpif.h'
#endif
#else
        ! dummy world and null communicator
        INTEGER, PARAMETER :: MPI_COMM_WORLD =  0
        INTEGER, PARAMETER :: MPI_COMM_NULL  = -1
        INTEGER, PARAMETER :: MPI_COMM_SELF  = -2
#endif
END MODULE xclib_parallel_include
