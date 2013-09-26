!
! Copyright (C) 2003-2004 Carlo Cavazzoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
!   SISSA Code Interface -- Carlo Cavazzoni
!------------------------------------------------------------------------------C
      MODULE parallel_include

#if defined __MPI 
!
!     Include file for MPI
!
         INCLUDE 'mpif.h'
#else
         ! dummy world and null communicator
         INTEGER, PARAMETER :: MPI_COMM_WORLD =  0
         INTEGER, PARAMETER :: MPI_COMM_NULL  = -1
         INTEGER, PARAMETER :: MPI_COMM_SELF  = -2
#endif

      END MODULE parallel_include
