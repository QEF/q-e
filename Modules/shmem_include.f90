!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
! Author :  Carlo Cavazzoni (CINECA)
! Update :  October 1999
!------------------------------------------------------------------------------!
!
!     This is a fixed format file
!
!------------------------------------------------------------------------------C
!
!     Holds External information for Message Passing Systems
!
!------------------------------------------------------------------------------C

      MODULE shmem_include

        USE kinds
        IMPLICIT NONE
        SAVE

         LOGICAL TSHMEM

#if defined __SHMEM
!
!     Include file for SHMEM Library
!
         INCLUDE 'mpp/shmem.fh'
         INTEGER, PARAMETER :: mp_shmem_bufsize = &
           MAX(524288,SHMEM_REDUCE_MIN_WRKDATA_SIZE)
         INTEGER PSYNCB(SHMEM_BARRIER_SYNC_SIZE)
         INTEGER PSYNCC(SHMEM_COLLECT_SYNC_SIZE)
         INTEGER PSYNC_STA(SHMEM_REDUCE_SYNC_SIZE)
         REAL(dbl), SAVE :: mp_shmem_buffer(mp_shmem_bufsize)
         REAL(dbl), SAVE :: mp_shmem_work(mp_shmem_bufsize)

         DATA PSYNC_STA /SHMEM_REDUCE_SYNC_SIZE*SHMEM_SYNC_VALUE/
         DATA PSYNCB /SHMEM_BARRIER_SYNC_SIZE*SHMEM_SYNC_VALUE/
         DATA PSYNCC /SHMEM_COLLECT_SYNC_SIZE*SHMEM_SYNC_VALUE/
         DATA TSHMEM /.TRUE./

#else
!
         DATA TSHMEM /.FALSE./
#endif
!
      END MODULE shmem_include
