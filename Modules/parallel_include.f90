!------------------------------------------------------------------------------!
!   SISSA Code Interface -- Carlo Cavazzoni
!------------------------------------------------------------------------------C
      MODULE parallel_include

         USE kinds
         LOGICAL tparallel

#if defined __MPI || defined SHMEM
!
!     Include file for MPI
!
         INCLUDE 'mpif.h'
         DATA tparallel /.true./
#else
         DATA tparallel /.false./
#endif

      END MODULE parallel_include
