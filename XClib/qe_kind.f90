!
!
MODULE kind_l
!------------------------------------------------------------------------------!
! #if defined(__MPI)
! #if defined(__MPI_MODULE)
!   USE mpi
! #else
!   INCLUDE 'mpif.h'
  ! #endif
  ! #endif
  !
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  PUBLIC :: DP
  !
END MODULE      
