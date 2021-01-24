!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------

! ... Several routines which write an error message to stderr:
! ... if ierr <= 0 it does nothing,
! ... if ierr  > 0 it stops.
!
! ...          **** Important note for parallel execution ***
! ... In QE parallel execution unit 6 is written only by the first node;
! ... all other nodes have unit 6 redirected to nothing (/dev/null).
! ... As a consequence an error not occurring on the first node
! ... will be invisible. In addition, "first node" depends on the parallel
! ... scopes like image, pool, bandgroup. It is better to avoid.
! ...
! ... For this reason, all the fft errors are printed to stderr defined in
! ... module fft_param using iso_fortran_env from F2003.
! ... fftx_error_uniform__ must be called by all the the ranks in an MPI
! ... communicator. The uniform error will be printed only by rank 0 secured
! ... by an MPI barrier before abort. Without the barrier, the error message
! ... may not be seen if a non rank 0 node calls abort first.
! ... In cases that errors don't occur uniformly, fftx_error__ prints the
! ... error message and aborts the code. If it is used for uniform errors,
! ... repeated error messages may be printed.
! ... Use fftx_error_uniform__ wherever errors can be handled cleanly.
! ... fftx_error__ is the last choice you may consider.

SUBROUTINE fftx_error_uniform__( calling_routine, message, ierr, comm )
  !----------------------------------------------------------------------------
  !
  ! ... Writes an uniform error via rank 0 of the given communicator
  !
  USE fft_param
  IMPLICIT NONE
  !
  ! the name of the calling calling_routine
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine
  ! the output message
  CHARACTER(LEN=*), INTENT(IN) :: message
  ! the error number
  INTEGER, INTENT(IN) :: ierr
  ! error scope MPI communicator
  INTEGER, INTENT(IN) :: comm
  !
  CHARACTER(LEN=6) :: cerr
  INTEGER          :: info
  INTEGER          :: my_rank
  !
  IF( ierr <= 0 ) THEN
     RETURN
  END IF
  !
  my_rank = 0
#if defined(__MPI)
  CALL mpi_comm_rank(comm, my_rank, info)
#endif
  if (my_rank == 0) then
    ! ... the error message is written on the "stderr" unit
    !
    WRITE( cerr, FMT = '(I6)' ) ierr
    WRITE( stderr, FMT = '(/,1X,78("%"))' )
    WRITE( stderr, FMT = '(5X,"Error in routine ",A," (",A,"):")' ) &
          TRIM(calling_routine), TRIM(ADJUSTL(cerr))
    WRITE( stderr, FMT = '(1X,A)' ) TRIM(message)
    WRITE( stderr, FMT = '(1X,78("%"),/)' )
    !
    WRITE( stderr, '("     stopping ...")' )
  endif
  !
#if defined(__MPI)
  !
  CALL mpi_barrier(MPI_COMM_WORLD, info)
  CALL mpi_abort(MPI_COMM_WORLD, ierr, info)
  !
#endif
  !
  STOP 1
  !
  RETURN
  !
END SUBROUTINE fftx_error_uniform__

SUBROUTINE fftx_error__( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... This is a simple routine which writes an error message to output:
  !
  USE fft_param
  IMPLICIT NONE
  !
  ! the name of the calling calling_routine
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine
  ! the output message
  CHARACTER(LEN=*), INTENT(IN) :: message
  ! the error number
  INTEGER, INTENT(IN) :: ierr
  !
  CHARACTER(LEN=6) :: cerr
  INTEGER          :: info
  !
  IF( ierr <= 0 ) THEN
     RETURN
  END IF
  !
  ! ... the error message is written on the "stderr" unit
  !
  WRITE( cerr, FMT = '(I6)' ) ierr
  WRITE( stderr, FMT = '(/,1X,78("%"))' )
  WRITE( stderr, FMT = '(5X,"Error in routine ",A," (",A,"):")' ) &
        TRIM(calling_routine), TRIM(ADJUSTL(cerr))
  WRITE( stderr, FMT = '(1X,A)' ) TRIM(message)
  WRITE( stderr, FMT = '(1X,78("%"),/)' )
  !
  WRITE( stderr, '("     stopping ...")' )
  !
#if defined(__MPI)
  !
  CALL mpi_abort(MPI_COMM_WORLD, ierr, info)
  !
#endif
  !
  STOP 1
  !
  RETURN
  !
END SUBROUTINE fftx_error__
