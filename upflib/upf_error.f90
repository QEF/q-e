!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE upf_error( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !
  ! ... Writes an error message to output (unit "*") if ierr != 0 and stops
  !
  USE upf_parallel_include
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine, message
    ! the name of the calling calling_routine
    ! the output message
  INTEGER,          INTENT(IN) :: ierr
  !
  CHARACTER(LEN=6) :: cerr
  INTEGER          :: info
  !
  IF( ierr /= 0 ) THEN
     !
     WRITE( cerr, FMT = '(I6)' ) ierr
     WRITE( UNIT = *, FMT = '(/,1X,78("%"))' )
     WRITE( UNIT = *, FMT = '(5X,"Error in routine ",A," (",A,"):")' ) &
          TRIM(calling_routine), TRIM(ADJUSTL(cerr))
     WRITE( UNIT = *, FMT = '(5X,A)' ) TRIM(message)
     WRITE( UNIT = *, FMT = '(1X,78("%"),/)' )
     !
     WRITE( *, '("     stopping ...")' )
     !
#if defined(__MPI)
     CALL mpi_abort(MPI_COMM_WORLD,ierr,info)
#endif
     STOP 1
  END IF
  !
END SUBROUTINE upf_error
