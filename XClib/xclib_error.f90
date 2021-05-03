!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE xclib_error( calling_routine, message, ierr )
  !----------------------------------------------------------------------------
  !! This is a simple routine which writes an error message to output (copied
  !! from laxlib).
  !
  USE xclib_utils_and_para
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: calling_routine
  !! the name of the calling calling_routine
  CHARACTER(LEN=*), INTENT(IN) :: message
  !! the output message
  INTEGER, INTENT(IN) :: ierr
  !! error code
  !
  CHARACTER(LEN=6) :: cerr
  INTEGER :: info
  !
  IF( ierr <= 0 ) THEN
     RETURN
  END IF
  !
  ! ... the error message is written un the "*" unit
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
  !
  CALL mpi_abort(MPI_COMM_WORLD,ierr,info)
  !
#endif
  !
  STOP 1
  !
  RETURN
  !
END SUBROUTINE xclib_error
!
!----------------------------------------------------------------------------
SUBROUTINE xclib_infomsg( calling_routine, message )
  !----------------------------------------------------------------------------
  !! This is a simple routine which writes an info/warning message to output.
  !
  USE xclib_utils_and_para,  ONLY: stdout
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*), INTENT(IN) :: calling_routine
  !! the name of the calling routine
  CHARACTER (LEN=*), INTENT(IN) :: message
  !! the output message
  !
  WRITE( UNIT=stdout ,FMT = '(5X,"Message from routine ",A,":")' ) calling_routine
  WRITE( UNIT=stdout ,FMT = '(5X,A)' ) message
  !
  RETURN
  !
END SUBROUTINE xclib_infomsg
