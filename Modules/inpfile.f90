!
! Copyright (C) 2002-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE input_from_file( )
  !
  !! Check command-line arguments for -i[nput] "file name"
  !! if "file name" is present, attach input unit 5 to the specified file
  !! In parallel execution, must be called by a single processor 
  !
  USE open_close_input_file, ONLY : get_file_name
  !
  IMPLICIT NONE
  !
  INTEGER              :: stdin = 5, stderr = 6, ierr = 0
  CHARACTER(LEN = 256) :: input_file
  !
  input_file = get_file_name ( )
  !
  IF ( TRIM ( input_file ) /= ' ' ) THEN
    !
    OPEN ( UNIT = stdin, FILE = input_file, FORM = 'FORMATTED', &
           STATUS = 'OLD', IOSTAT = ierr )
    !
    ! do not call "errore" here: it may hang in parallel execution
    ! if this routine is called by a single processor
    !
    IF ( ierr > 0 ) WRITE (stderr, &
    '(" *** Fatal error: input file ",A," not found ***")' ) TRIM( input_file )
    !
  ELSE
    ierr = -1
  ENDIF
  !
END SUBROUTINE input_from_file

!----------------------------------------------------------------------------
!
SUBROUTINE get_file( input_file )
  !
  !! This subroutine reads, either from command line or from terminal,
  !! the name of a file to be opened. To be used for serial codes only.
  !! Expected syntax: "code [filename]"  (one command-line option, or none)
  !
  USE open_close_input_file, ONLY : get_file_name
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*),INTENT(OUT)  :: input_file
  !! On output contains the path to the input file
  INTEGER :: stdin = 5, stdout = 6, stderr = 6
  LOGICAL :: exst
  !
  input_file = get_file_name ( )
  !
  IF ( TRIM ( input_file ) == ' ' ) THEN
10   WRITE(stdout,'(5x,"Input file > ")', advance="NO")
     READ (stdin,'(a)', end = 20, err=20) input_file
     IF ( TRIM(input_file) == ' ') GO TO 10
     INQUIRE ( FILE = input_file, EXIST = exst )
     IF ( .NOT. exst ) THEN
        WRITE(stderr,'(A,": file not found")') TRIM(input_file)
        GO TO 10
     END IF
  END IF
  RETURN
20 WRITE(stdout,'("Fatal error reading file name ",A)') TRIM(input_file)
  !
END SUBROUTINE get_file

