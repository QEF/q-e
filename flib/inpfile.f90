!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__ABSOFT)
#  define getenv getenv_
#  define getarg getarg_
#  define iargc  iargc_
#endif
!
SUBROUTINE get_env ( variable_name, variable_value )
  !
  ! Wrapper for intrinsic getenv - all machine-dependent stuff here
  !
  CHARACTER (LEN=*)  :: variable_name, variable_value
  !
  CALL getenv ( variable_name, variable_value)
  !
END SUBROUTINE get_env
!----------------------------------------------------------------------------
SUBROUTINE input_from_file( )
  !
  ! This subroutine checks program arguments and, if input file is present,
  ! attach input unit ( 5 ) to the specified file
  !
  IMPLICIT NONE
  !
  INTEGER             :: stdin = 5, stderr = 6, ierr = 0
  CHARACTER (LEN=256)  :: input_file
  LOGICAL, EXTERNAL :: input_file_name_getarg
  !
  IF ( input_file_name_getarg ( input_file ) ) THEN 
     !
     OPEN ( UNIT = stdin, FILE = input_file, FORM = 'FORMATTED', &
            STATUS = 'OLD', IOSTAT = ierr )
     !
     ! TODO: return error code ierr (-1 no file, 0 file opened, > 1 error)
     ! do not call "errore" here: it may hang in parallel execution
     ! if this routine is called by a single processor
     !
     IF ( ierr > 0 ) WRITE (stderr, &
            '(" *** input file ",A," not found ***")' ) TRIM( input_file )
     !
  ELSE
     ierr = -1
  END IF
  !
END SUBROUTINE input_from_file

FUNCTION input_file_name_getarg ( input_file ) RESULT ( found )
  !
  ! checks for presence of command-line option "-i" (or "-in", "-inp","-input")
  ! Returns true if found and the following name file in variable "input_ file"
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=256), INTENT (OUT) :: input_file
  LOGICAL             :: found
  !
  INTEGER :: iiarg, nargs
  ! Do not define iargc as external: gfortran doesn't like it
  INTEGER :: iargc
  !
  ! ... Input from file ?
  !
  nargs = iargc()
  found = .FALSE.
  input_file = ' '
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )
     !
     IF ( TRIM( input_file ) == '-i'     .OR. &
          TRIM( input_file ) == '-in'    .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-input' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )
        found =.TRUE.
        RETURN
        !
     END IF
     !
  END DO
  !
  RETURN 
  !
END FUNCTION input_file_name_getarg
!
!----------------------------------------------------------------------------
!
SUBROUTINE get_file( input_file )
  !
  ! This subroutine reads, either from command line or from terminal,
  ! the name of a file to be opened. To be used for serial codes only.
  ! Expected syntax: "code [filename]"  (one command-line option, or none)
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*)  :: input_file
  !
  CHARACTER (LEN=256) :: prgname
  INTEGER             :: nargs
  INTEGER             :: iargc
  LOGICAL             :: exst
  !
  nargs = iargc()
  CALL getarg (0,prgname)
  !
  IF ( nargs == 0 ) THEN
10   PRINT  '("Input file > ",$)'
     READ (5,'(a)', end = 20, err=20) input_file
     IF ( input_file == ' ') GO TO 10
     INQUIRE ( FILE = input_file, EXIST = exst )
     IF ( .NOT. exst) THEN
        PRINT  '(A,": file not found")', TRIM(input_file)
        GO TO 10
     END IF
  ELSE IF ( nargs == 1 ) then
     CALL getarg (1,input_file)
  ELSE
     PRINT  '(A,": too many arguments ",i4)', TRIM(prgname), nargs
  END IF
  RETURN
20 PRINT  '(A,": reading file name ",A)', TRIM(prgname), TRIM(input_file)
  !
END SUBROUTINE get_file

