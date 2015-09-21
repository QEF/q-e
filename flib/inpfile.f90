!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Wrappers for intrinsic iargc, getarg, getenv - machine-dependent stuff here
!
INTEGER FUNCTION i_argc ( )
#if defined(__NAG)
  USE F90_UNIX_ENV, ONLY : iargc
#else
  ! do not declare it external: gfortran doesn't like it
  INTEGER :: iargc
#endif
  i_argc = iargc ( )
END FUNCTION i_argc
  !
SUBROUTINE get_env ( variable_name, variable_value )
#if defined(__NAG)
  USE F90_UNIX_ENV, ONLY : getenv
#endif
  CHARACTER (LEN=*)  :: variable_name, variable_value
  CALL getenv ( variable_name, variable_value)
END SUBROUTINE get_env
!
SUBROUTINE get_arg ( iarg, arg )
#if defined(__NAG)
  USE F90_UNIX_ENV, ONLY : getarg
#endif
  INTEGER, INTENT(IN) :: iarg
  CHARACTER (LEN=*), INTENT(OUT) :: arg
  CALL getarg ( iarg, arg )
END SUBROUTINE get_arg
!
!----------------------------------------------------------------------------
SUBROUTINE input_from_file( )
  !
  ! This subroutine checks command-line arguments for -i[nput] "file name"
  ! if "file name" is present, attach input unit 5 to the specified file
  !
  IMPLICIT NONE
  !
  INTEGER             :: stdin = 5, stderr = 6, ierr = 0
  CHARACTER (LEN=256) :: input_file
  LOGICAL             :: found
  !
  INTEGER, EXTERNAL   :: i_argc
  INTEGER :: iiarg, nargs
  !
  nargs = i_argc()
  found = .FALSE.
  input_file = ' '
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL get_arg( iiarg, input_file )
     !
     IF ( TRIM( input_file ) == '-i'     .OR. &
          TRIM( input_file ) == '-in'    .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-input' ) THEN
        !
        CALL get_arg( ( iiarg + 1 ) , input_file )
        found =.TRUE.
        EXIT
        !
     END IF
     !
  END DO
  !
  IF ( found ) THEN
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
  RETURN 
  !
END SUBROUTINE input_from_file

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
  LOGICAL             :: exst
  INTEGER, EXTERNAL   :: i_argc
  !
  nargs = i_argc()
  CALL get_arg (0,prgname)
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
     CALL get_arg (1,input_file)
  ELSE
     PRINT  '(A,": too many arguments ",i4)', TRIM(prgname), nargs
  END IF
  RETURN
20 PRINT  '(A,": reading file name ",A)', TRIM(prgname), TRIM(input_file)
  !
END SUBROUTINE get_file

