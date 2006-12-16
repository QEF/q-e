!
! Copyright (C) 2002-2005 PWSCF group
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
!----------------------------------------------------------------------------
SUBROUTINE input_from_file( )
  !
  ! This subroutine checks program arguments and, if input file is present,
  ! attach input unit ( 5 ) to the specified file
  !
  !
  IMPLICIT NONE
  !
  INTEGER             :: unit = 5, &
                         ilen, iiarg, nargs, ierr
  ! do not define iargc as external: g95 does not like it
  INTEGER             :: iargc
  CHARACTER (LEN=80)  :: input_file
  !
  ! ... Input from file ?
  !
  nargs = iargc()
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )
     !
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )
        !
        OPEN ( UNIT = unit, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        !
        ! TODO: return error code instead
        !CALL errore( 'input_from_file', 'input file ' // TRIM( input_file ) &
        !     & // ' not found' , ierr )
        !
     END IF
     !
  END DO

END SUBROUTINE input_from_file
!
!----------------------------------------------------------------------------
!
SUBROUTINE get_file( input_file )
  !
  ! This subroutine reads, either from command line or from terminal,
  ! the name of a file to be opened
  ! TODO: return error code if an error occurs
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=*)  :: input_file
  !
  CHARACTER (LEN=256)  :: prgname
  ! do not define iargc as external: g95 does not like it
  INTEGER             :: nargs, iargc
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
!
!----------------------------------------------------------------------------
!
SUBROUTINE get_arg_nimage( nimage )
   !
   IMPLICIT NONE
   !
   INTEGER :: nimage
   !
   INTEGER :: nargs, iiarg
   CHARACTER(LEN=10) :: np
   INTEGER :: iargc
   !
   nimage = 1
   nargs = iargc()
   !
   DO iiarg = 1, ( nargs - 1 )
      !
      CALL getarg( iiarg, np )
      !
      IF ( TRIM( np ) == '-nimage' .OR. TRIM( np ) == '-nimages' ) THEN
         !
         CALL getarg( ( iiarg + 1 ), np )
         READ( np, * ) nimage
         !
      END IF
      !
   END DO
   !
   RETURN
END SUBROUTINE get_arg_nimage
!
!----------------------------------------------------------------------------
!
SUBROUTINE get_arg_ntg( ntask_groups )
   !
   IMPLICIT NONE
   !
   INTEGER :: ntask_groups
   !
   INTEGER :: nargs, iiarg
   CHARACTER(LEN=20) :: np
   INTEGER :: iargc
   !
   ntask_groups = 0
   nargs = iargc()
   !
   DO iiarg = 1, ( nargs - 1 )
      !
      CALL getarg( iiarg, np )
      !
      IF ( TRIM( np ) == '-ntg' .OR. TRIM( np ) == '-ntask_groups' ) THEN
         !
         CALL getarg( ( iiarg + 1 ), np )
         READ( np, * ) ntask_groups
         !
      END IF
      !
   END DO
   !
   RETURN
END SUBROUTINE get_arg_ntg

SUBROUTINE get_env ( variable_name, variable_value )
  !
  ! Wrapper for intrinsic getenv - all machine-dependent stuff here
  !
  CHARACTER (LEN=*)  :: variable_name, variable_value
  !
  CALL getenv ( variable_name, variable_value)
  !
END SUBROUTINE get_env
