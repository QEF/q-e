!
! Copyright (C) 2011-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE open_close_input_file
  !
  USE io_global,     ONLY : stdout, xmlinputunit
  ! ... BEWARE: xmlinputunit modified by open_input_file
  USE iotk_module,   ONLY : iotk_open_read, iotk_close_read,iotk_attlenx
  !
  INTEGER, SAVE :: unit_loc = 5
  LOGICAL, SAVE :: lxmlinput_loc = .false.
  CHARACTER(LEN=256), SAVE :: input_file = ' '
  PRIVATE
  PUBLIC :: open_input_file, close_input_file
  !
CONTAINS
  !----------------------------------------------------------------------------
  INTEGER FUNCTION open_input_file ( lxmlinput, attr, unit )
  !-----------------------------------------------------------------------------
  !
  ! ...  opens "unit" (optional, unit_loc=5 if unspecified) for input read
  ! ...  if lxmlinput=.true. and "attr" are present, tests for a xml file
  ! ...  to be opened with attribute "attr"
  ! ...  The code detects an input filename via options:
  ! ...     "-in failename", "-inp filename", "-input filename"
  ! ...  If no filename is not specified, the standard input is dumped to
  ! ...  temporary file "input_tmp.in" and this is opened for read
  ! ...  Returns -1 if standard input is dumped to file
  ! ...  Returns  0 if input file is successfully opened
  ! ...  Returns  1 if called with wrong arguments
  ! ...  Returns  2 if there was an error opening file
  ! ...  ---------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  LOGICAL, intent(inout), optional :: lxmlinput
  CHARACTER (len=*), intent(inout), optional :: attr
  INTEGER, intent(in), optional :: unit
  !
  LOGICAL :: lfound, lxmlinput_loc,lcheckxml
  INTEGER  :: iiarg, nargs, iargc, ierr
  CHARACTER (len=50) :: arg
  !
  INTEGER :: stdin=5, stdtmp
  CHARACTER(LEN=512) :: dummy
  !
#if defined(__ABSOFT)
#   define getarg getarg_
#   define iargc  iargc_
#endif
  !
  IF (     PRESENT(attr) .AND. .NOT.PRESENT(lxmlinput) .OR. &
      .NOT.PRESENT(attr) .AND.      PRESENT(lxmlinput) ) THEN
     !
     open_input_file = 1
     RETURN
     !
  ENDIF
  !
  lcheckxml = .false.
  IF ( PRESENT(attr) .AND. (PRESENT(lxmlinput)) ) lcheckxml = .true.
  !
  unit_loc = 5
  IF(PRESENT(unit)) unit_loc = unit
  !
  xmlinputunit = unit_loc
  lxmlinput_loc = .false.
  !
  ! ... Input from file ?
  !
  input_file=' '
  lfound=.false.
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
        lfound=.true.
        GO TO 10
     END IF
     !
  END DO
  !
10 stdtmp = find_free_unit()
  !
  IF ( .NOT. lfound ) THEN
     !
     ! if no file specified then copy from standard input
     !
     input_file="input_tmp.in"
     OPEN(UNIT = stdtmp, FILE=trim(input_file), FORM='formatted', &
          STATUS='unknown', IOSTAT = ierr )
     IF ( ierr > 0 ) GO TO 30
     !
     dummy=' '
     WRITE(stdout, '(5x,a)') "Waiting for input..."
     DO WHILE ( TRIM(dummy) .NE. "MAGICALME" )
        READ (stdin,fmt='(A512)',END=20) dummy
        WRITE (stdtmp,'(A)') trim(dummy)
     END DO
     !
20   CLOSE ( UNIT=stdtmp, STATUS='keep' )
  ENDIF
  !
  IF (lcheckxml) THEN
    !
    OPEN ( UNIT = stdtmp, FILE = TRIM(input_file) , FORM = 'FORMATTED', &
          STATUS = 'OLD', IOSTAT = ierr )
    IF ( ierr > 0 ) GO TO 30
    CALL test_input_xml (stdtmp, lxmlinput_loc )
    CLOSE ( UNIT=stdtmp, status='keep')
    !
    lxmlinput = lxmlinput_loc
    !
  ENDIF
  !
  IF (lxmlinput_loc) then
     IF ( input_file .NE. "input_tmp.in") THEN
        WRITE(stdout, '(5x,a)') "Reading xml input from "//TRIM(input_file)
     ELSE
        WRITE(stdout, '(5x,a)') "Reading xml input from standard input"
     END IF
     CALL iotk_open_read( unit_loc, TRIM(input_file), attr = attr, &
                          qe_syntax = .true., ierr = ierr)
  ELSE 
     IF ( input_file .NE. "input_tmp.in") THEN
         WRITE(stdout, '(5x,a)') "Reading input from "//TRIM(input_file)
     ELSE
         WRITE(stdout, '(5x,a)') "Reading input from standard input"
     END IF
     OPEN ( UNIT = unit_loc, FILE = TRIM(input_file), FORM = 'FORMATTED', &
            STATUS = 'OLD', IOSTAT = ierr )
  ENDIF
  IF ( ierr > 0 ) GO TO 30
  !
  open_input_file = 0
  RETURN
30 open_input_file = 2
  RETURN
  !
END FUNCTION open_input_file

INTEGER FUNCTION close_input_file ( )
  !
  ! ...  this subroutine closes the input file opened by open_input_file
  ! ...  removes temporary file if data was read from stdin
  ! ...  (not in the xml case, though, because it is not clear how to do it)
  ! ...  returns -1 if unit is not opened, 0 if no problem, > 0 if problems
  ! ...  ---------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  LOGICAL :: opened
  INTEGER :: ierr
  !
  INQUIRE ( unit_loc, opened = opened )
  IF (opened) THEN
     !
     IF (lxmlinput_loc) THEN
        CALL iotk_close_read(unit=unit_loc, ierr = ierr)
     ELSE
        IF ( TRIM(input_file) == "input_tmp.in") THEN
           CLOSE (UNIT=unit_loc, STATUS='delete', IOSTAT=ierr )
        ELSE
           CLOSE (UNIT=unit_loc, STATUS='keep', IOSTAT=ierr )
        ENDIF
     ENDIF
     !
  ELSE
     ierr = -1
  ENDIF 
  !
  close_input_file = ierr
  !
  END FUNCTION close_input_file
!
ENDMODULE open_close_input_file
