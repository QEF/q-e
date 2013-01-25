!
! Copyright (C) 2011-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE open_close_input_file
  !
  USE io_global,     ONLY : stdin, stdout, xmlstdin
  USE iotk_module,   ONLY : iotk_open_read, iotk_close_read,iotk_attlenx
  !
  LOGICAL, SAVE :: lxmlinput_loc = .false.
  CHARACTER(LEN=256), SAVE :: input_file = ' '
  PRIVATE
  PUBLIC :: open_input_file, close_input_file
  !
CONTAINS
  !----------------------------------------------------------------------------
  INTEGER FUNCTION open_input_file ( lxmlinput, attr )
  !-----------------------------------------------------------------------------
  !
  ! ...  Opens file for input read, connecting it to unit stdin (text input)
  ! ...  or xmlstdin (xml file).
  ! ...  On entry:
  ! ...    if optional variable lxmlinput is present and set to .true., 
  ! ...    test if the file is a valid xml file. In this case, optional
  ! ...    variable attr must be present and used to open the file.
  ! ...  The code detects an input filename via option "-i filename"
  ! ...  (also "-in", "-inp", "-input" are accepted)
  ! ...  If no filename is specified, the standard input is dumped to
  ! ...  temporary file "input_tmp.in" and this is opened for read
  ! ...  On exit:
  ! ...    Returns -1 if standard input is dumped to file
  ! ...    Returns  0 if input file is successfully opened
  ! ...    Returns  1 if called with wrong arguments
  ! ...    Returns  2 if there was an error opening file
  ! ...    lxmlinput=.true. if the file has extension '.xml' or '.XML'
  ! ...       or if either <xml...> or <?xml...> is found as first token
  ! ...  ---------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  LOGICAL, intent(inout), optional :: lxmlinput
  CHARACTER (len=*), intent(inout), optional :: attr
  !
  LOGICAL :: lfound, lxmlinput_loc,lcheckxml
  INTEGER :: ierr, len
  INTEGER :: stdtmp
  CHARACTER(LEN=512) :: dummy
  LOGICAL, EXTERNAL :: test_input_xml, input_file_name_getarg
  INTEGER, EXTERNAL :: find_free_unit
  !
  !
  lcheckxml = .false.
  IF ( PRESENT(lxmlinput) ) lcheckxml = lxmlinput
  !
  IF ( lcheckxml .AND. .NOT.PRESENT(attr) ) THEN
     open_input_file = 1
     RETURN
  ENDIF
  !
  ! ... Input from file ?
  !
  lfound = input_file_name_getarg ( input_file ) 
  !
  stdtmp = find_free_unit()
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
  lxmlinput_loc = .false.
  IF (lcheckxml) THEN
    !
    len = LEN_TRIM(input_file)
    IF ( len > 4) THEN
       lxmlinput_loc = ( input_file(len-3:len) == '.xml' .OR. &
                         input_file(len-3:len) == '.XML' )
    END IF
    IF ( .NOT. lxmlinput_loc ) THEN
       OPEN ( UNIT = stdtmp, FILE = TRIM(input_file) , FORM = 'FORMATTED', &
              STATUS = 'OLD', IOSTAT = ierr )
       IF ( ierr > 0 ) GO TO 30
       lxmlinput_loc = test_input_xml (stdtmp )
       CLOSE ( UNIT=stdtmp, status='keep')
    END IF
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
     CALL iotk_open_read( xmlstdin, TRIM(input_file), attr = attr, &
                          qe_syntax = .true., ierr = ierr)
  ELSE 
     IF ( input_file .NE. "input_tmp.in") THEN
         WRITE(stdout, '(5x,a)') "Reading input from "//TRIM(input_file)
     ELSE
         WRITE(stdout, '(5x,a)') "Reading input from standard input"
     END IF
     OPEN ( UNIT = stdin, FILE = TRIM(input_file), FORM = 'FORMATTED', &
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
  ! ...  removes temporary file if data was read from stdin (text file)
  ! ...  returns -1 if unit is not opened, 0 if no problem, > 0 if problems
  ! ...  ---------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  LOGICAL :: opened
  INTEGER :: ierr
  !
  INQUIRE ( stdin, opened = opened )
  IF (.NOT.opened .AND. stdin /= xmlstdin) INQUIRE ( xmlstdin, opened = opened )
  IF (opened) THEN
     !
     IF (lxmlinput_loc) THEN
        CALL iotk_close_read(unit=xmlstdin, ierr = ierr)
     ELSE
        IF ( TRIM(input_file) == "input_tmp.in") THEN
           CLOSE (UNIT=stdin, STATUS='delete', IOSTAT=ierr )
        ELSE
           CLOSE (UNIT=stdin, STATUS='keep', IOSTAT=ierr )
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
