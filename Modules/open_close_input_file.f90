!
! Copyright (C) 2011-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Aug 2020 (PG): streamlined and simplified
!                open_input_file() may read the file name either from
!                passed argument or directly from the command line
! Aug 2018 (PG): reading of old xml input file using iotk deleted

MODULE open_close_input_file
  !
  !! Input file management.
  !
  USE io_global,     ONLY : stdin, stdout, qestdin
  !
  CHARACTER(LEN=256), SAVE :: input_file = ' '
  PRIVATE
  PUBLIC :: get_file_name, open_input_file, close_input_file
  !
CONTAINS
  !----------------------------------------------------------------------------
  FUNCTION get_file_name ( )
    !
    !! Get the file name provided in command line
    !
    IMPLICIT NONE
    !
    CHARACTER (len=256) :: get_file_name
    !
    LOGICAL :: found
    INTEGER :: iiarg, nargs
    ! 
    nargs = command_argument_count()
    get_file_name = ' '
    found = .false.
    !
    DO iiarg = 1, ( nargs - 1 )
       !
       CALL get_command_argument( iiarg, get_file_name )
       !
       IF ( TRIM( get_file_name ) == '-i'     .OR. &
            TRIM( get_file_name ) == '-in'    .OR. &
            TRIM( get_file_name ) == '-inp'   .OR. &
            TRIM( get_file_name ) == '-input' ) THEN
          !
          CALL get_command_argument( ( iiarg + 1 ) , get_file_name )
          found = .true.
          EXIT
          !
       END IF
       !
    END DO
    !
    IF ( .NOT. found ) get_file_name = ' '
    !
  END FUNCTION get_file_name
  !
  !----------------------------------------------------------------------------
  FUNCTION open_input_file ( input_file_, is_xml) RESULT ( ierr )
  !-----------------------------------------------------------------------------
  !! Open file for input read, connecting it to unit \(\text{qestdin}\).  
  !! If optional variable \(\text{is_xml}\) is present, test if the file is a
  !! valid xml file.  
  !! In parallel execution, must be called by a single processor.  
  !! Module variable input_file is set to the file name actually read.
  !---------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CHARACTER (len=*), intent(in), OPTIONAL :: input_file_
  !! If \(\text{input_file_}\) is not present, read it from command line. 
  !! If \(\text{input_file_}\) is empty, the standard input is dumped to
  !! temporary file "input_tmp.in" and this is opened for read
  LOGICAL, intent(out), OPTIONAL :: is_xml
  !! \(\text{is_xml}=\text{TRUE}\) if the file has extension '.xml' or '.XML'
  !! or if either <xml...> or <?xml...> is found as first token.
  INTEGER :: ierr
  !! On exit:  
  !! \(\text{ierr} = -1\) if standard input is successfuly dumped to file;  
  !! \(\text{ierr} = 0\)  if input file is successfully opened;  
  !! \(\text{ierr} = 1\)  if there was an error opening file.
  !
  LOGICAL :: is_xml_, is_tmp
  INTEGER :: length
  CHARACTER(LEN=512) :: dummy
  LOGICAL, EXTERNAL :: test_input_xml
  !
  ! copy file to be opened into input_file
  !
  IF ( PRESENT(input_file_) ) THEN
     input_file = input_file_
  ELSE
     input_file = get_file_name ( )
  END IF
  !
  ! is_tmp: no file name read or provided, dump to "input_tmp.in"
  !
  is_tmp = ( TRIM(input_file) == ' ' )
  IF ( is_tmp ) THEN
     !
     input_file="input_tmp.in"
     OPEN(UNIT = qestdin, FILE=input_file, FORM='formatted', &
          STATUS='unknown', IOSTAT = ierr )
     IF ( ierr > 0 ) GO TO 30
     !
     dummy=' '
     WRITE(stdout, '(5x,a)') "Waiting for input..."
     DO
        READ (stdin,fmt='(A512)',END=20, ERR=30) dummy
        WRITE (qestdin,'(A)') trim(dummy)
     END DO
     !
20   CLOSE ( UNIT=qestdin, STATUS='keep' )
     !
  ENDIF
  !
  is_xml_ = PRESENT(is_xml)
  IF (is_xml_) THEN
    !
    length = LEN_TRIM(input_file)
    IF ( length > 4 ) THEN
       is_xml = ( input_file(length-3:length) == '.xml' .OR. &
                  input_file(length-3:length) == '.XML' )
    ELSE
       is_xml = .false.
    END IF
    IF ( .NOT. is_xml ) THEN
       OPEN ( UNIT = qestdin, FILE = input_file , FORM = 'FORMATTED', &
              STATUS = 'OLD', IOSTAT = ierr )
       IF ( ierr > 0 ) GO TO 30
       is_xml = test_input_xml (qestdin )
       CLOSE ( UNIT=qestdin, status='keep')
    END IF
    is_xml_ = is_xml
    !
  ENDIF
  !
  IF ( is_xml_ ) then
     IF ( is_tmp ) THEN
        WRITE(stdout, '(5x,a)') "Reading xml input from standard input"
     ELSE
        WRITE(stdout, '(5x,a)') "Reading xml input from "//TRIM(input_file)
     END IF
  ELSE 
     IF ( is_tmp ) THEN
        WRITE(stdout, '(5x,a)') "Reading input from standard input"
     ELSE
        WRITE(stdout, '(5x,a)') "Reading input from "//TRIM(input_file)
     END IF
  ENDIF
  !
  OPEN ( UNIT = qestdin, FILE = input_file, FORM = 'FORMATTED', &
          STATUS = 'OLD', IOSTAT = ierr )
  !
  IF ( ierr > 0 ) GO TO 30
  IF ( is_tmp ) ierr = -1
  RETURN
  !
30 ierr = 1
  !
  ! Do not call "errore" here: may hang in parallel execution
  !
  WRITE(stdout, "('open_input_file: fatal error opening ',A)") TRIM(input_file)
  !
END FUNCTION open_input_file

FUNCTION close_input_file ( ) RESULT ( ierr )
  !
  !! This subroutine closes the input file opened by \(\texttt{open_input_file}\),
  !! also removes temporary file if data was dumped from \(\text{stdin}\) and
  !! returns -1 if unit is not opened, 0 if no problem, > 0 if problems.
  !
  IMPLICIT NONE
  !
  INTEGER :: ierr
  LOGICAL :: opnd
  !
  INQUIRE ( qestdin, opened = opnd )
  IF (opnd) THEN
     IF ( TRIM(input_file) == "input_tmp.in") THEN
        CLOSE (UNIT=qestdin, STATUS='delete', IOSTAT=ierr )
     ELSE
        CLOSE (UNIT=qestdin, STATUS='keep', IOSTAT=ierr )
     ENDIF
  ELSE
     ierr = -1
  ENDIF 
  !
  END FUNCTION close_input_file
!
END MODULE open_close_input_file
