!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE open_input_file_x(lxmlinput,attr,unit)
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine opens the input file standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  or xml input. xml input can be opened on a different unit than 5.
  ! ...  ---------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  USE io_global,     ONLY : stdout, xmlinputunit
  !
  USE iotk_module,   ONLY : iotk_open_read, iotk_close_read,iotk_attlenx
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  LOGICAL, intent(inout), optional :: lxmlinput
  CHARACTER (len=*), intent(inout), optional :: attr
  INTEGER, intent(in), optional :: unit
  !
  LOGICAL :: lxmlinput_loc,lcheckxml
  INTEGER :: unit_loc
  !
  INTEGER  :: iiarg, nargs, iargc, ierr
  CHARACTER (len=50) :: arg
  !
  INTEGER :: stderr=6, stdin=5
  CHARACTER(LEN=256) :: dummy, input_file, temp_file
  INTEGER :: stdtmp
  LOGICAL :: lfound
  !
#if defined(__ABSOFT)
#   define getarg getarg_
#   define iargc  iargc_
#endif
  !
  unit_loc = 5
  lcheckxml = .false.
  !
  IF ( present(attr).and.(.not.present(lxmlinput)) ) THEN
     !
     CALL errore('open_input_file', 'xmlinput not present in routine call',1)
     !
  ELSE IF ( present(lxmlinput).and.(.not.present(attr)) ) THEN
     !
     CALL errore('open_input_file', 'attr not present in routine call',1)
     !
  ENDIF
  !
  IF (present(attr).and.(present(lxmlinput))) lcheckxml = .true.
  !
  IF(PRESENT(unit)) unit_loc = unit
  !
  xmlinputunit = unit_loc
  !
  ! ... check if use xml input or not
  !
  lxmlinput_loc = .false.
  !
  ! ... Input from file ?
  !
  temp_file="input_tmp.in"
  !
  input_file=""
  !
  lfound=.false.
  stdtmp = find_free_unit()
  OPEN(UNIT = stdtmp, FILE = trim(temp_file))
  !
  nargs = iargc()
  !
  ierr = -1
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
        OPEN ( UNIT = unit_loc, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        !
        ! copy file
        !
        IF ( ierr == 0 ) THEN
          lfound=.true.
          dummy=""
          do while (TRIM(dummy).ne."MAGICALME")
            read(unit_loc,fmt='(A256)',END=10) dummy
            write(stdtmp,*) trim(dummy)
          enddo
        ELSE IF ( ierr > 0 ) THEN
           !
           ! do not call "errore" here: it may hang in parallel execution
           ! since this routine is called by ionode only
           WRITE (stderr, '(" *** input file ",A," not found ***")' ) TRIM( input_file )
           !
        ENDIF
     END IF
     !
  END DO
10 CONTINUE
  !
  IF (lfound) THEN
     CLOSE(unit_loc)
  ELSE
     ! if no file specified then copy from standard input
     dummy=""
     ! also input_file which is used later must now be set to
     ! something sensible
     input_file="stdin"
     WRITE(stdout, '(5x,a)') "Waiting for input..."
     do while (TRIM(dummy).ne."MAGICALME")
       read(stdin,fmt='(A256)',END=20) dummy
       write(stdtmp,*) trim(dummy)
     enddo
  ENDIF
  !
20 CONTINUE
  !
  CLOSE(stdtmp)
  OPEN ( UNIT = unit_loc, FILE = trim(temp_file) , FORM = 'FORMATTED', &
          STATUS = 'OLD', IOSTAT = ierr )
  !
  IF (lcheckxml) THEN
    !
    CALL test_input_xml(unit_loc,lxmlinput_loc)
    !
    lxmlinput = lxmlinput_loc
    !
    IF(lxmlinput_loc) then
       CLOSE(unit_loc)
       WRITE(stdout, '(5x,a)') "Reading xml input from "//TRIM(input_file)
       CALL iotk_open_read( unit_loc, "input_tmp.in", attr = attr, qe_syntax = .true., ierr = ierr)
       IF (ierr /= 0) CALL errore('open_input_file','error opening xml file', abs(ierr))
    ENDIF
    !
  ENDIF
  !
  IF(.not.lxmlinput_loc) THEN
     CLOSE(unit_loc)
     WRITE(stdout, '(5x,a)') "Reading input from "//TRIM(input_file)
     OPEN ( UNIT = unit_loc, FILE = "input_tmp.in" , FORM = 'FORMATTED', &
            STATUS = 'OLD', IOSTAT = ierr )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE open_input_file_x

SUBROUTINE close_input_file_x(lxmlinput,unit)
  !
  ! ...  this subroutine close the input file for the specified unit
  ! ...  ( default is unit 5 )
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  or xml input
  ! ...  ---------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  !
  USE iotk_module,           ONLY : iotk_close_read
  !
  IMPLICIT NONE
  !
  LOGICAL, intent(inout), optional :: lxmlinput
  INTEGER, intent(in), optional :: unit
  !
  LOGICAL :: lxmlinput_loc
  LOGICAL :: opened
  INTEGER :: unit_loc, ierr
  !
  !
  unit_loc = 5
  lxmlinput_loc = .false.
  !
  IF (present(lxmlinput)) lxmlinput_loc = lxmlinput
  IF(present(unit)) unit_loc = unit
  !
  IF (lxmlinput_loc) THEN
     !
     CALL iotk_close_read(unit=unit_loc, ierr = ierr)
     IF (ierr /= 0) CALL errore('close_input_file','error closing xml file', abs(ierr) )
     !
  ELSE
     inquire( unit_loc, opened = opened )
     IF (opened) THEN
        close(unit_loc)
     ENDIF
  ENDIF 
  !
  !
END SUBROUTINE close_input_file_x
