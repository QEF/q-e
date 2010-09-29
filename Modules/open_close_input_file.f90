!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE open_input_file_x(xmlinput,attr,unit)
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine opens the input file standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  or xml input
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
  LOGICAL, intent(inout), optional :: xmlinput
  CHARACTER (len=*), intent(inout), optional :: attr
  INTEGER, intent(in), optional :: unit
  !
  LOGICAL :: xmlinput_loc,checkxml
  INTEGER :: unit_loc
  !
  INTEGER  :: iiarg, nargs, iargc, ierr
  CHARACTER (len=50) :: arg
  !
  !
#if defined(__ABSOFT)
#   define getarg getarg_
#   define iargc  iargc_
#endif
  !
  xmlinput_loc = .false.
  unit_loc = 5
  checkxml = .false.
  !
  IF(present(attr).and.(.not.present(xmlinput))) THEN
     !
     CALL errore('open_input_file', 'xmlinput not present in routine call')
     !
  ELSEIF(present(xmlinput).and.(.not.present(attr))) THEN
     !
     CALL errore('open_input_file', 'attr not present in routine call')
     !
  ENDIF
  !
  IF (present(attr).and.(present(xmlinput))) checkxml = .true.
  !
  IF(PRESENT(unit)) unit_loc = unit
  xmlinputunit = unit_loc
  !
  ! ... check if use xml input or not
  !
  xmlinput_loc = .false.
  nargs = iargc()
  !
  IF (checkxml) THEN
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, arg )
        !
        IF ( trim( arg ) == '-xmlinput') THEN
           CALL getarg( ( iiarg + 1 ) , arg )
           xmlinput_loc = .true.
           WRITE(stdout, '(5x,a)') "Waiting for xml input..."
           CALL iotk_open_read( unit_loc, arg, attr = attr, qe_syntax = .true., ierr = ierr)
           IF (ierr /= 0) CALL errore('open_input_file','error opening xml file', abs(ierr))
           EXIT
        ENDIF
        !
     ENDDO
     !
     xmlinput = xmlinput_loc
     !
  ENDIF
  !
  IF (.not.xmlinput_loc) THEN
     CALL input_from_file(unit_loc)
     WRITE(stdout, '(5x,a)') "Waiting for input..."
  ENDIF
  !
  RETURN
  !
END SUBROUTINE open_input_file_x

SUBROUTINE close_input_file_x(xmlinput,unit)
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
  LOGICAL, intent(inout), optional :: xmlinput
  INTEGER, intent(in), optional :: unit
  !
  LOGICAL :: xmlinput_loc
  LOGICAL :: opened
  INTEGER :: unit_loc, ierr
  !
  !
  unit_loc = 5
  xmlinput_loc = .false.
  !
  IF (present(xmlinput)) xmlinput_loc = xmlinput
  IF(present(unit)) unit_loc = unit
  !
  IF (xmlinput_loc) THEN
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
