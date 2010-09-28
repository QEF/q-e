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
  USE io_global,     ONLY : stdout
  !
  USE read_xml_module,       ONLY : read_xml
  USE iotk_module,           ONLY : iotk_open_read, iotk_close_read,iotk_attlenx
  !
  IMPLICIT NONE
  !
  LOGICAL, intent(inout), optional :: xmlinput
  CHARACTER (len=iotk_attlenx), intent(inout), optional :: attr
  INTEGER, intent(in), optional :: unit
  !
  LOGICAL :: xmlinput_loc
  CHARACTER (len=iotk_attlenx) :: attr_loc
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
  attr_loc = " "
  unit_loc = 5
  !
  IF(PRESENT(attr).and.(.not.PRESENT(xmlinput))) then
    CALL errore('open_input_file', 'xmlinput not present in routine call')
  ELSEIF(PRESENT(xmlinput).and.(.not.PRESENT(attr))) then
    CALL errore('open_input_file', 'attr not present in routine call')
  ENDIF
  !
  IF(PRESENT(unit)) unit_loc = unit
     !
     ! ... check if use xml input or not
     !
     xmlinput_loc = .false.
     nargs = iargc()
     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, arg )
        !
        IF ( trim( arg ) == '-xmlinput') THEN
           CALL getarg( ( iiarg + 1 ) , arg )
           xmlinput_loc = .true.
           WRITE(stdout, '(5x,a)') "Waiting for xml input..."
           CALL iotk_open_read( unit_loc, arg, attr = attr_loc, qe_syntax = .true., ierr = ierr)
           IF (ierr /= 0) CALL errore('iosys','error opening xml file', 1)
           EXIT
        ENDIF
        !
     ENDDO
     !
     xmlinput = xmlinput_loc
     !
     IF (.not.xmlinput) THEN
        CALL input_from_file(unit_loc)
        WRITE(stdout, '(5x,a)') "Waiting for input..."
     ENDIF
     !
  !
  RETURN
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
  INTEGER :: unit_loc
  !
  !
  unit_loc = 5
  xmlinput_loc = .false.
  !
  IF(PRESENT(xmlinput)) xmlinput_loc = xmlinput
  IF(PRESENT(unit)) unit_loc = unit
  !
  if(xmlinput_loc) then
    call iotk_close_read(unit=unit_loc)
  else
   close(unit_loc)
  endif 
  !
  !
END SUBROUTINE close_input_file_x

