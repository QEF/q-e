
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE close_io_units(myunit)
  !-----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: myunit
  !
  LOGICAL :: opnd
  !
  INQUIRE( UNIT = myunit, OPENED = opnd )
  IF ( opnd ) CLOSE( UNIT = myunit )
  !
END SUBROUTINE close_io_units
!
!----------------------------------------------------------------------------
SUBROUTINE open_io_units(myunit,file_name,lappend)
  !-----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, intent(in) :: myunit
  CHARACTER(LEN=256), intent(in) :: file_name
  LOGICAL, intent(in) :: lappend
  !
  LOGICAL :: opnd
  !
  INQUIRE( UNIT = myunit, OPENED = opnd )
  IF ( opnd ) CLOSE( UNIT = myunit )
  OPEN( UNIT = myunit, FILE = TRIM(file_name), &
  STATUS = 'UNKNOWN', POSITION = 'APPEND' )
  !
END SUBROUTINE open_io_units
