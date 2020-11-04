
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!
FUNCTION input_images_getarg( ) RESULT(input_images)
  !-----------------------------------------------------------------------------
  !
  ! check for command-line option "-input_images N" or "--input_images N",
  ! return N (0 if not found)
  !
  USE kinds,         ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER :: input_images
  CHARACTER(len=256) ::  myname
  INTEGER :: iiarg, nargs, i, i0
  !
  nargs = command_argument_count()
  input_images = 0
  !
  DO iiarg = 1, nargs
     !
     CALL get_command_argument( iiarg, myname)
     !
     IF ( TRIM( myname ) == '-input_images' .OR. &
          TRIM( myname ) == '--input_images' ) THEN
        !
        CALL get_command_argument( ( iiarg + 1 ) , myname )
        !
        READ(myname,*) input_images
        RETURN
        !
     END IF
     !
  ENDDO
  !
  RETURN
  !
END FUNCTION input_images_getarg

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
