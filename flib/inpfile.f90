!
! Copyright (C) 2002-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  INTEGER, EXTERNAL   :: iargc
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
        CALL errore( 'input_from_file', 'input file ' // TRIM( input_file ) &
             & // ' not found' , ierr )
        !
     END IF
     !
  END DO

END SUBROUTINE input_from_file
