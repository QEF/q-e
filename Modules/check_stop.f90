!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This module contains functions to check if the code should
! be stopped smootly.
! In particular th function check_stop_now return .TRUE. if
! either the user has created a given file or if the 
! elapsed time is larger than max_seconds

!------------------------------------------------------------------------------!
  MODULE check_stop
!------------------------------------------------------------------------------!

    USE kinds

    IMPLICIT NONE
    SAVE

    REAL(dbl)        :: max_seconds = 1.d+7
    LOGICAL, PRIVATE :: tinit = .FALSE.

!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!


    SUBROUTINE check_stop_init( val )

      USE io_global, ONLY: ionode, ionode_id, stdout
      USE io_files, ONLY: exit_file, stopunit

      IMPLICIT NONE
      REAL(dbl), INTENT(IN) :: val
      LOGICAL :: tex
      REAL(dbl) :: seconds
      REAL(dbl) :: elapsed_seconds
      EXTERNAL  :: elapsed_seconds

      IF( tinit ) THEN
        WRITE( stdout, fmt='("WARNING: check_stop already initialized *** ")' )
      END IF

      IF ( val > 0.0d0 ) max_seconds = val

      IF ( ionode ) THEN
        INQUIRE( FILE = exit_file, EXIST = tex )
        IF ( tex ) THEN
          OPEN(  stopunit, FILE = exit_file, STATUS = 'OLD' )
          CLOSE( stopunit, STATUS = 'DELETE' )
        END IF
      END IF

      seconds = elapsed_seconds()
      tinit   = .TRUE.

      RETURN
    END SUBROUTINE
      

!------------------------------------------------------------------------------!


    FUNCTION check_stop_now()

      USE mp, ONLY: mp_bcast
      USE io_global, ONLY: ionode, ionode_id, stdout
      USE io_files, ONLY: exit_file, stopunit

      IMPLICIT NONE

      LOGICAL :: check_stop_now, tex
      REAL(dbl) :: seconds
      REAL(dbl) :: elapsed_seconds
      EXTERNAL  :: elapsed_seconds

! ...     elapsed_seconds is a C function returning the elapsed solar time in
! ...     seconds since the first call to the function itself

      IF( .NOT. tinit ) THEN
        CALL errore( ' check_stop_now ', ' check_stop not initialized ', 1 )
      END IF

      check_stop_now = .FALSE.

      IF ( ionode ) THEN

         INQUIRE( FILE = TRIM( exit_file ), EXIST = tex )
         !!!! WRITE( *, * ) 'DEBUG from check_stop: ', exit_file, tex
         IF ( tex ) THEN
           check_stop_now = .TRUE.
           WRITE( stdout, fmt='(" *** Program stopped by user request *** ")' )
         END IF

         seconds = elapsed_seconds()
         IF( seconds  >  max_seconds ) THEN
           check_stop_now = .TRUE.
           WRITE( stdout, fmt='(" *** Maximum CPU time exceeded *** ")' )
           WRITE( stdout, fmt='(" *** max_seconds     = ", D10.2, " *** ")' ) max_seconds
           WRITE( stdout, fmt='(" *** elapsed seconds = ", D10.2, " *** ")' ) seconds
         END IF

      END IF

      CALL mp_bcast( check_stop_now, ionode_id )

      RETURN
    END FUNCTION check_stop_now

!------------------------------------------------------------------------------!
  END MODULE check_stop
!------------------------------------------------------------------------------!
