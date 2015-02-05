!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! ... This module contains functions nd variables used to check if the code
! ... should be smoothly stopped. In order to use this module, function
! ... check_stop_init must be called (only once) at the beginning of the calc.
! ... Function check_stop_now returns .TRUE. if either the user has created
! ... an "exit" file, or if the elapsed wall time is larger than max_seconds,
! ... or if these conditions have been met in a provious call of check_stop_now.
! ... Moreover, function check_stop_now removes the exit file and sets variable
! ... stopped_by_user to .true..
!
!------------------------------------------------------------------------------!
MODULE check_stop
!------------------------------------------------------------------------------!
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  SAVE
  !
  REAL(DP) :: max_seconds = 1.E+7_DP
  REAL(DP) :: init_second
  LOGICAL :: stopped_by_user = .FALSE.
  LOGICAL, PRIVATE :: tinit = .FALSE.
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_stop_init()
       !-----------------------------------------------------------------------
       !
       USE input_parameters, ONLY : max_seconds_ => max_seconds
       USE io_global,        ONLY : stdout
       USE io_files,         ONLY : prefix, exit_file
#if defined(__TRAP_SIGUSR1) || defined(__TERMINATE_GRACEFULLY)
       USE set_signal,       ONLY : signal_trap_init
#endif
       !
       IMPLICIT NONE
       !
       REAL(DP), EXTERNAL :: cclock
       !
       IF ( tinit ) &
          WRITE( UNIT = stdout, &
                 FMT = '(/,5X,"WARNING: check_stop already initialized")' )
       !
       ! ... the exit_file name is set here
       !
       exit_file = TRIM( prefix ) // '.EXIT'
       !
       IF ( max_seconds_ > 0.0_DP ) max_seconds = max_seconds_
       !
       init_second = cclock()
       tinit   = .TRUE.
       !
#if defined(__TRAP_SIGUSR1) || defined(__TERMINATE_GRACEFULLY)
       CALL signal_trap_init ( )
#endif
       !
       RETURN
       !
     END SUBROUTINE check_stop_init
     !
     !-----------------------------------------------------------------------
     FUNCTION check_stop_now( inunit )
       !-----------------------------------------------------------------------
       !
       USE mp,         ONLY : mp_bcast
       USE mp_images,  ONLY : intra_image_comm
       USE io_global,  ONLY : ionode, ionode_id, meta_ionode, stdout
       USE io_files,   ONLY : tmp_dir, exit_file, iunexit
#if defined(__TRAP_SIGUSR1) || defined(__TERMINATE_GRACEFULLY)
       USE set_signal, ONLY : signal_detected
#endif
       !
       IMPLICIT NONE
       !
       INTEGER, OPTIONAL, INTENT(IN) :: inunit
       !
       INTEGER            :: unit
       LOGICAL            :: check_stop_now, tex=.false.
       LOGICAL            :: signaled
       REAL(DP)           :: seconds
       REAL(DP), EXTERNAL :: cclock
       !
       IF ( stopped_by_user ) THEN
          check_stop_now = .TRUE.
          RETURN
       END IF
       !
       ! ... cclock is a C function returning the elapsed solar
       ! ... time in seconds since the Epoch ( 00:00:00 1/1/1970 )
       !
       IF ( .NOT. tinit ) &
          CALL errore( 'check_stop_now', 'check_stop not initialized', 1 )
       !
       unit = stdout
       IF ( PRESENT( inunit ) ) unit = inunit
       !
       check_stop_now = .FALSE.
       !
       signaled = .FALSE.
       !
       IF ( ionode ) THEN
          !
          ! ... Check first if exit file exists in current directory
          !
          INQUIRE( FILE = TRIM( exit_file ), EXIST = tex )
          !
          IF ( tex ) THEN
             !
             check_stop_now = .TRUE.
             OPEN( UNIT = iunexit, FILE = TRIM( exit_file ) )
             CLOSE( UNIT = iunexit, STATUS = 'DELETE' )
             !
          ELSE
             !
             ! ... Check if exit file exists in scratch directory
             !
             INQUIRE( FILE = TRIM(tmp_dir) // TRIM( exit_file ), EXIST = tex )
             !
             IF ( tex ) THEN
                !
                check_stop_now = .TRUE.
                OPEN( UNIT = iunexit, FILE = TRIM(tmp_dir) // TRIM(exit_file) )
                CLOSE( UNIT = iunexit, STATUS = 'DELETE' )
                !
             ELSE
                seconds = cclock() - init_second
                check_stop_now = ( seconds  >  max_seconds )
             END IF
             !
          END IF
          !
       END IF
       !
#if defined(__TRAP_SIGUSR1) || defined(__TERMINATE_GRACEFULLY)
       signaled = signal_detected()
       check_stop_now = check_stop_now .OR. signaled
       tex = tex .OR. signaled
#endif
       !
       CALL mp_bcast( check_stop_now, ionode_id, intra_image_comm )
       !
       IF ( check_stop_now .AND. meta_ionode ) THEN
          !
          IF ( tex ) THEN
             !
             WRITE( UNIT = unit, &
                    FMT = '(/,5X,"Program stopped by user request")' )
             !
          ELSE
             !
             WRITE( UNIT = unit, &
                    FMT = '(/,5X,"Maximum CPU time exceeded")' )
             WRITE( UNIT = unit, &
                    FMT = '(/,5X,"max_seconds     = ",F10.2)' ) max_seconds
             WRITE( UNIT = unit, &
                    FMT = '(5X,"elapsed seconds = ",F10.2)' ) seconds
             !
          END IF
          !
       END IF
       !
       stopped_by_user = check_stop_now
       !
       RETURN
       !
     END FUNCTION check_stop_now
     !
END MODULE check_stop
