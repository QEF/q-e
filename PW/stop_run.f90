!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run( lflag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  ! ... Called at the end of the run with flag = .TRUE. (removes 'restart')
  ! ... or during execution with flag = .FALSE. (does not remove 'restart')
  !
  USE io_global,          ONLY : ionode
  USE mp_global,          ONLY : nimage, mp_global_end
  USE environment,        ONLY : environment_end
  USE control_flags,      ONLY : lpath, twfcollect, lconstrain, &
                                 io_level, llondon
  USE io_files,           ONLY : iunwfc, iunigk, iunefield, iunefieldm,&
                                 iunefieldp, iuntmp
  USE buffers,            ONLY : close_buffer
  USE image_io_routines,   ONLY : io_image_stop
  USE london_module,      ONLY : dealloca_london
  USE constraints_module, ONLY : deallocate_constraint
  USE input_parameters,   ONLY : deallocate_input_parameters
  USE bp,                 ONLY : lelfield
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
  LOGICAL             :: exst, opnd, lflag2
  !
  !
#if defined (EXX)
  lflag2 = lpath .or. nimage > 1
#else
  lflag2 = lpath
#endif
!  IF ( lflag2 ) THEN
     !
!     CALL io_image_stop()
     !
!  ELSE
     !
     ! ... here we write all the data required to restart
     !
!     CALL punch( 'all' )
     !
!  END IF
  !
  ! ... iunwfc contains wavefunctions and is kept open during
  ! ... the execution - close the file and save it (or delete it 
  ! ... if the wavefunctions are already stored in the .save file)
  !
  IF (lflag .and. .not. lflag2 ) THEN
     CALL seqopn( iuntmp, 'restart', 'UNFORMATTED', exst )
     CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
  ENDIF

  IF ( lflag .AND. ionode ) THEN
     !
     ! ... all other files must be reopened and removed
     !
     CALL seqopn( iuntmp, 'update', 'FORMATTED', exst )
     CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     !
     CALL seqopn( iuntmp, 'para', 'FORMATTED', exst )
     CLOSE( UNIT = iuntmp, STATUS = 'DELETE' )
     !
  END IF
  !
  CALL close_files(lflag)
  !
  CALL print_clock_pw()
  !
  CALL environment_end( 'PWSCF' )
  !
  CALL mp_global_end ()
  !
  CALL clean_pw( .TRUE. )
  !
  IF ( lflag ) THEN
     !
     STOP
     !
  ELSE
     !
     STOP 1
     !
  END IF
  !
END SUBROUTINE stop_run
!
!----------------------------------------------------------------------------
SUBROUTINE closefile()
  !----------------------------------------------------------------------------
  !
  USE io_global,  ONLY :  stdout
  !
  ! ... Close all files and synchronize processes before stopping
  ! ... Called by "sigcatch" when it receives a signal
  !
  WRITE( stdout,'(5X,"Signal Received, stopping ... ")')
  !
  CALL stop_run( .FALSE. )
  !
  RETURN
  !
END SUBROUTINE closefile
