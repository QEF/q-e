!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stop_pw( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  ! ... Called at the end of the run with flag = .TRUE. (removes 'restart')
  ! ... or during execution with flag = .FALSE. (does not remove 'restart')
  !
  USE io_global,        ONLY : stdout, ionode
  USE control_flags,    ONLY : lpath, lneb, lsmd, twfcollect
  USE io_files,         ONLY : prefix, iunwfc, iunigk, iunres
  USE input_parameters, ONLY : deallocate_input_parameters
  USE path_variables,   ONLY : path_deallocation
  USE path_io_routines, ONLY : io_path_stop
  USE mp,               ONLY : mp_barrier, mp_end
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  LOGICAL             :: exst
  !
  !
  IF ( lpath ) CALL io_path_stop()
  !
  ! ... iunwfc contains wavefunctions and is kept open during
  ! ... the execution - close the file and save it (or delete it 
  ! ... if the wavefunctions are already stored in the .save file)
  !
  IF ( twfcollect .AND. flag ) THEN
     !
     CLOSE( UNIT = iunwfc, STATUS = 'DELETE' )
     !
  ELSE
     !
     CLOSE( UNIT = iunwfc, STATUS = 'KEEP' )
     !
  END IF
  !      
  IF ( flag .AND. ionode ) THEN
     !
     ! ... all other files must be reopened and removed
     !
     CALL seqopn( iunres, 'restart', 'UNFORMATTED', exst )
     CLOSE( UNIT = iunres, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, TRIM( prefix ) // '.bfgs', 'UNFORMATTED', exst )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, TRIM( prefix ) // '.md', 'FORMATTED', exst )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     !CALL seqopn( 4, TRIM( prefix ) // '.para', 'FORMATTED', exst )
     !CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, TRIM( prefix ) // '.BLOCK', 'FORMATTED', exst )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
  END IF
  !
  ! ... iunigk is kept open during the execution - close and remove
  !
  CLOSE( UNIT = iunigk, STATUS = 'DELETE' )
  !
  CALL print_clock_pw()
  !
  CALL mp_barrier()
  !
  CALL mp_end()
  !
#ifdef __T3E
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
#endif
  !
  CALL clean_pw(.true.)
  !
  CALL deallocate_input_parameters()
  !
  ! ... deallocation of variables specific of "path" optimizations
  !
  IF ( lneb ) THEN
     !
     CALL path_deallocation( 'neb' )
     !
  ELSE IF ( lsmd ) THEN
     !
     CALL path_deallocation( 'smd' )
     !
  END IF
  !
  IF ( flag ) THEN
     !
     STOP
     !
  ELSE
     !
     STOP 1
     !
  END IF
  !
END SUBROUTINE stop_pw
!
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
  CALL stop_pw( .FALSE. )
  !
  RETURN
  !
END SUBROUTINE closefile
