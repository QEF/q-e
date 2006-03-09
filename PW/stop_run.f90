!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_run( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes before stopping.
  ! ... Called at the end of the run with flag = .TRUE. (removes 'restart')
  ! ... or during execution with flag = .FALSE. (does not remove 'restart')
  !
  USE io_global,          ONLY : stdout, ionode
  USE control_flags,      ONLY : lpath, lneb, lsmd, twfcollect, lconstrain, &
                                 lcoarsegrained
  USE io_files,           ONLY : prefix, iunwfc, iunigk, iunres, iunefield
  USE path_variables,     ONLY : path_deallocation
  USE path_io_routines,   ONLY : io_path_stop
  USE constraints_module, ONLY : deallocate_constraint
  USE metadyn_vars,       ONLY : deallocate_metadyn_vars
  USE mp,                 ONLY : mp_barrier, mp_end
  USE bp,                 ONLY : lelfield
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: flag
  LOGICAL             :: exst
  !
  !
  ! ... here we write all the data required to restart
  !
  CALL punch()
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
     CALL seqopn( 4, 'bfgs', 'UNFORMATTED', exst )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, 'md', 'FORMATTED', exst )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, 'para', 'FORMATTED', exst )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     CALL seqopn( 4, 'BLOCK', 'FORMATTED', exst )
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
  END IF
  !
  ! ... close unit for electric field if needed
  !
  IF ( lelfield ) CLOSE( UNIT = iunefield, STATUS = 'KEEP' )
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
  CALL clean_pw( .TRUE. )
  !
  IF ( lconstrain ) CALL deallocate_constraint()
  !
  IF ( lcoarsegrained ) CALL deallocate_metadyn_vars()
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
