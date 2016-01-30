!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE punch( what )
  !----------------------------------------------------------------------------
  !
  ! ... This routine is called at the end of the run to save to a file
  ! ... the information needed for further processing (phonon etc.)
  !
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : prefix, iunpun, iunwfc, nwordwfc, diropn
  USE control_flags,        ONLY : io_level, twfcollect
  USE klist,                ONLY : nks
  USE pw_restart,           ONLY : pw_writefile
#ifdef __XSD
  USE pw_restart,           ONLY : pw_write_schema
#endif
  USE a2F,                  ONLY : la2F, a2Fsave
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: what
  LOGICAL :: exst
  !
  !
  IF (io_level < 0 ) RETURN
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing output data file ",A)' ) &
      TRIM( prefix ) // '.save'
  !
  ! ... if wavefunctions are stored in "distributed" format,
  ! ... save here wavefunctions to file if never saved before
  !
  IF ( .NOT. twfcollect .AND. nks == 1 ) THEN
     CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
     CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, 1 )
     CLOSE ( UNIT=iunwfc, STATUS='keep' )
  END IF
  iunpun = 4
  !
  CALL pw_writefile( TRIM( what ) )
  !
#ifdef __XSD
  CALL pw_write_schema( TRIM( what ) )
#endif
  !
  IF ( la2F ) CALL a2Fsave()
  !
  RETURN
  !
END SUBROUTINE punch
