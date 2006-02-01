!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE restart_from_file
  !----------------------------------------------------------------------------
  !
  USE io_global,     ONLY : stdout, ionode
  USE io_files,      ONLY : iunres, tmp_dir, prefix, delete_if_present
  USE control_flags, ONLY : restart
  USE mp_global,     ONLY : mpime
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=20) :: where_restart  
    ! parameter indicating from where to restart
  LOGICAL           :: exst
  !
  !
  ! ... check if restart file is present
  !
  IF ( ionode ) THEN
     !
     iunres = 1
     !
     IF ( .NOT. restart ) THEN
        !
        !WRITE( UNIT = stdout, &
        !     & FMT = '(/5X,"RECOVER from restart file has been", &
        !     &             " switched off on input")' )
        !
        CALL delete_if_present( TRIM(tmp_dir) // TRIM(prefix) // '.restart' )
        !
        RETURN
        !
     END IF
     !
     CALL seqopn( iunres, 'restart', 'UNFORMATTED', restart )
     !
     IF ( .NOT. restart ) THEN
        !
        WRITE( UNIT = stdout, &
             & FMT = '(/5X,"RECOVER from restart file failed:", &
             &             " file not found")')
        !
        CLOSE( UNIT = iunres, STATUS = 'DELETE' )
        !
        RETURN
        !
     END IF
     !
     WRITE( UNIT = stdout, FMT = '(/5X,"read information from restart file")' )
     !
     READ( iunres, ERR = 10, END = 10 ) where_restart 
     !
     WRITE( UNIT = stdout, FMT = '(5X,"Restarting in ",A)' ) where_restart 
     !
     IF ( where_restart /= 'ELECTRONS' .AND. where_restart /= 'IONS' ) THEN
        !
        WRITE( UNIT = stdout, FMT = * ) where_restart , '......?'
        !
        CALL errore( 'restart_from_file', ' wrong recover file ', 1 )
        !
     END IF
     !
     ! ... close the file for later use
     !
     CLOSE( UNIT = iunres, STATUS = 'KEEP' )
     !
  END IF
  !
  RETURN
  !
10 CALL errore( 'restart_from_file', 'problems in reading recover file', 1 )
  !
END SUBROUTINE restart_from_file
