!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE restart_from_file
  !----------------------------------------------------------------------------
  !
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE io_files,      ONLY : iunres, tmp_dir, prefix, delete_if_present, seqopn
  USE control_flags, ONLY : restart
  USE mp,            ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=20) :: where_restart  
    ! parameter indicating from where to restart
  INTEGER           :: ios
  !
  ! ... restart not required: delete restart file if present, return
  !
  IF ( .NOT. restart ) THEN
     !
     !WRITE( UNIT = stdout, &
     !     & FMT = '(/5X,"RECOVER from restart file has been", &
     !     &             " switched off on input")' )
     !
     IF ( ionode ) THEN
        !
        CALL delete_if_present( TRIM(tmp_dir) // TRIM(prefix) // '.restart' )
        !
     END IF
     !
     RETURN
     !
  END IF
  !
  ! ... restart required: check if restart file is present
  ! ... report the result of the check into variable "restart"
  !
  iunres = 1
  !
  IF ( ionode ) THEN
     !
     CALL seqopn( iunres, 'restart', 'UNFORMATTED', restart )
     !
  END IF
  !
  CALL mp_bcast ( restart, ionode_id )
  !
  IF ( .NOT. restart ) THEN
     !
     WRITE( UNIT = stdout, &
          & FMT = '(/5X,"RECOVER from restart file failed:", &
          &             " file not found")')
     !
     IF ( ionode ) THEN
        !
        CLOSE( UNIT = iunres, STATUS = 'DELETE' )
        !
     END IF
     !
     RETURN
     !
  END IF
  !
  IF ( ionode ) THEN
     !
     WRITE( UNIT = stdout, FMT = '(/5X,"read information from restart file")' )
     !
     READ( iunres, IOSTAT = ios ) where_restart 
     !
     IF ( where_restart /= 'ELECTRONS' .AND. where_restart /= 'IONS' ) THEN
        !
        ios = 1001
        !
     END IF
     ! ... close the file for later use
     !
     CLOSE( UNIT = iunres, STATUS = 'KEEP' )
     !
  END IF
  !
  CALL mp_bcast ( ios, ionode_id )
  CALL mp_bcast ( where_restart, ionode_id )
  !
  IF ( ios == 0 ) THEN
     !
     WRITE( UNIT = stdout, FMT = '(5X,"Restarting in ",A)' ) where_restart 
     !
  ELSE
     !
     CALL errore( 'restart_from_file', 'Cannot restart from here: '//TRIM(where_restart), ios)
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE restart_from_file
