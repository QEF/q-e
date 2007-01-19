!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE close_phq( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Close all files.
  ! ... Called at the end of the run with flag=.TRUE. (removes 'recover')
  ! ... or during execution with flag=.FALSE. (does not remove 'recover')
  !
  USE io_files,      ONLY : iunigk
  USE control_flags, ONLY : twfcollect
  USE mp_global,     ONLY : me_pool
  USE uspp,          ONLY : okvan
  USE units_ph,      ONLY : iuwfc, iudwf, iubar, iudrhous, iuebar, iudrho, &
                            iunrec, iudvscf
  USE control_ph,    ONLY : zue, epsil
  USE output,        ONLY : fildrho, fildvscf
  USE ramanm,        ONLY : lraman, elop, iuchf, iud2w, iuba2
  ! 
  IMPLICIT NONE
  !
  LOGICAL :: flag
  LOGICAL :: exst
  !
  !
  IF ( twfcollect ) THEN
     !
     CLOSE( UNIT = iuwfc, STATUS = 'DELETE' )
     !
  ELSE
     !
     CLOSE( UNIT = iuwfc, STATUS = 'KEEP' )
     !
  END IF
  !
  CLOSE( UNIT = iudwf, STATUS = 'KEEP' )
  CLOSE( UNIT = iubar, STATUS = 'KEEP' )
  !
  IF ( okvan ) CLOSE( UNIT = iudrhous, STATUS = 'KEEP' )
  !
  IF ( epsil .OR. zue ) CLOSE( UNIT = iuebar, STATUS = 'KEEP' )
  !
  IF ( me_pool == 0 .AND. &
       fildrho /= ' ') CLOSE( UNIT = iudrho, STATUS = 'KEEP' )
  !
  IF ( flag ) THEN
     !
     CALL seqopn( iunrec, 'recover', 'UNFORMATTED', exst )
     !
     CLOSE( UNIT = iunrec, STATUS = 'DELETE' )
     !
  END IF
  !
  IF ( fildvscf /= ' ' ) CLOSE( UNIT = iudvscf, STATUS = 'KEEP' )
  !
  IF (lraman .OR.elop) THEN
     CLOSE ( UNIT=iuchf, STATUS = 'keep' )
     CLOSE ( UNIT=iud2w, STATUS = 'keep' )
     CLOSE ( UNIT=iuba2, STATUS = 'keep' )
  ENDIF
  !
  CLOSE( UNIT = iunigk, STATUS = 'DELETE' )
  !
  RETURN
  !
END SUBROUTINE close_phq
