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
  USE paw_variables, ONLY : okpaw
  USE io_global,     ONLY : ionode, stdout
  USE uspp,          ONLY : okvan
  USE units_ph,      ONLY : iuwfc, iudwf, iubar, iudrhous, iuebar, iudrho, &
                            iudvscf, iucom, iudvkb3, iuint3paw
  USE control_ph,    ONLY : zue, epsil
  USE recover_mod,   ONLY : clean_recover
  USE output,        ONLY : fildrho, fildvscf
  USE ramanm,        ONLY : lraman, elop, iuchf, iud2w, iuba2
  USE el_phon,       ONLY : elph_mat,iunwfcwann
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  LOGICAL :: exst, opnd
  !
  !
  IF ( twfcollect ) THEN
     !
     INQUIRE( UNIT=iuwfc, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iuwfc, STATUS = 'DELETE' )
     !
  ELSE
     !
     INQUIRE( UNIT=iuwfc, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iuwfc, STATUS = 'KEEP' )
     !
  END IF
  !
  IF (flag) THEN
     INQUIRE( UNIT=iudwf, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iudwf, STATUS = 'DELETE' )
     INQUIRE( UNIT=iubar, OPENED=opnd ) 
     IF (opnd)  CLOSE( UNIT = iubar, STATUS = 'DELETE' )
     !
     IF ( okvan ) THEN
        INQUIRE( UNIT=iudrhous, OPENED=opnd ) 
        IF (opnd) CLOSE( UNIT = iudrhous, STATUS = 'DELETE' )
     ENDIF
     !
     IF ( epsil .OR. zue ) THEN
        INQUIRE( UNIT=iuebar, OPENED=opnd ) 
        IF (opnd) CLOSE( UNIT = iuebar, STATUS = 'DELETE' )
        IF (okvan) THEN
           INQUIRE( UNIT=iucom, OPENED=opnd ) 
           IF (opnd) CLOSE( UNIT = iucom, STATUS = 'DELETE' )
           INQUIRE( UNIT=iudvkb3, OPENED=opnd ) 
           IF (opnd) CLOSE( UNIT = iudvkb3, STATUS = 'DELETE' )
        ENDIF
     ENDIF
  ELSE
     INQUIRE( UNIT=iudwf, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iudwf, STATUS = 'KEEP' )
     INQUIRE( UNIT=iubar, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iubar, STATUS = 'KEEP' )
     !
     IF ( okvan ) THEN
        INQUIRE( UNIT=iudrhous, OPENED=opnd ) 
        IF (opnd) CLOSE( UNIT = iudrhous, STATUS = 'KEEP' )
     ENDIF
     !
     IF ( epsil .OR. zue ) THEN
        INQUIRE( UNIT=iuebar, OPENED=opnd ) 
        IF (opnd) CLOSE( UNIT = iuebar, STATUS = 'KEEP' )
        IF (okvan) THEN
           INQUIRE( UNIT=iucom, OPENED=opnd ) 
           IF (opnd) CLOSE( UNIT = iucom, STATUS = 'KEEP' )
           INQUIRE( UNIT=iudvkb3, OPENED=opnd ) 
           IF (opnd) CLOSE( UNIT = iudvkb3, STATUS = 'KEEP' )
        ENDIF
     ENDIF
  ENDIF
  !
  IF ( ionode .AND. fildrho /= ' ') THEN
     INQUIRE( UNIT=iudrho, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iudrho, STATUS = 'KEEP' )
  ENDIF
  !
  IF ( flag ) CALL clean_recover()
  !
  IF ( fildvscf /= ' ' ) THEN
     INQUIRE( UNIT=iudvscf, OPENED=opnd ) 
     IF (opnd) CLOSE( UNIT = iudvscf, STATUS = 'KEEP' )
     IF (okpaw) THEN
        INQUIRE( UNIT=iuint3paw, OPENED=opnd ) 
        IF (opnd) CLOSE( UNIT = iuint3paw, STATUS = 'KEEP' )
     ENDIF
  ENDIF
  !
  IF (lraman .OR.elop) THEN
     INQUIRE( UNIT=iuchf, OPENED=opnd ) 
     IF (opnd) CLOSE ( UNIT=iuchf, STATUS = 'KEEP' )
     INQUIRE( UNIT=iud2w, OPENED=opnd ) 
     IF (opnd) CLOSE ( UNIT=iud2w, STATUS = 'KEEP' )
     INQUIRE( UNIT=iuba2, OPENED=opnd ) 
     IF (opnd) CLOSE ( UNIT=iuba2, STATUS = 'KEEP' )
  ENDIF
  !
  INQUIRE( UNIT=iunigk, OPENED=opnd ) 
  IF (opnd) CLOSE( UNIT = iunigk, STATUS = 'DELETE' )

  IF (elph_mat) THEN
    INQUIRE( UNIT=iunwfcwann, OPENED=opnd ) 
    IF (opnd) CLOSE( UNIT = iunwfcwann, STATUS = 'KEEP' ) 
  ENDIF
  !
  RETURN
  !
END SUBROUTINE close_phq
