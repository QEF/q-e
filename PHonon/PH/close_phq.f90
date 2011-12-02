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
  USE mp_global,     ONLY : me_pool
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
  IF (flag) THEN
     CLOSE( UNIT = iudwf, STATUS = 'DELETE' )
     CLOSE( UNIT = iubar, STATUS = 'DELETE' )
     !
     IF ( okvan ) CLOSE( UNIT = iudrhous, STATUS = 'DELETE' )
     !
     IF ( epsil .OR. zue ) THEN
        CLOSE( UNIT = iuebar, STATUS = 'DELETE' )
        IF (okvan) CLOSE( UNIT = iucom, STATUS = 'DELETE' )
        IF (okvan) CLOSE( UNIT = iudvkb3, STATUS = 'DELETE' )
     ENDIF
  ELSE
     CLOSE( UNIT = iudwf, STATUS = 'KEEP' )
     CLOSE( UNIT = iubar, STATUS = 'KEEP' )
     !
     IF ( okvan ) CLOSE( UNIT = iudrhous, STATUS = 'KEEP' )
     !
     IF ( epsil .OR. zue ) THEN
        CLOSE( UNIT = iuebar, STATUS = 'KEEP' )
        IF (okvan) CLOSE( UNIT = iucom, STATUS = 'KEEP' )
        IF (okvan) CLOSE( UNIT = iudvkb3, STATUS = 'KEEP' )
     ENDIF
  ENDIF
  !
  IF ( ionode .AND. &
       fildrho /= ' ') CLOSE( UNIT = iudrho, STATUS = 'KEEP' )
  !
  IF ( flag ) CALL clean_recover()
  !
  IF ( fildvscf /= ' ' ) THEN
     CLOSE( UNIT = iudvscf, STATUS = 'KEEP' )
     IF (okpaw) CLOSE( UNIT = iuint3paw, STATUS = 'KEEP' )
  ENDIF
  !
  IF (lraman .OR.elop) THEN
     CLOSE ( UNIT=iuchf, STATUS = 'keep' )
     CLOSE ( UNIT=iud2w, STATUS = 'keep' )
     CLOSE ( UNIT=iuba2, STATUS = 'keep' )
  ENDIF
  !
  CLOSE( UNIT = iunigk, STATUS = 'DELETE' )

  if(elph_mat) CLOSE( UNIT = iunwfcwann, STATUS = 'keep' ) 

  !
  RETURN
  !
END SUBROUTINE close_phq
