!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE close_phq( flag )
  !----------------------------------------------------------------------------
  !! Close all files.
  !! Called at the end of the run with \(\text{flag}\)=TRUE (removes 'recover')
  !! or during execution with \(\text{flag}\)=FALSE (does not remove 'recover').
  !
  USE mp_pools,      ONLY : me_pool, root_pool
  USE paw_variables, ONLY : okpaw
  USE io_global,     ONLY : ionode, stdout
  USE io_files,      ONLY : iunhub, iunhub_noS
  USE buffers,       ONLY : close_buffer
  USE uspp,          ONLY : okvan
  USE units_ph,      ONLY : iubar, iudrhous, iuebar, iudrho, &
                            iudvscf, iucom, iudvkb3, iuint3paw, iudyn, &
                            iundnsscf, iudvpsi, iugauge
  USE units_lr,      ONLY : iuwfc, iuatwfc, iuatswfc, iudwf
  USE control_ph,    ONLY : zue, epsil, only_wfc
  USE recover_mod,   ONLY : clean_recover
  USE output,        ONLY : fildrho, fildvscf
  USE ramanm,        ONLY : lraman, elop, iuchf, iud2w, iuba2
  USE el_phon,       ONLY : elph_mat,iunwfcwann
  USE ldaU,          ONLY : lda_plus_u
  USE control_lr,    ONLY : lgamma
  USE dvscf_interpolate, ONLY : ldvscf_interpolate, dvscf_interpol_close
  USE ahc,           ONLY : elph_ahc
  USE control_flags, ONLY : io_level
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  LOGICAL :: exst, opnd
  !
  IF (only_wfc) RETURN
  !
  IF (io_level > 0) THEN
     CALL close_buffer(iuwfc,'keep')
  ELSE
     CALL close_buffer(iuwfc,'delete')
  ENDIF
  !
  IF (flag) THEN
     CALL close_buffer(iudwf,'delete')
     CALL close_buffer(iubar,'delete')
     !
     IF ( okvan ) CALL close_buffer(iudrhous,'delete')
     !
     IF ( epsil .OR. zue ) THEN
        CALL close_buffer(iuebar,'delete')
        IF (okvan) THEN
           CALL close_buffer(iucom,'delete')
           INQUIRE( UNIT=iudvkb3, OPENED=opnd )
           IF (opnd) CLOSE( UNIT = iudvkb3, STATUS = 'DELETE' )
        ENDIF
     ENDIF
  ELSE
     CALL close_buffer(iudwf,'keep')
     CALL close_buffer(iubar,'keep')
     !
     IF ( okvan ) CALL close_buffer(iudrhous,'keep')
     !
     IF ( epsil .OR. zue ) THEN
        CALL close_buffer(iuebar,'keep')
        IF (okvan) THEN
           CALL close_buffer(iucom,'keep')
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
  IF ( fildvscf /= ' ' .AND. ionode ) THEN
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
  IF (elph_mat) THEN
    INQUIRE( UNIT=iunwfcwann, OPENED=opnd )
    IF (opnd) CLOSE( UNIT = iunwfcwann, STATUS = 'KEEP' )
  ENDIF

  IF (ionode) THEN
     INQUIRE( UNIT=iudyn, OPENED=opnd )
     IF (opnd) CLOSE( UNIT = iudyn, STATUS = 'KEEP' )
  END IF
  !
  ! DFPT+U
  IF (lda_plus_u) THEN
     CALL close_buffer(iuatwfc, 'delete')
     CALL close_buffer(iuatswfc, 'delete')
     CLOSE( UNIT = iundnsscf, STATUS = 'KEEP' )
     IF (lgamma) THEN
        CALL close_buffer(iunhub, 'delete')
        CALL close_buffer(iunhub_noS, 'delete')
     ENDIF
  ENDIF
  !
  ! dVscf Fourier interpolation
  !
  IF (flag .AND. ldvscf_interpolate) THEN
    CALL dvscf_interpol_close()
  ENDIF
  !
  ! AHC e-ph
  !
  IF (elph_ahc) THEN
    CALL close_buffer(iudvpsi, 'DELETE')
    IF (me_pool == root_pool) THEN
      CLOSE(UNIT=iugauge, STATUS='KEEP')
    ENDIF
  ENDIF
  !
  RETURN
  !
END SUBROUTINE close_phq
