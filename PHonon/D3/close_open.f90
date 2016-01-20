!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE close_open (isw)
  !-----------------------------------------------------------------------
  !
  ! Close and open some units. It is useful in case of interrupted run
  !
  !
  USE pwcom,     ONLY : degauss
  USE phcom,     ONLY : iudwf, lrdwf
  USE io_files,  ONLY : prefix, diropn, seqopn
  USE d3com
  USE io_global, ONLY : ionode

  use control_lr, ONLY : lgamma
  !
  IMPLICIT NONE
  !
  INTEGER :: isw
  CHARACTER (len=256) :: file_extension
  ! the name of the file
  LOGICAL :: exst
  ! logical variable to check file existence

  IF (LEN_TRIM(prefix) == 0) CALL errore ('close_open', 'wrong prefix', 1)
  !
  IF (isw.EQ.3) THEN
     !
     ! This is to be used after gen_dwf(3)
     !
     IF ( ionode ) THEN
        !
        IF (degauss.NE.0.d0) THEN
           CLOSE (unit = iuef, status = 'keep')
           file_extension = 'efs'
           CALL seqopn (iuef, file_extension, 'unformatted', exst)
        ENDIF
        !
     END IF
     CLOSE (unit = iupd0vp, status = 'keep')
     file_extension = 'p0p'
     IF (lgamma) file_extension = 'pdp'

     CALL diropn (iupd0vp, file_extension, lrpdqvp, exst)
     CLOSE (unit = iudwf, status = 'keep')
     file_extension = 'dwf'

     CALL diropn (iudwf, file_extension, lrdwf, exst)
     !
  ELSE IF (isw.EQ.1) THEN
     !
     ! This is to be used after gen_dwf(1)
     !
     IF (lgamma) CALL errore (' close_open ', ' isw=1 ; lgamma', 1)
     CLOSE (unit = iupdqvp, status = 'keep')
     file_extension = 'pdp'

     CALL diropn (iupdqvp, file_extension, lrpdqvp, exst)
     CLOSE (unit = iudqwf, status = 'keep')
     file_extension = 'dqwf'

     CALL diropn (iudqwf, file_extension, lrdwf, exst)
  ELSEIF (isw.EQ.2) THEN
     !
     ! This is to be used after gen_dwf(2)
     !
     IF (lgamma) CALL errore (' close_open ', ' isw=2 ; lgamma', 1)
     CLOSE (unit = iud0qwf, status = 'keep')
     file_extension = 'd0wf'
     CALL diropn (iud0qwf, file_extension, lrdwf, exst)
  ELSEIF (isw.EQ.4) THEN
     !
     ! This is to be used after gen_dpdvp
     !
     IF (degauss.EQ.0.d0) RETURN
     CLOSE (unit = iudpdvp_1, status = 'keep')
     file_extension = 'pv1'

     CALL diropn (iudpdvp_1, file_extension, lrdpdvp, exst)
     IF (.NOT.lgamma) THEN
        CLOSE (unit = iudpdvp_2, status = 'keep')
        file_extension = 'pv2'

        CALL diropn (iudpdvp_2, file_extension, lrdpdvp, exst)
        CLOSE (unit = iudpdvp_3, status = 'keep')
        file_extension = 'pv3'
        CALL diropn (iudpdvp_3, file_extension, lrdpdvp, exst)
     ENDIF
  ENDIF
  RETURN
END SUBROUTINE close_open
