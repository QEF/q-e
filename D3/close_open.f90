!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE close_open (isw)
  !-----------------------------------------------------------------------
  !
  ! Close and open some units. It is useful in case of interrupted run
  !
  !
  USE pwcom,     ONLY : degauss
  USE phcom,     ONLY : iudwf, lrdwf, lgamma
  USE io_files,  ONLY : prefix
  USE d3com
  USE io_global, ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER :: isw
  CHARACTER (len=256) :: filint
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
           filint = TRIM(prefix) //'.efs'
           CALL seqopn (iuef, filint, 'unformatted', exst)
        ENDIF
        !
     END IF
     CLOSE (unit = iupd0vp, status = 'keep')
     filint = TRIM(prefix) //'.p0p'
     IF (lgamma) filint = TRIM(prefix) //'.pdp'

     CALL diropn (iupd0vp, filint, lrpdqvp, exst)
     CLOSE (unit = iudwf, status = 'keep')
     filint = TRIM(prefix) //'.dwf'

     CALL diropn (iudwf, filint, lrdwf, exst)
     !
  ELSE IF (isw.EQ.1) THEN
     !
     ! This is to be used after gen_dwf(1)
     !
     IF (lgamma) CALL errore (' close_open ', ' isw=1 ; lgamma', 1)
     CLOSE (unit = iupdqvp, status = 'keep')
     filint = TRIM(prefix) //'.pdp'

     CALL diropn (iupdqvp, filint, lrpdqvp, exst)
     CLOSE (unit = iudqwf, status = 'keep')
     filint = TRIM(prefix) //'.dqwf'

     CALL diropn (iudqwf, filint, lrdwf, exst)
  ELSEIF (isw.EQ.2) THEN
     !
     ! This is to be used after gen_dwf(2)
     !
     IF (lgamma) CALL errore (' close_open ', ' isw=2 ; lgamma', 1)
     CLOSE (unit = iud0qwf, status = 'keep')
     filint = TRIM(prefix) //'.d0wf'
     CALL diropn (iud0qwf, filint, lrdwf, exst)
  ELSEIF (isw.EQ.4) THEN
     !
     ! This is to be used after gen_dpdvp
     !
     IF (degauss.EQ.0.d0) RETURN
     CLOSE (unit = iudpdvp_1, status = 'keep')
     filint = TRIM(prefix) //'.pv1'

     CALL diropn (iudpdvp_1, filint, lrdpdvp, exst)
     IF (.NOT.lgamma) THEN
        CLOSE (unit = iudpdvp_2, status = 'keep')
        filint = TRIM(prefix) //'.pv2'

        CALL diropn (iudpdvp_2, filint, lrdpdvp, exst)
        CLOSE (unit = iudpdvp_3, status = 'keep')
        filint = TRIM(prefix) //'.pv3'
        CALL diropn (iudpdvp_3, filint, lrdpdvp, exst)
     ENDIF
  ENDIF
  RETURN
END SUBROUTINE close_open
