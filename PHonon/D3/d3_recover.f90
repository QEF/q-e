!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE d3_recover (ilab, isw)
  !-----------------------------------------------------------------------
  !
  !  isw = +1 Writes d3dyn in a file for possible recover
  !  isw = -1 Starts a recover run
  !
  USE pwcom
  USE phcom
  USE d3com
  USE io_global, ONLY : ionode
  USE mp,        ONLY: mp_bcast
  USE mp_world,  ONLY: world_comm
  USE io_files,  ONLY : seqopn
  !
  IMPLICIT NONE
  !
  INTEGER :: ilab, isw
  INTEGER :: root = 0
  LOGICAL :: exst

  iunrec = 98
  IF (isw.EQ.1) THEN
     !
     IF ( .NOT. ionode ) RETURN

     CALL seqopn (iunrec, 'recv_d3', 'unformatted', exst)
     IF (ilab.LE.4) THEN
        WRITE (iunrec) ilab
     ELSE
        WRITE (iunrec) ilab, d3dyn

     ENDIF

     CLOSE (unit = iunrec, status = 'keep')
  ELSEIF (isw.EQ. - 1) THEN
     !
     IF ( ionode ) THEN
        !
        CALL seqopn (iunrec, 'recv_d3', 'unformatted', exst)
        READ (iunrec) ilab
        IF (ilab.GE.5) THEN
           REWIND (iunrec)
           READ (iunrec) ilab, d3dyn

        ENDIF
        !
        CLOSE (unit = iunrec, status = 'keep')
        !
     END IF
     !
     CALL mp_bcast (d3dyn, root, world_comm)
     CALL mp_bcast (ilab, root, world_comm)
     !
  ENDIF
  RETURN
END SUBROUTINE d3_recover
