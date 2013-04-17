!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE check_restart_recover(exst_recover, exst_restart)
USE io_files, ONLY : seqopn
IMPLICIT NONE
INTEGER :: iunrec, iunres
LOGICAL :: exst_recover, exst_restart

iunrec = 99
iunres = 98
CALL seqopn (iunrec, 'recover', 'unformatted', exst_recover)
CALL seqopn( iunres, 'restart_k', 'UNFORMATTED', exst_restart )
IF (exst_recover) THEN
   close (unit = iunrec, status = 'keep')
ELSE
   close (unit = iunrec, status = 'delete')
ENDIF
IF (exst_restart) THEN
   close (unit = iunres, status = 'keep')
ELSE
   close (unit = iunres, status = 'delete')
ENDIF
RETURN
END SUBROUTINE check_restart_recover
