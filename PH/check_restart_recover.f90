!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE check_restart_recover(iq_start,start_q,current_iq)
IMPLICIT NONE
INTEGER, INTENT(IN) :: start_q, current_iq
INTEGER, INTENT(OUT) :: iq_start
INTEGER :: iunrec, iunres
LOGICAL :: exst, exst1

iunrec = 99
iunres = 98
CALL seqopn (iunrec, 'recover', 'unformatted', exst)
CALL seqopn( iunres, 'restart', 'UNFORMATTED', exst1 )
IF (.not.exst.and..not.exst1) THEN
   close (unit = iunrec, status = 'delete')
   close (unit = iunres, status = 'delete')
   iq_start=start_q
ELSE
   IF (exst) THEN
      close (unit = iunrec, status = 'keep')
   ELSE
      close (unit = iunrec, status = 'delete')
   ENDIF
   IF (exst1) THEN
      close (unit = iunres, status = 'keep')
   ELSE
      close (unit = iunres, status = 'delete')
   ENDIF
   iq_start=current_iq
ENDIF
RETURN
END SUBROUTINE check_restart_recover
