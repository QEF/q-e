!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_rec(where, irr, dr2, iter, convt, dvscfin, npe)
!-----------------------------------------------------------------------
!
!  This routine saves the information needed to recover the phonon 
!
USE kinds, ONLY : DP
USE lsda_mod,  ONLY : nspin
USE units_ph, ONLY : iunrec
USE gvect, ONLY : nrxx
USE uspp, ONLY : okvan
USE phus, ONLY : int1, int2, int3
USE control_ph, ONLY : where_rec, rec_code, reduce_io
USE ph_restart, ONLY : ph_writefile


IMPLICIT NONE
CHARACTER(LEN=10), INTENT(IN) :: where
INTEGER, INTENT(IN) :: irr, iter, npe
LOGICAL, INTENT(IN) :: convt
REAL(DP), INTENT(IN) :: dr2
COMPLEX(DP), INTENT(IN) :: dvscfin(nrxx,nspin,npe)

LOGICAL :: exst
CALL start_clock ('write_rec')
where_rec=where
CALL ph_writefile('data',0)
IF (where_rec=='done_drhod') CALL ph_writefile('data_dyn',irr)
CALL seqopn (iunrec, 'recover', 'unformatted', exst)
!
! info on current iteration (iter=0 potential mixing not available)
!
IF (reduce_io.or.convt) THEN
   WRITE (iunrec) 0, dr2
ELSE
   WRITE (iunrec) iter, dr2
ENDIF
WRITE (iunrec) dvscfin
IF (okvan) WRITE (iunrec) int1, int2, int3

CLOSE (UNIT = iunrec, STATUS = 'keep')

rec_code = 0
CALL stop_clock ('write_rec')

RETURN
END SUBROUTINE write_rec
