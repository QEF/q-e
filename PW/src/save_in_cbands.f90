!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE save_in_cbands (ik, ethr, avg_iter, et)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY: dp
  USE io_global,     ONLY: stdout
  USE io_files,      ONLY: iunres, seqopn
  USE klist,         ONLY: nks
  USE wvfct,         ONLY: nbnd
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: ik
  REAL(dp), INTENT(in) :: ethr, avg_iter, et(nbnd,nks)
  !
  LOGICAL :: exst
  !
  WRITE(stdout,'(5x,"Calculation stopped in k-point loop, point #",i6)') ik
  CALL seqopn (iunres, 'restart_k', 'formatted', exst)
  WRITE (iunres, *) ik, ethr, avg_iter
  WRITE (iunres, *) et(1:nbnd,1:nks)
  CLOSE ( unit=iunres, status='keep')
  !
END SUBROUTINE save_in_cbands
!
!-----------------------------------------------------------------------
SUBROUTINE restart_in_cbands (ik, ethr, avg_iter, et)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY: dp
  USE io_global,     ONLY: stdout
  USE io_files,      ONLY: iunres, seqopn
  USE klist,         ONLY: nks
  USE wvfct,         ONLY: nbnd
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (inout) :: ik
  REAL(dp), INTENT(inout) :: ethr, avg_iter, et(nbnd,nks)
  !
  REAL(dp), ALLOCATABLE :: et_(:,:)
  REAL(dp):: ethr_, avg_iter_
  INTEGER :: ios
  LOGICAL :: exst
  !
  CALL seqopn (iunres, 'restart_k', 'formatted', exst)
  IF ( exst ) THEN
     ios = 0
     READ (iunres, *, iostat=ios) ik, ethr_, avg_iter_
     IF ( ios /= 0 ) THEN
        ik = 0
     ELSE IF ( ik < 1 .OR. ik > nks ) THEN
        ik = 0
     ELSE
        ALLOCATE (et_(nbnd,nks))
        READ (iunres, *, iostat=ios) et_
        IF ( ios /= 0 ) THEN
           ik = 0
        ELSE
           IF ( ik == nks ) THEN
              WRITE( stdout, &
               '(5x,"Calculation restarted from end of k-point loop")' ) 
           ELSE
              WRITE( stdout, &
               '(5x,"Calculation restarted from kpoint #",i6)' ) ik + 1
           END IF
           ethr = ethr_
           avg_iter = avg_iter_
           et (:,:) = et_(:,:)
        END IF
        DEALLOCATE (et_)
     END IF
  ELSE
     ik = 0
  END IF
  CLOSE ( unit=iunres, status='delete')
  !
END SUBROUTINE restart_in_cbands
