!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE restart_in_electrons (iter, dr2, ethr, et)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY: dp
  USE io_global,     ONLY: stdout
  USE io_files,      ONLY: iunres, seqopn
  USE klist,         ONLY: nks
  USE wvfct,         ONLY: nbnd
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (inout) :: iter
  REAL(dp), INTENT(inout) :: dr2, ethr, et(nbnd,nks)
  !
  REAL(dp), ALLOCATABLE :: et_(:,:)
  REAL(dp):: dr2_, ethr_
  INTEGER :: ios
  LOGICAL :: exst
  !
  CALL seqopn (iunres, 'restart_scf', 'formatted', exst)
  IF ( exst ) THEN
     ios = 0
     READ (iunres, *, iostat=ios) iter, dr2_, ethr_
     IF ( ios /= 0 ) THEN
        iter = 0
     ELSE IF ( iter < 1 ) THEN
        iter = 0
     ELSE
        ALLOCATE (et_(nbnd,nks))
        READ (iunres, *, iostat=ios) et_
        IF ( ios /= 0 ) THEN
           iter = 0
        ELSE
           WRITE( stdout, &
           '(5x,"Calculation restarted from scf iteration #",i6)' ) iter + 1
           dr2 = dr2_
           ethr= ethr_
           et (:,:) = et_(:,:)
        END IF
        DEALLOCATE (et_)
     END IF
  ELSE
     iter = 0
  END IF
  CLOSE ( unit=iunres, status='delete')
  !
END SUBROUTINE restart_in_electrons
