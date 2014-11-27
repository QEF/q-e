!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE save_in_electrons (iter, dr2, ethr, et)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY: dp
  USE io_global,     ONLY: stdout
  USE io_files,      ONLY: iunres, seqopn
  USE klist,         ONLY: nks
  USE wvfct,         ONLY: nbnd
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: iter
  REAL(dp), INTENT(in) :: dr2, ethr, et(nbnd,nks)
  !
  LOGICAL :: exst
  !
  WRITE(stdout,'(5x,"Calculation stopped in scf loop at iteration #",i6)') iter
  CALL seqopn (iunres, 'restart_scf', 'formatted', exst)
  WRITE (iunres, *) iter, dr2, ethr
  WRITE (iunres, *) et(1:nbnd,1:nks)
  CLOSE ( unit=iunres, status='keep')
  !
END SUBROUTINE save_in_electrons
