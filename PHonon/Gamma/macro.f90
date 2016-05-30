!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE macro
  !----------------------------------------------------------------------
  !
  USE klist,  ONLY : nks
  USE cgcom,  ONLY : iubar, dvpsi
  USE io_files,   ONLY : seqopn
  !
  IMPLICIT NONE
  INTEGER:: ik, ipol
  CHARACTER(len=7) :: filbar
  LOGICAL :: here
  !
  DO ik=1,nks
     ! NB: this version works only for nks = 1 !
     DO ipol=1,3
        WRITE(filbar,'("filbar",i1)') ipol
        iubar=ipol
        CALL seqopn (iubar,filbar,'unformatted',here)
!!!            if (.not.here) then
        ! calculate x * psi  (if not already done)
        dvpsi(:,:) = (0.d0, 0.d0)
!!!            else
        ! otherwise restart from x * psi that is present on from file
!!!               read(iubar) dvpsi
!!!            end if
        CALL dvpsi_e(ik,ipol)
        ! write x * psi
        REWIND(iubar)
        WRITE(iubar) dvpsi
        CLOSE(unit=iubar,status='keep')
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE macro
