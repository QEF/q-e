!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE openfil_pp()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens all files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  !
  USE kinds,          ONLY : DP
  USE wvfct,          ONLY : nbnd, npwx
  USE control_flags,  ONLY:  twfcollect
  USE io_files,       ONLY : prefix, iunwfc, nwordwfc, diropn
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  LOGICAL       :: exst
  !
  twfcollect=.false.
  !
  ! ... nwordwfc is the record length for the direct-access file
  ! ... containing wavefunctions
  !
  nwordwfc = 2 * nbnd * npwx * npol
  !
  CALL diropn( iunwfc, 'wfc', nwordwfc, exst )
  !
  IF ( .not. exst ) THEN
     CALL errore ('openfil_pp','file '//trim( prefix )//'.wfc'//' not found',1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE openfil_pp
