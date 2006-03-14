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
  USE io_global,      ONLY : stdout
  USE basis,          ONLY : natomwfc
  USE wvfct,          ONLY : nbnd, npwx
  USE ldaU,           ONLY : lda_plus_U
  use control_flags,  only: twfcollect
  USE io_files,       ONLY : prefix, &
                             iunat, iunwfc, &
                             iunigk, nwordwfc, nwordatwfc
  USE restart_module,   ONLY : readfile_new
  USE mp_global,        ONLY : kunit
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  LOGICAL       :: exst
  INTEGER       :: ndr, kunittmp, ierr
  REAL(DP) :: edum(1,1), wdum(1,1)
  
  twfcollect=.false.
  !
  ! ... nwordwfc is the record length for the direct-access file
  ! ... containing wavefunctions
  !
  nwordwfc = 2 * nbnd * npwx * npol
  !
  CALL diropn( iunwfc, 'wfc', nwordwfc, exst )
  !
  IF ( .NOT. exst ) THEN
#if defined __OLDPUNCH
     ndr      = 4
#  ifdef __PARA
     kunittmp = kunit
# else
     kunittmp = 1
#  endif
     !
     CALL readfile_new( 'wave', ndr, edum, wdum, kunittmp, nwordwfc, &
                        iunwfc, ierr )
     IF ( ierr > 0 ) &
#endif
        call errore ('openfil_pp','file '//TRIM( prefix )//'.wfc'//' not found',1)     
     twfcollect=.not.exst
  END IF
  !
  RETURN
  !
END SUBROUTINE openfil_pp
