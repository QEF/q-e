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
  REAL(KIND=DP) :: edum(1,1), wdum(1,1)
  
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
     ndr      = 4
     kunittmp = 1
#  ifdef __PARA
     kunittmp = kunit
#  endif
     !
     CALL readfile_new( 'wave', ndr, edum, wdum, kunittmp, nwordwfc, &
                        iunwfc, ierr )
     IF ( ierr > 0 ) &
        call errore ('openfil_pp','file '//TRIM( prefix )//'.wfc'//' not found',1)     
     twfcollect=.not.exst
  END IF
  !
  ! ... Needed for LDA+U
  !
  ! ... iunat contains the orthogonalized wfcs
  !
  iunat = 13
  nwordatwfc = 2 * npwx * natomwfc
  !
  IF ( lda_plus_u ) &
     CALL diropn( iunat, 'atwfc', nwordatwfc, exst )
  !
  ! ... iunigk contains the number of PW and the indices igk
  ! ... Note that unit 15 is reserved for error messages 
  !
  iunigk = 16
  CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )
  !
  RETURN
  !
END SUBROUTINE openfil_pp
