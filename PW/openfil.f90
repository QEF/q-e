!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE openfil()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens some files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  ! ... All units are set in Modules/io_files.f90
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE basis,            ONLY : natomwfc, startingwfc
  USE wvfct,            ONLY : nbnd, npwx
  USE control_flags,    ONLY : order, lneb
  USE ldaU,             ONLY : lda_plus_U
  USE io_files,         ONLY : prefix, iunpun, iunat, iunwfc, iunigk, &
                               nwordwfc, nwordatwfc
  USE restart_module,   ONLY : readfile_new
  USE noncollin_module, ONLY : npol
  USE mp_global,        ONLY : kunit
  !
  IMPLICIT NONE
  !
  LOGICAL       :: exst
  INTEGER       :: ierr
  REAL(KIND=DP) :: edum(1,1), wdum(1,1)
  !
  !
  ! ... nwordwfc is the record length for the direct-access file
  ! ... containing wavefunctions
  !
  nwordwfc = 2 * nbnd * npwx * npol
  !
  CALL diropn( iunwfc, TRIM( prefix )//'.wfc', nwordwfc, exst )
  !
  IF ( startingwfc == 'file' .AND. .NOT. exst ) THEN
     !
     CALL readfile_new( 'wave', iunpun, edum, wdum, kunit, nwordwfc, &
                        iunwfc, ierr )
     !                   
     IF ( ierr > 0 ) THEN
        !
        WRITE( stdout, '(5X,"Cannot read wfc file: not found")' )
        !
        startingwfc = 'atomic'
        !
     END IF
     !
  END IF
  !
  ! ... Needed for LDA+U
  !
  ! ... iunat contains the orthogonalized wfcs
  ! ... iunocc contains the atomic occupations computed in new_ns
  ! ... it is opened and closed for each reading-writing operation  
  !
  nwordatwfc = 2 * npwx * natomwfc
  !
  IF ( lda_plus_u ) &
     CALL diropn( iunat, TRIM( prefix )//'.atwfc', nwordatwfc, exst )
  !
  ! ... iunigk contains the number of PW and the indices igk
  ! ... Note that unit 15 is reserved for error messages 
  !
  CALL seqopn( iunigk, TRIM( prefix )//'.igk', 'UNFORMATTED', exst )
  !
  RETURN
  !
END SUBROUTINE openfil
