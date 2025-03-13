SUBROUTINE openfil_bse()
!
! ... This routine opens all files needed to the self consistent run,
! ... sets various file names, units, record lengths
  USE wvfct,          ONLY : nbnd, npwx
  USE io_files,       ONLY : prefix, iunwfc, nwordwfc, diropn
  USE noncollin_module, ONLY : npol
  USE ions_base,        ONLY : nat, ityp
  USE noncollin_module, ONLY : noncolin



  IMPLICIT NONE
  !
  LOGICAL       :: exst
  !
  !
  ! ... nwordwfc is the record length for the direct-access file
  ! ... containing wavefunctions
  !
  nwordwfc = nbnd * npwx * npol
  !
  CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
  !
  IF ( .NOT. exst ) THEN
     call errore ('openfil_pw4gww','file '//TRIM( prefix )//'.wfc'//' not found',1)     
  END IF

  RETURN
  !
END SUBROUTINE openfil_bse

