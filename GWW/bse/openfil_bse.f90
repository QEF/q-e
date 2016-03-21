SUBROUTINE openfil_bse()
!
! ... This routine opens all files needed to the self consistent run,
! ... sets various file names, units, record lengths
  USE wvfct,          ONLY : nbnd, npwx
  use control_flags,  ONLY:  twfcollect
  USE io_files,       ONLY : prefix, iunwfc, nwordwfc,nwordatwfc, diropn
  USE noncollin_module, ONLY : npol
  USE basis,            ONLY : natomwfc
  USE ions_base,        ONLY : nat, ityp
  USE noncollin_module,   ONLY : noncolin
  USE uspp_param,         ONLY : n_atom_wfc



  IMPLICIT NONE
  !
  LOGICAL       :: exst
  !
  !
  twfcollect=.false.
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
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  nwordatwfc = 2*npwx*natomwfc*npol

  RETURN
  !
END SUBROUTINE openfil_bse

