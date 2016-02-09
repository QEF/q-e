!
!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

! Author: L. Martin-Samos
!
!----------------------------------------------------------------------------
SUBROUTINE openfil_pw4gww()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens all files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  !
  USE kinds,          ONLY : DP
  USE wvfct,          ONLY : nbnd, npwx
  use control_flags,  ONLY:  twfcollect
  USE io_files,       ONLY : prefix, tmp_dir, iunwfc, nwordwfc, iunsat, nwordatwfc, diropn
  USE noncollin_module, ONLY : npol
  USE ldaU,             ONLY : lda_plus_u
  USE basis,            ONLY : natomwfc
  USE ions_base,        ONLY : nat, ityp
  USE noncollin_module,   ONLY : noncolin
  USE uspp_param,         ONLY : n_atom_wfc
  !
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
  !
  ! ... Needed for LDA+U
  !
  ! ... iunat  contains the (orthogonalized) atomic wfcs 
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunocc contains the atomic occupations computed in new_ns
  ! ... it is opened and closed for each reading-writing operation  
  !
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  nwordatwfc = 2*npwx*natomwfc*npol
  !
  IF ( lda_plus_u ) then
     !CALL diropn( iunat,  'atwfc',  nwordatwfc, exst )
     IF ( .NOT. exst ) THEN
        call errore ('openfil_pw4gww','file '//TRIM( prefix )//'.atwfc'//' not found',1)
     END IF

     CALL diropn( iunsat, 'satwfc', nwordatwfc, exst )
     IF ( .NOT. exst ) THEN
        call errore ('openfil_pw4gww','file '//TRIM( prefix )//'.satwfc'//' not found',1)
     END IF
  END IF
  !

  RETURN
  !
END SUBROUTINE openfil_pw4gww
