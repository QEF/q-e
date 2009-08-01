!
! Copyright (C) 2001-2008 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE allocate_wfc()
  !----------------------------------------------------------------------------
  !
  ! ... dynamical allocation of arrays: wavefunctions
  ! ... must be called after allocate_nlpot 
  !
  USE io_global, ONLY : stdout
  USE wvfct,     ONLY : npwx, nbnd
  USE basis,     ONLY : natomwfc
  USE fixed_occ, ONLY : one_atom_occupations
  USE ldaU,      ONLY : swfcatom, lda_plus_u
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc
  USE wannier_new, ONLY : use_wannier
  !
  IMPLICIT NONE
  !
  !
  IF (noncolin) THEN
     ALLOCATE( evc( npwx*npol, nbnd ) )    
     IF ( lda_plus_u .OR. one_atom_occupations ) ALLOCATE( swfcatom( npwx*npol, natomwfc) )    
  ELSE
     ALLOCATE( evc( npwx, nbnd ) )    
     IF ( lda_plus_u .OR. use_wannier .OR. one_atom_occupations ) ALLOCATE( swfcatom( npwx, natomwfc) )    
  ENDIF
  !
  RETURN
  !
END subroutine allocate_wfc
