!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE allocate_wfc()
  !----------------------------------------------------------------------------
  !
  ! ... dynamical allocation of arrays: wavefunctions and eigenvectors
  ! ... must be called after allocate_nlpot 
  !
  USE io_global, ONLY : stdout
  USE wvfct,     ONLY : gamma_only, npwx, nbnd, nbndx, et, wg
  USE klist,     ONLY : nelec, nelup, neldw, two_fermi_energies
  USE basis,     ONLY : natomwfc
  USE gvect,     ONLY : ngl
  USE uspp,      ONLY : nkb
  USE ldaU,                 ONLY : swfcatom, lda_plus_u
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  !
  IF (noncolin) THEN
     ALLOCATE( evc( npwx*npol, nbnd ) )    
  ELSE
     ALLOCATE( evc( npwx, nbnd ) )    
  ENDIF
  !
  ! ... needed for LDA+U
  !
  IF ( lda_plus_u ) ALLOCATE( swfcatom( npwx, natomwfc) )    
  !  
  WRITE( stdout, 100) nbndx, nbnd, natomwfc, npwx, nelec, nkb, ngl
  IF (two_fermi_energies) WRITE( stdout, 101) nelup, neldw
  !

100 FORMAT(/5X,'nbndx  = ',I5,'  nbnd   = ',I5,'  natomwfc = ',I5, &
       &       '  npwx   = ',I7, &
       &   /5X,'nelec  =',F7.2,'  nkb   = ',I5,'  ngl    = ',I7)
101 FORMAT( 5X,'nelup  =',F7.2,' neldw  =',F7.2)
  !       
  RETURN
  !
END subroutine allocate_wfc
