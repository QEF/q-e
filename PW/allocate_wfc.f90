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
  USE wvfct,     ONLY : npwx, nbnd, nbndx
  USE basis,     ONLY : natomwfc
  USE ldaU,      ONLY : swfcatom, lda_plus_u
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  !
  IF (noncolin) THEN
     ALLOCATE( evc( npwx*npol, nbnd ) )    
     IF ( lda_plus_u ) ALLOCATE( swfcatom( npwx*npol, natomwfc) )    
  ELSE
     ALLOCATE( evc( npwx, nbnd ) )    
     IF ( lda_plus_u ) ALLOCATE( swfcatom( npwx, natomwfc) )    
  ENDIF
  !
  RETURN
  !
END subroutine allocate_wfc
