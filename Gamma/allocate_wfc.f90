!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_wfc
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: wavefunctions and eigenvectors
  !
#include "machine.h"
  USE io_global,      ONLY : stdout
  USE wvfct,          ONLY : npwx, nbnd, nbndx, et, wg
  USE klist,          ONLY : nkstot, nelec
  USE basis,          ONLY : natomwfc
  USE ldaU,           ONLY : swfcatom, lda_plus_u
  USE gvect,          ONLY : ngl
  USE us,             ONLY : nkb
  USE wavefunctions,  ONLY : evc
  implicit none
  !
  ! Allocate memory
  !
  allocate (et( nbnd, nkstot))    
  allocate (wg( nbnd, nkstot))    
  allocate (evc(npwx, nbnd))    
  !
  ! Needed for LDA+U
  !
  if (lda_plus_u) allocate (swfcatom( npwx, natomwfc))    

  et(:,:) = 0.d0
  WRITE( stdout, 100) nbndx, nbnd, natomwfc, npwx, nelec, nkb, ngl

100 format (/5x,'nbndx  = ',i5,'  nbnd   = ',i5,'  natomwfc = ',i5, &
       &     '  npwx   = ',i7, &
       &     /5x,'nelec  = ',f7.2,' nkb   = ',i5,'  ngl    = ',i7)
  return
end subroutine allocate_wfc

