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
  use pwcom
  USE wavefunctions,  ONLY: evc
  use rbecmod
  implicit none
  !
  ! Allocate memory
  !
  allocate (et( nbnd, nkstot))    
  allocate (wg( nbnd, nkstot))    
  allocate (evc(npwx, nbnd))    
  allocate(becp(nkb, nbndx))
  !
  ! Needed for LDA+U
  !
  if (lda_plus_u) allocate (swfcatom( npwx, natomwfc))    

  et(:,:) = 0.d0
  write (6, 100) nbndx, nbnd, natomwfc, npwx, nelec, nkb, ngl

100 format (/5x,'nbndx  = ',i5,'  nbnd   = ',i5,'  natomwfc = ',i5, &
       &     '  npwx   = ',i7, &
       &     /5x,'nelec  = ',f7.2,' nkb   = ',i5,'  ngl    = ',i7)
  return
end subroutine allocate_wfc

