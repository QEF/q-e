!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_fft  
  !-----------------------------------------------------------------------
  !     This routine computes the data structure associated to the FFT
  !     grid and allocate memory for all the arrays which depend upon
  !     these dimensions
  !
#include "machine.h"
  use pwcom
  use gamma
  use allocate 
  implicit none  
  !
  !     determines the data structure for fft arrays
  !
#ifdef PARA
  call data_structure_para
#else
  call data_structure_scal
#endif
  !
  if (nrxx.lt.ngm) then  
     write (6, '(/,4x," nr1=",i4," nr2= ", i4, " nr3=",i4, &
          &" nrxx = ",i8," ngm=",i8)') nr1, nr2, nr3, nrxx, ngm
     call error ('allocate_fft', 'the nr"s are too small!', 1)  

  endif
  if (nrxxs.lt.ngms) then  
     write (6, '(/,4x," nr1s=",i4," nr2s= ", i4, " nr3s=",i4, &
          &" nrxxs = ",i8," ngms=",i8)') nr1s, nr2s, nr3s, nrxxs, ngms
     call error ('allocate_fft', 'the nrs"s are too small!', 1)  

  endif
  if (ngm.le.0) call error ('allocate_fft', 'wrong ngm', 1)  
  if (ngms.le.0) call error ('allocate_fft', 'wrong ngms', 1)  
  if (nrxx.le.0) call error ('allocate_fft', 'wrong nrxx', 1)  
  if (nrxxs.le.0) call error ('allocate_fft', 'wrong nrxxs', 1)  
  if (nspin.le.0) call error ('allocate_fft', 'wrong nspin', 1)  
  !
  !     Allocate memory for all kind of stuff.
  !
  call mallocate(g, 3, ngm)  
  call mallocate(gg, ngm)  
  call mallocate(nl,  ngm)  
  call mallocate(nlm, ngm)  
  call mallocate(igtongl,  ngm)  
  call mallocate(ig1,  ngm)  
  call mallocate(ig2,  ngm)  
  call mallocate(ig3,  ngm)  
  call mallocate(rho, nrxx, nspin)  
  call mallocate(rho_save, nrxx, nspin) 
  call mallocate(vr, nrxx,nspin)  
  call mallocate(vltot, nrxx)  
  call mallocate(vnew  , nrxx, nspin)  
  call mallocate(rho_core, nrxx)  
  call mallocate(psic, nrxx)  
  call mallocate(vrs, nrxx, nspin)  
  if (doublegrid) then  
     call mallocate(nls, ngms)  
     call mallocate(nlsm,ngms)  
  else  
     nls => nl
     nlsm=> nlm
  endif
  return  
end subroutine allocate_fft

