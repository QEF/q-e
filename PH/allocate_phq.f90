!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_phq  
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
#include "machine.h"

  use pwcom 
  use allocate 
  use parameters, only : DP 
  use phcom
  use el_phon
  implicit none 
  !
  !  allocate space for the quantities needed in the phonon program
  !
  if (lgamma) then  
     !
     !  q=0  : evq and igkq are pointers to evc and igk
     !
     evq  => evc  
     igkq => igk  
  else  
     !
     !  q!=0 : evq, igkq are allocated and calculated at point k+q
     !
     call mallocate(evq , npwx , nbndx)  
     call mallocate(igkq , npwx)  

  endif
  call mallocate(dvpsi , npwx , nbndx)  
  call mallocate(dpsi , npwx , nbndx)  

  call mallocate(vlocq , ngm , ntyp)  
  call mallocate(dmuxc , nrxx , nspin , nspin)  

  call mallocate(eigqts , nat)  
  call mallocate(rtau , 3, 48, nat)  
  call mallocate(u , 3 * nat, 3 * nat)  
  call mallocate(ubar , 3 * nat)  
  call mallocate(dyn , 3 * nat, 3 * nat)  
  call mallocate(dyn00 , 3 * nat, 3 * nat)  
  call mallocate(w2 , 3 * nat)  
  call mallocate(t ,3, 3, 48,3 * nat)  
  call mallocate(tmq ,3, 3, 3 * nat)  
  call mallocate(npert , 3 * nat)  
  call mallocate(zstareu ,3, 3,  nat)  
  call mallocate(zstareu0 ,3, 3 * nat)  
  call mallocate(zstarue ,3 , nat, 3)  
  call mallocate(zstarue0 ,3 * nat, 3)  
  if (okvan) then  
     call mallocate(int1 , nhm, nhm, 3, nat, nspin)  
     call mallocate(int2 , nhm , nhm , 3 , nat , nat)  
     call mallocate(int3 , nhm , nhm , 3 , nat , nspin)  
     call mallocate(int4 , nhm * (nhm + 1)/2,  3 , 3 , nat, nspin)
     call mallocate(int5 , nhm * (nhm + 1)/2 , 3 , 3 , nat , nat)
     call mallocate(dpqq, nhm, nhm, 3, ntyp)
     call mallocate(alphasum , nhm * (nhm + 1)/2 , 3 , nat , nspin)
  endif
  call mallocate(alphap , nkb , nbnd , 3 , nksq)  
  call mallocate(becp1 , nkb , nbnd , nksq)  

  if (elph) call mallocate(el_ph_mat, nbnd, nbnd, nks, 3*nat)
  return  
end subroutine allocate_phq
