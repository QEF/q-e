!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_d3  
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the third
  ! derivative of the total energy
  !
#include "machine.h"
  use pwcom
  use phcom
  use d3com
  use allocate
  implicit none

  call allocate_phq  
  if (lgamma) then  
     vlocg0 => vlocq  
     npertg0=> npert  
     vkb0 => vkb
     ug0 => u  
     tg0 => t  
  else  
     call mallocate(vlocg0, ngm, ntyp)  
     call mallocate(ug0, 3*nat, 3*nat)  
     call mallocate(tg0, 3, 3, 48, 3*nat)  
     call mallocate(npertg0, 3*nat) 
     call mallocate(vkb0, npwx , nkb)  
  endif
  call mallocate(psidqvpsi, nbnd, nbnd)  
  call mallocate(d3dyn, 3 * nat, 3 * nat, 3 * nat)  

  if (degauss.ne.0.d0) call mallocate(ef_sh, 3 * nat)  
  call mallocate(d3dyn_aux1 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux2 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux3 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux4 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux5 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux6 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux7 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux8 , 3 * nat, 3 * nat, 3 * nat)
  call mallocate(d3dyn_aux9 , 3 * nat, 3 * nat, 3 * nat)
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux1, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux2, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux3, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux4, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux5, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux6, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux7, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux8, 1)  
  call setv (54 * nat * nat * nat, 0.d0, d3dyn_aux9, 1)  

  return  
end subroutine allocate_d3
