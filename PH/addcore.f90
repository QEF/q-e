!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine addcore (mode, drhoc)  
!
!    This routine computes the change of the core charge
!    when the atoms moves along the given mode
!
!
#include "machine.h"

use pwcom 
use allocate 
use parameters, only : DP 
use phcom
implicit none 

!
!   The dummy variables
!

integer :: mode  
                             ! input: the mode

complex(kind=DP) :: drhoc (nrxx)  
                             ! output: the change of the core along the


integer :: nt, ig, mu, na  
                             ! counter on types
                             ! counter on G vectors
                             ! counter on modes
                             ! counter on atoms


complex(kind=DP) :: fact, gu, gu0, u1, u2, u3, gtau  
                             ! a factor
                             ! auxiliary variables
                             ! auxiliary variable for the phase
if (.not.nlcc_any) return  
!
! compute the derivative of the core charge  along the given mode
!
call setv (2 * nrxx, 0.d0, drhoc, 1)  
do na = 1, nat  
nt = ityp (na)  
if (nlcc (nt) ) then  
   fact = tpiba * (0.d0, - 1.d0) * eigqts (na)  
   mu = 3 * (na - 1)  
   if (abs (u (mu + 1, mode) ) + abs (u (mu + 2, mode) ) + abs (u &
    (mu + 3, mode) ) .gt.1.0d-12) then
      u1 = u (mu + 1, mode)  
      u2 = u (mu + 2, mode)  
      u3 = u (mu + 3, mode)  
      gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3  
      do ig = 1, ngm  
      gtau = eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) &
       * eigts3 (ig3 (ig), na)
      gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) &
       * u3
      drhoc (nl (ig) ) = drhoc (nl (ig) ) + drc (ig, nt) * gu * &
       fact * gtau
      enddo  
   endif  
endif  
enddo  
!
!   transform to real space
!

call cft3 (drhoc, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)  
return  

end subroutine addcore
