!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_part
!-----------------------------------------------------------------------
!
! dynamical allocation of arrays for the control of partial computation
! of the dynamical matrix
!
#include "machine.h"


use pwcom
use parameters, only : DP
use phcom
implicit none

!
!  allocate space for several arrays which control the run
!
allocate (comp_irr ( 3 * nat))    
allocate (ifat     ( nat))    
allocate (done_irr (  3 * nat))    
allocate (list     (  3 * nat))    
allocate (atomo    (  nat))    
return
end subroutine allocate_part
