!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine psymrho (rho, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, &
     ftau)
  !-----------------------------------------------------------------------
  !  p-symmetrize the charge density.
  !
#include "machine.h"
#ifdef PARA
  use para
  use parameters, only : DP
  use allocate
  implicit none  
  include 'mpif.h'  
  integer :: nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau  

  real (kind=DP) :: rho (nxx)  
  real (kind=DP), pointer :: rrho (:)  
  call mallocate(rrho, nrx1 * nrx2 * nrx3)  

  call gather (rho, rrho)  
  if (me.eq.1) call symrho (rrho, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
       nsym, s, ftau)

  call scatter (rrho, rho)  

  call mfree (rrho)  
#endif
  return  
end subroutine psymrho

