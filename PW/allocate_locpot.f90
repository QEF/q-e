
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_locpot
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays:
  ! local potential for each kind of atom, structure factor
  !
#include "machine.h"
  use pwcom
  implicit none

  allocate (vloc( ngl, ntyp))    
  allocate (strf( ngm, ntyp))    

  allocate (eigts1( (2 * nr1 + 1), nat))    
  deallocate(eigts1)
  allocate( eigts1(-nr1:nr1,nat) )
  allocate (eigts2( (2 * nr2 + 1), nat))    
  deallocate(eigts2)
  allocate( eigts2(-nr2:nr2,nat) )
  allocate (eigts3( (2 * nr3 + 1), nat))    
  deallocate(eigts3)
  allocate( eigts3(-nr3:nr3,nat) )
  ! The above workaround is needed in order to cast the shape of eigts*
  !  (from (2*n+1,nat) --> (-n:n,nat))

  return
end subroutine allocate_locpot

