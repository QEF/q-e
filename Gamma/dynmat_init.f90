!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dynmat_init
  !-----------------------------------------------------------------------
  !
  !  Calculate part of the terms appearing in the dynamical matrix
  !
#include "machine.h"
  use pwcom
  use cgcom
  implicit none
  real(kind=DP), allocatable:: dyn0(:,:),dyn1(:,:), dyncc(:,:)
  integer :: i,j, na,nb
  !
  call start_clock('dynmat_init')
  !
  allocate  ( dyn0 ( 3*nat, nmodes))    
  allocate  ( dyn1 ( 3*nat, nmodes))    
  allocate  ( dyncc( 3*nat, nmodes))    
  !
  !  first electronic contribution arising from the term  <psi|d2v|psi>
  !
  call rhod2vkb(dyn0)
  !
  !  ionic contribution
  !
  call d2ion (nat,ntyp,ityp,zv,tau,alat,omega,                      &
       at,bg,g,gg,ngm,nmodes,u,has_equivalent,dyn1)
  !
  !  core-correction contribution
  !
  call dynmatcc(dyncc)
  !
  do j=1,nmodes
     do i=1,3*nat
        dyn(i,j)=dyn0(i,j)+dyn1(i,j)+dyncc(i,j)
     end do
  end do
  !
  deallocate(dyncc)
  deallocate(dyn1 )
  deallocate(dyn0 )
  !
  call stop_clock('dynmat_init')
  !
  return
end subroutine dynmat_init
