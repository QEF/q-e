!
!-----------------------------------------------------------------------
subroutine dynmat_init
  !-----------------------------------------------------------------------
  !
  !  Calculate part of the terms appearing in the dynamical matrix
  !
#include "machine.h"
  use allocate
  use pwcom
  use cgcom
  implicit none
  real(kind=DP), pointer:: dyn0(:,:),dyn1(:,:), dyncc(:,:)
  integer :: i,j, na,nb
  !
  call start_clock('dynmat_init')
  !
  call mallocate ( dyn0 , 3*nat, nmodes)
  call mallocate ( dyn1 , 3*nat, nmodes)
  call mallocate ( dyncc, 3*nat, nmodes)
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
  call mfree(dyncc)
  call mfree(dyn1 )
  call mfree(dyn0 )
  !
  call stop_clock('dynmat_init')
  !
  return
end subroutine dynmat_init
