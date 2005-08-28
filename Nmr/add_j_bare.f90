!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!************************************************************************

subroutine add_j_bare (phi1, phi2, weight, rho)

!************************************************************************

  use  kinds, only : DP
  use gvect,  only : nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx
  use wvfct ,only : nbnd, npwx
  
  implicit none

  complex(DP), intent(in):: phi1(npwx,nbnd), phi2(npwx,nbnd) 
  real (DP) :: weight
  complex(DP), intent(inout):: rho(nrxx)


  !local variable
  complex(DP), allocatable :: aux1(:),aux2(:)

  integer :: i,ibnd

  allocate (aux1(nrxx))
  allocate (aux2(nrxx))

  do ibnd=1,nbnd
     
     aux1=(0.d0,0.d0)
     aux2=(0.d0,0.d0)
     aux1(1:npwx)=phi1(1:npwx,ibnd)
     aux2(1:npwx)=phi2(1:npwx,ibnd)

     call cft3(aux1,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     call cft3(aux2,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
     
     do i = 1, nrxx
        
        rho(i) = rho(i) + CONJG(aux1(i)) * aux2(i) * weight
        
     enddo
  enddo

  return

end subroutine add_j_bare

