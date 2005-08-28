!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*********************************************************************

subroutine grad(ik,phi,dir,dphi)

!*********************************************************************

  use kinds,     only: dp
  use wvfct,     only: npw, npwx, igk, nbnd
  use klist,     only: xk
  use gvect,     only: g
  use cell_base, only: tpiba

  implicit none

! input
  integer, intent(in) :: ik,  & ! kpoint
                         dir    ! direction of the gradient
  complex(DP), intent(in) :: phi(npwx,nbnd)
  
! output

  complex(DP), intent(out):: dphi(npwx,nbnd)

! local variable
  integer :: ig
  real(DP), allocatable  :: gk (:,:)

  allocate (gk(3,npwx))
!  dphi=(0.d0,0.d0)

  do ig=1,npw
     gk(:,ig) = ( xk(:,ik) + g(:,igk(ig)) ) * tpiba
     dphi(ig,:) = dphi(ig,:) + phi(ig,:) * gk(dir,ig)
  enddo
  
  deallocate(gk)
end subroutine grad
