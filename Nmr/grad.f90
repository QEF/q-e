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
  complex(kind=dp), intent(in) :: phi(npwx,nbnd)
  
! output

  complex(kind=dp), intent(out):: dphi(npwx,nbnd)

! local variable
  integer :: ig
  real(kind=DP), allocatable  :: gk (:,:)

  allocate (gk(3,npwx))
!  dphi=(0.d0,0.d0)

  do ig=1,npw
     gk(:,ig) = ( xk(:,ik) + g(:,igk(ig)) ) * tpiba
     dphi(ig,:) = dphi(ig,:) + phi(ig,:) * gk(dir,ig)
  enddo
  
  deallocate(gk)
end subroutine grad
