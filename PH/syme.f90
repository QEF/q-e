!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

subroutine syme (dvsym)
  !---------------------------------------------------------------------
  !
  !     This routine symmetrize the change of the potential due to an
  !     electric field perturbation. It is assumed that the perturbations
  !     are on the basis of the crystal
  !
  !
#include "f_defs.h"

  use pwcom
  USE kinds, only : DP
  use phcom
  implicit none

  complex(DP) :: dvsym (nrx1, nrx2, nrx3, nspin, 3)
  complex(DP), allocatable ::  aux (:,:,:,:)
  ! the potential to symmetrize
  ! auxiliary quantity

  integer :: is, ri, rj, rk, i, j, k, irot, ipol, jpol, nspin0
  ! counter on spin polarization
  ! the rotated points
  ! the point
  ! counter on symmetries
  ! counter on polarizations

  nspin0=nspin
  if (nspin==4) nspin0=1
 
  do is = 1, nspin0
     do ipol = 1, 3
        dvsym(:,:,:,is,ipol) = CMPLX(DBLE(dvsym(:,:,:,is,ipol)),0.d0)
     end do
  end do
  if (nsym == 1) return
  allocate (aux(nrx1 , nrx2 , nrx3 , 3))    
  do is = 1, nspin0
     do ipol = 1, 3
        aux(:,:,:,ipol) = dvsym(:,:,:,is,ipol)
        dvsym(:,:,:,is,ipol) = (0.d0, 0.d0)
     enddo
     !
     !  symmmetrize
     !
     do k = 1, nr3
        do j = 1, nr2
           do i = 1, nr1
              do irot = 1, nsym
                 call ruotaijk (s(1,1,irot), ftau(1,irot), i, j, k, &
                                nr1, nr2, nr3, ri, rj, rk)
                 !
                 !    ruotaijk find the rotated of i,j,k with the inverse of S
                 !
                 do ipol = 1, 3
                    do jpol = 1, 3
                       dvsym(i,j,k,is,ipol) = dvsym(i,j,k,is,ipol) + &
                            s(ipol,jpol,irot) * aux(ri,rj,rk,jpol)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     do ipol = 1, 3
        dvsym(:,:,:,is,ipol) = dvsym(:,:,:,is,ipol) / DBLE(nsym)
     enddo
  enddo
  deallocate (aux)
  return
end subroutine syme
