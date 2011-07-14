!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
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

  USE fft_base,  only : dfftp
  USE symm_base, only : nsym, s, ftau
  USE noncollin_module, only : nspin_lsda, nspin_mag
  USE kinds, only : DP
  implicit none

  complex(DP) :: dvsym (dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin_mag, 3)
  complex(DP), allocatable ::  aux (:,:,:,:)
  ! the potential to symmetrize
  ! auxiliary quantity

  integer :: is, ri, rj, rk, i, j, k, irot, ipol, jpol
  ! counter on spin polarization
  ! the rotated points
  ! the point
  ! counter on symmetries
  ! counter on polarizations

  do is = 1, nspin_lsda
     do ipol = 1, 3
        dvsym(:,:,:,is,ipol) = CMPLX(DBLE(dvsym(:,:,:,is,ipol)),0.d0,kind=DP)
     end do
  end do
  if (nsym == 1) return
  allocate (aux(dfftp%nr1x , dfftp%nr2x , dfftp%nr3x , 3))
  do is = 1, nspin_lsda
     do ipol = 1, 3
        aux(:,:,:,ipol) = dvsym(:,:,:,is,ipol)
        dvsym(:,:,:,is,ipol) = (0.d0, 0.d0)
     enddo
     !
     !  symmmetrize
     !
     do k = 1, dfftp%nr3
        do j = 1, dfftp%nr2
           do i = 1, dfftp%nr1
              do irot = 1, nsym
                 call ruotaijk (s(1,1,irot), ftau(1,irot), i, j, k, &
                                dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
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
