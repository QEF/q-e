!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dv_of_drho_vdw (mode, dvscf, flag)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the self consistent potential
  !     due to the perturbation.
  !
#include "f_defs.h"
  use funct, only : dft_is_gradient, dmxc
  use pwcom
  USE kinds, only : DP
  use phcom
  use eff_v,  only : rho_veff
  implicit none

  integer :: mode
  ! input: the mode to do

  complex(kind=DP) :: dvscf (nrxx, nspin)
  ! input: the change of the charge,
  ! output: change of the potential

  logical :: flag
  ! input: if true add core charge

  integer :: ir, is, is1, ig
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors

  real(kind=DP) :: qg2, fac, ttd
  ! the modulus of (q+G)^2
  ! the structure factor
  ! constant 2/3

  complex(kind=DP), allocatable :: dvaux (:,:), drhoc (:),&
                                   dv_tfvw (:,:)
  ! auxiliary variable for potential
  !  the change of the core charge
  ! auxiliary variable for derivation of charge density
  !
  real(kind=dp) :: rhotot
  !
  call start_clock ('dv_of_drho')
  allocate (dvaux( nrxx,  nspin), dv_tfvw( nrxx, nspin))    
  allocate (drhoc( nrxx))    
  !
  dv_tfvw = dvscf
  !
  ! the exchange-correlation contribution is computed in real space
  !
  dvaux (:,:) = (0.d0, 0.d0)

  fac = 1.d0 / DBLE (nspin)
  if (nlcc_any.and.flag) then
     call addcore (mode, drhoc)
     do is = 1, nspin
        rho(:, is) = rho(:, is) + fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) + fac * drhoc (:)
     enddo
  endif
!  allocate ( dmuxc(nrxx, nspin, nspin) )
  dmuxc(:,:,:) = 0.d0
  do ir = 1, nrxx
     rhotot = rho (ir, nspin) + rho_core (ir)
     if (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
     if (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
  enddo
  do is = 1, nspin
     do is1 = 1, nspin
        do ir = 1, nrxx
           dvaux(ir,is) = dvaux(ir,is) + dmuxc(ir,is,is1) * dvscf(ir,is1)
        enddo
     enddo
  enddo
!  deallocate ( dmuxc )
  !
  ! add gradient correction to xc, NB: if nlcc is true we need to add here
  ! its contribution. grho contains already the core charge
  !
!  if (igcx /= 0 .or. igcc /= 0) call dgradcorr &
!       (rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
!       dvscf, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, nl, ngm, g, &
!       alat, omega, dvaux)
  if (nlcc_any.and.flag) then
     do is = 1, nspin
        rho(:, is) = rho(:, is) - fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) - fac * drhoc (:)
     enddo
  endif
  !
  ! copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  if (nspin == 2) then
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2) 
  end if
  !
  call cft3 (dvscf, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
  !
  ! hartree contribution is computed in reciprocal space
  !
  do is = 1, nspin
     call cft3 (dvaux (1, is), nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     do ig = 1, ngm
        qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
        if (qg2 > 1.d-8) then
           dvaux(nl(ig),is) = dvaux(nl(ig),is) + &
                              e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
!           dvaux(nl(ig),is) = e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
        endif
     enddo
     !
     !  and transformed back to real space
     !
     call cft3 (dvaux (1, is), nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
  enddo
  !
  ! TFvW contribution is computed in real space
  !
  ttd = 2.d0/3.d0
  is = 1
  dv_tfvw(:, is) = ttd * (0.125d0/ttd*fpi**2)**ttd * dv_tfvw (:, is) & 
                       * ( abs(rho_veff(:,is))** (-ttd/2.d0) )  
!  dv_tfvw(:, is) = ttd * (0.125d0/ttd*fpi**2)**ttd * dv_tfvw (:, is) & 
!                       * ( abs(rho(:,is)+rho_core(:))** (-ttd/2.d0) )  
  !
  ! at the end the three contributes are added
  !
  dvscf (:,:) = dvaux (:,:) + dv_tfvw (:,:)
  !
  deallocate (drhoc)
  deallocate (dvaux)

  call stop_clock ('dv_of_drho')
  return
end subroutine dv_of_drho_vdw
