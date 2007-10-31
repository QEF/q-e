!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dv_of_drho (mode, dvscf, flag)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the self consistent potential
  !     due to the perturbation.
  !
#include "f_defs.h"
  use funct, only : dft_is_gradient
  use pwcom
  use scf, only : rho, rho_core
  USE kinds, only : DP
  use phcom
  implicit none

  integer :: mode
  ! input: the mode to do

  complex(DP) :: dvscf (nrxx, nspin)
  ! input: the change of the charge,
  ! output: change of the potential

  logical :: flag
  ! input: if true add core charge

  integer :: ir, is, is1, ig, nspin0, nspin1
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors

  real(DP) :: qg2, fac
  ! the modulus of (q+G)^2
  ! the structure factor

  complex(DP), allocatable :: dvaux (:,:), drhoc (:)
  ! auxiliary variable for potential
  !  the change of the core charge

  call start_clock ('dv_of_drho')

  nspin0=nspin
  nspin1=nspin
  if (nspin==4) then
     nspin0=1
     nspin1=1
     if (domag) nspin1=2
  endif

  allocate (dvaux( nrxx,  nspin))    
  allocate (drhoc( nrxx))    
  !
  ! the exchange-correlation contribution is computed in real space
  !
  dvaux (:,:) = (0.d0, 0.d0)
  if (lrpa) goto 111

  fac = 1.d0 / DBLE (nspin0)
  if (nlcc_any.and.flag) then
     call addcore (mode, drhoc)
     do is = 1, nspin0
        rho(:, is) = rho(:, is) + fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) + fac * drhoc (:)
     enddo
  endif
  do is = 1, nspin
     do is1 = 1, nspin
        do ir = 1, nrxx
           dvaux(ir,is) = dvaux(ir,is) + dmuxc(ir,is,is1) * dvscf(ir,is1)
        enddo
     enddo
  enddo
  !
  ! add gradient correction to xc, NB: if nlcc is true we need to add here
  ! its contribution. grho contains already the core charge
  !
  if ( dft_is_gradient() ) call dgradcorr &
       (rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
       dvscf, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, nspin1, &
       nl, ngm, g, alat, omega, dvaux)
  if (nlcc_any.and.flag) then
     do is = 1, nspin0
        rho(:, is) = rho(:, is) - fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) - fac * drhoc (:)
     enddo
  endif
111 continue
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
  do is = 1, nspin0
     call cft3 (dvaux (1, is), nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     do ig = 1, ngm
        qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
        if (qg2 > 1.d-8) then
           dvaux(nl(ig),is) = dvaux(nl(ig),is) + &
                              e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
        endif
     enddo
     !
     !  and transformed back to real space
     !
     call cft3 (dvaux (1, is), nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
  enddo
  !
  ! at the end the two contributes are added
  !
  dvscf (:,:) = dvaux (:,:)
  !
  deallocate (drhoc)
  deallocate (dvaux)

  call stop_clock ('dv_of_drho')
  return
end subroutine dv_of_drho
