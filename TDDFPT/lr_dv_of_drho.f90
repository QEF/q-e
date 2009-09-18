!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
! OBM
!  050608 gamma_only correction
!  150608 rho -> rho%of_r

subroutine lr_dv_of_drho ( dvscf )
  !-----------------------------------------------------------------------
  !
  !   This routine computes the change of the self consistent potential
  !   due to the perturbation.
  !
#include "f_defs.h"
  use funct,         only : dft_is_gradient
  use pwcom
  USE kinds,         only : DP
  use lr_variables
  use control_flags, only : gamma_only !this was completely absent, added
  !very very strange, rho is again used, but not declared. Adding
  use scf, only : rho,rho_core,rhog_core
   USE io_global,      ONLY : stdout
 implicit none
  !
  real(DP) :: dvscf (nrxx, nspin)
  ! input: the change of the charge
  ! output: change of the potential
  complex(DP), allocatable :: dvscf_c(:,:)
  !
  integer :: ir, is, is1, ig
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors
  !
  real(DP) :: qg2, fac
  real(kind=dp) :: charge
  ! the structure factor
  !
  complex(DP), allocatable :: dvaux (:,:)
  complex(DP), allocatable :: dvhart(:,:)
  REAL (DP) :: xq(3)
  ! auxiliary variables for potential
  ! the coordinates of the q point
  !
  If (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_dv_of_drho>")')
  endif
  call start_clock ('lr_dv')
  !
  allocate (dvaux( nrxx,  nspin))    
  dvaux (:,:) = (0.d0, 0.d0)
  charge = 0.d0
  xq = 0.d0
  !
  !   The exchange-correlation contribution is computed in real space
  !
  do is = 1, nspin
     do is1 = 1, nspin
        do ir = 1, nrxx
           dvaux(ir,is) = dvaux(ir,is) + cmplx(dmuxc(ir,is,is1) * dvscf(ir,is1),0.d0,dp)
        enddo
     enddo
  enddo
  !
  !   Add gradient correction to xc, NB: if nlcc is true we need to add here
  !   its contribution. grho contains already the core charge
  !
  fac = 1.d0 / DBLE (nspin)
  !
  if ( dft_is_gradient() ) then
     !
     allocate (dvscf_c( nrxx,  nspin))
     dvscf_c = cmplx( dvscf, 0.d0, dp)
     !
     if ( nlcc_any ) then
        !
        do is = 1, nspin
           rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core (:)
        enddo
        !
     endif
     !
     call lr_dgradcorr (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
       dvscf_c, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, nl, ngm, g, &
       alat, omega, dvaux)
     !
     if ( nlcc_any ) then
        !
        do is = 1, nspin
           rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core (:)
        enddo
        !
     endif
     !
     deallocate (dvscf_c)
     !
  end if
  !
  !   Copy the total (up+down) delta rho in dvscf(*,1) 
  !
  if (nspin == 2) then
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2) 
  end if
  !
  allocate(dvhart(nrxx,nspin))
  dvhart(:,:) = (0.d0,0.d0)
  allocate (dvscf_c( nrxx,  1))
  dvscf_c(:,1) = cmplx( dvscf(:,1), 0.d0, dp)
  !
  call cft3 (dvscf_c, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
  !
  ! hartree contribution is computed in reciprocal space
  !
  do is = 1, nspin
     !
     do ig = 1, ngm
        !
        qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
        !
        if (qg2 > 1.d-8) then
           !
           dvhart(nl(ig),is) = e2 * fpi * dvscf_c(nl(ig),1) / (tpiba2 * qg2)
           if (gamma_only) dvhart(nlm(ig),is)=conjg(dvhart(nl(ig),is))
           !
        endif
        !
     enddo
     !
     !  and transformed back to real space
     !
     call cft3 (dvhart (1, is), nr1, nr2, nr3, nrx1, nrx2, nrx3, +1)
     !
  enddo
  !
  dvscf = dble( dvaux + dvhart )
  !
  deallocate (dvaux)
  deallocate (dvhart)
  deallocate (dvscf_c)
  !
  if (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dv>")') 
  call stop_clock ('lr_dv')
  !
  return
end subroutine lr_dv_of_drho
