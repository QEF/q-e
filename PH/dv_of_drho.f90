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
#include "machine.h"
  use funct
  use pwcom
  use parameters, only : DP
  use phcom
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

  real(kind=DP) :: qg2, fac
  ! the modulus of (q+G)^2
  ! the structure factor

  complex(kind=DP), allocatable :: dvaux (:,:), drhoc (:)
  ! auxiliary variable for potential
  !  the change of the core charge

  call start_clock ('dv_of_drho')
  allocate (dvaux( nrxx,  nspin))    
  allocate (drhoc( nrxx))    
  !
  ! the exchange-correlation contribution is computed in real space
  !
  call setv (2 * nrxx * nspin, 0.d0, dvaux, 1)

  fac = 1.d0 / float (nspin)
  if (nlcc_any.and.flag) then
     call addcore (mode, drhoc)
     do is = 1, nspin
        call DAXPY (nrxx,   fac, rho_core, 1, rho (1, is),   1)
        call DAXPY (2*nrxx, fac, drhoc,    1, dvscf (1, is), 1)
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
  if (igcx.ne.0.or.igcc.ne.0) call dgradcorr &
       (rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
       dvscf, nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, nspin, nl, ngm, g, &
       alat, omega, dvaux)
  if (nlcc_any.and.flag) then
     do is = 1, nspin
        call DAXPY (nrxx,   -fac, rho_core, 1, rho (1, is),   1)
        call DAXPY (2*nrxx, -fac, drhoc,    1, dvscf (1, is), 1)
     enddo
  endif
  !
  ! copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  do is = 2, nspin
     call DAXPY (2 * nrxx, 1.d0, dvscf(1,is), 1, dvscf(1,1), 1)
  enddo
  call cft3 (dvscf, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
  !
  ! hartree contribution is computed in reciprocal space
  !
  do is = 1, nspin
     call cft3 (dvaux (1, is), nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
     do ig = 1, ngm
        qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
        if (qg2.gt.1.d-8) then
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

  call ZCOPY (nrxx * nspin, dvaux, 1, dvscf, 1)
  !
  deallocate (drhoc)
  deallocate (dvaux)

  call stop_clock ('dv_of_drho')
  return
end subroutine dv_of_drho
