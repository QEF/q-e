!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhoscf (drhoscf, weight, ik, dbecsum)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  !
#include "f_defs.h"
  USE kinds, only : DP
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat
  USE gsmooth,   ONLY : nrxxs, nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
  USE wvfct,     ONLY : npw, igk
  USE uspp_param,ONLY: nhm
  USE wavefunctions_module,  ONLY: evc
  USE eqv,       ONLY : dpsi
  USE qpoint,    ONLY : npwq, igkq, ikks
  USE control_ph, ONLY : nbnd_occ

  implicit none

  integer :: ik
  ! input: the k point

  real(DP) :: weight
  ! input: the weight of the k point
  complex(DP) :: drhoscf (nrxxs), dbecsum (nhm*(nhm+1)/2,nat)
  ! output: the change of the charge densit
  ! inp/out: the accumulated dbec
  !
  !   here the local variable
  !

  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable  :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  integer :: ibnd, jbnd, ikk, ir, ig
  ! counters

  call start_clock ('incdrhoscf')
  allocate (dpsic(  nrxxs))    
  allocate (psi  (  nrxxs))    
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  do ibnd = 1, nbnd_occ (ikk)
     psi (:) = (0.d0, 0.d0)
     do ig = 1, npw
        psi (nls (igk (ig) ) ) = evc (ig, ibnd)
     enddo

     call cft3s (psi, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)

     dpsic(:) = (0.d0, 0.d0)
     do ig = 1, npwq
        dpsic (nls (igkq (ig) ) ) = dpsi (ig, ibnd)
     enddo

     call cft3s (dpsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
     do ir = 1, nrxxs
        drhoscf (ir) = drhoscf (ir) + wgt * CONJG(psi (ir) ) * dpsic (ir)
     enddo
  enddo


  call addusdbec (ik, weight, dpsi, dbecsum)
  deallocate (psi)
  deallocate (dpsic)

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf
