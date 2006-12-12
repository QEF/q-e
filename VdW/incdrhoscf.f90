!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhoscf_vdw (drhoscf, weight, ik, mode)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  !
#include "f_defs.h"
  USE kinds,                  ONLY : DP
  USE ions_base,              ONLY : nat
!  USE wavefunctions_module,  ONLY: evc
  USE uspp_param,             ONLY : nhm
  USE eff_v,                  ONLY : evc => evc_veff      
  use pwcom
  use phcom
  implicit none

  integer :: ik
  ! input: the k point

  real(kind=DP) :: weight
  ! input: the weight of the k point
  complex(kind=DP) :: drhoscf (nrxxs) , dbecsum (nhm*(nhm+1)/2,nat)
  ! output: the change of the charge densit
  ! inp/out: the accumulated dbec
  integer :: mode
  !
  !   here the local variable
  !

  real(kind=DP) :: wgt
  ! the effective weight of the k point

  complex(kind=DP), allocatable  :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  integer :: ibnd, jbnd, ikk, ir, ig
  ! counters

  call start_clock ('incdrhoscf')
  allocate (dpsic(  nrxxs))    
  allocate (psi  (  nrxxs))    
  wgt = 2.d0 * weight / omega
  if (lgamma) then
     ikk = ik
  else
     ikk = 2 * ik - 1
  endif
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
!        drhoscf (ir) = drhoscf (ir) + wgt * CONJG (psi (ir) ) * dpsic (ir)
        drhoscf (ir) = drhoscf (ir) + wgt * REAL( CONJG(psi(ir)) * dpsic (ir) )
     enddo
  enddo
!print*,'drhoscf=', sum(abs(REAL(drhoscf))), SUM(ABS(AIMAG(drhoscf)))
!  call addusdbec (ik, wgt, dpsi, dbecsum)
  deallocate (psi)
  deallocate (dpsic)
  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf_vdw
