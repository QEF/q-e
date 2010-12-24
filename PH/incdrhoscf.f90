!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhoscf (drhoscf, weight, ik, dbecsum, dpsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  !
  USE kinds, only : DP
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat
  USE fft_base,  ONLY: dffts
  USE fft_interfaces, ONLY: invfft
  USE gvecs,   ONLY : nls
  USE wvfct,     ONLY : npw, igk, npwx, nbnd
  USE uspp_param,ONLY: nhm
  USE wavefunctions_module,  ONLY: evc
  USE qpoint,    ONLY : npwq, igkq, ikks
  USE control_ph, ONLY : nbnd_occ

  implicit none
  ! I/O variables
  integer, INTENT (IN) :: ik
  ! input: the k point
  real(DP), INTENT (IN) :: weight
  ! input: the weight of the k point
  complex(DP), INTENT (IN) :: dpsi (npwx,nbnd)
  ! input: the perturbed wfc for the given k point
  complex(DP), INTENT (INOUT) :: drhoscf (dffts%nnr), dbecsum (nhm*(nhm+1)/2,nat)
  ! input/output: the accumulated change to the charge density and dbecsum
  !
  !
  !   here the local variable
  !
  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable  :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  integer :: ibnd, ikk, ir, ig
  ! counters

  call start_clock ('incdrhoscf')
  allocate (dpsic(  dffts%nnr))
  allocate (psi  (  dffts%nnr))
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

     CALL invfft ('Wave', psi, dffts)

     dpsic(:) = (0.d0, 0.d0)
     do ig = 1, npwq
        dpsic (nls (igkq (ig) ) ) = dpsi (ig, ibnd)
     enddo

     CALL invfft ('Wave', dpsic, dffts)
     do ir = 1, dffts%nnr
        drhoscf (ir) = drhoscf (ir) + wgt * CONJG(psi (ir) ) * dpsic (ir)
     enddo
  enddo


  call addusdbec (ik, weight, dpsi, dbecsum)
  deallocate (psi)
  deallocate (dpsic)

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf
