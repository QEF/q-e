!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhoscf_nc (drhoscf, weight, ik, dbecsum, dpsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat
  USE cell_base, ONLY : omega
  USE fft_base,  ONLY : dffts, dfftp
  USE fft_interfaces, ONLY: invfft
  USE gvecs,   ONLY : nls
  USE lsda_mod,  ONLY : nspin
  USE spin_orb,  ONLY : domag
  USE noncollin_module, ONLY : npol, nspin_mag
  USE uspp_param,ONLY : nhm
  USE wvfct,     ONLY : npw, npwx, igk, nbnd
  USE wavefunctions_module,  ONLY: evc
  USE qpoint,    ONLY : npwq, igkq, ikks
  USE control_ph, ONLY : nbnd_occ

  implicit none

  ! I/O variables
  INTEGER, INTENT(IN) :: ik
  ! input: the k point
  REAL(DP), INTENT(IN) :: weight
  ! input: the weight of the k point
  COMPLEX(DP), INTENT(IN) :: dpsi(npwx*npol,nbnd)
  ! input: the perturbed wfcs at the given k point
  COMPLEX(DP), INTENT(INOUT) :: drhoscf (dfftp%nnr,nspin_mag), dbecsum (nhm,nhm,nat,nspin)
  ! input/output: the accumulated change of the charge density and dbecsum
  !
  !
  !   here the local variable
  !

  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable  :: psi (:,:), dpsic (:,:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  integer :: ibnd, jbnd, ikk, ir, ig
  ! counters

  call start_clock ('incdrhoscf')
  allocate (dpsic(dffts%nnr, npol))
  allocate (psi  (dffts%nnr, npol))
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  do ibnd = 1, nbnd_occ (ikk)
     psi = (0.d0, 0.d0)
     do ig = 1, npw
        psi (nls (igk (ig) ), 1) = evc (ig, ibnd)
        psi (nls (igk (ig) ), 2) = evc (ig+npwx, ibnd)
     enddo
     CALL invfft ('Wave', psi(:,1), dffts)
     CALL invfft ('Wave', psi(:,2), dffts)

     dpsic = (0.d0, 0.d0)
     do ig = 1, npwq
        dpsic (nls (igkq (ig)), 1 ) = dpsi (ig, ibnd)
        dpsic (nls (igkq (ig)), 2 ) = dpsi (ig+npwx, ibnd)
     enddo

     CALL invfft ('Wave', dpsic(:,1), dffts)
     CALL invfft ('Wave', dpsic(:,2), dffts)
     do ir = 1, dffts%nnr
        drhoscf(ir,1)=drhoscf(ir,1)+wgt*(CONJG(psi(ir,1))*dpsic(ir,1)  +  &
                                     CONJG(psi(ir,2))*dpsic(ir,2) )

     enddo
     IF (domag) THEN
        do ir = 1, dffts%nnr
           drhoscf(ir,2)=drhoscf (ir,2) + wgt *(CONJG(psi(ir,1))*dpsic(ir,2)+ &
                                             CONJG(psi(ir,2))*dpsic(ir,1) )
           drhoscf(ir,3)=drhoscf (ir,3) + wgt *(CONJG(psi(ir,1))*dpsic(ir,2)- &
                          CONJG(psi(ir,2))*dpsic(ir,1) ) * (0.d0,-1.d0)
           drhoscf(ir,4)=drhoscf (ir,4) + wgt *(CONJG(psi(ir,1))*dpsic(ir,1)- &
                                            CONJG(psi(ir,2))*dpsic(ir,2) )
        enddo
     END IF
  enddo


  call addusdbec_nc (ik, weight, dpsi, dbecsum)
  deallocate (psi)
  deallocate (dpsic)

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf_nc
