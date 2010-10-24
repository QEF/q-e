!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE incdrhoscf_vdw (drhoscf, weight, ik, mode)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  !
  USE kinds,                  ONLY : DP
  USE ions_base,              ONLY : nat
!  USE wavefunctions_module,  ONLY: evc
  USE uspp_param,             ONLY : nhm
  USE eff_v,                  ONLY : evc => evc_veff
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE pwcom
  USE phcom
  IMPLICIT NONE

  INTEGER :: ik
  ! input: the k point

  real(kind=DP) :: weight
  ! input: the weight of the k point
  COMPLEX(kind=DP) :: drhoscf (dffts%nnr) , dbecsum (nhm*(nhm+1)/2,nat)
  ! output: the change of the charge densit
  ! inp/out: the accumulated dbec
  INTEGER :: mode
  !
  !   here the local variable
  !

  real(kind=DP) :: wgt
  ! the effective weight of the k point

  COMPLEX(kind=DP), ALLOCATABLE  :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  INTEGER :: ibnd, jbnd, ikk, ir, ig
  ! counters

  CALL start_clock ('incdrhoscf')
  ALLOCATE (dpsic(  dffts%nnr))
  ALLOCATE (psi  (  dffts%nnr))
  wgt = 2.d0 * weight / omega
  IF (lgamma) THEN
     ikk = ik
  ELSE
     ikk = 2 * ik - 1
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  DO ibnd = 1, nbnd_occ (ikk)
     psi (:) = (0.d0, 0.d0)
     DO ig = 1, npw
        psi (nls (igk (ig) ) ) = evc (ig, ibnd)
     ENDDO

     CALL invfft ('Wave', psi, dffts)

     dpsic(:) = (0.d0, 0.d0)
     DO ig = 1, npwq
        dpsic (nls (igkq (ig) ) ) = dpsi (ig, ibnd)
     ENDDO

     CALL invfft ('Wave', dpsic, dffts)
     DO ir = 1, dffts%nnr
!        drhoscf (ir) = drhoscf (ir) + wgt * CONJG (psi (ir) ) * dpsic (ir)
        drhoscf (ir) = drhoscf (ir) + wgt * REAL( conjg(psi(ir)) * dpsic (ir) )
     ENDDO
  ENDDO
!print*,'drhoscf=', sum(abs(REAL(drhoscf))), SUM(ABS(AIMAG(drhoscf)))
!  call addusdbec (ik, wgt, dpsi, dbecsum)
  DEALLOCATE (psi)
  DEALLOCATE (dpsic)
  CALL stop_clock ('incdrhoscf')
  RETURN
END SUBROUTINE incdrhoscf_vdw
