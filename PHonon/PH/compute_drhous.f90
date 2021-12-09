!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine compute_drhous (drhous, dbecsum, wgg, becq, alpq)
  !-----------------------------------------------------------------------
  !! This routine computes the part of the change of the charge density
  !! which is due to the orthogonalization constraint on wavefunctions.
  !
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  USE wavefunctions,  ONLY: evc
  USE buffers,    ONLY : get_buffer
  USE uspp,       ONLY : okvan, nkb, vkb
  USE uspp_param, ONLY : nhm
  USE lsda_mod,   ONLY : lsda, nspin, current_spin, isk
  USE klist,      ONLY : xk, wk, ngk, igk_k
  USE fft_base,   ONLY: dffts, dfftp
  USE fft_interfaces, ONLY: invfft
  USE wvfct,      ONLY : nbnd

  USE qpoint,     ONLY : nksq, ikks, ikqs
  USE eqv,        ONLY : evq
  USE control_lr, ONLY : lgamma
  USE units_lr,   ONLY : iuwfc, lrwfc
  USE becmod,     ONLY : bec_type
  USE uspp_init,        ONLY : init_us_2

  implicit none
  !
  !     the dummy variables
  !

  complex(DP) :: dbecsum (nhm * (nhm + 1) / 2, nat, nspin, 3 * nat)
  !! output:the derivative of becsum
  complex(DP) :: drhous (dfftp%nnr, nspin, 3 * nat)
  !! output: add the orthogonality term
  type (bec_type) :: becq(nksq)  ! (nkb, nbnd)
  !! input: the becp with \(\text{psi}_{k+q}\)
  type (bec_type) :: alpq (3, nksq)
  !! input: the alphap with \(\text{psi}_{k+q}\)
  real(DP) :: wgg (nbnd, nbnd, nksq)
  !! input: the weights

  integer :: npw, npwq, ik, ikq, ikk, ig, nu_i, ibnd, ios
  ! counter on k points
  ! the point k+q
  ! record for wfcs at k point
  ! counter on spin
  ! counter on g vectors
  ! counter on modes
  ! counter on the bands
  ! integer variable for I/O control

  real(DP) :: weight
  ! the weight of the k point

  complex(DP), allocatable :: evcr (:,:)
  ! the wavefunctions in real space

  if (.not.okvan) return

  call start_clock ('com_drhous')
  allocate (evcr( dffts%nnr, nbnd))
  !
  drhous(:,:,:) = (0.d0, 0.d0)
  dbecsum (:,:,:,:) = (0.d0, 0.d0)

  do ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw = ngk(ikk)
     npwq= ngk(ikq)
     weight = wk (ikk)
     if (lsda) current_spin = isk (ikk)
     !
     !   For each k point we construct the beta functions
     !
     call init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
     !
     !   Read the wavefunctions at k and transform to real space
     !
     call get_buffer (evc, lrwfc, iuwfc, ikk)
     evcr(:,:) = (0.d0, 0.d0)
     do ibnd = 1, nbnd
        do ig = 1, npw
           evcr (dffts%nl (igk_k(ig,ikk) ), ibnd) = evc (ig, ibnd)
        enddo
        CALL invfft ('Wave', evcr (:, ibnd), dffts)
     enddo
     !
     !   Read the wavefunctions at k+q
     !
     if (.not.lgamma.and.nksq.gt.1) call get_buffer (evq, lrwfc, iuwfc, ikq)
     !
     !   And compute the contribution of this k point to the change of
     !   the charge density
     !
     do nu_i = 1, 3 * nat
        call incdrhous (drhous (1, current_spin, nu_i), weight, ik, &
             dbecsum (1, 1, current_spin, nu_i), evcr, wgg, becq, alpq, nu_i)
     enddo

  enddo

  deallocate(evcr)
  call stop_clock ('com_drhous')
  return

end subroutine compute_drhous
