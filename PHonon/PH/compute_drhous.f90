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
  !
  !    This routine computes the part of the change of the charge density
  !    which is due to the orthogonalization constraint on wavefunctions
  !
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat
  USE wavefunctions_module,  ONLY: evc
  USE io_files,   ONLY : iunigk
  USE buffers,    ONLY : get_buffer
  USE uspp,       ONLY : okvan, nkb, vkb
  USE uspp_param, ONLY : nhm
  USE lsda_mod,   ONLY : lsda, nspin, current_spin, isk
  USE klist,      ONLY : xk, wk
  USE fft_base,   ONLY: dffts, dfftp
  USE fft_interfaces, ONLY: invfft
  USE gvecs,    ONLY : nls
  USE wvfct,      ONLY : npw, nbnd, igk

  USE qpoint,     ONLY : nksq, igkq, npwq, ikks, ikqs
  USE eqv,        ONLY : evq
  USE control_lr, ONLY : lgamma

  USE units_ph,   ONLY : iuwfc, lrwfc
  USE becmod,     ONLY : bec_type

  implicit none
  !
  !     the dummy variables
  !

  complex(DP) :: dbecsum (nhm * (nhm + 1) / 2, nat, nspin, 3 * nat) &
       , drhous (dfftp%nnr, nspin, 3 * nat)
  !output:the derivative of becsum
  ! output: add the orthogonality term
  type (bec_type) :: becq(nksq), & ! (nkb, nbnd)
                     alpq (3, nksq)
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k+q}

  real(DP) :: wgg (nbnd, nbnd, nksq)
  ! input: the weights

  integer :: ik, ikq, ikk, ig, nu_i, ibnd, ios
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

  if (nksq.gt.1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) then
        read (iunigk, err = 110, iostat = ios) npw, igk
110     call errore ('compute_drhous', 'reading igk', abs (ios) )
     endif
     if (lgamma) npwq = npw
     ikk = ikks(ik)
     ikq = ikqs(ik)
     weight = wk (ikk)

     if (lsda) current_spin = isk (ikk)
     if (.not.lgamma.and.nksq.gt.1) then
        read (iunigk, err = 210, iostat = ios) npwq, igkq
210     call errore ('compute_drhous', 'reading igkq', abs (ios) )
     endif
     !
     !   For each k point we construct the beta functions
     !
     call init_us_2 (npwq, igkq, xk (1, ikq), vkb)
     !
     !   Read the wavefunctions at k and transform to real space
     !

     call get_buffer (evc, lrwfc, iuwfc, ikk)
     evcr(:,:) = (0.d0, 0.d0)
     do ibnd = 1, nbnd
        do ig = 1, npw
           evcr (nls (igk (ig) ), ibnd) = evc (ig, ibnd)
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
