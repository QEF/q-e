!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine hp_allocate_q
  !-----------------------------------------------------------------------
  !
  ! Dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
  USE ions_base,            ONLY : nat
  USE wvfct,                ONLY : nbnd, npwx
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE fft_base,             ONLY : dfftp
  USE wavefunctions,        ONLY : evc
  USE becmod,               ONLY : allocate_bec_type
  USE uspp,                 ONLY : nkb, okvan
  USE qpoint,               ONLY : nksq, eigqts
  USE lrus,                 ONLY : becp1
  USE eqv,                  ONLY : dpsi, evq, dmuxc, dvpsi, &
                                   swfcatomk, swfcatomkpq
  USE control_lr,           ONLY : lgamma
  USE basis,                ONLY : natomwfc
  USE ldaU,                 ONLY : Hubbard_lmax
  USE ldaU_hp,              ONLY : this_pert_is_on_file
  !
  IMPLICIT NONE
  INTEGER :: ik
  !
  IF (lgamma) THEN
     ! q=0 : evq is a pointer to evc
     evq  => evc
  ELSE
     ! q/=0 : evq is allocated and calculated at point k+q
     ALLOCATE (evq(npwx*npol,nbnd))
  ENDIF
  !
  ALLOCATE (dvpsi(npwx*npol,nbnd))
  ALLOCATE (dpsi(npwx*npol,nbnd))
  ALLOCATE (dmuxc(dfftp%nnr,nspin_mag,nspin_mag))
  !
  ALLOCATE (this_pert_is_on_file(nksq))
  this_pert_is_on_file(:) = .FALSE.
  !
  IF (okvan) THEN
     ALLOCATE (eigqts(nat))
     ALLOCATE (becp1(nksq))
     DO ik = 1,nksq
        CALL allocate_bec_type ( nkb, nbnd, becp1(ik) )
     ENDDO
  ENDIF
  !
  ALLOCATE (swfcatomk(npwx, natomwfc))     
  IF (lgamma) THEN
     swfcatomkpq  => swfcatomk
  ELSE
     ALLOCATE (swfcatomkpq(npwx, natomwfc))
  ENDIF
  !
  RETURN
  !
end subroutine hp_allocate_q
