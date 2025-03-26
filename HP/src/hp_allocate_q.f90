!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
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
  USE noncollin_module,     ONLY : npol, nspin_mag, noncolin, domag
  USE fft_base,             ONLY : dfftp
  USE wavefunctions,        ONLY : evc
  USE becmod,               ONLY : allocate_bec_type
  USE uspp,                 ONLY : nkb, okvan
  USE qpoint,               ONLY : nksq, eigqts
  USE lrus,                 ONLY : becp1
  USE eqv,                  ONLY : dpsi, evq, dmuxc, dvpsi
  USE control_lr,           ONLY : lgamma
  USE ldaU,                 ONLY : Hubbard_lmax, nwfcU, lda_plus_u_kind, max_num_neighbors
  USE ldaU_lr,              ONLY : swfcatomk, swfcatomkpq, vh_u_save, vh_uv_save
  USE qpoint_aux,           ONLY : becpt
  USE lr_nc_mag,            ONLY : deeq_nc_save
  USE uspp_param,           ONLY : nhm 
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
     !$acc enter data create(evq)
  ENDIF
  !
  ALLOCATE (dvpsi(npwx*npol,nbnd))
  ALLOCATE (dpsi(npwx*npol,nbnd))
  ALLOCATE (dmuxc(dfftp%nnr,nspin_mag,nspin_mag))
  !
  IF (noncolin.AND.domag) THEN
    ALLOCATE (becpt(nksq))
    DO ik=1,nksq
       CALL allocate_bec_type ( nkb, nbnd, becpt(ik) )
    ENDDO
    IF (okvan) THEN
       ALLOCATE (deeq_nc_save( nhm, nhm, nat, nspin, 2))
    ENDIF
  ENDIF
  !
  ALLOCATE (eigqts(nat))
  ALLOCATE (becp1(nksq))
  DO ik = 1,nksq
     CALL allocate_bec_type ( nkb, nbnd, becp1(ik) )
  ENDDO
  !
  ALLOCATE (swfcatomk(npwx*npol,nwfcU))     
  IF (lgamma) THEN
     swfcatomkpq  => swfcatomk
  ELSE
     ALLOCATE (swfcatomkpq(npwx*npol,nwfcU))
  ENDIF
  !
  IF (lda_plus_u_kind == 0) THEN
     ALLOCATE (vh_u_save(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, 2))
  ELSEIF (lda_plus_u_kind == 2) THEN
     ALLOCATE (vh_uv_save(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, max_num_neighbors, nat, nspin, 2))
  ENDIF
  !
  RETURN
  !
end subroutine hp_allocate_q
