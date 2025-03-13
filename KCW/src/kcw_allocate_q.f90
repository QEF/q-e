!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine kcw_allocate_q
  !-----------------------------------------------------------------------
  !
  ! Dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
  USE ions_base,            ONLY : nat
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol, nspin_mag
  USE fft_base,             ONLY : dfftp
  USE wavefunctions,        ONLY : evc
  USE becmod,               ONLY : allocate_bec_type, becp
  USE uspp,                 ONLY : nkb, okvan
  USE qpoint,               ONLY : nksq, eigqts
  USE lrus,                 ONLY : becp1
  USE eqv,                  ONLY : dpsi, evq, dmuxc, dvpsi
  USE control_lr,           ONLY : lgamma
  USE qpoint_aux,           ONLY : becpt, alphapt
  USE noncollin_module,     ONLY : domag, noncolin
  !
  IMPLICIT NONE
  INTEGER :: ik, ipol
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
  IF (noncolin.AND.domag) THEN
     ALLOCATE (becpt(nksq))
     ALLOCATE (alphapt(3,nksq))
     DO ik=1,nksq
        CALL allocate_bec_type ( nkb, nbnd, becpt(ik) )
        DO ipol=1,3
           CALL allocate_bec_type ( nkb, nbnd, alphapt(ipol,ik) )
        ENDDO
     ENDDO
  ENDIF
  CALL allocate_bec_type ( nkb, nbnd, becp )
  ALLOCATE (eigqts(nat))
  ALLOCATE (becp1(nksq))
  DO ik = 1,nksq
     CALL allocate_bec_type ( nkb, nbnd, becp1(ik) )
  ENDDO
  !
  RETURN
  !
end subroutine kcw_allocate_q
