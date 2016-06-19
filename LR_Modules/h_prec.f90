!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine h_prec (ik, h_diag)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the precondition vector h_diag used in the solution of the
  ! ... linear system - requires h_diag to be already allocated;
  ! ... g2kin to contain kinetic energy (k+q+G)^2 for the current k-point;
  ! ... eprec to be initialized to 1.35*<Ekin> for all bands and k-points
  !
  USE kinds,      ONLY : dp
  USE klist,      ONLY : ngk
  USE qpoint,     ONLY : ikqs, ikks
  USE wvfct,      ONLY : g2kin, npwx, nbnd
  USE eqv,        ONLY : eprec
  USE control_lr, ONLY : nbnd_occ
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ik
  REAL(dp), INTENT(out) :: h_diag(npwx*npol, nbnd)
  !
  INTEGER :: ibnd, ig, ikk, ikq

  ikk = ikks(ik)
  ikq = ikqs(ik)
  
  h_diag=0.0_dp
  do ibnd = 1, nbnd_occ (ikk)
     do ig = 1, ngk(ikq)
        h_diag(ig,ibnd)=1.0_dp/max(1.0_dp,g2kin(ig)/eprec(ibnd,ik))
        !
        ! I think the following version is more correct:
        ! h_diag(ig,ibnd)=1.0_dp/max(1.0_dp,g2kin(ig)/eprec(ibnd,ikq))
        !
     enddo
     IF (noncolin) THEN
        h_diag(ig,ibnd)=1.0_dp/max(1.0_dp,g2kin(ig)/eprec(ibnd,ik))
        do ig = 1, ngk(ikq)
           h_diag(ig+npwx,ibnd)=h_diag(ig,ibnd)
        enddo
     END IF
  enddo
end subroutine h_prec
