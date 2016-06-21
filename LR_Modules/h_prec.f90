!
! Copyright (C) 2016 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine h_prec (ik, evq, h_diag)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the precondition vector h_diag used in the solution of the
  ! ... linear system - On input:
  ! ... ik     index of k-point
  ! ... evq    wavefunction at k+q point
  !...  h_diag must be allocated
  ! ... g2kin  contains kinetic energy (k+q+G)^2 for the current k+q point
  !
  USE kinds,      ONLY : dp
  USE klist,      ONLY : ngk
  USE qpoint,     ONLY : ikqs, ikks
  USE wvfct,      ONLY : g2kin, npwx, nbnd
  USE gvect,      ONLY : gstart
  USE control_lr, ONLY : nbnd_occ
  USE mp,         ONLY : mp_sum
  USE mp_global,  ONLY : intra_bgrp_comm
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ik
  COMPLEX(dp), INTENT(in) :: evq(npwx*npol, nbnd)
  REAL(dp), INTENT(out) :: h_diag(npwx*npol, nbnd)
  !
  REAL(dp), ALLOCATABLE :: eprec(:)
  COMPLEX(dp), ALLOCATABLE :: aux(:)
  INTEGER :: ibnd, nbnd_, ig, ikk, ikq, npwq
  REAL(dp), EXTERNAL :: DDOT
  !
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npwq = ngk(ikq)
  nbnd_=nbnd_occ(ikk)
  
  ALLOCATE ( eprec(nbnd_) )
  ALLOCATE ( aux(npol*npwx) )
  DO ibnd = 1, nbnd_
     aux=(0.d0,0.d0)
     DO ig = 1, npwq
        aux (ig) = g2kin (ig) * evq (ig, ibnd)
     END DO
     ! NOTE: eprec(i) = 1.35*<\psi_i|Ek|\psi_i> is always real
     IF (noncolin) THEN
        DO ig = 1, npwq
           aux (ig+npwx) = g2kin (ig)* evq (ig+npwx, ibnd)
        END DO
        eprec(ibnd) = DDOT(2*npwx*npol,evq(1,ibnd),1,aux(1),1)
     ELSE IF ( gamma_only) THEN
        eprec(ibnd) = 2.0_dp*DDOT(2*npwq,evq(1,ibnd),1,aux(1),1)
        ! the following line is actually not needed
        ! because q=0 in gamma-only case, so |k+q+G|=0 for G=0
        IF (gstart==2) eprec(ibnd) = eprec(ibnd)-DBLE(evq(1,ibnd))*DBLE(aux(1))
     ELSE
        eprec(ibnd) = DDOT(2*npwq,evq(1,ibnd),1,aux(1),1)
     END IF
     eprec(ibnd) = 1.35_dp * eprec(ibnd)
     !
  END DO
  DEALLOCATE (aux)
  CALL mp_sum(eprec, intra_bgrp_comm)
  !
  h_diag=0.0_dp
  DO ibnd = 1, nbnd_
     DO ig = 1, npwq
        ! Diagonal preconditining:
        ! h_diag(G) = <Ek>/|k+q+G|^2 if |k+q+G|^2 .gt. <Ek>
        ! h_diag(G) = 1 otherwise
        ! written in this funny way because g2kin may be zero
        h_diag(ig,ibnd)=1.0_dp/max(1.0_dp,g2kin(ig)/eprec(ibnd))
     END DO
     IF (noncolin) THEN
        DO ig = 1, npwq
           h_diag(ig+npwx,ibnd)=h_diag(ig,ibnd)
        END DO
     END IF
  END DO
  DEALLOCATE (eprec)
  
END SUBROUTINE h_prec
