!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE ks_hamiltonian (evc, ik, h_dim)
  !---------------------------------------------------------------------
  !
  !! This routine compute and diagonalize the KS Hamiltonian 
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : npwx, npw, et
  USE uspp,                 ONLY : nkb
  USE becmod,               ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
  USE klist,                ONLY : init_igk, xk, nkstot
  USE mp,                   ONLY : mp_sum
  USE constants,            ONLY : rytoev
  USE control_kcw,          ONLY : Hamlt, calculation, spin_component, check_ks
  USE lsda_mod,             ONLY : nspin
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)    :: ik, h_dim
  ! 
  COMPLEX(DP), INTENT(IN) :: evc(npwx,h_dim)
  !
  COMPLEX(DP) :: hpsi(npwx,h_dim), ham(h_dim,h_dim), hij, eigvc(npwx,h_dim)
  !
  REAL(DP) :: eigvl(h_dim), check
  !
  INTEGER :: iband, jband, ig, ik_eff
  !
  IF (check_ks ) WRITE(stdout,'(/,8x, "KS Hamiltonian calculation at k=", 3f12.4, 2x, " ... ")', advance="no" )  xk(:,ik)
  !
  CALL allocate_bec_type ( nkb, h_dim, becp, intra_bgrp_comm )
  !
  CALL init_igk ( npwx, ngm, g, gcutw )
  !
  call g2_kin( ik )
  hpsi(:,:) = (0.D0, 0.D0)
  !
  CALL h_psi( npwx, npw, h_dim, evc, hpsi )
  !
  ! ##### Build up the KI Hamiltonian 
  !
  ham(:,:)= (0.D0,0D0)
  !
  DO iband = 1, h_dim
     ! 
     DO jband = iband, h_dim
        !
        hij = 0.D0
        DO ig = 1, npw
           hij = hij + CONJG(evc(ig,iband)) * hpsi(ig,jband)
        ENDDO
        CALL mp_sum (hij, intra_bgrp_comm)
        !
        ham(iband,jband) = hij
        ham(jband,iband) = CONJG(ham(iband,jband))
        !
        !IF (iband==jband) WRITE(*,*) iband,jband, REAL(hij)*rytoev
     ENDDO
     !
  ENDDO
  !
  ! Store the hamiltonian in the Wannier Gauge
  !
  IF (calculation == 'ham') then 
    ik_eff = ik - (spin_component -1)*nkstot/nspin
    !WRITE(*,*) ik, ik_eff
    Hamlt(ik_eff,1:h_dim,1:h_dim) = ham(1:h_dim,1:h_dim)
  ENDIF
  !
  ! Check the eigenvalue are consistent with the PWSCF calculation
  IF (check_ks) THEN
    CALL cdiagh( h_dim, ham, h_dim, eigvl, eigvc )
    WRITE(stdout,'(2x, " DONE " ,/)')
  ENDIF
  !
  check = 0.D0
  DO iband = 1, h_dim
    check = check + (eigvl(iband)-et(iband,ik))/h_dim
  ENDDO 
  !
  IF ( check_ks ) THEN 
     !WRITE(stdout,'(/,8x, "WARNING: Eig DIFFERS! k=", 3f12.4, 3x)' )  xk(:,ik)
     WRITE( stdout, '(8X, "WANN  ",8F11.4)' ) (eigvl(iband)*rytoev, iband=1,h_dim)
     WRITE( stdout, '(8X, "PWSCF ",8F11.4)' ) (et(iband,ik)*rytoev, iband=1,h_dim)
  ENDIF
  !
  CALL deallocate_bec_type (becp)
  !
END SUBROUTINE ks_hamiltonian
