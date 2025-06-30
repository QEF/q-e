!
! Copyright (C) 2007-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE h_psi_meta( ldap, np, mp, psip, hpsi )
  !-----------------------------------------------------------------------
  !! This routine computes the specific contribution from the meta-GGA
  !! potential to H*psi; the result is added to hpsi.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba
  USE lsda_mod,             ONLY : nspin, current_spin
  USE wvfct,                ONLY : npwx, current_k
  USE gvect,                ONLY : g
  USE scf,                  ONLY : kedtau
  USE klist,                ONLY : xk, igk_k
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : psic
  USE fft_base,             ONLY : dffts
  USE fft_wave,             ONLY : wave_r2g, wave_g2r
  !
  IMPLICIT NONE
  !
  INTEGER :: ldap
  !! leading dimension of arrays psip, hpsip
  INTEGER :: np
  !! true dimension of psip, hpsip
  INTEGER :: mp
  !! number of states psi
  COMPLEX(DP) :: psip(ldap,mp)
  !! the wavefunction
  COMPLEX(DP) :: hpsi(ldap,mp)
  !! Hamiltonian dot psip
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: psi_g(:,:)
  INTEGER :: im, i, j, nrxxs, ebnd, brange, dim_g
  REAL(DP) :: kplusgi, fac
  COMPLEX(DP), PARAMETER :: ci=(0.d0,1.d0)
  !
  CALL start_clock( 'h_psi_meta' )
  !
  nrxxs = dffts%nnr
  dim_g = 1
  IF (gamma_only) dim_g = 2
  !
  ALLOCATE( psi_g(npwx,dim_g) )
  !$acc enter data create(psi_g, psic) copyin(kedtau) 
  !$acc data present_or_copyin(psip) present_or_copyout(hpsi) 
  !
  IF (gamma_only) THEN
     !
     ! ... Gamma algorithm
     !
     DO im = 1, mp, 2
        !
        fac = 1.d0
        IF ( im < mp ) fac = 0.5d0
        !
        DO j = 1, 3
           !
           !$acc parallel loop  
           DO i = 1, np
              kplusgi = (xk(j,current_k)+g(j,i)) * tpiba
              psi_g(i,1) = CMPLX(0._DP,kplusgi,KIND=DP) * psip(i,im)
              IF ( im < mp ) psi_g(i,2) = CMPLX(0._DP,kplusgi,KIND=DP) * psip(i,im+1)
           ENDDO
           !
           ebnd = im
           IF ( im < mp ) ebnd = ebnd + 1
           brange = ebnd-im+1
           !
           CALL wave_g2r( psi_g(1:np,1:brange), psic, dffts )
           !
           !$acc kernels 
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs)
           !$acc end kernels
           !
           CALL wave_r2g( psic(1:dffts%nnr), psi_g(:,1:brange), dffts )
           !
           !$acc parallel loop 
           DO i = 1, np
              kplusgi = (xk(j,current_k)+g(j,i)) * tpiba
              hpsi(i,im) = hpsi(i,im) - ci * kplusgi * fac * psi_g(i,1)
              IF ( im < mp ) hpsi(i,im+1) = hpsi(i,im+1) - ci * kplusgi * fac * psi_g(i,2)
           ENDDO
           !
        ENDDO
     ENDDO
     !
  ELSE
     !
     ! ... generic k algorithm
     !
     DO im = 1, mp
        DO j = 1, 3
           !
           !$acc parallel loop 
           DO i = 1, np
              kplusgi = (xk(j,current_k)+g(j,igk_k(i,current_k)))*tpiba
              psi_g(i,1) = CMPLX(0._DP,kplusgi,KIND=DP) * psip(i,im)
           ENDDO
           !
           CALL wave_g2r( psi_g(1:np,1:1), psic, dffts, igk=igk_k(:,current_k) )
           !
           !$acc kernels 
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs)
           !$acc end kernels 
           !
           CALL wave_r2g( psic(1:dffts%nnr), psi_g(1:np,1:1), dffts, igk=igk_k(:,current_k) )
           !
           !$acc parallel loop 
           DO i = 1, np
              kplusgi = (xk(j,current_k)+g(j,igk_k(i,current_k)))*tpiba
              hpsi(i,im) = hpsi(i,im) - CMPLX(0._DP,kplusgi,KIND=DP) * psi_g(i,1)
           ENDDO
           !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  !$acc end data
  !$acc exit data delete(psi_g,psic,kedtau)
  DEALLOCATE( psi_g )
  !
  CALL stop_clock( 'h_psi_meta' )
  !
  RETURN
  !
END SUBROUTINE h_psi_meta
