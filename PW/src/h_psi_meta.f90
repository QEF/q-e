!
! Copyright (C) 2007 Quantum ESPRESSO group
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
  USE wvfct,                ONLY : current_k
  USE gvect,                ONLY : g
  USE scf,                  ONLY : kedtau
  USE klist,                ONLY : xk, igk_k
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : psic
  USE fft_base,             ONLY : dffts
  USE fft_wave,             ONLY : wave_g2r
  USE fft_interfaces,       ONLY : fwfft
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
  REAL(DP), ALLOCATABLE :: kplusg(:)
  INTEGER :: im, j, nrxxs
  
  INTEGER :: i, ebnd, brange
  REAL(DP) :: kplusgi
  
  COMPLEX(DP), ALLOCATABLE :: kplusg_evc(:,:)
  COMPLEX(DP), PARAMETER :: ci=(0.d0,1.d0)
  !
  CALL start_clock( 'h_psi_meta' )
  !
  nrxxs = dffts%nnr
  ALLOCATE( kplusg(np) )
  
  ALLOCATE( kplusg_evc(np,2) )
  
  !
  IF (gamma_only) THEN
     !
     ! ... gamma algorithm
     !
     DO im = 1, mp, 2
        DO j = 1, 3
           !
           DO i = 1, np
              kplusgi = (xk(j,current_k)+g(j,i)) * tpiba
              kplusg_evc(i,1) = CMPLX(0.D0,kplusgi) * psip(i,im)
              IF ( im < mp ) kplusg_evc(i,2) = CMPLX(0.d0,kplusgi) * psip(i,im+1)
           ENDDO
           !
           ebnd = im
           IF ( im < mp ) ebnd = ebnd + 1
           brange = ebnd-im+1
           !
           CALL wave_g2r( kplusg_evc(1:np,1:brange), psic, dffts )
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           CALL fwfft( 'Wave', psic, dffts )
           !
           
           kplusg (1:np) = (xk(j,current_k)+g(j,1:np)) * tpiba
           
           IF ( im < mp ) THEN
              hpsi(1:np,im) = hpsi(1:np,im)   - ci * kplusg(1:np) * 0.5d0 * &
                       ( psic(dffts%nl(1:np)) + CONJG(psic(dffts%nlm(1:np))) )
              hpsi(1:np,im+1) = hpsi(1:np,im+1) - kplusg(1:np) * 0.5d0 * &
                       ( psic(dffts%nl(1:np)) - CONJG(psic(dffts%nlm(1:np))) )
           ELSE
              hpsi(1:np,im) = hpsi(1:np,im) - ci * kplusg(1:np) * &
                              psic(dffts%nl(1:np))
           ENDIF
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
           DO i = 1, np
              kplusgi = (xk(j,current_k)+g(j,igk_k(i,current_k))) * tpiba
              kplusg_evc(i,1) = CMPLX(0.D0,kplusgi,kind=DP) * psip(i,im)
           ENDDO
           !
           CALL wave_g2r( kplusg_evc(1:np,1:1), psic, dffts, igk=igk_k(:,current_k) )
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           
           kplusg (1:np) = (xk(j,current_k)+g(j,igk_k(1:np,current_k)))*tpiba
           
           CALL fwfft( 'Wave', psic, dffts )
           !
           hpsi(1:np,im) = hpsi(1:np,im) - CMPLX(0d0, kplusg(1:np), KIND=DP) &
                                         * psic(dffts%nl(igk_k(1:np,current_k)))
        ENDDO
     ENDDO
     !
  ENDIF
  !
  DEALLOCATE( kplusg_evc )
  
  DEALLOCATE( kplusg )
  !
  CALL stop_clock( 'h_psi_meta' )
  !
  RETURN
  !
END SUBROUTINE h_psi_meta
