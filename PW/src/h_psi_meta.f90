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
  USE fft_interfaces,       ONLY : fwfft, invfft
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
  COMPLEX(DP), PARAMETER :: ci=(0.d0,1.d0)
  !
  CALL start_clock( 'h_psi_meta' )
  !
  nrxxs = dffts%nnr
  ALLOCATE( kplusg(np) )
  !
  IF (gamma_only) THEN
     !
     ! gamma algorithm
     !
     DO im = 1, mp, 2
        DO j = 1, 3
           !
           psic(1:nrxxs) = ( 0.D0, 0.D0 )
           !
           kplusg (1:np) = (xk(j,current_k)+g(j,1:np)) * tpiba
           IF (im < mp ) THEN
              psic(dffts%nl (1:np)) =  ci * kplusg(1:np) * &
                             ( psip (1:np,im) + ci * psip(1:np,im+1) )
              psic(dffts%nlm(1:np)) = -ci * kplusg(1:np) * &
                        CONJG( psip (1:np,im) - ci * psip(1:np,im+1) )
           ELSE
              psic(dffts%nl (1:np)) =  ci * kplusg(1:np) *       psip(1:np,im) 
              psic(dffts%nlm(1:np)) = -ci * kplusg(1:np) * CONJG(psip(1:np,im))
           ENDIF
           !
           CALL invfft( 'Wave', psic, dffts )
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           CALL fwfft( 'Wave', psic, dffts )
           !
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
           psic(1:nrxxs) = ( 0.D0, 0.D0 )
           !
           kplusg (1:np) = (xk(j,current_k)+g(j,igk_k(1:np,current_k)))*tpiba
           psic(dffts%nl(igk_k(1:np,current_k))) = CMPLX(0d0, kplusg(1:np), KIND=DP) &
                                                   * psip(1:np,im)
           !
           CALL invfft( 'Wave', psic, dffts )
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           CALL fwfft( 'Wave', psic, dffts )
           !
           hpsi(1:np,im) = hpsi(1:np,im) - CMPLX(0d0, kplusg(1:np), KIND=DP) &
                                         * psic(dffts%nl(igk_k(1:np,current_k)))
        ENDDO
     ENDDO
     !
  ENDIF
  !
  DEALLOCATE( kplusg )
  !
  CALL stop_clock( 'h_psi_meta' )
  !
  RETURN
  !
END SUBROUTINE h_psi_meta
