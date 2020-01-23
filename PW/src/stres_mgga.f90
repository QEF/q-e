!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_mgga( sigmaxc )
  !----------------------------------------------------------------------------
  !
  ! Analytic stress tensor contribution from metagga is added to sigmaxc
  !
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin
  USE cell_base,              ONLY : alat, at, bg, omega, tpiba
  USE gvect,                  ONLY : g
  USE scf,                    ONLY : rho, v
  USE wavefunctions,          ONLY : evc
  USE funct,                  ONLY : dft_is_meta
  USE klist,                  ONLY : nks, xk, ngk
  USE buffers,                ONLY : get_buffer
  USE io_files,               ONLY : iunwfc, nwordwfc
  USE wvfct,                  ONLY : nbnd, npwx, wg 
  USE lsda_mod,               ONLY : lsda, nspin, current_spin, isk
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE fft_base,               ONLY : dfftp, dffts
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : inter_pool_comm
  USE mp_bands,               ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT)  :: sigmaxc(3,3)
  !
  ! Internal variables
  !
  INTEGER                   :: ix, iy, ir, iss, ipol, incr, ibnd, ik, npw
  INTEGER                   :: ipol2xy(3,3) 
  !! ipol2xy(i,j) = ipol2x(j,i) is a collapsed symmetric index
  DATA ipol2xy / 1, 2, 3, 2, 4, 5, 3, 5, 6/
  REAL(DP), PARAMETER       :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  COMPLEX(DP), ALLOCATABLE  :: gradwfc (:,:), crosstaus(:,:,:)
  REAL(DP)                  :: w1, w2, delta, sigma_mgga(3,3)
  !
  if ( .not. dft_is_meta() ) return
  !
  current_spin=1
  !
  !
  ! Stop if something is not yet implemented
  !
  if (noncolin) call errore('stres_mgga', &
                    'noncollinear stress + meta-GGA not implemented',1)
  !
  ! Initialization of a set of variables
  !
  allocate (gradwfc( dffts%nnr, 3))    
  allocate (crosstaus( dffts%nnr,6,nspin))    
  !
  ! For gamma_only efficiency
  !
  incr=1
  IF ( gamma_only ) incr=2 
  !
  crosstaus(:,:,:) = 0.d0
  gradwfc(:,:) = 0.d0
  !
  ! Loop over the k points
  !
  k_loop: DO ik = 1, nks
  !
    !
    IF ( lsda ) current_spin = isk(ik)
    !
    npw = ngk(ik)
    !
    ! Read the wavefunctions 
    !
    IF ( nks > 1 ) THEN
       !
       CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
       !
    END IF
    !
    do ibnd = 1, nbnd, incr 
       !
       ! w1, w2: weights for each k point and band
       !
       w1 = wg(ibnd,ik) / omega 
       !
       IF ( (ibnd < nbnd) .and. (gamma_only) ) THEN
          !
          ! ... two ffts at the same time
          !
          w2 = wg(ibnd+1,ik) / omega
          !
       ELSE
          !
          w2 = w1
          !
       END IF
       !
       ! Gradient of the wavefunction in real space
       ! 
       CALL wfc_gradient( ibnd, ik, npw, gradwfc )
       !
       ! Cross terms of kinetic energy density
       !
       do ix=1,3
          !
          do iy=1,ix
             !
             ipol = ipol2xy(iy,ix)
             !
             do ir=1,dffts%nnr
                !
                crosstaus(ir,ipol,current_spin) = crosstaus(ir,ipol,current_spin) +&
                                        2.0_DP*w1*DBLE(gradwfc(ir,ix))*DBLE(gradwfc(ir,iy)) +&
                                        2.0_DP*w2*AIMAG(gradwfc(ir,ix))*AIMAG(gradwfc(ir,iy))
                !
             end do
             !
          end do
          !
       end do
       !
    end do
    !
  END DO k_loop 
  !
  call mp_sum(  crosstaus, inter_pool_comm )
  !
  ! gradwfc not used anymore
  !
  deallocate (gradwfc)    
  !
  sigma_mgga(:,:) = 0.D0
  !
  ! metagga contribution to the stress tensor
  !
  do iss=1,nspin
     !
     do ix=1,3
        !
        do iy=1,3
           !
           delta=0.
           if (ix==iy) delta=1.
           !
           do ir=1,dffts%nnr
              !
              sigma_mgga(ix,iy) = sigma_mgga(ix,iy) + v%kin_r(ir,iss) &
                                * ( rho%kin_r(ir,iss) * delta &
                                + crosstaus(ir,ipol2xy(ix,iy),iss) )
              !
           end do
           !
           !
        end do
        !
     end do
     !
  end do
  deallocate( crosstaus )
  !
  call mp_sum(  sigma_mgga, intra_bgrp_comm )
  sigmaxc(:,:) = sigmaxc(:,:) + sigma_mgga(:,:) / &
                 (dffts%nr1 * dffts%nr2 * dffts%nr3)
  !
 return
 !  
END SUBROUTINE stres_mgga

SUBROUTINE wfc_gradient ( ibnd, ik, npw, gradpsi )
  !
  ! Returns the gradient of the wavefunction in real space
  ! Slightly adapted from sum_bands.f90 
  !
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : gamma_only
  USE wavefunctions,          ONLY : psic, evc
  USE wvfct,                  ONLY : npwx, nbnd
  USE cell_base,              ONLY : omega, tpiba
  USE klist,                  ONLY : xk, igk_k
  USE gvect,                  ONLY : g
  USE fft_base,               ONLY : dffts
  USE fft_interfaces,         ONLY : invfft
  !
  IMPLICIT NONE 
  !
  INTEGER                   :: ibnd, ik, npw 
  COMPLEX(DP)               :: gradpsi(dffts%nnr,3)  
  !
  ! Internal variables
  !
  REAL(DP)               :: kplusg(npwx)
  INTEGER                :: ipol
  !
  ! Compute the gradient of the wavefunction in reciprocal space
  !
  IF ( gamma_only ) THEN
     !
     DO ipol=1,3
        !
        psic(:) = ( 0.D0, 0.D0 )
        !
        kplusg (1:npw) = (xk(ipol,ik)+g(ipol,igk_k(1:npw,ik))) * tpiba
        !
        IF ( ibnd < nbnd ) THEN
           !
           ! ... two ffts at the same time
           !
           psic(dffts%nl(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP)* &
                                   ( evc(1:npw,ibnd) + &
                                   ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
           !
           psic(dffts%nlm(1:npw)) = CMPLX(0d0,-kplusg(1:npw),kind=DP) * &
                                    CONJG( evc(1:npw,ibnd) - &
                                    ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
           !
        ELSE
           !
           psic(dffts%nl(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP)* &
                                   evc(1:npw,ibnd)
           !
           psic(dffts%nlm(1:npw)) = CMPLX(0d0,-kplusg(1:npw),kind=DP) * &
                                    CONJG( evc(1:npw,ibnd) )
           !
        END IF
        !
        ! Gradient of the wavefunction in real space
        !
        CALL invfft ('Wave', psic, dffts)
        !
        gradpsi(:,ipol) = psic
        !
     END DO
     !
  ELSE
     !
     DO ipol=1,3
         !
         psic(:) = ( 0.D0, 0.D0 )
         !
         kplusg (1:npw) = (xk(ipol,ik)+g(ipol,igk_k(1:npw,ik))) * tpiba
         psic(dffts%nl(igk_k(1:npw,ik))) = CMPLX(0d0, kplusg(1:npw),kind=DP)* &
                                 evc(1:npw,ibnd)
         !
         ! Gradient of the wavefunction in real space
         !
         CALL invfft ('Wave', psic, dffts)
         !
         gradpsi(:,ipol) = psic
         !
     END DO 
     !
  END IF
  !
END SUBROUTINE wfc_gradient
