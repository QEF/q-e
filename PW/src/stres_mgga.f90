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
  !! Analytic stress tensor contribution from metagga is added to sigmaxc.
  !
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin
  USE cell_base,              ONLY : omega
  USE gvect,                  ONLY : g
  USE scf,                    ONLY : rho, v
  USE wavefunctions,          ONLY : evc
  USE wavefunctions_gpum,     ONLY : using_evc
  USE xc_lib,                 ONLY : xclib_dft_is
  USE klist,                  ONLY : nks, xk, ngk
  USE buffers,                ONLY : get_buffer
  USE io_files,               ONLY : iunwfc, nwordwfc
  USE wvfct,                  ONLY : nbnd, npwx, wg 
  USE lsda_mod,               ONLY : lsda, nspin, current_spin, isk
  USE fft_interfaces,         ONLY : fwfft, invfft
  USE fft_base,               ONLY : dffts
  USE mp,                     ONLY : mp_sum
  USE mp_pools,               ONLY : inter_pool_comm
  USE mp_bands,               ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: sigmaxc(3,3)
  !
  ! ... local variables
  !
  INTEGER :: ix, iy, K, ir, ipol, iss, incr, ibnd, ik, npw
  !
  REAL(DP), PARAMETER :: epsr = 1.E-6_DP, epsg = 1.E-10_DP, e2 = 2._DP
  !
  COMPLEX(DP), ALLOCATABLE :: gradwfc(:,:), crosstaus(:,:,:)
  REAL(DP), ALLOCATABLE :: vkin(:), rhokin(:)
  REAL(DP) :: sigma1, sigma2, sigma3, &
              sigma4, sigma5, sigma6
  !
  REAL(DP) :: w1, w2, delta, sigma_mgga(3,3)
  !
  !
  IF ( .NOT. xclib_dft_is('meta') ) RETURN
  !
  CALL using_evc(1)
  !
  current_spin = 1
  !
  IF ( noncolin ) CALL errore( 'stres_mgga', 'noncollinear stress + meta-GGA &
                               &not implemented', 1 )
  !
  ALLOCATE( gradwfc(dffts%nnr,3) )
  ALLOCATE( crosstaus(dffts%nnr,6,nspin) )
  !$acc data create(crosstaus)
  !$acc data create(gradwfc)
  !
  !$acc kernels
  crosstaus = 0.d0
  !$acc end kernels
  !
  ! ... For gamma_only efficiency
  !
  incr = 1
  IF ( gamma_only ) incr = 2
  !
  ! ... Loop over the k points
  !
  k_loop: DO ik = 1, nks
    !
    IF ( lsda ) current_spin = isk(ik)
    !
    npw = ngk(ik)
    !
    ! ... Read the wavefunctions 
    !
    IF ( nks > 1 ) THEN
      CALL get_buffer( evc, nwordwfc, iunwfc, ik )
      !$acc update device(evc)
    ENDIF  
    !
    !
    DO ibnd = 1, nbnd, incr
       !
       ! ... w1, w2: weights for each k point and band
       !
       w1 = wg(ibnd,ik) / omega 
       !
       IF ( (ibnd < nbnd) .AND. (gamma_only) ) THEN
          !
          ! ... two ffts at the same time
          !
          w2 = wg(ibnd+1,ik) / omega
          !
       ELSE
          !
          w2 = w1
          !
       ENDIF
       !
       ! ... Gradient of the wavefunction in real space
       !
       CALL wfc_gradient( ibnd, ik, npw, gradwfc )
       !
       ! ... Cross terms of kinetic energy density
       !
       !$acc parallel loop collapse(2)
       DO ir = 1, dffts%nnr
         DO ipol = 1, 6
            !
            ! ... explanation here: https://stackoverflow.com/a/244550
            !
            !      M*(M+1)/ 2
            K = (3.0*4.0)/2.0 - 1 - (ipol - 1)
            K = FLOOR((SQRT(FLOAT(8*K+1))-1)/2)
            ix = (ipol-1) - (3.0*4.0)/2.0 + (K+1)*(K+2)/2.0 + 1 + (2-K)
            iy = 3 - K
            !
            crosstaus(ir,ipol,current_spin) = crosstaus(ir,ipol,current_spin) + &
                          2.0_DP*w1*DBLE(gradwfc(ir,ix))*DBLE(gradwfc(ir,iy)) + &
                          2.0_DP*w2*AIMAG(gradwfc(ir,ix))*AIMAG(gradwfc(ir,iy))
         ENDDO
       ENDDO
       !
    ENDDO !ibnd
    !
  ENDDO k_loop 
  !
  !$acc end data
  DEALLOCATE( gradwfc )
  !
  !$acc host_data use_device(crosstaus)
  CALL mp_sum( crosstaus, inter_pool_comm )
  !$acc end host_data
  !
  ALLOCATE( vkin(dffts%nnr) )
  ALLOCATE( rhokin(dffts%nnr) )
  !$acc data create(vkin,rhokin)
  !
  ! ... metagga contribution to the stress tensor
  sigma_mgga(:,:) = 0._DP
  !
  sigma1 = 0.d0 ;  sigma4 = 0.d0
  sigma2 = 0.d0 ;  sigma5 = 0.d0
  sigma3 = 0.d0 ;  sigma6 = 0.d0
  !
  DO iss = 1, nspin 
     !
     vkin = v%kin_r(:,iss)
     rhokin = rho%kin_r(:,iss)
     !$acc update device(vkin,rhokin)
     !
     !$acc parallel loop reduction(+:sigma1,sigma2,sigma3,sigma4,sigma5,sigma6)
     DO ir = 1, dffts%nnr
        !
        sigma1 = sigma1 + vkin(ir) * ( rhokin(ir) &
                        + DBLE(crosstaus(ir,1,iss)) )
        sigma2 = sigma2 + vkin(ir) * DBLE(crosstaus(ir,2,iss))
        sigma3 = sigma3 + vkin(ir) * DBLE(crosstaus(ir,3,iss))
        sigma4 = sigma4 + vkin(ir) * ( rhokin(ir) &
                        + DBLE(crosstaus(ir,4,iss)) )
        sigma5 = sigma5 + vkin(ir) * DBLE(crosstaus(ir,5,iss))
        sigma6 = sigma6 + vkin(ir) * ( rhokin(ir) &
                        + DBLE(crosstaus(ir,6,iss)) )
        !
     ENDDO
     !
  ENDDO
  !
  sigma_mgga(1,1) = sigma1  ;  sigma_mgga(2,3) = sigma5
  sigma_mgga(1,2) = sigma2  ;  sigma_mgga(3,1) = sigma3
  sigma_mgga(1,3) = sigma3  ;  sigma_mgga(3,2) = sigma5
  sigma_mgga(2,1) = sigma2  ;  sigma_mgga(3,3) = sigma6
  sigma_mgga(2,2) = sigma4
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( vkin, rhokin )
  DEALLOCATE( crosstaus )
  !
  CALL mp_sum( sigma_mgga, intra_bgrp_comm )
  !
  sigmaxc(:,:) = sigmaxc(:,:) + sigma_mgga(:,:) / &
                 (dffts%nr1 * dffts%nr2 * dffts%nr3)
  !
  RETURN
  !
END SUBROUTINE stres_mgga
!
!
!----------------------------------------------------------
SUBROUTINE wfc_gradient( ibnd, ik, npw, gradpsi )
  !----------------------------------------------------------
  !! Returns the gradient of the wavefunction in real space.
  !! Slightly adapted from sum_bands.f90 
  !
  USE kinds,                  ONLY: DP
  USE control_flags,          ONLY: gamma_only
  USE wvfct,                  ONLY: npwx, nbnd
  USE cell_base,              ONLY: omega, tpiba
  USE klist,                  ONLY: xk, igk_k
  USE wavefunctions,          ONLY: psic, evc
  USE fft_base,               ONLY: dffts
  USE fft_interfaces,         ONLY: invfft
  USE gvect,                  ONLY: ngm, g
  !
  IMPLICIT NONE 
  !
  INTEGER, INTENT(IN) :: ibnd, ik, npw
  COMPLEX(DP), INTENT(OUT) :: gradpsi(dffts%nnr,3)
  !
  ! ... local variables
  !
  INTEGER  :: ipol, j
  REAL(DP) :: kplusg
  INTEGER, ALLOCATABLE :: nld(:), nlmd(:)
  REAL(DP) :: xki(3)
  !
  !$acc data present(gradpsi) copyin(evc) create(psic)
  !
  ALLOCATE( nld(npw) )
  nld = dffts%nl
  !
  xki(1:3) = xk(1:3,ik)
  !
  ! ... Compute the gradient of the wavefunction in reciprocal space
  !
  IF ( gamma_only ) THEN
     !
     ALLOCATE( nlmd(npw) )
     nlmd = dffts%nlm
     !$acc data copyin(xki,nld,nlmd)
     !
     DO ipol = 1, 3
        !
        !$acc kernels
        psic(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        IF ( ibnd < nbnd ) THEN
           !
           ! ... two ffts at the same time
           !
           !$acc parallel loop
           DO j = 1, npw
              kplusg = (xki(ipol)+g(ipol,igk_k(j,ik))) * tpiba
              !
              psic(nld(j))  = CMPLX(0._DP, kplusg, kind=DP) * &
                              ( evc(j,ibnd) + (0._DP,1._DP) * evc(j,ibnd+1) )
              !
              psic(nlmd(j)) = CMPLX(0._DP,-kplusg, kind=DP) * &
                              CONJG( evc(j,ibnd) - (0._DP,1._DP) * evc(j,ibnd+1) )
           ENDDO
           !
        ELSE
           !
           !$acc parallel loop
           DO j = 1, npw
              kplusg = (xki(ipol)+g(ipol,igk_k(j,ik))) * tpiba
              !
              psic(nld(j))  = CMPLX(0._DP, kplusg, kind=DP) * evc(j,ibnd)
              !
              psic(nlmd(j)) = CMPLX(0._DP,-kplusg,kind=DP) * CONJG(evc(j,ibnd))
           ENDDO
           !
        ENDIF
        !
        ! ... Gradient of the wavefunction in real space
        !
        !$acc host_data use_device(psic)
        CALL invfft( 'Wave', psic, dffts )
        !$acc end host_data
        !
        !$acc kernels
        gradpsi(:,ipol) = psic
        !$acc end kernels
        !
     ENDDO
     !$acc end data
     DEALLOCATE( nlmd )
     !
  ELSE
     !
     !$acc data copyin(xki,nld)
     DO ipol = 1, 3
        !
        !$acc kernels
        psic(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        !$acc parallel loop
        DO j = 1, npw
           kplusg = (xki(ipol)+g(ipol,igk_k(j,ik))) * tpiba
           !
           psic(nld(j)) = CMPLX(0._DP,kplusg,kind=DP)*evc(j,ibnd)
        ENDDO
        !
        ! ... Gradient of the wavefunction in real space
        !
        !$acc host_data use_device(psic)
        CALL invfft( 'Wave', psic, dffts )
        !$acc end host_data
        !
        !$acc kernels
        gradpsi(:,ipol) = psic
        !$acc end kernels
        !
     ENDDO
     !$acc end data
  ENDIF 
  !
  DEALLOCATE( nld )
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE wfc_gradient
