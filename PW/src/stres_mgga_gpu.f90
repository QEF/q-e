!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_mgga_gpu( sigmaxc )
  !----------------------------------------------------------------------------
  !! Analytic stress tensor contribution from metagga is added to sigmaxc
  !
  USE kinds,                  ONLY : DP
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin
  USE cell_base,              ONLY : alat, at, bg, omega, tpiba
  USE gvect,                  ONLY : g
  USE scf,                    ONLY : rho, v
  USE wavefunctions,          ONLY : evc
  USE xc_lib,                 ONLY : xclib_dft_is
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
  USE wavefunctions_gpum,     ONLY : using_evc
  !
  USE device_fbuff_m ,              ONLY : dev_buf
  USE device_memcpy_m,          ONLY : dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: sigmaxc(3,3)
  !
  ! Internal variables
  !
  INTEGER :: ix, iy, K, ir, ipol, iss, incr, ibnd, ik, npw
  !
  INTEGER :: ierrs(4)
  !
  REAL(DP), PARAMETER :: epsr = 1.E-6_DP, epsg = 1.E-10_DP, e2 = 2._DP
  !
  COMPLEX(DP), POINTER :: gradwfc_d(:,:), crosstaus_d(:,:,:)
  INTEGER :: ix_d(6), iy_d(6)
  REAL(DP), POINTER :: vkin_d(:), rhokin_d(:)
  REAL(DP) :: sigma1_d, sigma2_d, sigma3_d, &
              sigma4_d, sigma5_d, sigma6_d
  !
  REAL(DP) :: w1, w2, delta, sigma_mgga(3,3)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: gradwfc_d, crosstaus_d, vkin_d, rhokin_d, ix_d, iy_d
#endif
  !
  if ( .not. xclib_dft_is('meta') ) return
  !
  current_spin = 1
  !
  ! Stop if something is not yet implemented
  !
  IF ( noncolin ) CALL errore( 'stres_mgga', &
                    'noncollinear stress + meta-GGA not implemented', 1 )
  !
  ! Initialization of a set of variables
  !
  CALL dev_buf%lock_buffer( gradwfc_d, (/ dffts%nnr, 3 /), ierrs(1) )
  CALL dev_buf%lock_buffer( crosstaus_d, (/ dffts%nnr, 6, nspin /), ierrs(2) )
  IF (ANY(ierrs(1:2) /= 0)) CALL errore( 'stres_mgga_gpu', 'cannot allocate buffers', -1 )
  !
  CALL dev_memset(crosstaus_d , (0._DP,0._DP) )
  !
  ! For gamma_only efficiency
  !
  incr = 1
  IF ( gamma_only ) incr = 2
  !
  !Polarization indexes
  ! ix_d(1) = 1  ;  iy_d(1) = 1
  ! ix_d(2) = 2  ;  iy_d(2) = 1
  ! ix_d(3) = 3  ;  iy_d(3) = 1
  ! ix_d(4) = 2  ;  iy_d(4) = 2
  ! ix_d(5) = 3  ;  iy_d(5) = 2
  ! ix_d(6) = 3  ;  iy_d(6) = 3
  !
  ! Loop over the k points
  !
  k_loop: DO ik = 1, nks
    !
    IF ( lsda ) current_spin = isk(ik)
    !
    npw = ngk(ik)
    !
    ! Read the wavefunctions 
    !
    IF ( nks > 1 ) THEN
       CALL get_buffer( evc, nwordwfc, iunwfc, ik )
       CALL using_evc(2)
    ENDIF
    !
    DO ibnd = 1, nbnd, incr
       !
       ! w1, w2: weights for each k point and band
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
       ! Gradient of the wavefunction in real space
       ! 
       CALL wfc_gradient_gpu( ibnd, ik, npw, gradwfc_d )
       !
       ! Cross terms of kinetic energy density
       !
       ! FIX ME: PGI complains here if I set do(2)
       !$cuf kernel do (1) <<<*,*>>>
       DO ir = 1, dffts%nnr
         DO ipol = 1, 6
            !
            ! explenation here: https://stackoverflow.com/a/244550
            !
            !      M*(M+1)/ 2
            K = (3.0*4.0)/2.0 - 1 - (ipol - 1)
            K = floor((SQRT(FLOAT(8*K+1))-1)/2)
            ix = (ipol-1) - (3.0*4.0)/2.0 + (K+1)*(K+2)/2.0 + 1 + (2-K)
            iy = 3 - K
            !
            crosstaus_d(ir,ipol,current_spin) = crosstaus_d(ir,ipol,current_spin) + &
                       2.0_DP*w1*DBLE(gradwfc_d(ir,ix))*DBLE(gradwfc_d(ir,iy)) + &
                       2.0_DP*w2*AIMAG(gradwfc_d(ir,ix))*AIMAG(gradwfc_d(ir,iy))
         ENDDO
       ENDDO
       !
       !crosstaus(:,:,current_spin) = crosstaus(:,:,current_spin) + &
       !                              crosstaus_d(:,:,current_spin)
       !
    ENDDO !ibnd
    !
  ENDDO k_loop 
  !
  !
  CALL dev_buf%release_buffer( gradwfc_d, ierrs(1) )
  !
  CALL mp_sum( crosstaus_d, inter_pool_comm )
  !
  CALL dev_buf%lock_buffer( vkin_d, dffts%nnr, ierrs(3) )
  CALL dev_buf%lock_buffer( rhokin_d, dffts%nnr, ierrs(4) )
  IF (ANY(ierrs(3:4) /= 0)) CALL errore( 'stres_mgga_gpu', 'cannot allocate buffers', -1 )
  !
  ! metagga contribution to the stress tensor
  sigma_mgga(:,:) = 0._DP
  !
  sigma1_d = 0.d0 ;  sigma4_d = 0.d0
  sigma2_d = 0.d0 ;  sigma5_d = 0.d0
  sigma3_d = 0.d0 ;  sigma6_d = 0.d0
  !
  !
  DO iss = 1, nspin 
     !
     CALL dev_memcpy( vkin_d, v%kin_r(:,iss) )
     CALL dev_memcpy( rhokin_d, rho%kin_r(:,iss) )
     !
     !$cuf kernel do (1) <<<*,*>>>
     DO ir = 1, dffts%nnr
         !
         sigma1_d = sigma1_d + vkin_d(ir) * (rhokin_d(ir) &
                             + DBLE(crosstaus_d(ir,1,iss)) )
         sigma2_d = sigma2_d + vkin_d(ir) * DBLE(crosstaus_d(ir,2,iss))
         sigma3_d = sigma3_d + vkin_d(ir) * DBLE(crosstaus_d(ir,3,iss))
         sigma4_d = sigma4_d + vkin_d(ir) * (rhokin_d(ir) &
                             + DBLE(crosstaus_d(ir,4,iss)) )
         sigma5_d = sigma5_d + vkin_d(ir) * DBLE(crosstaus_d(ir,5,iss))
         sigma6_d = sigma6_d + vkin_d(ir) * (rhokin_d(ir) &
                             + DBLE(crosstaus_d(ir,6,iss)) )
         !
     ENDDO
     !
  ENDDO
  !
  sigma_mgga(1,1) = sigma1_d  ;  sigma_mgga(2,3) = sigma5_d
  sigma_mgga(1,2) = sigma2_d  ;  sigma_mgga(3,1) = sigma3_d
  sigma_mgga(1,3) = sigma3_d  ;  sigma_mgga(3,2) = sigma5_d
  sigma_mgga(2,1) = sigma2_d  ;  sigma_mgga(3,3) = sigma6_d
  sigma_mgga(2,2) = sigma4_d
  !
  CALL dev_buf%release_buffer( vkin_d, ierrs(3) )
  CALL dev_buf%release_buffer( rhokin_d, ierrs(4) )
  CALL dev_buf%release_buffer( crosstaus_d, ierrs(2) )
  !
  CALL mp_sum( sigma_mgga, intra_bgrp_comm )
  !
  sigmaxc(:,:) = sigmaxc(:,:) + sigma_mgga(:,:) / &
                 (dffts%nr1 * dffts%nr2 * dffts%nr3)
  !
  RETURN
  !
END SUBROUTINE stres_mgga_gpu
!
!
!----------------------------------------------------------
SUBROUTINE wfc_gradient_gpu( ibnd, ik, npw, gradpsi_d )
  !----------------------------------------------------------
  !! Returns the gradient of the wavefunction in real space.
  !! Slightly adapted from sum_bands.f90 
  !
  USE kinds,                  ONLY: DP
  USE control_flags,          ONLY: gamma_only
  USE wvfct,                  ONLY: npwx, nbnd
  USE cell_base,              ONLY: omega, tpiba
  USE klist,                  ONLY: xk, igk_k
  
  USE fft_base,               ONLY: dffts
  USE fft_interfaces,         ONLY: invfft
  !
  USE gvect_gpum,             ONLY: g_d
  USE wavefunctions_gpum,     ONLY: using_evc, using_evc_d, evc_d, &
                                    using_psic, using_psic_d, psic_d
  USE device_fbuff_m,               ONLY: dev_buf
  USE device_memcpy_m,          ONLY: dev_memcpy
  !
  IMPLICIT NONE 
  !
  INTEGER, INTENT(IN) :: ibnd, ik, npw
  COMPLEX(DP), INTENT(OUT) :: gradpsi_d(dffts%nnr,3)
  !
  ! ... local variables
  !
  INTEGER  :: ipol, j
  REAL(DP) :: kplusg
  INTEGER, POINTER :: nl_d(:), nlm_d(:)
  INTEGER, ALLOCATABLE :: igk_k_d(:)
  REAL(DP) :: xk_d(3)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: gradpsi_d, nl_d, nlm_d, xk_d, igk_k_d
#endif  
  !
  CALL using_evc_d(0)
  CALL using_psic_d(2)
  !
  nl_d    => dffts%nl_d
  nlm_d   => dffts%nlm_d
  !
  ALLOCATE( igk_k_d(npwx) )
  igk_k_d = igk_k(:,ik)
  !
  xk_d(1:3) = xk(1:3,ik)
  !
  ! Compute the gradient of the wavefunction in reciprocal space
  !
  IF ( gamma_only ) THEN
     !
     DO ipol = 1, 3
        !
        psic_d(:) = (0._DP,0._DP)
        !
        IF ( ibnd < nbnd ) THEN
           !
           ! ... two ffts at the same time
           !
           !$cuf kernel do (1) <<<*,*>>>
           DO j = 1, npw
              kplusg = (xk_d(ipol)+g_d(ipol,igk_k_d(j))) * tpiba
              !
              psic_d(nl_d(j))  = CMPLX(0._DP, kplusg, kind=DP) * &
                                   ( evc_d(j,ibnd) + &
                                   ( 0._DP, 1._DP ) * evc_d(j,ibnd+1) )
              !
              psic_d(nlm_d(j)) = CMPLX(0._DP,-kplusg, kind=DP) * &
                                    CONJG( evc_d(j,ibnd) - &
                                    ( 0._DP, 1._DP ) * evc_d(j,ibnd+1) )
           ENDDO
           !
        ELSE
           !
           !$cuf kernel do (1) <<<*,*>>>
           DO j = 1, npw
              kplusg = (xk_d(ipol)+g_d(ipol,igk_k_d(j))) * tpiba
              !
              psic_d(nl_d(j))  = CMPLX(0._DP, kplusg, kind=DP) * &
                                       evc_d(j,ibnd)
              !
              psic_d(nlm_d(j)) = CMPLX(0._DP,-kplusg,kind=DP) * &
                                       CONJG( evc_d(j,ibnd) )
           ENDDO
           !
        ENDIF
        !
        ! Gradient of the wavefunction in real space
        !
        CALL invfft( 'Wave', psic_d, dffts )
        !
        CALL dev_memcpy( gradpsi_d(:,ipol), psic_d )
        !
     ENDDO
     !
  ELSE
     !     
     DO ipol = 1, 3
         !
         psic_d(:) = (0._DP,0._DP)
         !
         !$cuf kernel do(1) <<<*,*>>>
         DO j = 1, npw
            kplusg = (xk_d(ipol)+g_d(ipol,igk_k_d(j))) * tpiba
            !
            psic_d(nl_d(j)) = CMPLX(0._DP,kplusg,kind=DP)*evc_d(j,ibnd)
         ENDDO
         !
         ! Gradient of the wavefunction in real space
         !
         CALL invfft( 'Wave', psic_d, dffts )
         !
         CALL dev_memcpy( gradpsi_d(:,ipol), psic_d )
         !
     ENDDO 
     !
  ENDIF
  !
  DEALLOCATE( igk_k_d )
  !
END SUBROUTINE wfc_gradient_gpu
