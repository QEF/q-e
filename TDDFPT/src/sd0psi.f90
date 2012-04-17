!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE sd0psi()
    !
    !  S * d0psi for US case
    !
    ! Modified by Osman Baris Malcioglu (2009)
    !
    USE klist,                ONLY : nks,xk
    USE lr_variables,         ONLY : n_ipol!, real_space
    USE lr_variables,         ONLY : d0psi
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE wvfct,                ONLY : nbnd, npwx
    USE control_flags,        ONLY : gamma_only
    USE becmod,               ONLY : bec_type, becp, calbec
    !use real_beta,            only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
    USE realus,              ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                           bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, v_loc_psir, &
                           s_psir_gamma, igk_k, npw_k, real_space_debug
   USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout


    !
    IMPLICIT NONE
    !
    INTEGER :: ik, ip,ibnd
    !
    IF (lr_verbosity > 5) THEN
      WRITE(stdout,'("<sd0psi>")')
    ENDIF
    IF ( nkb==0 .or. (.not.okvan) ) RETURN
    !
    DO ip=1,n_ipol
       !
       IF (gamma_only) THEN
          !
          IF (real_space_debug>4) THEN
           DO ibnd=1,nbnd,2
            CALL fft_orbital_gamma(d0psi(:,:,1,ip),ibnd,nbnd)
            CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
            CALL s_psir_gamma(ibnd,nbnd)
            CALL bfft_orbital_gamma(d0psi(:,:,1,ip),ibnd,nbnd)
           ENDDO
           ! makedo part until spsi is in place
           !call s_psi(npwx,npw_k(1),nbnd,d0psi(:,:,:,ip),d0psi(:,:,:,ip))
          ELSE
           !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,d0psi(:,:,:,ip),npwx,rbecp,nkb)
           CALL calbec(npw_k(1),vkb,d0psi(:,:,1,ip),becp)
           !notice the third index given as :, whereas in the above routine it is 1. Inquire.
           ! I think it is spin index, spin is not considered yet, leave it for later
           !call s_psi(npwx,npw_k(1),nbnd,d0psi(:,:,:,ip),d0psi(:,:,:,ip))
           CALL s_psi(npwx,npw_k(1),nbnd,d0psi(:,:,1,ip),d0psi(:,:,1,ip))
          ENDIF
          !
       ELSE
         !
         DO ik=1,nks
            !
            CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
            !call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp,vkb,d0psi(:,:,ik,ip))
            CALL calbec(npw_k(ik),vkb,d0psi(:,:,ik,ip),becp)
            CALL s_psi(npwx,npw_k(ik),nbnd,d0psi(:,:,ik,ip),d0psi(:,:,ik,ip))
            !
         ENDDO
         !
       ENDIF
       !
    ENDDO
    !
  END SUBROUTINE sd0psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
