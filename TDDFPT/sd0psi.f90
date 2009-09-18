!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sd0psi()
    !
    !  S * d0psi for US case 
    !
    ! OBM :
    ! 050608 Modified for calbec interface in v4.0

#include "f_defs.h"
    !
    use klist,                only : nks,xk
    use lr_variables,         only : n_ipol!, real_space
    use lr_variables,         only : d0psi
    use uspp,                 only : vkb, nkb, okvan
    use wvfct,                only : nbnd, npwx
    use control_flags,        only : gamma_only
    use becmod,               only : becp, rbecp, calbec
    !use real_beta,            only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
    USE realus,              ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                           bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, v_loc_psir, &
                           s_psir_gamma, igk_k, npw_k, real_space_debug
   USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout
                          

    !
    implicit none
    !
    integer :: ik, ip,ibnd
    !
    If (lr_verbosity > 5) THEN
      WRITE(stdout,'("<sd0psi>")')
    endif
    if ( nkb==0 .or. (.not.okvan) ) return
    !
    do ip=1,n_ipol
       !
       if (gamma_only) then
          !
          if (real_space_debug>4) then
           do ibnd=1,nbnd,2
            call fft_orbital_gamma(d0psi(:,:,1,ip),ibnd,nbnd)
            call calbec_rs_gamma(ibnd,nbnd,rbecp)
            call s_psir_gamma(ibnd,nbnd)
            call bfft_orbital_gamma(d0psi(:,:,1,ip),ibnd,nbnd)
           enddo
           ! makedo part until spsi is in place
           !call s_psi(npwx,npw_k(1),nbnd,d0psi(:,:,:,ip),d0psi(:,:,:,ip))
          else
           !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,d0psi(:,:,:,ip),npwx,rbecp,nkb)
           call calbec(npw_k(1),vkb,d0psi(:,:,1,ip),rbecp)
           !notice the third index given as :, whereas in the above routine it is 1. Inquire.
           ! I think it is spin index, spin is not considered yet, leave it for later
           !call s_psi(npwx,npw_k(1),nbnd,d0psi(:,:,:,ip),d0psi(:,:,:,ip))
           call s_psi(npwx,npw_k(1),nbnd,d0psi(:,:,1,ip),d0psi(:,:,1,ip))
          endif
          !
       else
         !
         do ik=1,nks
            !
            call init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
            !call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp,vkb,d0psi(:,:,ik,ip))
            call calbec(npw_k(ik),vkb,d0psi(:,:,ik,ip),becp)
            call s_psi(npwx,npw_k(ik),nbnd,d0psi(:,:,ik,ip),d0psi(:,:,ik,ip))
            !
         enddo
         !
       end if
       !
    end do
    !
  end subroutine sd0psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
