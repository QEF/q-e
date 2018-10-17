!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE s_1psi_gpu( npwx, n, psi_d, spsi_d )
  !----------------------------------------------------------------------------
  !
  ! ... spsi = S*psi for one wavefunction
  ! ... Wrapper routine - calls calbec and s_psi
  !
  USE kinds,  ONLY : DP
  USE uspp,   ONLY : nkb
  USE becmod, ONLY : bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only 
  USE noncollin_module, ONLY : noncolin, npol 
  USE realus, ONLY : real_space, &
              invfft_orbital_gamma, fwfft_orbital_gamma, calbec_rs_gamma, s_psir_gamma, &
              invfft_orbital_k, fwfft_orbital_k, calbec_rs_k, s_psir_k
  USE wvfct,  ONLY: nbnd
  !
  USE uspp_gpum,        ONLY : vkb_d, using_vkb_d
  USE becmod_gpum,      ONLY : becp_d, using_becp_r
  USE becmod_subs_gpum, ONLY : using_becp_d_auto, calbec_gpu
  !
  IMPLICIT NONE
  !
  INTEGER     :: npwx, n, ibnd
  COMPLEX(DP) :: psi_d(npwx*npol,1), spsi_d(npwx*npol,1)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, spsi_d
#endif
  COMPLEX(DP), ALLOCATABLE :: psi_h(:,:), spsi_h(:,:)
  !
  CALL start_clock( 's_1psi' )
  !
  IF ( real_space) then
     ALLOCATE(psi_h(npwx*npol,1), spsi_h(npwx*npol,1))
     psi_h(1:npwx*npol,1) = psi_d(1:npwx*npol,1)
     spsi_h(1:npwx*npol,1) = spsi_d(1:npwx*npol,1)
     IF ( gamma_only ) then
        do ibnd=1,nbnd,2
           ! transform the orbital to real space
           call invfft_orbital_gamma(psi_h,ibnd,nbnd) 
           ! global becp%r is updated
           CALL using_becp_r(2)
           call calbec_rs_gamma(ibnd,nbnd,becp%r) 
        enddo
        call s_psir_gamma(1,1)
        call fwfft_orbital_gamma(spsi_h,1,1)
        !
     ELSE
        do ibnd=1,nbnd
           ! transform the orbital to real space
           call invfft_orbital_k(psi_h,ibnd,nbnd) 
           ! global becp%r is updated
           call calbec_rs_k(ibnd,nbnd) 
        enddo
        call s_psir_k(1,1)
        call fwfft_orbital_k(spsi_h,1,1)
     END IF
     spsi_d(1:npwx*npol,1) = spsi_h(1:npwx*npol,1)
     DEALLOCATE(psi_h, spsi_h)
  ELSE
     !
     CALL using_vkb_d(0); CALL using_becp_d_auto(1)
     CALL calbec_gpu( n, vkb_d, psi_d, becp_d )
     CALL s_psi_gpu( npwx, n, 1, psi_d, spsi_d )
     !
  END IF
  !
  CALL stop_clock( 's_1psi' )
  !
  RETURN
  !
END SUBROUTINE s_1psi_gpu
