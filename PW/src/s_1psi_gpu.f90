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
  !! spsi = S*psi for one wavefunction. Wrapper routine - calls \(texttt{calbec}
  !! and \texttt{s_psi}.
  !
  USE kinds,              ONLY: DP
  USE uspp,               ONLY: nkb, vkb
  USE becmod,             ONLY: bec_type, becp, calbec
  USE control_flags,      ONLY: gamma_only 
  USE noncollin_module,   ONLY: noncolin, npol 
  USE realus,             ONLY: real_space, invfft_orbital_gamma,     &
                                fwfft_orbital_gamma, calbec_rs_gamma, &
                                s_psir_gamma, invfft_orbital_k,       &
                                fwfft_orbital_k, calbec_rs_k, s_psir_k
  USE wvfct,              ONLY: nbnd
  USE becmod_gpum,        ONLY: becp_d, using_becp_r
  USE becmod_subs_gpum,   ONLY: using_becp_d_auto, calbec_gpu
  IMPLICIT NONE
  !
  INTEGER :: npwx
  !! maximum number of PW for wavefunctions
  INTEGER :: n
  !! the number of plane waves
  COMPLEX(DP) :: psi_d(npwx*npol,1)
  !! input vector
  COMPLEX(DP) :: spsi_d(npwx*npol,1)
  !! S*psi
  !
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, spsi_d
#endif
  COMPLEX(DP), ALLOCATABLE :: psi_h(:,:), spsi_h(:,:)
  ! ... local variables
  !
  INTEGER :: ibnd
  !
  !
  CALL start_clock_gpu( 's_1psi' )
  !
  IF ( real_space) THEN
     !
     ALLOCATE(psi_h(npwx*npol,1), spsi_h(npwx*npol,1))
     psi_h(1:npwx*npol,1) = psi_d(1:npwx*npol,1)
     spsi_h(1:npwx*npol,1) = spsi_d(1:npwx*npol,1)
     IF ( gamma_only ) THEN
        !
        DO ibnd = 1, nbnd, 2
           ! transform the orbital to real space
           CALL invfft_orbital_gamma(psi_h,ibnd,nbnd) 
           ! global becp%r is updated
           CALL using_becp_r(2)
           CALL calbec_rs_gamma(ibnd,nbnd,becp%r) 
        ENDDO
        !
        CALL s_psir_gamma(1,1)
        CALL fwfft_orbital_gamma(spsi_h,1,1)
        !
     ELSE
        !
        DO ibnd = 1, nbnd
           ! transform the orbital to real space
           CALL invfft_orbital_k(psi_h,ibnd,nbnd) 
           ! global becp%r is updated
           CALL calbec_rs_k( ibnd, nbnd )
        ENDDO
        !
        CALL s_psir_k( 1, 1 )
        CALL fwfft_orbital_k( spsi_h, 1, 1 )
        !
     ENDIF
     !
     spsi_d(1:npwx*npol,1) = spsi_h(1:npwx*npol,1)
     DEALLOCATE(psi_h, spsi_h)
  ELSE
     !
     CALL using_becp_d_auto(1)
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
     CALL calbec_gpu( n, vkb, psi_d, becp_d )
!$acc end host_data 
!$acc end data 
     CALL s_psi_gpu( npwx, n, 1, psi_d, spsi_d )
     !
  ENDIF
  !
  CALL stop_clock_gpu( 's_1psi' )
  !
  RETURN
  !
END SUBROUTINE s_1psi_gpu
