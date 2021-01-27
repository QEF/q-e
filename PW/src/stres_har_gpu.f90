!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE stres_har_gpu( sigmahar )
  !--------------------------------------------------------------------
  !! Calculates the Hartree contribution to the stress
  !
  USE kinds,              ONLY: DP
  USE constants,          ONLY: e2, fpi
  USE cell_base,          ONLY: omega, tpiba2
  USE ener,               ONLY: ehart
  USE fft_base,           ONLY: dfftp
  USE fft_interfaces,     ONLY: fwfft
  USE gvect,              ONLY: ngm, gstart
  USE scf,                ONLY: rho
  USE control_flags,      ONLY: gamma_only
  USE wavefunctions,      ONLY: psic
  USE mp_bands,           ONLY: intra_bgrp_comm
  USE mp,                 ONLY: mp_sum
  USE Coul_cut_2D,        ONLY: do_cutoff_2D, cutoff_stres_sigmahar_gpu
  !
  USE gvect_gpum,         ONLY: g_d, gg_d
  USE wavefunctions_gpum, ONLY: using_psic, using_psic_d, psic_d
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sigmahar(3,3)
  !! Hartree term of the stress tensor
  !
  ! ... local variables
  !
  REAL(DP) :: shart, g2
  REAL(DP), PARAMETER :: eps = 1.E-8_DP
  INTEGER :: ig, l, m
  !
  INTEGER, POINTER :: nl_d(:)
  REAL(DP) :: sigmahar11, sigmahar31, sigmahar21, &
              sigmahar32, sigmahar22, sigmahar33
  !
#if defined(__CUDA)
  attributes(DEVICE) :: nl_d
#endif
  !
  sigmahar(:,:) = 0.0_DP
  !
  CALL using_psic(2)
  !
  psic(:) = CMPLX( rho%of_r(:,1), KIND=DP )
  !
  CALL using_psic_d(1)
  !
  CALL fwfft( 'Rho', psic_d, dfftp )
  ! psic contains now the charge density in G space
  ! the  G=0 component is not computed
  !
  CALL using_psic_d(0) ; CALL using_psic(0)
  !
  IF (do_cutoff_2D) THEN
     !
     CALL cutoff_stres_sigmahar_gpu( psic_d, sigmahar )
     !
  ELSE
     !
     sigmahar11 = 0._DP  ;  sigmahar31 = 0._DP
     sigmahar21 = 0._DP  ;  sigmahar32 = 0._DP
     sigmahar22 = 0._DP  ;  sigmahar33 = 0._DP
     !
     nl_d => dfftp%nl_d
     !
     !$cuf kernel do (1) <<<*,*>>>
     DO ig = gstart, ngm
       !
       g2 = gg_d(ig)
       !
       shart = DBLE(psic_d(nl_d(ig))*CONJG(psic_d(nl_d(ig)))) / g2
       !
       sigmahar11 = sigmahar11 + shart *2._DP * &
                                 g_d(1,ig) * g_d(1,ig) / g2
       sigmahar21 = sigmahar21 + shart *2._DP * &
                                 g_d(2,ig) * g_d(1,ig) / g2
       sigmahar22 = sigmahar22 + shart *2._DP * &
                                 g_d(2,ig) * g_d(2,ig) / g2
       sigmahar31 = sigmahar31 + shart *2._DP * &
                                 g_d(3,ig) * g_d(1,ig) / g2
       sigmahar32 = sigmahar32 + shart *2._DP * &
                                 g_d(3,ig) * g_d(2,ig) / g2
       sigmahar33 = sigmahar33 + shart *2._DP * &
                                 g_d(3,ig) * g_d(3,ig) / g2
     ENDDO
     !
     sigmahar(1,1) = sigmahar(1,1) + sigmahar11 / tpiba2
     sigmahar(2,1) = sigmahar(2,1) + sigmahar21 / tpiba2
     sigmahar(2,2) = sigmahar(2,2) + sigmahar22 / tpiba2
     sigmahar(3,1) = sigmahar(3,1) + sigmahar31 / tpiba2
     sigmahar(3,2) = sigmahar(3,2) + sigmahar32 / tpiba2
     sigmahar(3,3) = sigmahar(3,3) + sigmahar33 / tpiba2
     !
  ENDIF 
  !
  CALL mp_sum( sigmahar, intra_bgrp_comm )
  !
  IF (gamma_only) THEN
     sigmahar(:,:) = fpi * e2 * sigmahar(:,:)
  ELSE
     sigmahar(:,:) = fpi * e2 * sigmahar(:,:) * 0.5_DP
  ENDIF
  !
  DO l = 1, 3
     sigmahar(l,l) = sigmahar(l,l) - ehart / omega
  ENDDO
  !
  DO l = 1, 3
     DO m = 1, l-1
        sigmahar(m,l) = sigmahar(l,m)
     ENDDO
  ENDDO
  !
  sigmahar(:,:) = -sigmahar(:,:)
  !
  !
  RETURN
  !
END SUBROUTINE stres_har_gpu

