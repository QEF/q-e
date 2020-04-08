!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE stres_har( sigmahar )
  !--------------------------------------------------------------------
  !! Calculates the Hartree contribution to the stress
  !
  USE kinds,           ONLY: DP
  USE constants,       ONLY: e2, fpi
  USE cell_base,       ONLY: omega, tpiba2
  USE ener,            ONLY: ehart
  USE fft_base,        ONLY: dfftp
  USE fft_interfaces,  ONLY: fwfft
  USE gvect,           ONLY: ngm, gstart, g, gg
  USE scf,             ONLY: rho
  USE control_flags,   ONLY: gamma_only
  USE wavefunctions,   ONLY: psic
  USE mp_bands,        ONLY: intra_bgrp_comm
  USE mp,              ONLY: mp_sum
  USE Coul_cut_2D,     ONLY: do_cutoff_2D, cutoff_stres_sigmahar
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sigmahar(3,3)
  !! Hartree term of the stress tensor
  !
  ! ... local variables
  !
  REAL(DP) :: shart, g2
  REAL(DP), PARAMETER :: eps = 1.d-8
  INTEGER :: ig, l, m
  !
  sigmahar(:,:) = 0.0_DP
  psic(:) = CMPLX( rho%of_r(:,1), KIND=DP )
  !
  CALL fwfft('Rho', psic, dfftp)
  ! psic contains now the charge density in G space
  ! the  G=0 component is not computed
  IF (do_cutoff_2D) THEN  
     CALL cutoff_stres_sigmahar( psic, sigmahar )
  ELSE
     DO ig = gstart, ngm
        g2 = gg(ig) * tpiba2
        shart = psic(dfftp%nl(ig)) * CONJG(psic(dfftp%nl(ig))) / g2
        DO l = 1, 3
           DO m = 1, l
              sigmahar(l,m) = sigmahar(l,m) + shart * tpiba2 * 2 * &
                              g(l,ig) * g(m,ig) / g2
           ENDDO
        ENDDO
     ENDDO
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
END SUBROUTINE stres_har

