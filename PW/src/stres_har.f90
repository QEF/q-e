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
  USE kinds,              ONLY: DP
  USE constants,          ONLY: e2, fpi
  USE cell_base,          ONLY: omega, tpiba2
  USE ener,               ONLY: ehart
  USE fft_base,           ONLY: dfftp
  USE fft_rho,            ONLY: rho_r2g
  USE gvect,              ONLY: ngm, gstart, g, gg
  USE scf,                ONLY: rho
  USE control_flags,      ONLY: gamma_only
  USE mp_bands,           ONLY: intra_bgrp_comm
  USE mp,                 ONLY: mp_sum
  USE Coul_cut_2D,        ONLY: do_cutoff_2D, cutoff_stres_sigmahar
  !
  IMPLICIT NONE
  !
  REAL(DP) :: sigmahar(3,3)
  !! Hartree term of the stress tensor
  !
  ! ... local variables
  !
  REAL(DP) :: shart, g2
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:)
  REAL(DP), PARAMETER :: eps = 1.E-8_DP
  INTEGER :: ig, l, m
  REAL(DP) :: sigmahar11, sigmahar31, sigmahar21, &
              sigmahar32, sigmahar22, sigmahar33
  !
  ALLOCATE( rhog(dfftp%nnr,1) )
  !$acc data create(rhog)
  !
  CALL rho_r2g( dfftp, rho%of_r(:,1), rhog )
  !
  ! ... the G=0 component is not computed
  !
  sigmahar(:,:) = 0.0_DP
  !
  IF (do_cutoff_2D) THEN
     !
     CALL cutoff_stres_sigmahar( rhog(:,1), sigmahar )
     !
  ELSE
     !
     sigmahar11 = 0._DP  ;  sigmahar31 = 0._DP
     sigmahar21 = 0._DP  ;  sigmahar32 = 0._DP
     sigmahar22 = 0._DP  ;  sigmahar33 = 0._DP
     !
     !$acc parallel loop reduction(+:sigmahar11,sigmahar21,sigmahar22,&
     !$acc&                          sigmahar31,sigmahar32,sigmahar33)
     DO ig = gstart, ngm
       !
       g2 = gg(ig)
       !
       shart = DBLE(rhog(ig,1)*CONJG(rhog(ig,1))) / g2
       !
       sigmahar11 = sigmahar11 + shart *2._DP * &
                                 g(1,ig) * g(1,ig) / g2
       sigmahar21 = sigmahar21 + shart *2._DP * &
                                 g(2,ig) * g(1,ig) / g2
       sigmahar22 = sigmahar22 + shart *2._DP * &
                                 g(2,ig) * g(2,ig) / g2
       sigmahar31 = sigmahar31 + shart *2._DP * &
                                 g(3,ig) * g(1,ig) / g2
       sigmahar32 = sigmahar32 + shart *2._DP * &
                                 g(3,ig) * g(2,ig) / g2
       sigmahar33 = sigmahar33 + shart *2._DP * &
                                 g(3,ig) * g(3,ig) / g2
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
  !$acc end data
  DEALLOCATE( rhog )
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
  RETURN
  !
END SUBROUTINE stres_har
