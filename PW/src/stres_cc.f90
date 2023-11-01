!
! Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE stres_cc( sigmaxcc )
  !-----------------------------------------------------------------------
  !! Core correction term of the stress.
  !
  USE kinds,                ONLY : DP
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, gstart, ngl, gl, igtongl, g, gg
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE rhoc_mod,             ONLY : interp_rhc, interp_drhc
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: sigmaxcc(3,3)
  !
  ! ... local variables
  !
  INTEGER :: nt, ng, l, m, ir, dfftp_nnr
  REAL(DP) :: fact
  REAL(DP), ALLOCATABLE :: rhocg(:), vxc(:,:)
  COMPLEX(DP), ALLOCATABLE :: vaux(:,:)
  !
  REAL(DP) :: rhocg1, sigma_rid, sigmadiag
  REAL(DP) :: sigma1, sigma2, sigma3, &
              sigma4, sigma5, sigma6
  !
  sigmaxcc(:,:) = 0._DP
  !
  IF ( .NOT. ANY( upf(1:ntyp)%nlcc ) ) RETURN
  !
  dfftp_nnr = dfftp%nnr !to avoid unnecessary copies in acc loop
  !
  ! ... recalculate the exchange-correlation potential
  !
  ALLOCATE( vxc(dfftp%nnr,nspin), vaux(dfftp%nnr,1) )
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
  !
  !$acc data create( vaux )
  !$acc data copyin( vxc )
  IF ( nspin==2 ) then
     !$acc parallel loop
     DO ir = 1, dfftp_nnr
        vxc(ir,1) = 0.5d0 * ( vxc(ir,1) + vxc(ir,2) )
     ENDDO
  ENDIF
  !
  CALL rho_r2g( dfftp, vxc(:,1), vaux(:,1:1) )
  !
  !$acc end data
  DEALLOCATE( vxc )
  !
  ! ... vaux contains now Vxc(G)
  !
  sigmadiag = 0._DP
  !
  fact = 1._DP
  IF (gamma_only) fact = 2._DP
  !
  !$acc data copyin( gl, strf ) present( igtongl )
  !
  ALLOCATE( rhocg(ngl) )
  !$acc data create( rhocg )
  !
  sigma1 = 0._DP ;  sigma4 = 0._DP
  sigma2 = 0._DP ;  sigma5 = 0._DP
  sigma3 = 0._DP ;  sigma6 = 0._DP
  !
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%nlcc ) THEN
        !
        CALL interp_rhc( nt, ngl, gl, tpiba2, rhocg )
        !
        ! ... diagonal term
        IF (gstart==2) THEN
          !$acc kernels
          rhocg1 = rhocg(igtongl(1))
          sigmadiag = sigmadiag + DBLE(CONJG(vaux(1,1)) * &
                                  strf(1,nt)) * rhocg1
          !$acc end kernels
        ENDIF
        !
        !$acc parallel loop reduction(+:sigmadiag)
        DO ng = gstart, ngm
           sigmadiag = sigmadiag + DBLE(CONJG(vaux(ng,1)) * &
                                   strf(ng,nt)) * rhocg(igtongl(ng)) * fact
        ENDDO
        !
        CALL interp_drhc( nt, ngl, gl, tpiba2, rhocg )
        !
        ! ... non diagonal term (g=0 contribution missing)
        !
        !$acc parallel loop reduction(+:sigma1,sigma2,sigma3,sigma4,sigma5,sigma6)
        DO ng = gstart, ngm
          !
          sigma_rid = DBLE(CONJG(vaux(ng,1)) * strf(ng,nt)) * &
                      rhocg(igtongl(ng)) * tpiba / SQRT(gg(ng)) * fact
          !
          sigma1 = sigma1 + sigma_rid * g(1,ng)*g(1,ng)
          sigma2 = sigma2 + sigma_rid * g(1,ng)*g(2,ng)
          sigma3 = sigma3 + sigma_rid * g(1,ng)*g(3,ng)
          sigma4 = sigma4 + sigma_rid * g(2,ng)*g(2,ng)
          sigma5 = sigma5 + sigma_rid * g(3,ng)*g(2,ng)
          sigma6 = sigma6 + sigma_rid * g(3,ng)*g(3,ng)
          !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  !$acc end data
  !
  sigmaxcc(1,1) = sigma1  ;  sigmaxcc(2,3) = sigma5
  sigmaxcc(1,2) = sigma2  ;  sigmaxcc(3,1) = sigma3
  sigmaxcc(1,3) = sigma3  ;  sigmaxcc(3,2) = sigma5
  sigmaxcc(2,1) = sigma2  ;  sigmaxcc(3,3) = sigma6
  sigmaxcc(2,2) = sigma4
  !
  DO l = 1, 3
     sigmaxcc(l,l) = sigmaxcc(l,l) + sigmadiag
  ENDDO
  !
  CALL mp_sum( sigmaxcc, intra_bgrp_comm )
  !
  !$acc end data
  !$acc end data
  !
  DEALLOCATE( rhocg )
  DEALLOCATE( vaux  )
  !
  RETURN
  !
END SUBROUTINE stres_cc

