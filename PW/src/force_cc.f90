!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE force_cc( forcecc )
  !----------------------------------------------------------------------
  !! Calculates the NLCC contribution to the force.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi, e2
  USE atom,                 ONLY : rgrid, msh
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, gstart, g, gg, ngl, gl, igtongl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core, tau_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE xc_lib,               ONLY : xclib_dft_is
  USE rhoc_mod,             ONLY : interp_rhc, interp_tac
  !
  IMPLICIT NONE
  !
  REAL(DP) :: forcecc(3,nat)
  !! output: the NLCC forces on atoms
  !
  ! ... local variables
  !
  INTEGER :: ir, nt
  INTEGER :: dfftp_nnr
  REAL(DP), ALLOCATABLE :: vxc(:,:), rhocg(:), kedtaur(:,:)
  COMPLEX(DP), ALLOCATABLE :: vaux(:,:)
  REAL(DP) :: prod, arg, fact
  REAL(DP) :: forcecc_x, forcecc_y, forcecc_z, tau1, tau2, tau3
  REAL(DP) :: etxc_loc, vtxc_loc
  !
  forcecc(:,:) = 0.d0
  !
  IF ( .NOT. ANY(upf(1:ntyp)%nlcc) ) RETURN
  !
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
  !
  dfftp_nnr = dfftp%nnr !to avoid unnecessary copies in acc loop
  !
  ! ... recalculate the exchange-correlation potential
  !
  ALLOCATE( vxc(dfftp%nnr,nspin), vaux(dfftp%nnr,1) )
  !
  IF ( xclib_dft_is('meta') ) THEN
     ALLOCATE( kedtaur(dfftp%nnr,nspin) )
     vxc = 0._DP
     kedtaur = 0._DP
     CALL v_xc_meta( rho, rho_core, rhog_core, tau_core, etxc_loc, vtxc_loc, vxc, kedtaur )
  ELSE
     CALL v_xc( rho, rho_core, rhog_core, etxc_loc, vtxc_loc, vxc )
  END IF
  !
  !$acc data copyin(vxc) create(vaux)
  !
  IF ( nspin==2 ) THEN
     !$acc parallel loop
     DO ir = 1, dfftp_nnr
        vxc(ir,1) = 0.5d0 * ( vxc(ir,1) + vxc(ir,2) )
     ENDDO
  ENDIF
  !
  CALL rho_r2g( dfftp, vxc(:,1:1), vaux(:,1:1) ) 
  !
  ! ... vaux contains now Vxc(G)
  !
  ALLOCATE( rhocg(ngl) )
  !$acc data create(rhocg) present(igtongl)
  !
  ! ... core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  !     g = 0 term gives no contribution
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%nlcc ) THEN
        !
        CALL interp_rhc( nt, ngl, gl, tpiba2, rhocg )
        CALL add_nlcc_force( nt, 1.0_DP )
     ENDIF
  ENDDO
  !
  !$acc end data
  !
  ! ... tau_core displacement term: v_kin(G) * d(tau_core)/d(R) (meta-GGA only)
  !
  IF ( xclib_dft_is('meta') ) THEN
     IF ( nspin==2 ) kedtaur(:,1) = 0.5_DP * ( kedtaur(:,1) + kedtaur(:,2) )
     CALL rho_r2g( dfftp, kedtaur(:,1), vaux(:,1:1) )
     DEALLOCATE( kedtaur )
     !
     ! ... vaux contains now v_kin(G)
     !
     !$acc data create(rhocg) present(igtongl)
     DO nt = 1, ntyp
        IF ( upf(nt)%nlcc .AND. ALLOCATED(upf(nt)%tau_core) ) THEN
           CALL interp_tac( nt, ngl, gl, tpiba2, rhocg )
           CALL add_nlcc_force( nt, e2 )
        ENDIF
     ENDDO
     !$acc end data
  END IF
  !
  CALL mp_sum( forcecc, intra_bgrp_comm )
  !
  !$acc end data
  DEALLOCATE( rhocg )
  DEALLOCATE( vxc, vaux )
  !
  RETURN
  !
CONTAINS

  SUBROUTINE add_nlcc_force( nt_, scale )
    !! Accumulates into forcecc the G-space sum: i*G * exp(-i*R*G) * rhocg(G) * scale * vaux(G)
    INTEGER,  INTENT(IN) :: nt_
    REAL(DP), INTENT(IN) :: scale
    INTEGER  :: na, ig
    REAL(DP) :: tau1, tau2, tau3, arg, prod, fx, fy, fz
    !
#if !defined(_OPENACC)
    !$omp parallel do private( tau1,tau2,tau3,fx,fy,fz,ig,arg,prod )
#endif
    DO na = 1, nat
       IF ( nt_ == ityp(na) ) THEN
          tau1 = tau(1,na)
          tau2 = tau(2,na)
          tau3 = tau(3,na)
          fx = 0.d0 ;  fy = 0.d0 ;  fz = 0.d0
          !$acc parallel loop reduction(+:fx,fy,fz)
          DO ig = gstart, ngm
             arg = (g(1,ig)*tau1 + g(2,ig)*tau2 + g(3,ig)*tau3) * tpi
             prod = tpiba * omega * rhocg(igtongl(ig)) * scale * &
                    DBLE( CONJG(vaux(ig,1)) * &
                    CMPLX(SIN(arg), COS(arg), KIND=DP) ) * fact
             fx = fx + g(1,ig) * prod
             fy = fy + g(2,ig) * prod
             fz = fz + g(3,ig) * prod
          ENDDO
          forcecc(1,na) = forcecc(1,na) + fx
          forcecc(2,na) = forcecc(2,na) + fy
          forcecc(3,na) = forcecc(3,na) + fz
       ENDIF
    ENDDO
#if !defined(_OPENACC)
    !$omp end parallel do
#endif
    !
  END SUBROUTINE add_nlcc_force

END SUBROUTINE force_cc
