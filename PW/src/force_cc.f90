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
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid, msh
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, gstart, g, gg, ngl, gl, igtongl
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE rhoc_mod,             ONLY : interp_rhc
  !
  IMPLICIT NONE
  !
  REAL(DP) :: forcecc(3,nat)
  !! output: the NLCC forces on atoms
  !
  ! ... local variables
  !
  INTEGER :: ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms
  INTEGER :: dfftp_nnr
  REAL(DP), ALLOCATABLE :: vxc(:,:), rhocg(:)
  ! exchange-correlation potential
  ! radial fourier transform of rho core
  COMPLEX(DP), ALLOCATABLE :: vaux(:,:)
  REAL(DP) :: prod, arg, fact
  REAL(DP) :: forcecc_x, forcecc_y, forcecc_z, tau1, tau2, tau3
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
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
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
        !
#if !defined(_OPENACC)
        !$omp parallel do private( tau1,tau2,tau3,forcecc_x,forcecc_y,forcecc_z,&
        !$omp                      ig,arg,prod )
#endif
        DO na = 1, nat
          IF (nt == ityp(na) ) THEN
             !
             tau1 = tau(1,na)
             tau2 = tau(2,na)
             tau3 = tau(3,na)
             forcecc_x = 0.d0
             forcecc_y = 0.d0
             forcecc_z = 0.d0
             !
             !$acc parallel loop reduction(+:forcecc_x,forcecc_y,forcecc_z)
             DO ig = gstart, ngm
                arg = (g(1,ig)*tau1 + g(2,ig)*tau2 + g(3,ig)*tau3) * tpi
                prod = tpiba * omega * rhocg(igtongl(ig)) * &
                       DBLE( CONJG(vaux(ig,1)) * &
                       CMPLX(SIN(arg), COS(arg), KIND=DP) ) * fact
                forcecc_x = forcecc_x + g(1,ig) * prod
                forcecc_y = forcecc_y + g(2,ig) * prod
                forcecc_z = forcecc_z + g(3,ig) * prod
             ENDDO
             !
             forcecc(1,na) = forcecc_x
             forcecc(2,na) = forcecc_y
             forcecc(3,na) = forcecc_z
             !
          ENDIF
        ENDDO
#if !defined(_OPENACC)
        !$omp end parallel do
#endif
     ENDIF
  ENDDO
  !
  CALL mp_sum( forcecc, intra_bgrp_comm )
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( rhocg )
  DEALLOCATE( vxc, vaux )
  !
  RETURN
  !
END SUBROUTINE force_cc
