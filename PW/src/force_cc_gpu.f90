!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE force_cc_gpu( forcecc )
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
  USE gvect,                ONLY : ngm, gstart, g, g_d, gg, ngl, gl, gl_d, &
                                   igtongl, igtongl_d
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
#if defined(__CUDA)
  USE device_fbuff_m,             ONLY : dev_buf
  USE device_memcpy_m,        ONLY : dev_memcpy
#endif
  !
  IMPLICIT NONE
  !
  REAL(DP) :: forcecc(3,nat)
  !! output: the NLCC forces on atoms
  !
  ! ... local variables
  !
  INTEGER :: ipol, ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms
  REAL(DP), ALLOCATABLE :: vxc(:,:), rhocg(:)
  ! exchange-correlation potential
  ! radial fourier transform of rho core
  COMPLEX(DP), ALLOCATABLE :: vaux(:,:)
  REAL(DP) ::  prod, arg, fact
  !
  real(DP), pointer :: rhocg_d (:),r_d(:), rab_d(:), rhoc_d(:)
  real(DP):: forcelc_x, forcelc_y, forcelc_z, tau1, tau2, tau3
  integer           :: ierrs(4)
  integer           :: maxmesh
#if defined(__CUDA)
  attributes(DEVICE) :: rhocg_d, r_d, rab_d, rhoc_d
  !
  forcecc(:,:) = 0.d0
  !
  IF ( ANY( upf(1:ntyp)%nlcc ) ) GOTO 15
  RETURN
  !
15 CONTINUE
  IF (gamma_only) THEN
     fact = 2.d0
  ELSE
     fact = 1.d0
  ENDIF
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
     DO ir = 1, dfftp%nnr
        vxc(ir,1) = 0.5d0 * ( vxc(ir,1) + vxc(ir,2) )
     ENDDO
  ENDIF
  !
  CALL rho_r2g( dfftp, vxc(:,1:1), vaux(:,1:1) ) 
  !
  ! ... vaux contains now Vxc(G)
  !
  ALLOCATE( rhocg(ngl) )
  !
  ! ... core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  maxmesh = MAXVAL(msh(1:ntyp)) 
  CALL dev_buf%lock_buffer(rhocg_d, ngl, ierrs(1) )
  CALL dev_buf%lock_buffer(r_d, maxmesh, ierrs(2) )
  CALL dev_buf%lock_buffer(rab_d, maxmesh, ierrs(3) )
  CALL dev_buf%lock_buffer(rhoc_d, maxmesh, ierrs(4) )
  IF (ANY(ierrs /= 0)) CALL errore('force_cc_gpu', 'cannot allocate buffers', ABS(MAXVAL(ierrs)) )
  !
  ! ... core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%nlcc ) THEN
        r_d(1:msh(nt)) = rgrid(nt)%r(1:msh(nt)) 
        rab_d(1:msh(nt)) = rgrid(nt)%rab(1:msh(nt))
        rhoc_d(1:msh(nt)) = upf(nt)%rho_atc 
        !
        CALL drhoc_gpu( ngl, gl_d, omega, tpiba2, msh(nt), r_d, &
                         rab_d, rhoc_d, rhocg_d)
        DO na = 1, nat
           IF (nt == ityp (na) ) THEN
              tau1 = tau(1, na)
              tau2 = tau(2, na)
              tau3 = tau(3, na)

              forcelc_x = 0.d0
              forcelc_y = 0.d0
              forcelc_z = 0.d0

              !$acc parallel loop reduction(+:forcelc_x,forcelc_y,forcelc_z)
              do ig = gstart, ngm
                 arg = (g_d (1, ig) * tau1 + g_d (2, ig) * tau2 &
                      + g_d (3, ig) * tau3 ) * tpi
                 prod = tpiba * omega * &
                      rhocg_d (igtongl_d (ig) ) * dble(CONJG(vaux(ig,1)) * &
                      CMPLX( sin (arg), cos (arg), kind=DP))* fact

                 forcelc_x = forcelc_x + g_d (1, ig) * prod
                 forcelc_y = forcelc_y + g_d (2, ig) * prod
                 forcelc_z = forcelc_z + g_d (3, ig) * prod
              ENDDO

              forcecc(1, na) = forcecc(1, na) + forcelc_x
              forcecc(2, na) = forcecc(2, na) + forcelc_y
              forcecc(3, na) = forcecc(3, na) + forcelc_z
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !
  CALL mp_sum( forcecc, intra_bgrp_comm )
  !
  !$acc end data
  DEALLOCATE( vxc, vaux )
  CALL dev_buf%release_buffer(rhocg_d, ierrs(1) )
  CALL dev_buf%release_buffer(r_d, ierrs(2) )
  CALL dev_buf%release_buffer(rab_d, ierrs(3) )
  CALL dev_buf%release_buffer(rhoc_d, ierrs(4) )
#endif
  !
  RETURN
  !
END SUBROUTINE force_cc_gpu
