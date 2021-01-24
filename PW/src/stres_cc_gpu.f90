!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE stres_cc_gpu( sigmaxcc )
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : rgrid, msh
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, ngl, gl, igtongl, igtongl_d
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  USE gvect_gpum,           ONLY : g_d, gg_d
  USE wavefunctions_gpum,   ONLY : using_psic, using_psic_d, psic_d
  USE device_fbuff_m,             ONLY : dev_buf
  USE device_memcpy_m,        ONLY : dev_memcpy
  !
  !
  IMPLICIT NONE
  !
  ! output
  REAL(DP) :: sigmaxcc(3,3)
  ! local variables

  INTEGER :: nt, ng, l, m, ir
  ! counters
  REAL(DP) :: fact
  REAL(DP), ALLOCATABLE :: vxc(:,:)
  !
  INTEGER,  POINTER :: nl_d(:)
  REAL(DP), POINTER :: rhocg_d(:), r_d(:), rab_d(:), rhoc_d(:), gl_d(:)
  COMPLEX(DP), POINTER :: strf_d(:)
  !
  INTEGER :: maxmesh, ierrs(6)
  REAL(DP) :: rhocg1(1), sigma_rid, sigmadiag
  REAL(DP) :: sigma1, sigma2, sigma3, &
              sigma4, sigma5, sigma6
  !
#if defined(__CUDA)
  attributes(DEVICE) :: rhocg_d, nl_d, r_d, rab_d, rhoc_d, &
                        gl_d, strf_d, nl_d
#endif
  !
  nl_d => dfftp%nl_d
  !
  sigmaxcc(:,:) = 0._DP
  IF ( ANY( upf(1:ntyp)%nlcc ) ) GOTO 15
  !
  RETURN
  !
15 CONTINUE
  !
  ! recalculate the exchange-correlation potential
  !
  ALLOCATE( vxc(dfftp%nnr,nspin) )
  !
  CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
  !
  CALL using_psic(2)
  !
  IF (nspin==1 .OR. nspin==4) THEN
     DO ir = 1, dfftp%nnr
        psic(ir) = CMPLX(vxc(ir,1))
     ENDDO
  ELSE
     DO ir = 1, dfftp%nnr
        psic(ir) = CMPLX(0.5_DP * (vxc(ir,1) + vxc(ir,2)))
     ENDDO
  ENDIF
  !
  DEALLOCATE( vxc )
  !
  CALL using_psic(0)
  CALL using_psic_d(1)
  !
  CALL fwfft( 'Rho', psic_d, dfftp )
  !
  CALL using_psic(0)
  ! psic contains now Vxc(G)
  !
  sigmadiag = 0._DP
  !
  fact = 1._DP
  IF (gamma_only) fact = 2._DP
  !
  maxmesh = MAXVAL(msh(1:ntyp)) 
  CALL dev_buf%lock_buffer( gl_d, ngl, ierrs(1) )
  CALL dev_memcpy( gl_d, gl, (/ 1, ngl /) )
  CALL dev_buf%lock_buffer( rhocg_d,   ngl, ierrs(2) )
  CALL dev_buf%lock_buffer( r_d,   maxmesh, ierrs(3) )
  CALL dev_buf%lock_buffer( rab_d, maxmesh, ierrs(4) )
  CALL dev_buf%lock_buffer( rhoc_d,maxmesh, ierrs(5) )
  CALL dev_buf%lock_buffer( strf_d,    ngm, ierrs(6) )
  IF (ANY(ierrs /= 0)) CALL errore( 'stres_cc_gpu', 'cannot allocate buffers', -1 )
  !
  sigma1 = 0._DP ;  sigma4 = 0._DP
  sigma2 = 0._DP ;  sigma5 = 0._DP
  sigma3 = 0._DP ;  sigma6 = 0._DP
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%nlcc ) THEN
        !
        CALL dev_memcpy( strf_d, strf(:,nt),  (/1, ngm/)     )
        CALL dev_memcpy( r_d,    rgrid(nt)%r, (/1, msh(nt)/) )
        CALL dev_memcpy( rab_d,  rgrid(nt)%rab,   (/1, msh(nt)/) )
        CALL dev_memcpy( rhoc_d, upf(nt)%rho_atc, (/1, msh(nt)/) )
        !
        CALL drhoc_gpu( ngl, gl_d, omega, tpiba2, msh(nt), r_d, &
                        rab_d, rhoc_d, rhocg_d )
        !
        ! diagonal term
        IF (gstart==2) THEN
          rhocg1=rhocg_d(igtongl(1:1))
          sigmadiag = sigmadiag + DBLE(CONJG(psic(dfftp%nl(1))) * &
                       strf(1,nt)) * rhocg1(1)
        ENDIF
        !
        !$cuf kernel do (1) <<<*,*>>>
        DO ng = gstart, ngm
           sigmadiag = sigmadiag + DBLE(CONJG(psic_d(nl_d(ng))) * &
                        strf_d(ng)) * rhocg_d(igtongl_d(ng) ) * fact
        ENDDO
        !
        
        !
        CALL deriv_drhoc_gpu( ngl, gl_d, omega, tpiba2, msh(nt), &
                              r_d, rab_d, rhoc_d, rhocg_d )
        !
        ! non diagonal term (g=0 contribution missing)
        !
        !$cuf kernel do (1) <<<*,*>>>
        DO ng = gstart, ngm
          !
          sigma_rid = DBLE(CONJG(psic_d(nl_d(ng))) &
                      * strf_d(ng)) * rhocg_d(igtongl_d(ng)) * tpiba &
                      / SQRT(gg_d(ng)) * fact
          !
          sigma1 = sigma1 + sigma_rid * g_d(1,ng)*g_d(1,ng)
          sigma2 = sigma2 + sigma_rid * g_d(1,ng)*g_d(2,ng)
          sigma3 = sigma3 + sigma_rid * g_d(1,ng)*g_d(3,ng)
          sigma4 = sigma4 + sigma_rid * g_d(2,ng)*g_d(2,ng)
          sigma5 = sigma5 + sigma_rid * g_d(3,ng)*g_d(2,ng)
          sigma6 = sigma6 + sigma_rid * g_d(3,ng)*g_d(3,ng)
          !
        ENDDO
        !
     ENDIF
     !
  ENDDO
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
  CALL dev_buf%release_buffer( gl_d,   ierrs(1) )
  CALL dev_buf%release_buffer( rhocg_d,ierrs(2) )
  CALL dev_buf%release_buffer( r_d,    ierrs(3) )
  CALL dev_buf%release_buffer( rab_d,  ierrs(4) )
  CALL dev_buf%release_buffer( rhoc_d, ierrs(5) )
  CALL dev_buf%release_buffer( strf_d, ierrs(6) )
  !
  RETURN
  !
END SUBROUTINE stres_cc_gpu

