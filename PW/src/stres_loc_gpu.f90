
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE stres_loc_gpu( sigmaloc )
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : msh, rgrid
  USE m_gth,                ONLY : dvloc_gth_gpu
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : omega, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, ngl, gl, igtongl, &
                                   gl_d, igtongl_d
  USE gvect_gpum,           ONLY : g_d
  USE scf,                  ONLY : rho
  USE vlocal,               ONLY : strf, vloc
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : psic
  USE uspp_param,           ONLY : upf
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE Coul_cut_2D,          ONLY : do_cutoff_2D, cutoff_stres_evloc_gpu, cutoff_stres_sigmaloc_gpu
  !
  USE wavefunctions_gpum,   ONLY : using_psic, using_psic_d, psic_d
  USE device_fbuff_m ,            ONLY : dev_buf
  USE device_memcpy_m,        ONLY : dev_memcpy
  !
  implicit none
  !
  REAL(DP) :: sigmaloc(3,3)
  REAL(DP) :: evloc, fact
  INTEGER :: ng, nt, l, m
  ! counter on g vectors
  ! counter on atomic type
  ! counter on angular momentum
  ! counter on spin components
  !
  INTEGER :: ierrs(5)
  REAL(DP) :: evloc_d, zp_d
  REAL(DP) :: spart, sigma11, sigma21, sigma22, sigma31, sigma32, sigma33
  !
  INTEGER :: mshd
  !
  INTEGER,  POINTER :: nl_d(:)
  REAL(DP), POINTER :: rab_d(:), r_d(:), vloc_d(:,:), dvloc_d(:), &
                       upfvloc_d(:)
  COMPLEX(DP), POINTER :: strf_d(:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: vloc_d, strf_d, nl_d, rab_d, r_d, dvloc_d, &
                        upfvloc_d
#endif
  !
  nl_d => dfftp%nl_d
  !
  sigmaloc(:,:) = 0._DP
  !
  CALL using_psic(2)
  !
  psic(:) = CMPLX(rho%of_r(:,1), KIND=dp)
  !
  CALL using_psic_d(1)
  !
  CALL fwfft( 'Rho', psic_d, dfftp )
  !
  !CALL using_psic_d(0) ; 
  CALL using_psic(0)
  !
  gl_d = gl ; igtongl_d = igtongl
  !
  ! psic contains now the charge density in G space
  fact = 1._DP
  IF (gamma_only) fact = 2._DP
  !
  evloc = 0.0_DP
  !
  CALL dev_buf%lock_buffer( strf_d, (/ ngm, ntyp /), ierrs(1) )
  CALL dev_memcpy( strf_d, strf )
  !
  CALL dev_buf%lock_buffer( vloc_d, (/ ngl, ntyp /), ierrs(2) )
  CALL dev_memcpy( vloc_d, vloc )
  !
  !
  DO nt = 1, ntyp
     IF (gstart==2) evloc = evloc + &
          psic(dfftp%nl(1)) * strf(1,nt) * vloc(igtongl(1),nt)
     !
     evloc_d = 0._DP
     !
     !$cuf kernel do(1)<<<*,*>>>
     DO ng = gstart, ngm
        evloc_d = evloc_d + DBLE(CONJG(psic_d(nl_d(ng))) * strf_d(ng,nt)) &
                            * vloc_d(igtongl_d(ng),nt)
     ENDDO
     !
     evloc = evloc + evloc_d * fact
     !
  ENDDO
  !  
  !
  CALL dev_buf%release_buffer( vloc_d, ierrs(2) )
  !
  CALL dev_buf%lock_buffer( dvloc_d, ngl,  ierrs(2) )
  mshd = MAXVAL(rgrid(1:ntyp)%mesh)
  CALL dev_buf%lock_buffer( rab_d,   mshd, ierrs(3) )
  CALL dev_buf%lock_buffer( r_d,     mshd, ierrs(4) )
  CALL dev_buf%lock_buffer( upfvloc_d, mshd, ierrs(5) )
  !
  !
  ! 2D:  add contribution from cutoff long-range part of Vloc
  IF (do_cutoff_2D) CALL cutoff_stres_evloc_gpu( psic_d, strf_d, evloc )
  !
  !      WRITE( 6,*) ' evloc ', evloc, evloc*omega   ! DEBUG
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%is_gth ) THEN
        !
        ! special case: GTH pseudopotential
        !
        CALL dvloc_gth_gpu( nt, upf(nt)%zp, tpiba2, ngl, gl_d, omega, dvloc_d )
        !
     ELSE IF ( upf(nt)%tcoulombp ) THEN
        !
        ! special case: pseudopotential is coulomb 1/r potential
        !
        CALL dvloc_coul_gpu( upf(nt)%zp, tpiba2, ngl, gl_d, omega, dvloc_d )
        !
     ELSE
        !
        ! normal case: dvloc contains dV_loc(G)/dG
        !
        ! the  G=0 component is not computed
        rab_d(1:msh(nt)) = rgrid(nt)%rab(1:msh(nt))
        r_d(1:msh(nt))   = rgrid(nt)%r(1:msh(nt))
        upfvloc_d(1:msh(nt)) = upf(nt)%vloc(1:msh(nt))
        zp_d                 = upf(nt)%zp
        !
        CALL dvloc_of_g_gpu( rgrid(nt)%mesh, msh(nt), rab_d(1:rgrid(nt)%mesh),   &
                             r_d(1:rgrid(nt)%mesh), upfvloc_d(1:rgrid(nt)%mesh), &
                             zp_d, tpiba2, ngl, gl_d, omega, dvloc_d )
        !
     END IF
     !
     sigma11 = 0._DP ; sigma21 = 0._DP ; sigma22 = 0._DP
     sigma31 = 0._DP ; sigma32 = 0._DP ; sigma33 = 0._DP
     !
     !$cuf kernel do (1) <<<*,*>>>
     DO ng = 1, ngm
       spart = DBLE(CONJG(psic_d(nl_d(ng))) * strf_d(ng,nt)) * 2.0_DP *&
               dvloc_d(igtongl_d(ng))
       sigma11 = sigma11 + spart * g_d(1,ng) * g_d(1,ng)
       sigma21 = sigma21 + spart * g_d(2,ng) * g_d(1,ng)
       sigma22 = sigma22 + spart * g_d(2,ng) * g_d(2,ng)
       sigma31 = sigma31 + spart * g_d(3,ng) * g_d(1,ng)
       sigma32 = sigma32 + spart * g_d(3,ng) * g_d(2,ng)
       sigma33 = sigma33 + spart * g_d(3,ng) * g_d(3,ng)
     ENDDO  
     !
     sigmaloc(1,1) = sigmaloc(1,1) + sigma11 * fact * tpiba2
     sigmaloc(2,1) = sigmaloc(2,1) + sigma21 * fact * tpiba2
     sigmaloc(2,2) = sigmaloc(2,2) + sigma22 * fact * tpiba2
     sigmaloc(3,1) = sigmaloc(3,1) + sigma31 * fact * tpiba2
     sigmaloc(3,2) = sigmaloc(3,2) + sigma32 * fact * tpiba2
     sigmaloc(3,3) = sigmaloc(3,3) + sigma33 * fact * tpiba2
     !
  ENDDO
  !
  IF (do_cutoff_2D) CALL cutoff_stres_sigmaloc_gpu( psic_d, strf_d, sigmaloc ) ! 2D: re-add LR Vloc to sigma here
  !
  do l = 1, 3
     sigmaloc (l, l) = sigmaloc (l, l) + evloc
     do m = 1, l - 1
        sigmaloc (m, l) = sigmaloc (l, m)
     enddo
  enddo
  !
  call mp_sum(  sigmaloc, intra_bgrp_comm )
  !
  CALL dev_buf%release_buffer( strf_d, ierrs(1) )
  CALL dev_buf%release_buffer( dvloc_d, ierrs(2) )
  CALL dev_buf%release_buffer( rab_d, ierrs(3) )
  CALL dev_buf%release_buffer( r_d, ierrs(4) )
  CALL dev_buf%release_buffer( upfvloc_d, ierrs(5) )
  !
  RETURN
  !
END SUBROUTINE stres_loc_gpu

