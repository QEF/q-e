!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_dvpsi_eels (ik, dvpsi1, dvpsi2)
  !--------------------------------------------------------------------------
  !
  ! This subroutine calculates the starting Lanczos vectors.
  ! EELS: The external perturbation is an electron (plane wave) 
  ! with a wavevector q.
  !
  ! In the case of metals, the projector P^+_c must not be present
  ! in the calculation of the zeta coefficients (see lr_lanczos),
  ! hence we distinguish between dvpsi1 and dvpsi2:
  !
  ! dvpsi1 = S^{-1} P^+_c \int K(r') exp(iq*r') dr' psi_n,k
  ! dvpsi2 =              \int K(r') exp(iq*r') dr' psi_n,k 
  ! 
  ! K(r;r1,r2) = \delta(r-r1) \delta(r-r2) + \sum Q(r-R) beta(r1-R) beta^*(r2-R)
  ! See Eq.(3) in Ref. J. Tobik and A. Dal Corso, JCP 120, 9934 (2004)
  !
  ! Writtent by Iurii Timrov (2013)
  !
  USE kinds,                 ONLY : DP
  USE wvfct,                 ONLY : npwx, nbnd
  USE fft_base,              ONLY : dffts, dtgs
  USE gvecw,                 ONLY : gcutw
  USE qpoint,                ONLY : ikks, ikqs, nksq 
  USE eqv,                   ONLY : evq, dpsi 
  USE wavefunctions_module,  ONLY : evc
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag
  use klist,                 only : xk, igk_k, ngk
  use gvect,                 only : ngm, g
  USE control_lr,            ONLY : nbnd_occ
  USE io_files,              ONLY : iunwfc, nwordwfc
  use uspp,                  only : vkb, okvan
  USE mp_bands,              ONLY : ntask_groups
  USE buffers,               ONLY : get_buffer
 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik 
  COMPLEX(DP), INTENT(OUT) :: dvpsi1(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(OUT) :: dvpsi2(npwx*npol,nbnd)
  INTEGER :: ibnd, ig, ios, incr, v_siz
  INTEGER :: ikk, & ! index of the point k
             ikq, & ! index of the point k+q
             npwq   ! number of the plane-waves at point k+q
  COMPLEX(DP), ALLOCATABLE :: revc(:,:), &  ! wavefunctions in R-space
                            & tg_psic(:,:)  ! wavefunctions in R-space (for task groups) 
  LOGICAL :: exst
  !
  CALL start_clock ('lr_dvpsi_eels') 
  !
  ALLOCATE (revc(dffts%nnr,npol))
  revc(:,:)   = (0.d0, 0.d0)
  dpsi(:,:)   = (0.d0, 0.d0)
  dvpsi1(:,:) = (0.d0, 0.d0)
  dvpsi2(:,:) = (0.d0, 0.d0)
  !
  incr = 1
  !
  IF ( dtgs%have_task_groups ) THEN
     !
     v_siz =  dtgs%tg_nnr * dtgs%nogrp
     !
     ALLOCATE( tg_psic(v_siz,npol) )
     !
     incr = dtgs%nogrp
     !
  ENDIF 
  !
  ikk  = ikks(ik)
  ikq  = ikqs(ik)
  npwq = ngk(ikq)
  !
  ! Read the ground-state wavefunctions evc(k) and evq(k+q).
  ! Note: even in the case when there is only one k point (nksq=1),
  ! the ground-state wavefunctions evc(k) and evq(k+q) are read here
  ! and kept in memory for the rest of the code.
  !
  CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
  CALL get_buffer (evq, nwordwfc, iunwfc, ikq)
  !
  ! Re-ordering of the G vectors.
  !
  DO ibnd = 1, nbnd_occ(ikk), incr
     !
     IF ( dtgs%have_task_groups ) THEN
        !
        ! FFT to R-space
        CALL cft_wave_tg(ik, evc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
        !
        ! back-FFT to G-space
        CALL cft_wave_tg(ik, dvpsi1, tg_psic, -1, v_siz, ibnd, nbnd_occ(ikk))
        !
     ELSE
        !
        ! FFT to R-space
        CALL cft_wave(ik, evc(1,ibnd), revc, +1)
        !
        ! back-FFT to G-space
        CALL cft_wave(ik, dvpsi1(1,ibnd), revc, -1) 
        !
     ENDIF
     !
  ENDDO
  !
  ! Compute and add the US term.
  !
  IF (okvan) THEN
     !
     ! Calculate beta-functions vkb at point k+q
     CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb)
     !
     CALL lr_addus_dvpsi (ik, npwx, npwq, nbnd_occ(ikk), dvpsi1, dpsi)
     !
     dvpsi1 = dpsi
     dpsi(:,:) = (0.d0, 0.d0)
     !
  ENDIF
  !
  dvpsi2 = dvpsi1
  !
  ! Ortogonalize dvpsi1 to valence states.
  ! Apply -P_c^+, and then change the sign, because we need P_c^+.
  !
  CALL orthogonalize(dvpsi1, evq, ikk, ikq, dpsi, npwq, .false.) 
  dvpsi1 = -dvpsi1
  !
  ! In the ultrasoft case apply the S^{-1} operator.
  !
  IF (okvan) THEN
     !
     dpsi(:,:) = (0.0d0, 0.0d0)
     !
     CALL lr_sm1_psiq (.TRUE., ik, npwx, npwq, nbnd_occ(ikk), dvpsi1, dpsi)
     !
     dvpsi1 = dpsi
     !
  ENDIF
  !
  DEALLOCATE (revc)
  !
  IF ( dtgs%have_task_groups ) THEN
     DEALLOCATE( tg_psic )
  ENDIF
  !
  CALL stop_clock ('lr_dvpsi_eels')
  !
  RETURN
  !
END SUBROUTINE lr_dvpsi_eels
