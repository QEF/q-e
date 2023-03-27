!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------------
SUBROUTINE lr_apply_liouvillian_magnons( evc1, evc1_new, L_dag )
  !------------------------------------------------------------------------------
  !
  ! This subroutine applies the linear response operator to response wavefunctions.
  ! S^{-1} { (H - E*S)*psi(k+q) + P_c^+(k+q) V_HXC(q)*psi0(k) }
  !
  ! Inspired by PH/solve_linter.f90
  !
  ! Written by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dfftp, dffts
  USE klist,                ONLY : xk, igk_k, ngk
  USE lr_variables,         ONLY : no_hxc, iunTwfc
  USE lsda_mod,             ONLY : current_spin
  USE wvfct,                ONLY : nbnd, npwx, et, current_k
  USE uspp,                 ONLY : vkb
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE wavefunctions,        ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, domag, npol, nspin_mag
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : ntask_groups, me_bgrp
  USE buffers,              ONLY : get_buffer
  USE qpoint,               ONLY : ikks, ikqs, nksq
  USE eqv,                  ONLY : evq, dpsi, dvpsi
  USE control_lr,           ONLY : nbnd_occ, nbnd_occx
  USE dv_of_drho_lr
  USE fft_helper_subroutines
  USE fft_interfaces,       ONLY : fft_interpolate
  USE scf,                  ONLY : vrs
 
  USE io_global,             ONLY : stdout
  USE uspp_init,             ONLY : init_us_2
  USE scf_gpum,              ONLY : vrs_d

  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: evc1(npwx*npol,nbnd_occx,nksq,2)
  COMPLEX(DP), INTENT(OUT) :: evc1_new(npwx*npol,nbnd_occx,nksq,2)
  LOGICAL,     INTENT(IN)  :: L_dag
  !
  ! Local variables
  !
  LOGICAL :: interaction1
  INTEGER :: i,j,ir, ibnd, ig, ia, ios, is, incr, v_siz, ipol
  INTEGER :: ik,  &
             ikk, & ! index of the point k
             ikq, & ! index of the point k+q
             npwq   ! number of the plane-waves at point k+q
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:), revc(:,:), &
                            & dvrsc(:,:), dvrssc(:,:), &
                            & sevc1_new(:,:,:,:), &
  ! Change of the Hartree and exchange-correlation (HXC) potential
                            & tg_psic(:,:), tg_dvrssc(:,:)
  ! Task groups: wfct in R-space and the response HXC potential
  !
  ! -k part
  !
  INTEGER :: imk, imkq
  !
  COMPLEX(DP)  :: Tevc(npol*npwx,nbnd)      ! T-rev op. applied to u_{-k}
  COMPLEX(DP)  :: Tevq(npol*npwx,nbnd)      ! T-rev op. applied to u_{-k-Q}
  !
  CALL start_clock('lr_apply')
  !
  ALLOCATE (hpsi(npwx*npol,nbnd))
  ALLOCATE (spsi(npwx*npol,nbnd))
  ALLOCATE (sevc1_new(npwx*npol,nbnd_occx,nksq,2))
  ALLOCATE (revc(dffts%nnr,npol))
  ALLOCATE (dvrsc(dfftp%nnr,nspin_mag))
  ALLOCATE (dvrssc(dffts%nnr,nspin_mag))
  IF (.not. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  hpsi(:,:)   = (0.d0,0.d0)
  spsi(:,:)   = (0.d0,0.d0)
  revc(:,:)   = (0.d0,0.d0)
  dvrsc(:,:)  = (0.0d0,0.0d0)
  dvrssc(:,:) = (0.0d0,0.0d0)
  sevc1_new(:,:,:,:) = (0.0d0,0.0d0)
  Tevc(:,:) = (0.0d0,0.0d0)
  Tevq(:,:) = (0.0d0,0.0d0)
  evc1_new(:,:,:,:) = (0.0d0,0.0d0)
  !
  incr = 1
  !
  IF ( dffts%has_task_groups ) THEN
     !
     v_siz =  dffts%nnr_tg 
     !
     ALLOCATE( tg_dvrssc(v_siz,nspin_mag) )
     ALLOCATE( tg_psic(v_siz,npol) )
     !
     incr = fftx_ntgrp(dffts)
     !
  ENDIF
  !
  IF (no_hxc) THEN
     write(stdout,*) 'no_hxc = .true., only D-part computed'
  ENDIF
  !
  IF (.not. no_hxc)      CALL start_clock('lr_apply_int')
  IF (no_hxc) CALL start_clock('lr_apply_no')
  !
  ! Calculation of the charge density response and the 
  ! corresponding change of the Hartree and exchange-correlation 
  ! potentials. If no_hxc=.true. this part is not needed and thus
  ! is not calculated.
  !
  IF ( .not. no_hxc ) THEN
     ! 
     ! Calculation of the response charge density 
     ! and its symmetrization.
     !
     call lr_calc_dens_magnons (dvrsc(:,:), evc1(:,:,:,:), L_dag)
     !
     ! Calculation of the response HXC potential
     ! from the response charge density.
     !
     CALL dv_of_drho (dvrsc, .false.)
     !
     ! Interpolation of the HXC potential from the thick mesh 
     ! to a smoother mesh (if doublegrid=.true.)
     ! dvrsc -> dvrssc
     !
     DO is = 1, nspin_mag
        CALL fft_interpolate (dfftp, dvrsc(:,is), dffts, dvrssc(:,is))
     ENDDO
     !
  ENDIF
  !
  ! Now we will calculate two terms:
  ! 1) HXC term : P_c^+(k+q) V_HXC(q)*psi0(k) 
  ! 2) (H - E)*psi(k+q)
  !
  DO ik = 1, nksq
     !
     !====================================================
     !   Sternheimer equation for X
     !====================================================
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npwq = ngk(ikq)
     !
     ! Calculate beta-functions vkb at k+q (Kleinman-Bylander projectors)
     ! The vkb's are needed for the non-local potential in h_psi,
     ! and for the ultrasoft term.
     !
     CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkb, .true.)
     !$acc update host(vkb)
     !
     ! Read unperturbed wavefuctions evc (wfct at k) 
     ! and evq (wfct at k+q)
     !
     CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
     CALL get_buffer (evq, nwordwfc, iunwfc, ikq)
     !
     dpsi(:,:)  = (0.d0,0.d0)
     dvpsi(:,:) = (0.d0,0.d0)
     !
     ! 1) Hartree and exchange-correlation term.  
     !    The multiplication of the HXC potential dvrssc with
     !    the unperturbed wavefunctions revc is done in R-space.
     !    If interaction1=.true. calculate HXC term
     !    If interaction1=.false. skip this step and go to 2)
     !
!     IF (interaction1) THEN
     IF ( .not. no_hxc) THEN
        !
        !$acc data copyin(evc, dvrssc) create(revc) copy (dvpsi)
        !
        ! The potential in dvrssc is distributed across all processors.
        ! We need to redistribute it so that it is completely contained in the
        ! processors of an orbital TASK-GROUP.
        !
        IF ( dffts%has_task_groups ) THEN
           !
           IF (noncolin) THEN
              !
              CALL tg_cgather( dffts, dvrssc(:,1), tg_dvrssc(:,1))
              !
              IF (domag) THEN
                 DO ipol = 2, 4
                    CALL tg_cgather( dffts, dvrssc(:,ipol), tg_dvrssc(:,ipol))
                 ENDDO
              ENDIF
              !
           ELSE
              !
              CALL tg_cgather( dffts, dvrssc(:,current_spin), tg_dvrssc(:,1))
              !
           ENDIF
           !
        ENDIF
        !
        DO ibnd = 1, nbnd_occ(ikk), incr
           !
           IF ( dffts%has_task_groups ) THEN
              !
              ! FFT to R-space
              !
              CALL cft_wave_tg(ik, evc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
              !
              ! Multiply the HXC potential with unperturbed wfct's
              !
              CALL apply_dpot(v_siz, tg_psic, tg_dvrssc, 1)
              !
              ! back-FFT to G-space
              !
              CALL cft_wave_tg(ik,dvpsi, tg_psic,-1, v_siz, ibnd, nbnd_occ(ikk))
              !
           ELSE
              !
              ! FFT to R-space
              !
              CALL cft_wave(ik, evc(1,ibnd), revc, +1)
              !
              ! Multiply the HXC potential with unperturbed wfct's 
              !
              CALL apply_dpot(dffts%nnr, revc, dvrssc, current_spin)
              !
              ! back-FFT to G-space
              !
              CALL cft_wave(ik, dvpsi(1,ibnd), revc, -1)
              !
           ENDIF
           !
        ENDDO
        !$acc end data
        !
        ! In the case of US pseudopotentials there is an additional term.
        ! See the second term in Eq.(11) in J. Chem. Phys. 127, 164106 (2007).
        !
!        IF (okvan) CALL adddvscf(1, ik) 
        !
        ! Ortogonalize dvpsi to valence states.
        ! Apply -P_c^+, and then change the sign, because we need P_c^+.
        !
        CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .false.)
        dvpsi = -dvpsi
        !
     ENDIF
     !
     ! 2) (H - E*S) * psi(k+q)
     !
     ! Compute the kinetic energy g2kin: (k+q+G)^2
     !
     CALL g2_kin(ikq)
     !
     IF (.NOT. ALLOCATED(psic_nc)) ALLOCATE(psic_nc(dffts%nnr,npol))
     !
     ! Apply the operator ( H - \epsilon S + alpha_pv P_v) to evc1
     ! where alpha_pv = 0
     !
     !$acc data copyin(evq) copy(evc1(1:npwx*npol,1:nbnd,ik,1),sevc1_new(1:npwx*npol,1:nbnd,ik), et(:,ikk))
     CALL ch_psi_all (npwq, evc1(:,:,ik,1), sevc1_new(:,:,ik,1), et(:,ikk), ik, nbnd_occ(ikk)) 
     !$acc end data
     !
     IF (ALLOCATED(psic_nc)) DEALLOCATE(psic_nc)
     !
     ! 3) Sum up the two terms : (H - E*S)*psi(k+q) + P_c^+(k+q) V_HXC(q)*psi0(k)
     !
!     IF (interaction1) THEN
!     IF ( .not. no_hxc) THEN
        !
        DO ibnd = 1, nbnd_occ(ikk)
           DO ig = 1, npwq
              sevc1_new(ig,ibnd,ik,1) = sevc1_new(ig,ibnd,ik,1) + dvpsi(ig,ibnd)
              !
              sevc1_new(ig+npwx,ibnd,ik,1) = sevc1_new(ig+npwx,ibnd,ik,1) &
                                                & + dvpsi(ig+npwx,ibnd)
           ENDDO
        ENDDO
        !
!     ENDIF
     !
     evc1_new(:,:,ik,1) = sevc1_new(:,:,ik,1)
     !
     !==========================================================
     !   Sternheimer equation for Y
     !==========================================================
     !
     dvrssc(:,2) = -dvrssc(:,2)
     dvrssc(:,3) = -dvrssc(:,3)
     dvrssc(:,4) = -dvrssc(:,4)
     !
     ! indexes of  -k and -k-Q, needed for Y-equation
     !
     IF ( mod(ik,2) == 0) THEN
        imk  = ikk - 3      ! position of -k
        imkq = ikk - 1      ! position of -k-Q
     ELSE
        imk  = ikk + 3      ! position of -k
        imkq = ikk + 5      ! position of -k-Q
     ENDIF
     !
     ! Read unperturbed wavefuctions Tu_{-k} and Tu_{-k-Q}
     !
     CALL get_buffer (Tevc, nwordwfc, iunTwfc, 2*ik-1)
     CALL get_buffer (Tevq, nwordwfc, iunTwfc, 2*ik  )
     !
     dpsi(:,:)  = (0.d0,0.d0)
     dvpsi(:,:) = (0.d0,0.d0)
     !
     ! 1) Hartree and exchange-correlation term.  
     !    The multiplication of the HXC potential dvrssc with
     !    the unperturbed wavefunctions revc is done in R-space.
     !    If interaction1=.true. calculate HXC term
     !    If interaction1=.false. skip this step and go to 2)
     !
!     IF (interaction1) THEN
     IF ( .not. no_hxc) THEN
        !
        !$acc data copyin(evc, dvrssc) create(revc) copy (dvpsi)
        !
        ! The potential in dvrssc is distributed across all processors.
        ! We need to redistribute it so that it is completely contained in the
        ! processors of an orbital TASK-GROUP.
        !
        IF ( dffts%has_task_groups ) THEN
           !
           IF (noncolin) THEN
              !
              CALL tg_cgather( dffts, dvrssc(:,1), tg_dvrssc(:,1))
              !
              IF (domag) THEN
                 DO ipol = 2, 4
                    CALL tg_cgather( dffts, dvrssc(:,ipol), tg_dvrssc(:,ipol))
                 ENDDO
              ENDIF
              !
           ELSE
              !
              CALL tg_cgather( dffts, dvrssc(:,current_spin), tg_dvrssc(:,1))
              !
           ENDIF
           !
        ENDIF
        !
        DO ibnd = 1, nbnd_occ(imk), incr
           !
           IF ( dffts%has_task_groups ) THEN
              !
              ! FFT to R-space
              !
              CALL cft_wave_tg(ik, Tevc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
              !
              ! Multiply the HXC potential with unperturbed wfct's
              !
              CALL apply_dpot(v_siz, tg_psic, tg_dvrssc, 1)
              !
              ! back-FFT to G-space
              !
              CALL cft_wave_tg(ik,dvpsi, tg_psic,-1, v_siz, ibnd, nbnd_occ(ikk))
              !
           ELSE
              !
              ! FFT to R-space
              !
              CALL cft_wave(ik, Tevc(1,ibnd), revc, +1)
              !
              ! Multiply the HXC potential with unperturbed wfct's 
              !
              CALL apply_dpot(dffts%nnr, revc, dvrssc, current_spin)
              !
              ! back-FFT to G-space
              !
              CALL cft_wave(ik, dvpsi(1,ibnd), revc, -1)
              !
           ENDIF
           !
        ENDDO
        !$acc end data
        !
        ! In the case of US pseudopotentials there is an additional term.
        ! See the second term in Eq.(11) in J. Chem. Phys. 127, 164106 (2007).
        !
!        IF (okvan) CALL adddvscf(1, ik)
        !
        ! Ortogonalize dvpsi to valence states.
        ! Apply -P_c^+, and then change the sign, because we need P_c^+.
        !
        CALL orthogonalize(dvpsi, Tevq, imk, imkq, dpsi, npwq, .false.)
        IF ( .not. L_dag ) dvpsi = -dvpsi
        !
     ENDIF
     !
     ! 2) (H - E*S) * psi(k+q)
     !
     IF (.NOT. ALLOCATED(psic_nc)) ALLOCATE(psic_nc(dffts%nnr,npol))
     !
     ! Change the sign of b_xc
     !
     vrs(:,2) = - vrs(:,2)
     vrs(:,3) = - vrs(:,3)
     vrs(:,4) = - vrs(:,4)     
     !
#if defined(__CUDA)
     vrs_d = vrs
#endif
     ! Apply the operator ( H - \epsilon S + alpha_pv P_v) to evc1
     ! where alpha_pv = 0
     !
     !$acc data copyin(evq) copy(evc1(1:npwx*npol,1:nbnd,ik,2),sevc1_new(1:npwx*npol,1:nbnd,ik,2), et(:,imk))
     CALL ch_psi_all (npwq, evc1(:,:,ik,2), sevc1_new(:,:,ik,2), et(:,imk), ik, nbnd_occ(imk))
     !$acc end data
     !
     ! Change the sign of b_xc back
     !
     vrs(:,2) = - vrs(:,2)
     vrs(:,3) = - vrs(:,3)
     vrs(:,4) = - vrs(:,4)
     !
#if defined(__CUDA)
     vrs_d = vrs
#endif
     !
     IF (ALLOCATED(psic_nc)) DEALLOCATE(psic_nc)
     !
     ! 3) Sum up the two terms : (H - E*S)*psi(-k-q) + P_c^+(-k-q)
     ! V_HXC(q)*psi0(k)
     !
!     IF (interaction1) THEN
!     IF (.not. no_hxc) THEN
        !
        DO ibnd = 1, nbnd_occ(imk)
           DO ig = 1, npwq
              sevc1_new(ig,ibnd,ik,2) = sevc1_new(ig,ibnd,ik,2) + dvpsi(ig,ibnd)
              !
              sevc1_new(ig+npwx,ibnd,ik,2) = sevc1_new(ig+npwx,ibnd,ik,2) &
                                                & + dvpsi(ig+npwx,ibnd)
           ENDDO
        ENDDO
        !
!     ENDIF
     !
     evc1_new(:,:,ik,2) = - sevc1_new(:,:,ik,2)
     !
     ! Change sign to b'_xc
     !
     dvrssc(:,2) = - dvrssc(:,2)
     dvrssc(:,3) = - dvrssc(:,3)
     dvrssc(:,4) = - dvrssc(:,4)
     !
  ENDDO ! loop on ik
  !
  DEALLOCATE (dvrsc)
  DEALLOCATE (dvrssc) 
  DEALLOCATE (hpsi)
  DEALLOCATE (spsi)
  DEALLOCATE (revc)
  DEALLOCATE (sevc1_new)
  IF (ALLOCATED(psic)) DEALLOCATE(psic)
  !
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE( tg_dvrssc )
     DEALLOCATE( tg_psic )
  ENDIF
  !
  IF (.not. no_hxc)      CALL stop_clock('lr_apply_int')
  IF (no_hxc) CALL stop_clock('lr_apply_no')
  CALL stop_clock('lr_apply')
  !
  RETURN
  !
END SUBROUTINE lr_apply_liouvillian_magnons
