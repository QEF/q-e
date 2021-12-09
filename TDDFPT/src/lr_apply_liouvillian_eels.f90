!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------------
SUBROUTINE lr_apply_liouvillian_eels ( evc1, evc1_new, interaction )
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
  USE lr_variables,         ONLY : no_hxc
  USE lsda_mod,             ONLY : current_spin
  USE wvfct,                ONLY : nbnd, npwx, et, current_k
  USE uspp,                 ONLY : vkb
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE wavefunctions, ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : ntask_groups, me_bgrp
  USE buffers,              ONLY : get_buffer
  USE qpoint,               ONLY : ikks, ikqs, nksq
  USE eqv,                  ONLY : evq, dpsi, dvpsi
  USE control_lr,           ONLY : nbnd_occ
  USE dv_of_drho_lr
  USE fft_helper_subroutines
  USE fft_interfaces,       ONLY : fft_interpolate
  USE apply_dpot_mod,       ONLY : apply_dpot_allocate, apply_dpot_deallocate, &
                                   apply_dpot_bands
  USE uspp_init,            ONLY : init_us_2
 
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: evc1(npwx*npol,nbnd,nksq)
  COMPLEX(DP), INTENT(OUT) :: evc1_new(npwx*npol,nbnd,nksq)
  LOGICAL,     INTENT(IN)  :: interaction
  !
  ! Local variables
  !
  LOGICAL :: interaction1
  INTEGER :: i,j,ir, ibnd, ig, ia, ios, is
  INTEGER :: ik,  &
             ikk, & ! index of the point k
             ikq, & ! index of the point k+q
             npwq   ! number of the plane-waves at point k+q
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:), &
                            & dvrsc(:,:), dvrssc(:,:), &
                            & sevc1_new(:,:,:)
  ! Change of the Hartree and exchange-correlation (HXC) potential
  !
  CALL start_clock('lr_apply')
  !
  ALLOCATE (hpsi(npwx*npol,nbnd))
  ALLOCATE (spsi(npwx*npol,nbnd))
  ALLOCATE (sevc1_new(npwx*npol,nbnd,nksq))
  ALLOCATE (dvrsc(dfftp%nnr,nspin_mag))
  ALLOCATE (dvrssc(dffts%nnr,nspin_mag))
  IF (.not. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  CALL apply_dpot_allocate()
  hpsi(:,:)   = (0.d0,0.d0)
  spsi(:,:)   = (0.d0,0.d0)
  dvrsc(:,:)  = (0.0d0,0.0d0)
  dvrssc(:,:) = (0.0d0,0.0d0)
  sevc1_new(:,:,:) = (0.0d0,0.0d0)
  !
  IF (no_hxc) THEN
     interaction1 = .false.
  ELSE
     interaction1 = interaction
  ENDIF
  !
  IF (interaction1)      CALL start_clock('lr_apply_int')
  IF (.not.interaction1) CALL start_clock('lr_apply_no')
  !
  ! Calculation of the charge density response and the 
  ! corresponding change of the Hartree and exchange-correlation 
  ! potentials. If no_hxc=.true. this part is not needed and thus
  ! is not calculated.
  !
  IF ( interaction1 ) THEN
     ! 
     ! Calculation of the response charge density 
     ! and its symmetrization.
     !
     !if (.not. allocated(psic)) allocate(psic(dfftp%nnr))   
     !
     IF (noncolin) THEN
        call lr_calc_dens_eels_nc (dvrsc(1,1), evc1(1,1,1))
     ELSE
        call lr_calc_dens_eels (dvrsc(1,current_spin), evc1(1,1,1))
     ENDIF
     !
     !if (allocated(psic)) deallocate(psic) 
     !
     ! Calculation of the response HXC potential
     ! from the response charge density.
     !
     CALL dv_of_drho (dvrsc, .false.)
     !
     ! USPP: Compute the integral of the HXC response potential with the Q function.
     ! Input : dvrsc = V_HXC(r), Output: int3 = \int V_HXC(r) * Q_nm(r) dr 
     ! See Eq.(B22) in Ref. A. Dal Corso, PRB 64, 235118 (2001)
     !
     IF (okvan) CALL newdq(dvrsc, 1)
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
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     npwq = ngk(ikq)
     !
     ! Calculate beta-functions vkb at k+q (Kleinman-Bylander projectors)
     ! The vkb's are needed for the non-local potential in h_psi,
     ! and for the ultrasoft term.
     !
     CALL init_us_2 (npwq, igk_k(1,ikq), xk(1,ikq), vkb)
     !
     ! Read unperturbed wavefuctions evc (wfct at k) 
     ! and evq (wfct at k+q)
     !
     IF (nksq > 1) THEN 
        CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
        CALL get_buffer (evq, nwordwfc, iunwfc, ikq)
     ENDIF
     !
     dpsi(:,:)  = (0.d0,0.d0)
     dvpsi(:,:) = (0.d0,0.d0)
     !
     ! 1) Hartree and exchange-correlation term.
     !    If interaction1=.true. calculate HXC term
     !    If interaction1=.false. skip this step and go to 2)
     !
     IF (interaction1) THEN
        !
        ! The potential in dvrssc is distributed across all processors.
        ! We need to redistribute it so that it is completely contained in the
        ! processors of an orbital TASK-GROUP.
        !
        CALL apply_dpot_bands(ik, nbnd_occ(ikk), dvrssc, evc, dvpsi)
        !
        ! In the case of US pseudopotentials there is an additional term.
        ! See the second term in Eq.(11) in J. Chem. Phys. 127, 164106 (2007).
        !
        IF (okvan) CALL adddvscf(1, ik) 
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
     IF (noncolin) THEN
        IF (.NOT. ALLOCATED(psic_nc)) ALLOCATE(psic_nc(dfftp%nnr,npol))
     !ELSE
     !  IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
     ENDIF
     !
     ! Apply the operator ( H - \epsilon S + alpha_pv P_v) to evc1
     ! where alpha_pv = 0
     !
     CALL ch_psi_all (npwq, evc1(:,:,ik), sevc1_new(:,:,ik), et(:,ikk), ik, nbnd_occ(ikk)) 
     !
     IF (noncolin) THEN
        IF (ALLOCATED(psic_nc)) DEALLOCATE(psic_nc)
     !ELSE
     !  IF (ALLOCATED(psic)) DEALLOCATE(psic)
     ENDIF 
     !
     ! 3) Sum up the two terms : (H - E*S)*psi(k+q) + P_c^+(k+q) V_HXC(q)*psi0(k)
     !
     IF (interaction1) THEN
        !
        DO ibnd = 1, nbnd_occ(ikk)
           DO ig = 1, npwq
              sevc1_new(ig,ibnd,ik) = sevc1_new(ig,ibnd,ik) + dvpsi(ig,ibnd)
           ENDDO
        ENDDO
        !
        IF (noncolin) THEN
           DO ibnd = 1, nbnd_occ(ikk)
              DO ig = 1, npwq
                 sevc1_new(ig+npwx,ibnd,ik) = sevc1_new(ig+npwx,ibnd,ik) &
                                                   & + dvpsi(ig+npwx,ibnd)
              ENDDO
           ENDDO
        ENDIF
        !
     ENDIF
     !
     ! 4) Ultrasoft case: apply the S^{-1} operator.
     !    evc1_new = S^{-1} * sevc1_new
     !    If not ultrasoft: evc1_new = sevc1_new
     !
     CALL lr_sm1_psi(ik, npwx, npwq, nbnd_occ(ikk), &
                          & sevc1_new(1,1,ik), evc1_new(1,1,ik))
     !
  ENDDO ! loop on ik
  !
  CALL apply_dpot_deallocate()
  DEALLOCATE (dvrsc)
  DEALLOCATE (dvrssc)
  DEALLOCATE (hpsi)
  DEALLOCATE (spsi)
  DEALLOCATE (sevc1_new)
  IF (ALLOCATED(psic)) DEALLOCATE(psic)
  !
  IF (interaction1)      CALL stop_clock('lr_apply_int')
  IF (.NOT.interaction1) CALL stop_clock('lr_apply_no')
  CALL stop_clock('lr_apply')
  !
  RETURN
  !
END SUBROUTINE lr_apply_liouvillian_eels
