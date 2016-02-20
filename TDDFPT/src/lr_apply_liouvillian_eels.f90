!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------------
SUBROUTINE lr_apply_liouvillian_eels ( evc1, evc1_new, sevc1_new, interaction )
  !------------------------------------------------------------------------------
  !
  ! This subroutine applies the linear response operator to response wavefunctions.
  ! (H - E)*psi(k+q) + P_c V_HXC(q)*psi0(k)
  !
  ! Inspired by PH/solve_linter.f90
  !
  ! Written by Iurii Timrov (2013)
  !
  USE ions_base,            ONLY : ityp, nat, ntyp=>nsp
  USE cell_base,            ONLY : tpiba2
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_parallel,         ONLY : tg_cgather
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvecs,                ONLY : nls, nlsm, ngms, doublegrid
  USE gvect,                ONLY : nl, nlm, ngm, g, gg
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE klist,                ONLY : nks, xk
  USE lr_variables,         ONLY : evc0, no_hxc
  USE lsda_mod,             ONLY : nspin, current_spin
  USE wvfct,                ONLY : nbnd, npwx, g2kin, et, npw, igk
  USE gvecw,                ONLY : gcutw
  USE io_global,            ONLY : stdout
  USE uspp,                 ONLY : vkb
  USE io_files,             ONLY : iunigk
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE units_ph,             ONLY : lrwfc, iuwfc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE uspp,                 ONLY : okvan
  USE mp_bands,             ONLY : ntask_groups, me_bgrp
  USE spin_orb,             ONLY : domag

  USE qpoint,               ONLY : npwq, igkq, ikks, ikqs, nksq
  USE eqv,                  ONLY : evq, dpsi, dvpsi
  USE control_lr,           ONLY : nbnd_occ
 
  IMPLICIT NONE
  !
  COMPLEX(kind=dp),INTENT(in)  :: evc1(npwx*npol,nbnd,nksq)
  COMPLEX(kind=dp),INTENT(out) :: evc1_new(npwx*npol,nbnd,nksq), &
                                  sevc1_new(npwx*npol,nbnd,nksq) 
  ! output : sevc1_new = S * evc1_new
  LOGICAL, INTENT(in) :: interaction
  LOGICAL :: interaction1
  INTEGER :: i,j,ir, ibnd, ik, ig, ia, ios, ikk, ikq, is, &
             incr, v_siz, ipol
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:), revc(:,:), &
                            & dvrsc(:,:), dvrssc(:,:), &
  ! Change of the Hartree and exchange-correlation (HXC) potential
                            & tg_psic(:,:), tg_dvrssc(:,:)
  ! Task groups: wfct in R-space
  ! Task groups: HXC potential
  INTEGER, ALLOCATABLE :: ibuf(:)
  !
  CALL start_clock('lr_apply')
  !
  IF ( ntask_groups > 1 ) dffts%have_task_groups = .TRUE.
  !
  ALLOCATE (hpsi(npwx*npol,nbnd))
  ALLOCATE (spsi(npwx*npol,nbnd))
  ALLOCATE (revc(dffts%nnr,npol))
  ALLOCATE (dvrsc(dfftp%nnr,nspin_mag))
  ALLOCATE (dvrssc(dffts%nnr,nspin_mag))
  IF (.not. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  hpsi(:,:)   = (0.d0,0.d0)
  spsi(:,:)   = (0.d0,0.d0)
  revc(:,:)   = (0.d0,0.d0)
  dvrsc(:,:)  = (0.0d0,0.0d0)
  dvrssc(:,:) = (0.0d0,0.0d0)
  !
  incr = 1
  !
  IF ( dffts%have_task_groups ) THEN
     !
     v_siz =  dffts%tg_nnr * dffts%nogrp
     !
     ALLOCATE( tg_dvrssc(v_siz,nspin_mag) )
     ALLOCATE( tg_psic(v_siz,npol) )
     !
     incr = dffts%nogrp
     !
  ENDIF
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
     ! Calculation of the charge density response, and symmetrization of it.
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
     ! Calculation of the HXC potential
     ! input:  the change of the charge density (dvrsc)
     ! output: the change of the HXC potential  (dvrsc)
     ! Note: check the implementation of the non-linear core correction.
     !
     CALL dv_of_drho(0, dvrsc, .false.)
     !
     ! Interpolation of the HXC potential from the thick mesh 
     ! to a smoother mesh (if doublegrid=.true.)
     ! dvrsc -> dvrssc
     !
     DO is = 1, nspin_mag
        CALL cinterpolate (dvrsc(1,is), dvrssc(1,is), -1)
     ENDDO
     !
  ENDIF
  !
  ! Now we will calculate two terms:
  ! 1) HXC term : P_c (delta V_HXC(q)) psi(k) 
  ! 2) (H - E) * psi(k+q)
  !
  ! rewind (unit = iunigk)
  ! 
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     !
     ! Determination of npw, igk, and npwq, igkq;
     ! g2kin is used here as a workspace.
     !
     CALL gk_sort( xk(1,ikk), ngm, g, gcutw, npw,  igk,  g2kin )
     CALL gk_sort( xk(1,ikq), ngm, g, gcutw, npwq, igkq, g2kin ) 
     !
     ! Calculate beta-functions vkb at k+q (Kleinman-Bylander projectors)
     ! The vks's are needed for the non-local potential in h_psiq,
     ! and for the ultrasoft term.
     !
     CALL init_us_2 (npwq, igkq, xk(1,ikq), vkb)
     !
!    IF (nksq > 1) THEN
!        read (iunigk, err = 100, iostat = ios) npw, igk
!100     call errore ('lr_apply_liouvillian', 'reading igk', abs (ios) )
!        read (iunigk, err = 200, iostat = ios) npwq, igkq
!200     call errore ('lr_apply_liouvillian', 'reading igkq', abs (ios) )
!    ENDIF
     !
     ! Read unperturbed wavefuctions psi(k) and psi(k+q)
     !
     IF (nksq > 1) THEN 
        CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
        CALL davcio (evq, lrwfc, iuwfc, ikq, - 1)
     ENDIF
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
     IF (interaction1) THEN
        !
        IF ( ntask_groups > 1 ) dffts%have_task_groups = .TRUE.
        !
        ! The potential in dvrssc is distributed across all processors.
        ! We need to redistribute it so that it is completely contained in the
        ! processors of an orbital TASK-GROUP.
        !
        IF ( dffts%have_task_groups ) THEN
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
           IF ( dffts%have_task_groups ) THEN
              !
              ! FFT to R-space
              !
              CALL cft_wave_tg(evc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
              !
              ! Multiply the HXC potential with unperturbed wfct's
              !
              CALL apply_dpot(v_siz, tg_psic, tg_dvrssc, 1)
              !
              ! back-FFT to G-space
              !
              CALL cft_wave_tg(dvpsi, tg_psic, -1, v_siz, ibnd, nbnd_occ(ikk))
              !
           ELSE
              !
              ! FFT to R-space
              !
              CALL cft_wave(evc(1,ibnd), revc, +1)
              !
              ! Multiply the HXC potential with unperturbed wfct's 
              !
              CALL apply_dpot(dffts%nnr, revc, dvrssc, current_spin)
              !
              ! back-FFT to G-space
              !
              CALL cft_wave(dvpsi(1,ibnd), revc, -1)
              !
           ENDIF
           !
        ENDDO
        !
        dffts%have_task_groups = .FALSE.
        !
        ! In the case of US pseudopotentials there is an additional term.
        ! See second part of Eq.(39) in J. Chem. Phys. 127, 164106 (2007)
        !
        IF (okvan) THEN
           !
           ! Compute the integral of the HXC response potential with the Q function.
           ! Input : dvrsc = V_HXC(r)
           ! Output: int3 = \int V_HXC(r) * Q^*_nm(r) dr 
           ! See Eq.(B22) in Ref. A. Dal Corso, PRB 64, 235118 (2001)
           !
           CALL newdq(dvrsc, 1)
           !
           CALL adddvscf(1, ik) 
           !
        ENDIF
        !
        ! Ortogonalize dvpsi to valence states.
        ! Apply -P_c, and then change the sign, because we need +P_c.
        !
        CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq)
        dvpsi = -dvpsi
        !
     ENDIF
     !
     ! 2) (H - E) * psi(k+q)
     !
     ! Compute the kinetic energy g2kin: (k+q+G)^2
     !
     DO ig = 1, npwq
        g2kin (ig) = ( (xk (1,ikq) + g (1, igkq(ig)) ) **2 + &
                       (xk (2,ikq) + g (2, igkq(ig)) ) **2 + &
                       (xk (3,ikq) + g (3, igkq(ig)) ) **2 ) * tpiba2
     ENDDO
     !
     ! Apply the operator ( H - \epsilon S + alpha_pv P_v) to evc1
     ! where alpha_pv = 0
     !
     !call ch_psi_all (npwq, evc1(1,1,ik), sevc1_new(1,1,ik), et(1,ikk), ik, nbnd_occ(ikk)) 
     !
     ! Compute H*psi
     !
     IF (noncolin) THEN
        IF (.NOT. ALLOCATED(psic_nc)) ALLOCATE(psic_nc(dfftp%nnr,npol))
     !ELSE
     !  IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
     ENDIF
     !
     IF (ntask_groups > 1) dffts%have_task_groups = .TRUE. 
     !
     IF (dffts%have_task_groups) THEN
        !
        ! With task groups we use the H*psi routine of PW parallelized on task groups
        ! (see PH/ch_psi_all.f90)
        !
        ALLOCATE(ibuf(npwx))
        ibuf = igk
        igk = igkq
        CALL h_psi (npwx, npwq, nbnd_occ(ikk), evc1(:,:,ik), hpsi)
        CALL s_psi (npwx, npwq, nbnd_occ(ikk), evc1(:,:,ik), spsi)
        igk = ibuf
        DEALLOCATE(ibuf)
        !
     ELSE
        !
        CALL h_psiq (npwx, npwq, nbnd_occ(ikk), evc1(:,:,ik), hpsi, spsi)
        !
     ENDIF
     !
     dffts%have_task_groups = .FALSE.
     !
     IF (noncolin) THEN
        IF (ALLOCATED(psic_nc)) DEALLOCATE(psic_nc)
     !ELSE
     !  IF (ALLOCATED(psic)) DEALLOCATE(psic)
     ENDIF
     !
     ! Subtract the eigenevalues H*psi(k+q) - et*psi(k+q)
     !
     DO ibnd = 1, nbnd_occ(ikk)
        DO ig = 1, npwq
           sevc1_new(ig,ibnd,ik) = hpsi(ig,ibnd) - &
                    &  cmplx(et(ibnd,ikk),0.0d0,dp) * spsi(ig,ibnd)
        ENDDO
     ENDDO
     !
     IF (noncolin) THEN
        DO ibnd = 1, nbnd_occ(ikk)
           DO ig = 1, npwq
              sevc1_new(ig+npwx,ibnd,ik) = hpsi(ig+npwx,ibnd) - &
                     cmplx(et(ibnd,ikk),0.0d0,dp) * spsi(ig+npwx,ibnd)
           ENDDO
        ENDDO
     ENDIF
     !
     ! 3) Sum up the two terms : (H - E)*psi(k+q) + HXC
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
                 sevc1_new(ig+npwx,ibnd,ik) = sevc1_new(ig+npwx,ibnd,ik) + dvpsi(ig+npwx,ibnd)
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
     CALL sm1_psi(.FALSE.,ik, npwx, npwq, nbnd_occ(ikk), sevc1_new(1,1,ik), evc1_new(1,1,ik))
     !
  ENDDO ! loop on ik
  !
  DEALLOCATE (dvrsc)
  DEALLOCATE (dvrssc) 
  DEALLOCATE (hpsi)
  DEALLOCATE (spsi)
  DEALLOCATE (revc)
  IF (ALLOCATED(psic)) DEALLOCATE(psic)
  !
  IF ( ntask_groups > 1) dffts%have_task_groups = .TRUE.
  !
  IF ( dffts%have_task_groups ) THEN
     DEALLOCATE( tg_dvrssc )
     DEALLOCATE( tg_psic )
  ENDIF
  !
  dffts%have_task_groups = .FALSE.
  !
  IF (interaction1)      CALL stop_clock('lr_apply_int')
  IF (.NOT.interaction1) CALL stop_clock('lr_apply_no')
  CALL stop_clock('lr_apply')
  !
  RETURN
  !
END SUBROUTINE lr_apply_liouvillian_eels
