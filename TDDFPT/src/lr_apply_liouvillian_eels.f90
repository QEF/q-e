!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
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
  ! (H - E)*psi(k+q) + V_HXC(q)*psi0(k)
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
  USE lr_variables,         ONLY : evc0, no_hxc, lr_periodic
  USE lsda_mod,             ONLY : nspin, current_spin
  USE wvfct,                ONLY : nbnd, npwx, g2kin, et, npw, igk
  USE gvecw,                ONLY : gcutw
  USE io_global,            ONLY : stdout
  USE uspp,                 ONLY : vkb
  USE qpoint,               ONLY : npwq, igkq, ikks, ikqs, nksq
  USE eqv,                  ONLY : evq, dpsi, dvpsi
  USE io_files,             ONLY : iunigk
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE units_ph,             ONLY : lrwfc, iuwfc
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE control_ph,           ONLY : nbnd_occ
  USE uspp,                 ONLY : okvan
  USE nlcc_ph,              ONLY : nlcc_any
  USE iso_c_binding,        ONLY : c_int
  USE mp_bands,             ONLY : ntask_groups, me_bgrp
  USE spin_orb,             ONLY : domag
 
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
  INTEGER(kind=c_int) :: kilobytes
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
  ! Memory usage
  !
  !CALL memstat( kilobytes )
  !IF ( kilobytes > 0 ) WRITE(stdout,'(5X,"lr_apply_liouvillian_eels, & 
  !         & per-process dynamical memory:",f7.1,"Mb")' ) kilobytes/1000.0
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
     !CALL dv_of_drho(0, dvrsc, .false.)
     CALL lr_dv_of_drho_eels(dvrsc)
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
     IF (lr_periodic) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = ikks(ik)
        ikq = ikqs(ik)
     ENDIF
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
     IF (lr_periodic) THEN
        evc(:,:) = evc0(:,:,ik)
        evq(:,:) = evc0(:,:,ik)
     ELSE
        IF (nksq > 1) THEN 
           CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
           CALL davcio (evq, lrwfc, iuwfc, ikq, - 1)
        ENDIF
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
CONTAINS

!-------------------------------------------------------------------------
SUBROUTINE lr_dv_of_drho_eels (dvscf)
  !-----------------------------------------------------------------------
  !
  ! This subroutine computes the change of the self consistent potential
  ! (Hartree and XC) due to the perturbation.
  ! Inspired by PH/dv_of_drho.f90
  !
  ! Written by I. Timrov, Feb 2015
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : e2, fpi
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE gvect,             ONLY : nl, ngm, g, nlm
  USE cell_base,         ONLY : alat, tpiba2, omega
  USE noncollin_module,  ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE funct,             ONLY : dft_is_gradient
  USE scf,               ONLY : rho, rho_core
  USE eqv,               ONLY : dmuxc
  USE nlcc_ph,           ONLY : nlcc_any
  USE qpoint,            ONLY : xq
  USE gc_ph,             ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  USE control_ph,        ONLY : lrpa
  USE control_flags,     ONLY : gamma_only
  USE lr_variables,      ONLY : clfe  !eps

  IMPLICIT NONE

  COMPLEX(DP), INTENT(inout):: dvscf(dfftp%nnr, nspin_mag)
  ! input:  the change of the charge density
  ! output: the change of the HXC potential
  INTEGER :: ir, is, is1, ig, ngm_, gstart
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors
  ! numver of G vectors to be considered
  ! initial starting G vector
  REAL(DP) :: qg2, fac
  ! the modulus of (q+G)^2
  ! the structure factor

  complex(DP), allocatable :: dvaux(:,:), drhoc(:)
  !  the change of the core charge
  complex(DP), allocatable :: dvhart(:,:) 
  complex(DP), allocatable :: dvaux_mt(:), rgtot(:)
  ! auxiliary array for Martyna-Tuckerman correction in TDDFPT
  ! total response density  
  real(DP) :: eh_corr
  ! Correction to response Hartree energy due to Martyna-Tuckerman correction 
  ! (only TDDFT). Not used.

  CALL start_clock ('lr_dv_of_drho_eels')
  !
  ALLOCATE (dvaux(dfftp%nnr,nspin_mag))
  dvaux(:,:) = (0.d0, 0.d0)
  !
  ! 1) The exchange-correlation contribution is computed in real space.
  !
  IF (lrpa) goto 111
  !
  DO is = 1, nspin_mag
     DO is1 = 1, nspin_mag
        DO ir = 1, dfftp%nnr
           dvaux(ir,is) = dvaux(ir,is) + dmuxc(ir,is,is1) * dvscf(ir,is1)
        ENDDO
     ENDDO
  ENDDO
  !
  ! Add a gradient correction to XC.
  ! If nlcc=.true. we need to add here its contribution.
  ! grho contains the core charge.
  !
  fac = 1.d0 / DBLE (nspin_lsda)
  !
  IF (nlcc_any) THEN
     DO is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core (:)
     ENDDO
  ENDIF
  !
  IF (dft_is_gradient()) CALL dgradcorr &
       (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
       dvscf, dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, alat, dvaux)
  !
  IF (nlcc_any) THEN
     DO is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core (:)
     ENDDO
  ENDIF
  !
111 CONTINUE
  !
  ! Copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  IF (nspin_mag == 2) THEN
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
  ENDIF
  !
  CALL fwfft ('Dense', dvscf(:,1), dfftp)
  !
  ! 2) The Hartree contribution is computed in reciprocal space.
  !
  ! An extension to gamma_ionly case can be done from PH/dv_of_drho.f90
  !
  IF (gamma_only)  CALL errore( 'lr_dv_of_drho_eels', 'gamma_only is not supported', 1 )  
  !
  !IF (eps) THEN
  !   ! No G=0 term
  !   gstart = 2 
  !ELSE
  !   ! With G=0 term
  !   gstart = 1
  !ENDIF 
  !
  gstart = 1
  !
  IF (clfe) THEN
     ! All G vectors are considered
     ngm_ = ngm
  ELSE
     ! Only G=0 is considered
     ngm_ = 1
  ENDIF
  !
  !IF (eps .AND. .NOT.clfe) CALL errore( 'lr_dv_of_drho_eels', &
  !             & 'No Hartree term, because eps=.true. and clfe=.false.', 1 )
  ! 
  DO is = 1, nspin_lsda
     !
     ! FFT from R-space to G-space.
     !
     CALL fwfft ('Dense', dvaux (:, is), dfftp)
     ! 
     DO ig = gstart, ngm_
        !
        qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
        !
        ! Hartree term: 4*pi*e2/|q+G|^2 * n'(G)
        !  
        IF (qg2 > 1.d-8) THEN
             dvaux(nl(ig),is) = dvaux(nl(ig),is) + &
                                e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
        ENDIF
        !
     ENDDO
     !
     ! back-FFT from G-space to R-space
     !
     CALL invfft ('Dense', dvaux (:, is), dfftp)
     !
  ENDDO
  !
  ! At the end the two contributes are added
  !  
  dvscf (:,:) = dvaux (:,:)
  !
  DEALLOCATE(dvaux)
  !
  CALL stop_clock ('lr_dv_of_drho_eels')
  !
  RETURN
  !
END SUBROUTINE lr_dv_of_drho_eels
  !
END SUBROUTINE lr_apply_liouvillian_eels
