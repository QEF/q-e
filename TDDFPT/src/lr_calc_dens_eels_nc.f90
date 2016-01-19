!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE lr_calc_dens_eels_nc (drhoscf, dpsi)
  !-------------------------------------------------------------------------
  !
  ! Calculates response charge-density from linear-response
  ! orbitals and ground-state orbitals (noncollinear case).
  ! See Eq.(36) in B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
  ! Inspired by PH/incdrhoscf_nc.f90 and PH/solve_linter.f90
  !
  ! It does:
  ! 1) Calculates the response charge-density
  ! 2) Calculates the contribution to the charge-density due to US PP's
  ! 3) Sums up the normal and ultrasoft terms
  ! 4) Symmetrizes the total response charge-density
  !
  ! TODO: Make a call to the routine PH/incdrhoscf_nc.f90 and remove
  ! all the duplicated stuff here.
  !
  ! Written by Iurii Timrov (2013)
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE ions_base,             ONLY : nat
  USE gvecs,                 ONLY : nls
  USE gvect,                 ONLY : ngm, g
  USE fft_base,              ONLY : dffts, dfftp
  USE fft_interfaces,        ONLY : invfft
  USE klist,                 ONLY : xk, wk
  USE lr_variables,          ONLY : evc0, lr_periodic
  USE wvfct,                 ONLY : nbnd,wg,npwx,npw,igk,g2kin
  USE gvecw,                 ONLY : gcutw
  USE qpoint,                ONLY : npwq, igkq, nksq, ikks, ikqs
  USE control_lr,            ONLY : nbnd_occ
  USE units_ph,              ONLY : lrwfc, iuwfc
  USE wavefunctions_module,  ONLY : evc
  USE eqv,                   ONLY : evq
  USE noncollin_module,      ONLY : npol, nspin_mag
  USE spin_orb,              ONLY : domag
  USE uspp_param,            ONLY : nhm
  USE uspp,                  ONLY : okvan, vkb
  USE lsda_mod,              ONLY : nspin
  USE mp_global,             ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp_bands,              ONLY : me_bgrp, ntask_groups
  USE mp,                    ONLY : mp_sum
  USE io_global,             ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: dpsi(npwx*npol,nbnd,nksq)
  ! input: the perturbed wavefunctions
  COMPLEX(DP), INTENT(out) :: drhoscf(dfftp%nnr,nspin_mag)
  ! input/output: the accumulated change of the charge density
  COMPLEX(DP), ALLOCATABLE :: drhoscfh(:,:), dbecsum(:,:,:), dbecsum_nc(:,:,:,:)
  ! the accumulated change to the charge density on the smooth mesh
  ! the ultrasoft term
  INTEGER :: ir, ik, ikk, ikq, ibnd, ig, is, incr, v_siz, idx, ioff, ipol
  REAL(DP) :: wgt
  ! the effective weight of the k point
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), dpsic(:,:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  COMPLEX(DP), ALLOCATABLE :: tg_psi (:,:), tg_dpsi (:,:), tg_drho(:,:)
  ! arrays for task groups parallelization
  !
  CALL start_clock ('lr_calc_dens')
  !
  IF (ntask_groups > 1 ) dffts%have_task_groups = .TRUE.
  !
  ALLOCATE (dpsic(dffts%nnr,npol))
  ALLOCATE (psi(dffts%nnr,npol))
  ALLOCATE (drhoscfh(dffts%nnr,nspin_mag))
  dpsic(:,:)    = (0.0d0, 0.0d0)
  psi(:,:)      = (0.0d0, 0.0d0)
  drhoscfh(:,:) = (0.0d0, 0.0d0)
  !
  ! Step in the loop for bands
  !
  incr = 1
  !
  IF (dffts%have_task_groups) THEN
     !
     v_siz = dffts%tg_nnr * dffts%nogrp
     !
     ALLOCATE ( tg_psi(v_siz,npol) )
     ALLOCATE ( tg_dpsi(v_siz,npol) )
     ALLOCATE ( tg_drho(v_siz,nspin_mag) )
     !
     incr  = dffts%nogrp
     !
  ENDIF
  !
  IF (okvan) THEN
     ALLOCATE (dbecsum((nhm*(nhm + 1))/2,nat,nspin_mag))
     ALLOCATE (dbecsum_nc(nhm,nhm,nat,nspin))
     dbecsum(:,:,:)      = (0.0d0, 0.0d0)
     dbecsum_nc(:,:,:,:) = (0.0d0, 0.0d0)
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
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
     ! g2kin is used here as work space.
     !
     CALL gk_sort( xk(1,ikk), ngm, g, gcutw, npw,  igk,  g2kin )
     CALL gk_sort( xk(1,ikq), ngm, g, gcutw, npwq, igkq, g2kin )
     !
     ! Read the unperturbed wavefuctions evc(k)
     !
     IF (lr_periodic) THEN
        evc(:,:) = evc0(:,:,ik)
     ELSE
        IF (nksq > 1) CALL davcio (evc, lrwfc, iuwfc, ikk, - 1)
     ENDIF
     !
     DO ibnd = 1, nbnd_occ(ikk), incr
        !
        ! The weight
        !
        wgt = 2.0d0*wk(ikk)/omega
        !
        IF (dffts%have_task_groups) THEN
           !
           ! Task groups case
           !
           tg_drho = (0.0_DP, 0.0_DP)
           tg_psi  = (0.0_DP, 0.0_DP)
           tg_dpsi = (0.0_DP, 0.0_DP)
           !
           ioff   = 0
           !
           DO idx = 1, dffts%nogrp
              !
              ! ... dffts%nogrp ffts at the same time. We prepare both
              ! evc (at k) and dpsi (at k+q)
              !
              IF ( idx + ibnd - 1 <= nbnd_occ(ikk) ) THEN
                 !
                 ! Unperturbed wavefunctions
                 !
                 DO ig = 1, npw
                    tg_psi( nls( igk( ig ) ) + ioff, 1 ) = evc( ig, idx+ibnd-1 )
                    tg_psi( nls( igk( ig ) ) + ioff, 2 ) = evc( npwx+ig, idx+ibnd-1)
                 ENDDO
                 !
                 ! Perturbed wavefunctions
                 !
                 DO ig = 1, npwq
                    tg_dpsi( nls( igkq( ig ) ) + ioff, 1 ) = dpsi( ig, idx+ibnd-1, ik )
                    tg_dpsi( nls( igkq( ig ) ) + ioff, 2 ) = dpsi( npwx+ig, idx+ibnd-1, ik )
                 ENDDO
                 !
              ENDIF
              !
              ioff = ioff + dffts%tg_nnr
              !
           ENDDO
           !
           CALL invfft ('Wave', tg_psi(:,1), dffts)
           CALL invfft ('Wave', tg_psi(:,2), dffts)
           CALL invfft ('Wave', tg_dpsi(:,1), dffts)
           CALL invfft ('Wave', tg_dpsi(:,2), dffts)
           !
           ! Calculation of the density
           !
           DO ir = 1, dffts%tg_npp( me_bgrp + 1 ) * dffts%nr1x * dffts%nr2x
              tg_drho (ir,1) = tg_drho (ir,1) + wgt * (CONJG(tg_psi (ir,1) )*tg_dpsi (ir,1) &
                                                     & + CONJG(tg_psi (ir,2) )*tg_dpsi (ir,2) )
           ENDDO
           !
           IF (domag) THEN
              DO ir = 1, dffts%tg_npp( me_bgrp + 1 ) * dffts%nr1x * dffts%nr2x
                   tg_drho(ir,2) = tg_drho(ir,2) + wgt *(CONJG(tg_psi(ir,1))*tg_dpsi(ir,2) &
                                                     & + CONJG(tg_psi(ir,2))*tg_dpsi(ir,1) )
                   tg_drho(ir,3) = tg_drho(ir,3) + wgt *(CONJG(tg_psi(ir,1))*tg_dpsi(ir,2) &
                                                     & - CONJG(tg_psi(ir,2))*tg_dpsi(ir,1) )*(0.d0,-1.d0)
                   tg_drho(ir,4) = tg_drho(ir,4) + wgt *(CONJG(tg_psi(ir,1))*tg_dpsi(ir,1) &
                                                     & - CONJG(tg_psi(ir,2))*tg_dpsi(ir,2) )
              ENDDO
           ENDIF
           !
           ! Reduce the group charge (equivalent to sum over the bands of the
           ! orbital group)
           !
           CALL mp_sum( tg_drho, gid = dffts%ogrp_comm )
           !
           ioff = 0
           !
           DO idx = 1, dffts%nogrp
              !
              IF ( me_bgrp == dffts%nolist( idx ) ) EXIT
              !
              ioff = ioff + dffts%nr1x * dffts%nr2x * &
                                 &  dffts%npp( dffts%nolist( idx ) + 1 )
           ENDDO
           !
           ! Copy the charge back to the proper processor location
           !
           DO ipol=1,nspin_mag
              DO ir = 1, dffts%nnr
                 drhoscfh(ir,ipol) = drhoscfh(ir,ipol) + tg_drho(ir+ioff,ipol)
              ENDDO
           ENDDO
           !
        ELSE
           !        
           ! Normal case: no task groups
           !
           ! FFT to R-space of the unperturbed wfct's evc
           !
           psi = (0.d0, 0.d0)
           DO ig = 1, npw
              psi(nls(igk(ig)),1) = evc(ig,ibnd)
              psi(nls(igk(ig)),2) = evc(ig+npwx,ibnd)
           ENDDO
           CALL invfft ('Wave', psi(:,1), dffts)
           CALL invfft ('Wave', psi(:,2), dffts)
           !
           ! FFT to R-space of the perturbed wfct's dpsi
           !
           dpsic = (0.d0, 0.d0)
           DO ig = 1, npwq
              dpsic(nls(igkq(ig)),1) = dpsi(ig,ibnd,ik)
              dpsic(nls(igkq(ig)),2) = dpsi(ig+npwx,ibnd,ik)
           ENDDO
           CALL invfft ('Wave', dpsic(:,1), dffts)
           CALL invfft ('Wave', dpsic(:,2), dffts)
           !
           ! Calculation of the response charge density
           !
           DO ir = 1, dffts%nnr
              drhoscfh(ir,1) = drhoscfh(ir,1) + wgt * (CONJG(psi(ir,1))*dpsic(ir,1) + &
                                                    &  CONJG(psi(ir,2))*dpsic(ir,2) )
           ENDDO
           !
           IF (domag) THEN
              DO ir = 1, dffts%nnr
                 drhoscfh(ir,2) = drhoscfh(ir,2) + wgt * (CONJG(psi(ir,1))*dpsic(ir,2) + &
                                                      &   CONJG(psi(ir,2))*dpsic(ir,1) )
                 drhoscfh(ir,3) = drhoscfh(ir,3) + wgt * (CONJG(psi(ir,1))*dpsic(ir,2) - &
                                                      &   CONJG(psi(ir,2))*dpsic(ir,1) ) * (0.d0,-1.d0)
                 drhoscfh(ir,4) = drhoscfh(ir,4) + wgt * (CONJG(psi(ir,1))*dpsic(ir,1)- &
                                                      &   CONJG(psi(ir,2))*dpsic(ir,2) )
              ENDDO
           ENDIF
           !
        ENDIF
        !
     ENDDO ! loop on bands
     !
     ! Ultrasoft contribution.
     ! Distribute the calculation of dbecsum_nc across processors.
     !
     IF (okvan) THEN
        !
        ! Calculate beta-functions vkb at point k+q
        ! 
        CALL init_us_2 (npwq, igkq, xk(1,ikq), vkb)
        !
        ! Calculate dbecsum_nc = <evc|vkb><vkb|evc1> 
        !
        CALL addusdbec_nc (ik, wk(ikk), dpsi(1,1,ik), dbecsum_nc(1,1,1,1))
        !
     ENDIF
     !
  ENDDO ! loop on k points
  !
  DEALLOCATE (psi)
  DEALLOCATE (dpsic)
  !
  IF (dffts%have_task_groups) THEN
     DEALLOCATE(tg_psi)
     DEALLOCATE(tg_dpsi)
     DEALLOCATE(tg_drho)
  ENDIF
  !
  dffts%have_task_groups = .FALSE. 
  !
  ! Interpolate the charge density from the smooth mesh
  ! to a thicker mesh (if doublegrid=.true.)
  ! drhoscfh -> drhoscf
  !
  DO is = 1, nspin_mag
     CALL cinterpolate(drhoscf(1,is), drhoscfh(1,is), 1)
  ENDDO
  !
  IF (okvan) THEN
     !
     ! The calculation of dbecsum is distributed across processors (see addusdbec_nc).
     ! Sum over processors the contributions coming from each slice of bands.
     !
#ifdef __MPI
     CALL mp_sum (dbecsum_nc, intra_bgrp_comm)
#endif
     !
     ! In the noncolinear, spin-orbit case rotate dbecsum
     !
     CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 1)
     !
     ! Calculate the total response charge-density
     ! (sum up the normal and ultrasoft terms)
     !
     CALL lr_addusddens (drhoscf, dbecsum)
     !
  ENDIF
  !
  ! Reduce the charge density response across pools.
  !
#ifdef __MPI
  CALL mp_sum (drhoscf, inter_pool_comm)
#endif
  !
  ! Symmetrization of the charge density response
  ! wrt the small group of q.
  !
#ifdef __MPI
  CALL lr_psym_eels(drhoscf)
#else
  CALL lr_sym_eels(drhoscf)
#endif
  !
  DEALLOCATE (drhoscfh)
  !
  IF (okvan) THEN
     DEALLOCATE (dbecsum)
     DEALLOCATE (dbecsum_nc)
  ENDIF
  !
  CALL stop_clock ('lr_calc_dens')
  !
  RETURN
  !
END SUBROUTINE lr_calc_dens_eels_nc
