!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE lr_calc_dens_eels (drhoscf, dpsi)
  !--------------------------------------------------------------------------
  !
  ! Calculates response charge-density from linear-response
  ! orbitals and ground-state orbitals.
  ! See Eq.(36) in B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007)
  ! Inspired by PH/incdrhoscf.f90 and PH/solve_linter.f90
  !
  ! It does:
  ! 1) Calculates the response charge-density
  ! 2) Calculates the contribution to the charge-density due to US PP's
  ! 3) Sums up the normal and ultrasoft terms
  ! 4) Symmetrizes the total response charge-density
  !
  ! TODO: Make a call to the routine PH/incdrhoscf.f90 and remove
  ! all the duplicated stuff here.
  !
  ! Written by Iurii Timrov (2013)
  !
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE cell_base,             ONLY : omega
  USE ions_base,             ONLY : nat
  USE gvecs,                 ONLY : nls
  USE gvect,                 ONLY : ngm, g
  USE fft_base,              ONLY : dfftp, dffts
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
  USE uspp_param,            ONLY : nhm 
  USE uspp,                  ONLY : okvan, vkb
  USE mp_global,             ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp_bands,              ONLY : me_bgrp, ntask_groups
  USE mp,                    ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: dpsi(npwx,nbnd,nksq)
  ! input: the perturbed wavefunctions
  COMPLEX(DP), INTENT(out) :: drhoscf(dfftp%nnr)
  ! input/output: the accumulated change to the charge density on the dense mesh
  COMPLEX(DP), ALLOCATABLE :: drhoscfh(:), dbecsum(:,:)
  ! the accumulated change to the charge density on the smooth mesh
  ! the ultrasoft term
  INTEGER :: ir, ik, ikk, ikq, ibnd, ig, incr, v_siz, idx, ioff
  REAL(DP) :: wgt  
  ! the effective weight of the k point
  COMPLEX(DP), ALLOCATABLE :: psi(:), dpsic(:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_dpsi(:), tg_drho(:)
  ! arrays for task groups parallelization
  !
  CALL start_clock('lr_calc_dens')
  !
  IF (ntask_groups > 1) dffts%have_task_groups = .TRUE.
  !
  ALLOCATE ( dpsic(dffts%nnr) )
  ALLOCATE ( psi(dffts%nnr) )
  ALLOCATE ( drhoscfh(dffts%nnr) )
  dpsic(:) = (0.0d0, 0.0d0) 
  psi(:) = (0.0d0, 0.0d0) 
  drhoscfh(:) = (0.0d0, 0.0d0)
  !
  ! Step in the loop for bands
  !
  incr = 1 
  !
  IF (dffts%have_task_groups) THEN
     !
     v_siz = dffts%tg_nnr * dffts%nogrp
     !
     ALLOCATE ( tg_psi(v_siz) )
     ALLOCATE ( tg_dpsi(v_siz) )
     ALLOCATE ( tg_drho(v_siz) )
     !
     incr  = dffts%nogrp
     !
  ENDIF
  !
  IF (okvan) THEN
     ALLOCATE (dbecsum(nhm*(nhm+1)/2,nat)) 
     dbecsum(:,:) = (0.0d0, 0.0d0)
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
     ! Read the unperturbed wavefuctions evc at k
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
           ioff = 0
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
                    tg_psi( nls( igk( ig ) ) + ioff ) = evc( ig, idx+ibnd-1 )
                 ENDDO
                 !
                 ! Perturbed wavefunctions
                 !
                 DO ig = 1, npwq
                    tg_dpsi( nls( igkq( ig ) ) + ioff ) = dpsi( ig, idx+ibnd-1, ik )
                 ENDDO
                 !
              ENDIF
              !
              ioff = ioff + dffts%tg_nnr
              !
           ENDDO
           !
           CALL invfft ('Wave', tg_psi, dffts)
           CALL invfft ('Wave', tg_dpsi, dffts)
           !
           ! Calculation of the density
           !
           DO ir = 1, dffts%tg_npp( me_bgrp + 1 ) * dffts%nr1x * dffts%nr2x
              tg_drho (ir) = tg_drho (ir) + wgt * CONJG(tg_psi (ir) ) * tg_dpsi (ir)
           ENDDO
           !
           ! Reduce the group charge (equivalent to sum over bands of 
           ! orbital group)
           !
           CALL mp_sum( tg_drho, gid = dffts%ogrp_comm )
           !
           ioff = 0
           !
           DO idx = 1, dffts%nogrp
              !
              IF ( me_bgrp == dffts%nolist(idx) ) EXIT
              !
              ioff = ioff + dffts%nr1x * dffts%nr2x * &
                                   dffts%npp( dffts%nolist(idx) + 1 )
              !
           ENDDO
           !
           ! Copy the charge back to the proper processor location
           !
           DO ir = 1, dffts%nnr
              drhoscfh(ir) = drhoscfh(ir) + tg_drho(ir+ioff)
           ENDDO
           !
        ELSE
           !
           ! Normal case: no task groups
           !
           ! FFT to R-space of the unperturbed wfct's evc
           !
           psi(:) = (0.d0, 0.d0)
           DO ig = 1, npw
              psi(nls(igk(ig))) = evc(ig,ibnd)
           ENDDO
           CALL invfft ('Wave', psi, dffts)
           !
           ! FFT to R-space of the perturbed wfct's dpsi
           !
           dpsic(:) = (0.d0, 0.d0)
           DO ig = 1, npwq
              dpsic(nls(igkq(ig))) = dpsi(ig,ibnd,ik)
           ENDDO
           CALL invfft ('Wave', dpsic, dffts)
           !
           ! Calculation of the response charge-density 
           !
           DO ir = 1, dffts%nnr
              drhoscfh(ir) = drhoscfh(ir) + wgt * CONJG(psi(ir)) * dpsic(ir)
           ENDDO
           !
        ENDIF
        !
     ENDDO ! loop on bands
     !
     ! Ultrasoft contribution.
     ! Distribute the calculation of dbecsum across processors.
     !
     IF (okvan) THEN
        !
        ! Calculate beta-functions vkb at point k+q
        !
        CALL init_us_2 (npwq, igkq, xk(1,ikq), vkb)
        !
        ! Calculate dbecsum = <evc|vkb><vkb|evc1>
        !
        CALL addusdbec (ik, wk(ikk), dpsi(1,1,ik), dbecsum(1,1))
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
  CALL cinterpolate(drhoscf, drhoscfh, 1)
  !
  IF (okvan) THEN
     !
     ! The calculation of dbecsum is distributed across processors (see addusdbec).
     ! Sum over processors the contributions coming from each slice of bands.
     !
#ifdef __MPI
     CALL mp_sum (dbecsum, intra_bgrp_comm)
#endif
     !
     ! Calculate the total charge density response
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
  IF (okvan) DEALLOCATE (dbecsum)
  !
  CALL stop_clock('lr_calc_dens')
  !
  RETURN
  !
END SUBROUTINE lr_calc_dens_eels
