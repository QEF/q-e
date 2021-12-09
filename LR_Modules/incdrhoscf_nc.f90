!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhoscf_nc (drhoscf, weight, ik, dbecsum, dpsi, rsign)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : invfft
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, domag, nspin_mag
  USE uspp_param,           ONLY : nhm
  USE wvfct,                ONLY : npwx, nbnd
  USE wavefunctions, ONLY : evc
  USE klist,                ONLY : ngk,igk_k
  USE qpoint,               ONLY : ikks, ikqs
  USE control_lr,           ONLY : nbnd_occ
  USE qpoint_aux,           ONLY : becpt
  USE lrus,                 ONLY : becp1
  USE mp_bands,             ONLY : me_bgrp, inter_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum
  USE fft_helper_subroutines

  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER, INTENT(IN) :: ik
  ! input: the k point
  REAL(DP), INTENT(IN) :: weight
  REAL(DP), INTENT(IN) :: rsign
  ! input: the weight of the k point
  ! the sign in front of the response of the magnetization density
  COMPLEX(DP), INTENT(IN) :: dpsi(npwx*npol,nbnd)
  ! input: the perturbed wfcs at the given k point
  COMPLEX(DP), INTENT(INOUT) :: drhoscf (dfftp%nnr,nspin_mag), dbecsum (nhm,nhm,nat,nspin)
  ! input/output: the accumulated change of the charge density and dbecsum
  !
  !   here the local variable
  !
  REAL(DP) :: wgt
  ! the effective weight of the k point
  !
  COMPLEX(DP), ALLOCATABLE :: psi (:,:), dpsic (:,:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  !
  COMPLEX(DP), ALLOCATABLE :: tg_psi (:,:), tg_dpsi (:,:), tg_drho(:,:)
  !
  INTEGER :: npw, npwq, ikk, ikq
  INTEGER :: ibnd, jbnd, ir, ir3, ig, incr, v_siz, idx, ioff, ioff_tg, nxyp
  INTEGER :: ntgrp, right_inc
  ! counters
  !
  !
  CALL start_clock ('incdrhoscf')
  !
  ALLOCATE (dpsic(dffts%nnr, npol))
  ALLOCATE (psi  (dffts%nnr, npol))
  !
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  incr = 1
  !
  IF (dffts%has_task_groups) THEN
     !
     v_siz = dffts%nnr_tg
     !
     ALLOCATE( tg_psi( v_siz, npol ) )
     ALLOCATE( tg_dpsi( v_siz, npol ) )
     ALLOCATE( tg_drho( v_siz, nspin_mag ) )
     !
     incr  = fftx_ntgrp(dffts)
     !
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  do ibnd = 1, nbnd_occ(ikk), incr

     IF (dffts%has_task_groups) THEN
        !
        tg_drho=(0.0_DP, 0.0_DP)
        tg_psi=(0.0_DP, 0.0_DP)
        tg_dpsi=(0.0_DP, 0.0_DP)
        !
        ioff   = 0
        CALL tg_get_recip_inc( dffts, right_inc )
        ntgrp = fftx_ntgrp( dffts )
        !
        DO idx = 1, ntgrp
           !
           ! ... dtgs%nogrp ffts at the same time. We prepare both
           ! evc (at k) and dpsi (at k+q)
           !
           IF( idx + ibnd - 1 <= nbnd_occ(ikk) ) THEN
              !
              DO ig = 1, npw
                 tg_psi( dffts%nl( igk_k( ig,ikk ) ) + ioff, 1 ) = evc( ig, idx+ibnd-1 )
                 tg_psi( dffts%nl( igk_k( ig,ikk ) ) + ioff, 2 ) = evc( npwx+ig, idx+ibnd-1 )
              END DO
              DO ig = 1, npwq
                 tg_dpsi( dffts%nl( igk_k( ig,ikq ) ) + ioff, 1 ) = dpsi( ig, idx+ibnd-1 )
                 tg_dpsi( dffts%nl( igk_k( ig,ikq ) ) + ioff, 2 ) = dpsi( npwx+ig, idx+ibnd-1 )
              END DO
              !
           END IF
           !
           ioff = ioff + right_inc
           !
        END DO
        CALL invfft ('tgWave', tg_psi(:,1), dffts)
        CALL invfft ('tgWave', tg_psi(:,2), dffts)
        CALL invfft ('tgWave', tg_dpsi(:,1), dffts)
        CALL invfft ('tgWave', tg_dpsi(:,2), dffts)

        do ir = 1, dffts%nr1x * dffts%nr2x * dffts%my_nr3p
           tg_drho (ir,1) = tg_drho (ir,1) + wgt * (CONJG(tg_psi (ir,1) )*tg_dpsi (ir,1) &
                                                  + CONJG(tg_psi (ir,2) )*tg_dpsi (ir,2) )
        enddo
        IF (domag) THEN
           do ir = 1, dffts%nr1x * dffts%nr2x * dffts%my_nr3p
              tg_drho(ir,2)= tg_drho(ir,2) + (rsign) *wgt * (CONJG(tg_psi(ir,1))*tg_dpsi(ir,2) &
                                                  + CONJG(tg_psi(ir,2))*tg_dpsi(ir,1) )
              tg_drho(ir,3)= tg_drho(ir,3) + (rsign) *wgt * (CONJG(tg_psi(ir,1))*tg_dpsi(ir,2) &
                                                  - CONJG(tg_psi(ir,2))*tg_dpsi(ir,1) ) * (0.d0,-1.d0)
              tg_drho(ir,4)= tg_drho(ir,4) + (rsign) *wgt * (CONJG(tg_psi(ir,1))*tg_dpsi(ir,1) &
                                                  - CONJG(tg_psi(ir,2))*tg_dpsi(ir,2) )
           enddo
        ENDIF
        !
        ! reduce the group charge (equivalent to sum over the bands of the
        ! orbital group)
        !
        CALL tg_reduce_rho( drhoscf, tg_drho, dffts )
        !
     ELSE
        !
        ! Normal case: no task groups
        !
        ! FFT to R-space of the unperturbed wfct's evc
        !
        psi = (0.d0, 0.d0)
        do ig = 1, npw
           psi (dffts%nl (igk_k(ig,ikk) ), 1) = evc (ig, ibnd)
           psi (dffts%nl (igk_k(ig,ikk) ), 2) = evc (ig+npwx, ibnd)
        enddo
        CALL invfft ('Wave', psi(:,1), dffts)
        CALL invfft ('Wave', psi(:,2), dffts)
        !
        ! FFT to R-space of the perturbed wfct's dpsi
        !
        dpsic = (0.d0, 0.d0)
        do ig = 1, npwq
           dpsic (dffts%nl (igk_k(ig,ikq)), 1 ) = dpsi (ig, ibnd)
           dpsic (dffts%nl (igk_k(ig,ikq)), 2 ) = dpsi (ig+npwx, ibnd)
        enddo
        CALL invfft ('Wave', dpsic(:,1), dffts)
        CALL invfft ('Wave', dpsic(:,2), dffts)
        !
        ! Calculation of the response charge density
        !
        do ir = 1, dffts%nnr
           drhoscf(ir,1)=drhoscf(ir,1)+wgt*(CONJG(psi(ir,1))*dpsic(ir,1)  +  &
                                            CONJG(psi(ir,2))*dpsic(ir,2) )
        enddo
        IF (domag) THEN
           do ir = 1, dffts%nnr
              drhoscf(ir,2)=drhoscf (ir,2) + (rsign) *wgt * (CONJG(psi(ir,1))*dpsic(ir,2) &
                                                  + CONJG(psi(ir,2))*dpsic(ir,1) )
              drhoscf(ir,3)=drhoscf (ir,3) + (rsign) *wgt * (CONJG(psi(ir,1))*dpsic(ir,2) &
                                                  - CONJG(psi(ir,2))*dpsic(ir,1) ) * (0.d0,-1.d0)
              drhoscf(ir,4)=drhoscf (ir,4) + (rsign) *wgt * (CONJG(psi(ir,1))*dpsic(ir,1) &
                                                  - CONJG(psi(ir,2))*dpsic(ir,2) )
           enddo
        END IF
        !
     END IF
     !
  enddo
  !
  ! Ultrasoft contribution
  ! Calculate dbecsum_nc = <evc|vkb><vkb|dpsi>
  !
  IF (rsign==1.0d0) THEN
     CALL addusdbec_nc (ik, weight, dpsi, dbecsum, becp1)
  ELSE
     CALL addusdbec_nc (ik, weight, dpsi, dbecsum, becpt)
  ENDIF
  !
  DEALLOCATE(psi)
  DEALLOCATE(dpsic)
  !
  IF (dffts%has_task_groups) THEN
     DEALLOCATE(tg_psi)
     DEALLOCATE(tg_dpsi)
     DEALLOCATE(tg_drho)
  END IF
  !
  CALL stop_clock ('incdrhoscf')
  !
  RETURN
  !
end subroutine incdrhoscf_nc
