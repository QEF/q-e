!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhoscf (drhoscf, weight, ik, dbecsum, dpsi)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE ions_base,            ONLY : nat
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE wvfct,                ONLY : npwx, nbnd
  USE uspp_param,           ONLY : nhm
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : ngk,igk_k
  USE qpoint,               ONLY : ikks, ikqs
  USE control_lr,           ONLY : nbnd_occ
  USE mp_bands,             ONLY : me_bgrp, inter_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum
  USE fft_helper_subroutines

  IMPLICIT NONE
  !
  ! I/O variables
  INTEGER, INTENT (IN) :: ik
  ! input: the k point
  REAL(DP), INTENT (IN) :: weight
  ! input: the weight of the k point
  COMPLEX(DP), INTENT (IN) :: dpsi (npwx,nbnd)
  ! input: the perturbed wfc for the given k point
  COMPLEX(DP), INTENT (INOUT) :: drhoscf (dffts%nnr), dbecsum (nhm*(nhm+1)/2,nat)
  ! input/output: the accumulated change to the charge density and dbecsum
  !
  !   here the local variables
  !
  REAL(DP) :: wgt
  ! the effective weight of the k point

  COMPLEX(DP), ALLOCATABLE :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  COMPLEX(DP), ALLOCATABLE :: tg_psi(:), tg_dpsi(:), tg_drho(:)

  INTEGER :: npw, npwq, ikk, ikq, itmp
  INTEGER :: ibnd, ir, ir3, ig, incr, v_siz, idx, ioff, ioff_tg, nxyp
  INTEGER :: right_inc, ntgrp
  ! counters

  ! For device buffer 
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nl_d(:)
  !
  nl_d  => dffts%nl_d
  !$acc update device(evc) 
#else
  INTEGER, ALLOCATABLE :: nl_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
#endif
  

  CALL start_clock_gpu ('incdrhoscf')
  !
  ALLOCATE(dpsic(dffts%nnr))
  ALLOCATE(psi(dffts%nnr))
  !
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  ikq = ikqs(ik)
  npw = ngk(ikk)
  npwq= ngk(ikq)
  incr = 1
  !
  IF ( dffts%has_task_groups ) THEN
     !
     v_siz = dffts%nnr_tg
     !
     ALLOCATE( tg_psi( v_siz ) )
     ALLOCATE( tg_dpsi( v_siz ) )
     ALLOCATE( tg_drho( v_siz ) )
     !
     incr = fftx_ntgrp(dffts)
     !
  ELSE
     v_siz = dffts%nnr
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  !$acc data copyin(dpsi(1:npwx,1:nbnd)) copy(drhoscf(1:v_siz)) create(psi(1:v_siz),dpsic(1:v_siz)) present(igk_k) deviceptr(nl_d) 
  do ibnd = 1, nbnd_occ(ikk), incr
     !
     IF ( dffts%has_task_groups ) THEN
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
                 tg_psi( dffts%nl( igk_k( ig,ikk ) ) + ioff ) = evc( ig, idx+ibnd-1 )
              END DO
              DO ig = 1, npwq
                 tg_dpsi( dffts%nl( igk_k( ig,ikq ) ) + ioff ) = dpsi( ig, idx+ibnd-1 )
              END DO
              !
           END IF
           !
           ioff = ioff + right_inc
           !
        END DO
        CALL invfft ('tgWave', tg_psi, dffts)
        CALL invfft ('tgWave', tg_dpsi, dffts)

        do ir = 1, dffts%nr1x * dffts%nr2x * dffts%my_nr3p
           tg_drho (ir) = tg_drho (ir) + wgt * CONJG(tg_psi (ir) ) *  tg_dpsi (ir)
        enddo
        !
        ! reduce the group charge (equivalent to sum over bands of 
        ! orbital group)
        !
        CALL tg_reduce_rho( drhoscf, tg_drho, dffts )
        !
     ELSE
        !
        ! Normal case: no task groups
        !
        ! Initialize psi and dpsic
        !
        !$acc kernels 
        psi (:) = (0.d0, 0.d0)
        dpsic(:) = (0.d0, 0.d0)
        !$acc end kernels
        !
        !$acc parallel loop 
        do ig = 1, npw
           itmp = nl_d (igk_k(ig,ikk) )
           psi (itmp ) = evc (ig, ibnd)
        enddo
        !$acc parallel loop
        do ig = 1, npwq
           itmp = nl_d (igk_k(ig,ikq) )
           dpsic ( itmp ) = dpsi (ig, ibnd)
        enddo
        !
        ! FFT to R-space of the unperturbed/perturbed wfcts psi/dpsi
        !
        !$acc host_data use_device(psi)
        CALL invfft ('Wave', psi, dffts)
        !$acc end host_data
        !$acc host_data use_device(dpsic)
        CALL invfft ('Wave', dpsic, dffts)
        !$acc end host_data
        !
        ! Calculation of the response charge-density
        !
        !$acc parallel loop 
        do ir = 1, v_siz
           drhoscf (ir) = drhoscf (ir) + wgt * CONJG(psi (ir) ) * dpsic (ir)
        enddo
        !
     ENDIF
     !
  enddo ! loop on bands
  !$acc end data
  !
  ! Ultrasoft contribution
  ! Calculate dbecsum = <evc|vkb><vkb|dpsi>
  ! 
  CALL addusdbec (ik, weight, dpsi, dbecsum)
  !
  DEALLOCATE(psi)
  DEALLOCATE(dpsic)
  !
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE(tg_psi)
     DEALLOCATE(tg_dpsi)
     DEALLOCATE(tg_drho)
  ENDIF
  !
  CALL stop_clock_gpu ('incdrhoscf')
  !
  RETURN
  !
end subroutine incdrhoscf
