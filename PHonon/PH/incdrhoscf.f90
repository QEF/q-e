!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
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
  !
  USE kinds, only : DP
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat
  USE fft_base,  ONLY: dffts
  USE fft_interfaces, ONLY: invfft
  USE gvecs,   ONLY : nls
  USE wvfct,     ONLY : npw, igk, npwx, nbnd
  USE uspp_param,ONLY: nhm
  USE wavefunctions_module,  ONLY: evc
  USE qpoint,    ONLY : npwq, igkq, ikks
  USE control_lr, ONLY : nbnd_occ
  USE mp_bands,   ONLY : me_bgrp, inter_bgrp_comm, ntask_groups
  USE mp, ONLY : mp_sum


  implicit none
  ! I/O variables
  integer, INTENT (IN) :: ik
  ! input: the k point
  real(DP), INTENT (IN) :: weight
  ! input: the weight of the k point
  complex(DP), INTENT (IN) :: dpsi (npwx,nbnd)
  ! input: the perturbed wfc for the given k point
  complex(DP), INTENT (INOUT) :: drhoscf (dffts%nnr), dbecsum (nhm*(nhm+1)/2,nat)
  ! input/output: the accumulated change to the charge density and dbecsum
  !
  !
  !   here the local variable
  !
  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable  :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space
  complex(DP), allocatable :: tg_psi(:), tg_dpsi(:), tg_drho(:)

  integer :: ibnd, ikk, ir, ig, incr, v_siz, idx, ioff
  ! counters

  call start_clock ('incdrhoscf')
  IF (ntask_groups > 1) dffts%have_task_groups=.TRUE.
  allocate (dpsic(  dffts%nnr))
  allocate (psi  (  dffts%nnr))
  wgt = 2.d0 * weight / omega
  ikk = ikks(ik)
  incr=1
  !
  IF (dffts%have_task_groups) THEN
     !
     v_siz = dffts%tg_nnr * dffts%nogrp
     !
     ALLOCATE( tg_psi( v_siz ) )
     ALLOCATE( tg_dpsi( v_siz ) )
     ALLOCATE( tg_drho( v_siz ) )
     !
     incr  = dffts%nogrp
     !
  ENDIF
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  do ibnd = 1, nbnd_occ(ikk), incr
     !
     IF (dffts%have_task_groups) THEN
        !
        tg_drho=(0.0_DP, 0.0_DP)
        tg_psi=(0.0_DP, 0.0_DP)
        tg_dpsi=(0.0_DP, 0.0_DP)
        !
        ioff   = 0
        !
        DO idx = 1, dffts%nogrp
           !
           ! ... dffts%nogrp ffts at the same time. We prepare both
           ! evc (at k) and dpsi (at k+q)
           !
           IF( idx + ibnd - 1 <= nbnd_occ(ikk) ) THEN
              !
              DO ig = 1, npw
                 tg_psi( nls( igk( ig ) ) + ioff ) = evc( ig, idx+ibnd-1 )
              END DO
              DO ig = 1, npwq
                 tg_dpsi( nls( igkq( ig ) ) + ioff ) = dpsi( ig, idx+ibnd-1 )
              END DO
              !
           END IF
           !
           ioff = ioff + dffts%tg_nnr
           !
        END DO
        CALL invfft ('Wave', tg_psi, dffts)
        CALL invfft ('Wave', tg_dpsi, dffts)

        do ir = 1, dffts%tg_npp( me_bgrp + 1 ) * dffts%nr1x * dffts%nr2x
           tg_drho (ir) = tg_drho (ir) + wgt * CONJG(tg_psi (ir) ) *  &
                                                     tg_dpsi (ir)
        enddo
        !
        ! reduce the group charge (equivalent to sum over bands of 
        ! orbital group)
        !
        CALL mp_sum( tg_drho, gid = dffts%ogrp_comm )
        !
        ioff = 0
        DO idx = 1, dffts%nogrp
           IF( me_bgrp == dffts%nolist( idx ) ) EXIT
           ioff = ioff + dffts%nr1x * dffts%nr2x * &
                                      dffts%npp( dffts%nolist( idx ) + 1 )
        END DO
        !
        ! copy the charge back to the proper processor location
        !
        DO ir = 1, dffts%nnr
           drhoscf(ir) = drhoscf(ir) + tg_drho(ir+ioff)
        END DO
     ELSE
        psi (:) = (0.d0, 0.d0)
        do ig = 1, npw
           psi (nls (igk (ig) ) ) = evc (ig, ibnd)
        enddo

        CALL invfft ('Wave', psi, dffts)

        dpsic(:) = (0.d0, 0.d0)
        do ig = 1, npwq
           dpsic (nls (igkq (ig) ) ) = dpsi (ig, ibnd)
        enddo

        CALL invfft ('Wave', dpsic, dffts)

        do ir = 1, dffts%nnr
           drhoscf (ir) = drhoscf (ir) + wgt * CONJG(psi (ir) ) * dpsic (ir)
        enddo
     ENDIF
  enddo
  

  call addusdbec (ik, weight, dpsi, dbecsum)
  deallocate (psi)
  deallocate (dpsic)
  IF (dffts%have_task_groups) THEN
     DEALLOCATE(tg_psi)
     DEALLOCATE(tg_dpsi)
     DEALLOCATE(tg_drho)
  ENDIF
  dffts%have_task_groups=.FALSE.

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf
