!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_Opsi_magnons (ik, ip, dOpsi)
  !--------------------------------------------------------------------------
  !
  ! This subroutine calculates the starting Lanczos vectors.
  ! Magnons case
  !
  ! V0psi_X = { P_C V_ext(q,\omega) \psi_nk }
  ! V0psi_Y = { T P_C \sigma V_ext(-q,-\omega) \psi_n-k }
  !
  ! X,Y standign for the resonant and anti-resonant batch respectively,
  ! T time-reversal operator,
  ! n,k running indexes.
  !
  ! Adapted from lr_dvspi_eels
  !
  ! OBS.: for the time being, working only for Q=q \in BZ
  !
  USE kinds,                 ONLY : DP
  USE wvfct,                 ONLY : npwx, nbnd
  USE fft_base,              ONLY : dffts
  USE gvecw,                 ONLY : gcutw
  USE qpoint,                ONLY : ikks, ikqs, nksq 
  USE eqv,                   ONLY : evq, dpsi 
  USE wavefunctions,         ONLY : evc
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag
  use klist,                 only : xk, igk_k, ngk
  use gvect,                 only : ngm, g
  USE control_lr,            ONLY : nbnd_occ, nbnd_occx
  USE io_files,              ONLY : iunwfc, nwordwfc
  use uspp,                  only : vkb, okvan
  USE mp_bands,              ONLY : ntask_groups
  USE buffers,               ONLY : get_buffer
  USE fft_helper_subroutines
  USE lr_variables,          ONLY : iunTwfc
  USE lr_magnons_routines,   ONLY : pauli
 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik 
  INTEGER, INTENT(IN) :: ip
  COMPLEX(DP), INTENT(OUT) :: dOpsi(npwx*npol,nbnd_occx,2)
  INTEGER :: ibnd, ig, ios, incr, v_siz
  INTEGER :: ikk, & ! index of the point k
             ikq, & ! index of the point k+q
             npwq, & ! number of the plane-waves at point k+q
             imk, &
             imkq
  COMPLEX(DP), ALLOCATABLE :: revc(:,:), &  ! wavefunctions in R-space
                            & tg_psic(:,:)  ! wavefunctions in R-space (for task groups) 
  LOGICAL :: exst
  !
  COMPLEX(DP) :: Tevc(npol*npwx,nbnd)       ! T-rev op. applied to u_{-k}
  COMPLEX(DP) :: Tevq(npol*npwx,nbnd)       ! T-rev op. applied to u_{-k-Q}
  !
  INTEGER :: ibnd1, ibnd2 
  !
  CALL start_clock ('lr_Opsi_magnons') 
  !
  ALLOCATE (revc(dffts%nnr,npol))
  revc(:,:)    = (0.d0, 0.d0)
  dpsi(:,:)    = (0.d0, 0.d0)
  dOpsi(:,:,:) = (0.d0, 0.d0)
  !
  incr = 1
  !
  IF ( dffts%has_task_groups ) THEN
     !
     v_siz =  dffts%nnr_tg
     !
     ALLOCATE( tg_psic(v_siz,npol) )
     !
     incr = fftx_ntgrp(dffts)
     !
  ENDIF 
  !
  ! Setting up positions of k-points
  !
  ikk  = ikks(ik)
  ikq  = ikqs(ik)
  npwq = ngk(ikq)
  !
  if ( mod(ik,2) == 0) then
     imk  = ikk -3      ! position of -k
     imkq = ikk -1      ! position of -k-Q
  else
     imk  = ikk +3      ! position of -k
     imkq = ikk +5      ! position of -k-Q
  endif  
  !
  ! Read the ground-state wavefunctions evc(k) and evq(k+q).
  ! Note: even in the case when there is only one k point (nksq=1),
  ! the ground-state wavefunctions evc(k) and evq(k+q) are read here
  ! and kept in memory for the rest of the code.
  !
  ! Resonant Batch
  !
  CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
  CALL get_buffer (evq, nwordwfc, iunwfc, ikq)
  !
  ! Re-ordering of the G vectors.
  !
  DO ibnd = 1, nbnd_occ(ikk), incr
     !
     IF ( dffts%has_task_groups ) THEN
        !
        ! FFT to R-space
        CALL cft_wave_tg(ik, evc, tg_psic, 1, v_siz, ibnd, nbnd_occ(ikk) )
        !
        ! back-FFT to G-space
        CALL cft_wave_tg(ik, dOpsi(:,:,1), tg_psic, -1, v_siz, ibnd, nbnd_occ(ikk))
        !
        CALL pauli(dOpsi(:,ibnd,1), ip)
        !
     ELSE
        !
        ! FFT to R-space
!$acc data copyin(evc(1:npwx*npol,ibnd)) copy(revc, dOpsi(1:npwx*npol,ibnd,1))        
        CALL cft_wave(ik, evc(1,ibnd), revc, +1)
        !
        ! back-FFT to G-space
        CALL cft_wave(ik, dOpsi(1,ibnd,1), revc, -1) 
!$acc end data
        !
        CALL pauli(dOpsi(:,ibnd,1), ip)
        !
     ENDIF
     !
  ENDDO
  !
  ! Anti-Resonant Batch
  !
  ! CHange k-point --> set to zero again all variables to avoid leftovers in
  ! wfc
  dpsi(:,:)    = (0.d0, 0.d0)
  revc(:,:)    = (0.d0, 0.d0)  
  !
  ! Re-ordering of the G vectors of Tevc.
  ! The routine cft_wave transforms according to +k and +k+Q indexes, which is
  ! correct also in this case since (i sigma_y) u^*_{-k} transforms with +k
  ! indexes, and e^{i Q r} (i sigma_y) u^*_{-k} according to +k+Q indexes.
  !
  ! The No. of occupied bands remains the same as of -k
  !
  ! Read unperturbed wavefuctions Tu_{-k} and Tu_{-k-Q}
  !
  Tevc(:,:) = (0.0d0, 0.0d0)
  Tevq(:,:) = (0.0d0, 0.0d0)
  !
  CALL get_buffer (Tevc, nwordwfc, iunTwfc, 2*ik-1)
  CALL get_buffer (Tevq, nwordwfc, iunTwfc, 2*ik)
  !
  ! Re-ordering of the G vectors.
  !
  DO ibnd = 1, nbnd_occ(imk), incr
     !
     IF ( dffts%has_task_groups ) THEN
        !
        ! FFT to R-space
        CALL cft_wave_tg(ik, Tevc, tg_psic, 1, v_siz, ibnd, nbnd_occ(imk) )
        !
        ! back-FFT to G-space
        CALL cft_wave_tg(ik, dOpsi(:,:,2), tg_psic, -1, v_siz, ibnd, nbnd_occ(imk))
        !
        CALL pauli(dOpsi(:,ibnd,2), ip)
        !
     ELSE
        !
        ! FFT to R-space
!$acc data copyin(Tevc(1:npwx*npol,ibnd)) copy(revc, dOpsi(1:npwx*npol,ibnd,2))        
        CALL cft_wave(ik, Tevc(1,ibnd), revc, +1)
        !
        ! back-FFT to G-space
        CALL cft_wave(ik, dOpsi(1,ibnd,2), revc, -1) 
!$acc end data
        !
        CALL pauli(dOpsi(:,ibnd,2), ip)
        !
     ENDIF
     !
  ENDDO
  !
  ! Put the correct sign
  !
  dOpsi(:,:,2) = - dOpsi(:,:,2)
  !
  DEALLOCATE (revc)
  !
  IF ( dffts%has_task_groups ) THEN
     DEALLOCATE( tg_psic )
  ENDIF
  !
  CALL stop_clock ('lr_Opsi_magnons')
  !
  RETURN
  !
END SUBROUTINE lr_Opsi_magnons
