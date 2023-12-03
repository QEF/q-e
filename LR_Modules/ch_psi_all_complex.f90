!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE ch_psi_all_complex (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in ah.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba
  USE wvfct,                ONLY : npwx, nbnd, current_k
  USE becmod,               ONLY : becp, calbec
  USE uspp,                 ONLY : nkb, vkb
  USE fft_base,             ONLY : dffts
  USE gvect,                ONLY : g
  USE klist,                ONLY : xk, igk_k
  USE noncollin_module,     ONLY : noncolin, npol
  USE eqv,                  ONLY : evq
  USE qpoint,               ONLY : ikqs
  USE mp_bands,             ONLY : use_bgrp_in_hpsi, inter_bgrp_comm, intra_bgrp_comm
  USE xc_lib,               ONLY : exx_is_active
  USE mp,                   ONLY : mp_sum
  USE control_lr,           ONLY : alpha_pv, nbnd_occ, lgamma
  USE control_flags,        ONLY : gamma_only, offload_type
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer
  USE io_files,             ONLY : nwordwfcU
  USE ldaU,                 ONLY : lda_plus_u, wfcU, lda_plus_u_kind
  USE units_lr,             ONLY : iuatswfc

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  COMPLEX(DP), INTENT(IN) :: e (m)
  ! input: the eigenvalue

  COMPLEX(DP), INTENT(IN)  :: h (npwx*npol, m)
  complex(DP), INTENT(OUT) :: ah (npwx*npol, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  INTEGER :: ibnd, ig
  ! counter on bands
  ! counter on G vetors

  COMPLEX(DP), ALLOCATABLE :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h

  CALL start_clock ('ch_psi')
  !
  !  This routine is task groups aware
  !
  ALLOCATE (ps  ( nbnd , m))
  ALLOCATE (hpsi( npwx*npol , m))
  ALLOCATE (spsi( npwx*npol , m))
  !
  current_k = ikqs(ik) ! k+q
  !
  ! DFT+U case
  !
  IF (lda_plus_u) THEN
     !
     ! Read the atomic orbitals (S*phi) at k+q from the file
     ! and put the result in wfcU, because it is used
     ! in vhpsi (via the call of h_psi) to compute the 
     ! Hubbard potential.
     !
     CALL get_buffer (wfcU, nwordwfcU, iuatswfc, current_k)
     !
     ! Compute the phase factor at k+q
     !IF (lda_plus_u_kind.EQ.2) CALL phase_factor(ikq)
     !
  ENDIF
  !
  ! Compute an action of the Hamiltonian and the S operator
  ! on the h vector (i.e. H*h and S*h, respectively).
  !
  !$acc enter data create(hpsi(1:npwx*npol, 1:m), spsi(1:npwx*npol, 1:m), ps(1:nbnd, 1:m))
  !$acc kernels present(hpsi,spsi)
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !$acc end kernels
#if defined(__CUDA)
  !$acc data present(h, hpsi, spsi)
  !$acc host_data use_device(h, hpsi, spsi)
  CALL h_psi_gpu (npwx, n, m, h, hpsi)
  CALL s_psi_acc (npwx, n, m, h, spsi)
  !$acc end host_data
  !$acc end data

#else

  CALL h_psi (npwx, n, m, h, hpsi)
  CALL s_psi (npwx, n, m, h, spsi)

#endif

  CALL start_clock ('last')
  !
  !   then we compute ( H - \epsilon S ) * h
  !   and put the result in ah
  !
  CALL start_clock ('Hesh')
  !$acc data copyin(e)

  !$acc kernels present(ah)
  ah=(0.d0,0.d0)
  !$acc end kernels

  !$acc parallel loop collapse(2) present(hpsi, spsi, ah, e)
  DO ibnd = 1, m
     DO ig = 1, n
        ah(ig, ibnd) = hpsi(ig, ibnd) - e(ibnd) * spsi(ig, ibnd)
     ENDDO
  ENDDO
  IF (noncolin) THEN
     CALL start_clock ('Hesh:noncolin')
     !$acc parallel loop collapse(2) present(hpsi, spsi, ah, e)
     DO ibnd = 1, m
        DO ig = 1, n
           ah(ig+npwx, ibnd) = hpsi(ig+npwx, ibnd) - e(ibnd) * spsi(ig+npwx, ibnd)
        ENDDO
     ENDDO
     CALL stop_clock ('Hesh:noncolin')
  ENDIF
  !$acc end data
  CALL stop_clock ('Hesh')

  !
  !   lastly we compute alpha_pv * P_v * h (if alpha_pv.NE.0.0d0)
  !   and add it to ah 
  !
  IF (alpha_pv.NE.0.0d0) THEN
     IF (gamma_only) THEN
        CALL ch_psi_all_gamma_complex()
     ELSE
        CALL ch_psi_all_k_complex()
     ENDIF
  ENDIF
!!!!  !$acc exit data copyout(ah) delete(ps)
  !$acc exit data delete(hpsi, spsi, ps)
  DEALLOCATE (spsi)
  DEALLOCATE (hpsi)
  DEALLOCATE (ps)

  CALL stop_clock ('last')
  CALL stop_clock ('ch_psi')
  RETURN
CONTAINS

  SUBROUTINE ch_psi_all_k_complex()
    !
    ! K-point part
    !
#if defined(__CUDA)
    USE cublas
#endif
    
    IMPLICIT NONE
    INTEGER :: m_start, m_end
    INTEGER :: k
    INTEGER :: ibnd, ig

    k = nbnd_occ (ikqs(ik))
    CALL start_clock_gpu ('ch_psi_all_k_complex')
    !$acc data present(evq, ps, hpsi, spsi)
    !
    !   Here we compute the projector in the valence band
    !
    !$acc kernels present(ps)
    ps (:,:) = (0.d0, 0.d0)
    !$acc end kernels
    !
    ! ikqs(ik) is the index of the point k+q if q\=0
    !          is the index of the point k   if q=0
    !
    !$acc host_data use_device(spsi, ps, evq)
    IF (noncolin) THEN
       CALL zgemm ('C', 'N', k, m, npwx*npol, (1.d0, 0.d0) , evq, &
            npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
    ELSE
       CALL zgemm ('C', 'N', k, m, n, (1.d0, 0.d0) , evq, &
            npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
    ENDIF
    !$acc end host_data
    !$acc kernels present(ps,hpsi) 
    ps (:,:) = ps(:,:) * alpha_pv
    hpsi (:,:) = (0.d0, 0.d0)
    !$acc end kernels
    !$acc host_data use_device(ps)
    CALL mp_sum ( ps, intra_bgrp_comm )
    !$acc end host_data
    !$acc host_data use_device(hpsi, ps, evq)
    IF (noncolin) THEN
       CALL zgemm ('N', 'N', npwx*npol, m, k, (1.d0, 0.d0) , evq, &
            npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
    ELSE
       CALL zgemm ('N', 'N', n, m, k, (1.d0, 0.d0) , evq, &
            npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
    END IF
    !$acc end host_data
    !$acc kernels present(spsi, hpsi)
    spsi(:,:) = hpsi(:,:)
    !$acc end  kernels
    !$acc end data
    !
    !    And apply S again
    !
    CALL start_clock_gpu ('ch_psi_calbec')
    if (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) then
       call divide (inter_bgrp_comm, m, m_start, m_end)
       if (m_end >= m_start) then
          CALL calbec (offload_type, n, vkb, hpsi(:,m_start:m_end), becp, m_end- m_start + 1)
       endif
    else
       CALL calbec (offload_type, n, vkb, hpsi, becp, m)
    endif
    CALL stop_clock_gpu ('ch_psi_calbec')
    !$acc host_data use_device(hpsi, spsi)
    CALL s_psi_acc (npwx, n, m, hpsi, spsi)
    !$acc end host_data
    !$acc parallel loop collapse(2) present(ah, spsi)
    DO ibnd = 1, m
       DO ig = 1, n
          ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
       ENDDO
    ENDDO
    !$acc end parallel loop 
    IF (noncolin) THEN
       !$acc parallel loop collapse(2) present(ah, spsi)
       DO ibnd = 1, m
          DO ig = 1, n
             ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
          ENDDO
       ENDDO
       !$acc end parallel loop
    END IF
    CALL stop_clock_gpu ('ch_psi_all_k_complex')
    return
  END SUBROUTINE ch_psi_all_k_complex

  SUBROUTINE ch_psi_all_gamma_complex()
    !
    ! gamma_only case
    !  
    USE realus, ONLY : real_space, invfft_orbital_gamma, &
                       fwfft_orbital_gamma, calbec_rs_gamma,  s_psir_gamma
    use gvect,  only : gstart

    IMPLICIT NONE
    INTEGER :: m_start, m_end

    ps (:,:) = 0.d0
    
    IF (noncolin) THEN
       CALL errore('ch_psi_all', 'non collin in gamma point not implemented',1)
    ELSE
       CALL DGEMM( 'C', 'N', nbnd, m, 2*n, 2.D0,evc, 2*npwx*npol, spsi, 2*npwx*npol, 0.D0, ps, nbnd )
       if(gstart==2) CALL DGER(nbnd, m, -1.0_DP, evc, 2*npwx, spsi, 2*npwx, ps, nbnd )
    ENDIF
    ps (:,:) = ps(:,:) * alpha_pv
    CALL mp_sum ( ps, intra_bgrp_comm )

    hpsi (:,:) = (0.d0, 0.d0)

    IF (noncolin) THEN
       CALL ZGEMM ('N', 'N', npwx*npol, m, nbnd_occ (ik) , (1.d0, 0.d0) , evc, &
            npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
    ELSE
       CALL DGEMM ('N', 'N', 2*n, m, nbnd_occ (ik) , 1.d0 , evc, &
            2*npwx, ps, nbnd, 1.d0 , hpsi, 2*npwx)
    ENDIF
    spsi(:,:) = hpsi(:,:)
    !
    !    And apply S again
    !
    IF (real_space ) THEN
       DO ibnd=1,m,2
          CALL invfft_orbital_gamma(hpsi,ibnd,m)
          CALL calbec_rs_gamma(ibnd,m,becp%r)
          CALL s_psir_gamma(ibnd,m)
          CALL fwfft_orbital_gamma(spsi,ibnd,m)
       ENDDO
    ELSE
       if (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) then
          call divide( inter_bgrp_comm, m, m_start, m_end)
          if (m_end >= m_start) CALL calbec (n, vkb, hpsi(:,m_start:m_end), becp, m_end- m_start + 1)
       else
          CALL calbec (n, vkb, hpsi, becp, m)
       end if
       CALL s_psi (npwx, n, m, hpsi, spsi)
    ENDIF
    DO ibnd = 1, m
       DO ig = 1, n
          ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
       ENDDO
    ENDDO
    IF (noncolin) THEN
       DO ibnd = 1, m
          DO ig = 1, n
             ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
          ENDDO
       ENDDO
    ENDIF
    return
  END SUBROUTINE ch_psi_all_gamma_complex
 
END SUBROUTINE ch_psi_all_complex
