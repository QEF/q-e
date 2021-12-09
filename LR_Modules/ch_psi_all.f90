!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE ch_psi_all (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in ah.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba
  USE wvfct,                ONLY : npwx, nbnd, current_k
  USE becmod,               ONLY : bec_type, becp, calbec
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
  USE control_flags,        ONLY : gamma_only
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

  REAL(DP), INTENT(IN) :: e (m)
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
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
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
     ! Compute the phase factor at k+q (needed for DFT+U+V)
     !
     IF (lda_plus_u_kind.EQ.2) CALL phase_factor (current_k)
     !
  ENDIF
  !
  ! Compute an action of the Hamiltonian and the S operator
  ! on the h vector (i.e. H*h and S*h, respectively).
  !
  CALL h_psi (npwx, n, m, h, hpsi)
  CALL s_psi (npwx, n, m, h, spsi)
  !
  CALL start_clock ('last')
  !
  !   then we compute ( H - \epsilon S ) * h
  !   and put the result in ah
  !
  ah=(0.d0,0.d0)
  DO ibnd = 1, m
     DO ig = 1, n
        ah(ig, ibnd) = hpsi(ig, ibnd) - e(ibnd) * spsi(ig, ibnd)
     ENDDO
  ENDDO
  IF (noncolin) THEN
     DO ibnd = 1, m
        DO ig = 1, n
           ah(ig+npwx, ibnd) = hpsi(ig+npwx, ibnd) - e(ibnd) * spsi(ig+npwx, ibnd)
        ENDDO
     ENDDO
  ENDIF
  !
  !   lastly we compute alpha_pv * P_v * h (if alpha_pv.NE.0.0d0)
  !   and add it to ah 
  !
  IF (alpha_pv.NE.0.0d0) THEN
     IF (gamma_only) THEN
        CALL ch_psi_all_gamma()
     ELSE
        CALL ch_psi_all_k()
     ENDIF
  ENDIF
  
  DEALLOCATE (spsi)
  DEALLOCATE (hpsi)
  DEALLOCATE (ps)

  CALL stop_clock ('last')
  CALL stop_clock ('ch_psi')
  RETURN
CONTAINS

  SUBROUTINE ch_psi_all_k()
    !
    ! K-point part
    !
    USE becmod, ONLY : becp, calbec
    
    IMPLICIT NONE
    INTEGER :: m_start, m_end
    !
    !   Here we compute the projector in the valence band
    !
    ps (:,:) = (0.d0, 0.d0)
    !
    ! ikqs(ik) is the index of the point k+q if q\=0
    !          is the index of the point k   if q=0
    ! 
    IF (noncolin) THEN
       CALL zgemm ('C', 'N', nbnd_occ (ikqs(ik)) , m, npwx*npol, (1.d0, 0.d0) , evq, &
            npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
    ELSE
       CALL zgemm ('C', 'N', nbnd_occ (ikqs(ik)) , m, n, (1.d0, 0.d0) , evq, &
            npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
    ENDIF
    ps (:,:) = ps(:,:) * alpha_pv
    CALL mp_sum ( ps, intra_bgrp_comm )
    
    hpsi (:,:) = (0.d0, 0.d0)
    IF (noncolin) THEN
       CALL zgemm ('N', 'N', npwx*npol, m, nbnd_occ (ikqs(ik)) , (1.d0, 0.d0) , evq, &
            npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
    ELSE
       CALL zgemm ('N', 'N', n, m, nbnd_occ (ikqs(ik)) , (1.d0, 0.d0) , evq, &
            npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
    END IF
    spsi(:,:) = hpsi(:,:)
    !
    !    And apply S again
    !
    if (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. m > 1) then
       call divide (inter_bgrp_comm, m, m_start, m_end)
       if (m_end >= m_start) CALL calbec (n, vkb, hpsi(:,m_start:m_end), becp, m_end- m_start + 1)
    else
       CALL calbec (n, vkb, hpsi, becp, m)
    end if
    CALL s_psi (npwx, n, m, hpsi, spsi)
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
    END IF
    return
  END SUBROUTINE ch_psi_all_k

  SUBROUTINE ch_psi_all_gamma()
    !
    ! gamma_only case
    !  
    USE becmod, ONLY : becp,  calbec
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
  END SUBROUTINE ch_psi_all_gamma
 
END SUBROUTINE ch_psi_all
