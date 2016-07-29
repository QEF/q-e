!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
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
  USE fft_base,             ONLY : dffts, dtgs
  USE gvect,                ONLY : g
  USE klist,                ONLY : xk, igk_k
  USE noncollin_module,     ONLY : noncolin, npol
  USE eqv,                  ONLY : evq
  USE qpoint,               ONLY : ikqs
  USE mp_bands,             ONLY : intra_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum
  USE control_lr,           ONLY : alpha_pv, nbnd_occ, lgamma
  !Needed only for TDDFPT
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc

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

  COMPLEX(DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h
  INTEGER, ALLOCATABLE :: ibuf(:)

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
  !   compute an action of the Hamiltonian and the S operator
  !   on the h vector (i.e. H*h and S*h, respectively).
  !
  current_k = ikqs(ik)
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
    CALL calbec (n, vkb, hpsi, becp, m)
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
    USE realus, ONLY : real_space, real_space_debug, invfft_orbital_gamma, &
                       fwfft_orbital_gamma, calbec_rs_gamma,  s_psir_gamma
    use gvect,  only : gstart

    IMPLICIT NONE

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
    IF (real_space_debug >6 ) THEN
       DO ibnd=1,m,2
          CALL invfft_orbital_gamma(hpsi,ibnd,m)
          CALL calbec_rs_gamma(ibnd,m,becp%r)
          CALL s_psir_gamma(ibnd,m)
          CALL fwfft_orbital_gamma(spsi,ibnd,m)
       ENDDO
    ELSE
       CALL calbec (n, vkb, hpsi, becp, m)
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
