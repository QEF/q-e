!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE ch_psi_all (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  ! Merged with lr_ch_psi_all June 2011. This function is now used both
  ! in ph.x and turbo_lanczos.x
  !

  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npwx, nbnd
  USE becmod,               ONLY : bec_type, becp, calbec
  USE uspp,                 ONLY : nkb, vkb
  USE fft_base,             ONLY : dffts
  USE wvfct,                ONLY : npwx, igk
  USE qpoint,               ONLY : igkq
  USE noncollin_module,     ONLY : noncolin, npol

  USE control_ph,           ONLY : alpha_pv, nbnd_occ, lgamma
  USE eqv,                  ONLY : evq
  USE qpoint,               ONLY : ikqs

  USE mp_bands,             ONLY : intra_bgrp_comm, ntask_groups
  USE mp,                   ONLY : mp_sum

  !Needed only for TDDFPT
  USE control_flags,        ONLY : gamma_only, tddfpt
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
  INTEGER :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
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
  IF (ntask_groups > 1) dffts%have_task_groups=.TRUE.

  ALLOCATE (ps  ( nbnd , m))
  ALLOCATE (hpsi( npwx*npol , m))
  ALLOCATE (spsi( npwx*npol , m))
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  IF (dffts%have_task_groups) THEN
  !
  !   With task groups we use the Hpsi routine of PW parallelized
  !   on task groups
  !
     IF (.NOT.lgamma) THEN
        ALLOCATE(ibuf(npwx))
        ibuf=igk
        igk=igkq 
     ENDIF   
     CALL h_psi (npwx, n, m, h, hpsi)
     CALL s_psi (npwx, n, m, h, spsi)
     IF (.NOT.lgamma) THEN
        igk=ibuf
        DEALLOCATE(ibuf)
     ENDIF
  ELSE
    CALL h_psiq (npwx, n, m, h, hpsi, spsi)
  ENDIF

  CALL start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah=(0.d0,0.d0)
  DO ibnd = 1, m
     DO ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     ENDDO
  ENDDO
  IF (noncolin) THEN
     DO ibnd = 1, m
        DO ig = 1, n
           ah (ig+npwx,ibnd)=hpsi(ig+npwx,ibnd)-e(ibnd)*spsi(ig+npwx,ibnd)
        ENDDO
     ENDDO
  ENDIF

  IF(gamma_only) THEN

       CALL ch_psi_all_gamma()
    ELSE

       IF(tddfpt) THEN
          ikq = ik
          evq => evc
       ELSE
          ikq = ikqs(ik)
       ENDIF
          
       CALL ch_psi_all_k()
  ENDIF
  
  DEALLOCATE (spsi)
  DEALLOCATE (hpsi)
  DEALLOCATE (ps)

  IF (tddfpt) NULLIFY(evq)
  dffts%have_task_groups=.FALSE.

  CALL stop_clock ('last')
  CALL stop_clock ('ch_psi')
  RETURN
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!K-point part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ch_psi_all_k()

    USE becmod, ONLY : becp, calbec
    
    IMPLICIT NONE
    !
    !   Here we compute the projector in the valence band
    !
    ps (:,:) = (0.d0, 0.d0)
    
    IF (noncolin) THEN
       CALL zgemm ('C', 'N', nbnd_occ (ikq) , m, npwx*npol, (1.d0, 0.d0) , evq, &
            npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
    ELSE
       CALL zgemm ('C', 'N', nbnd_occ (ikq) , m, n, (1.d0, 0.d0) , evq, &
            npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
    ENDIF
    ps (:,:) = ps(:,:) * alpha_pv
    CALL mp_sum ( ps, intra_bgrp_comm )
    
    hpsi (:,:) = (0.d0, 0.d0)
    IF (noncolin) THEN
       CALL zgemm ('N', 'N', npwx*npol, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
            npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
    ELSE
       CALL zgemm ('N', 'N', n, m, nbnd_occ (ikq) , (1.d0, 0.d0) , evq, &
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gamma part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ch_psi_all_gamma()
    
    USE becmod, ONLY : becp,  calbec
    USE realus, ONLY : real_space, real_space_debug, fft_orbital_gamma, &
         bfft_orbital_gamma, calbec_rs_gamma,  s_psir_gamma
    use gvect,                only : gstart

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
          CALL fft_orbital_gamma(hpsi,ibnd,m)
          CALL calbec_rs_gamma(ibnd,m,becp%r)
          CALL s_psir_gamma(ibnd,m)
          CALL bfft_orbital_gamma(spsi,ibnd,m)
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
