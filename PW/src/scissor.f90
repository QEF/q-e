! ----------------------------------------------
!
MODULE sci_mod
   !
   ! ... Module for scissor operator
   !     Written by Stefano Falletta
   !
   USE buffers,          ONLY : get_buffer, save_buffer
   USE ener,             ONLY : esci
   USE io_files,         ONLY : iunwfc, nwordwfc
   USE kinds,            ONLY : DP
   USE klist,            ONLY : nelec, nks, nkstot
   USE lsda_mod,         ONLY : isk, current_spin
   USE mp,               ONLY : mp_sum
   USE mp_bands,         ONLY : inter_bgrp_comm, intra_bgrp_comm, inter_bgrp_comm
   USE wavefunctions,    ONLY : evc
   USE wvfct,            ONLY : nbnd, current_k, wg, npwx
   USE io_global,        ONLY : stdout
   USE uspp,             ONLY : nkb, vkb, okvan
   USE sic_mod,          ONLY : isp, pol_type
   USE control_flags,    ONLY : sic
   ! 
   IMPLICIT NONE
   SAVE
   !
   REAL(DP) :: sci_vb, sci_cb             ! prefactor scissor operator 
   COMPLEX(DP), ALLOCATABLE :: evcc(:,:)  ! evc vector
   INTEGER  :: sci_iter = 0
   !
   CONTAINS
   !
   !---------------------------------
   SUBROUTINE allocate_scissor()
   !---------------------------------
      !
      IMPLICIT NONE
      !
      ALLOCATE(evcc(npwx,nbnd))
      !
   END SUBROUTINE allocate_scissor
   !
   !---------------------------------
   SUBROUTINE deallocate_scissor()
   !---------------------------------
      !
      IMPLICIT NONE
      !
      IF(ALLOCATED(evcc)) DEALLOCATE(evcc)
      !
   END SUBROUTINE deallocate_scissor
   !
   !---------------------------------
   SUBROUTINE p_psi(lda,n,m,psi,hpsi)
   !---------------------------------
     !  
     ! ... calculate the scissor operator with the current version of evc and
     !     apply it to psi. In scf calculations, skip the first iteration.
     !     in nscf calculations, no need for that
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN)        :: lda, n, m
     COMPLEX(DP), INTENT(IN)    :: psi(lda,m)
     COMPLEX(DP), INTENT(INOUT) :: hpsi(lda,m)
     COMPLEX(DP), ALLOCATABLE   :: coeff(:,:)
     INTEGER                    :: ibnd, ik, ibnd_p, ibnd_1, ibnd_2
     REAL(DP)                   :: fac
     REAL(DP), PARAMETER        :: ry2ev = 13.605698066
     !
     IF( sci_iter == 0 ) THEN
        CALL save_buffer ( evc, nwordwfc, iunwfc, current_k )
     ELSE
        ik = current_k
        esci = 0
        !
        !   ... move occupied states by sci_vb and unoccupied states by sci_cb
        !
        IF (.not. sic) THEN
           ALLOCATE(coeff(nbnd,m))
           CALL ZGEMM( 'C', 'N', nbnd, m, lda, (1.0_dp,0.0_dp), evcc, lda, psi, lda, (0.0_dp,0.0_dp), coeff, nbnd )
           DO ibnd = 1, nbnd
              fac = ( wg(ibnd,ik)*sci_vb + (1-wg(ibnd,ik))*sci_cb ) / ry2eV
              ibnd_p = nelec/2+1
              coeff(ibnd,:) = coeff(ibnd,:) * fac 
           END DO
           CALL mp_sum(coeff,inter_bgrp_comm) 
           CALL mp_sum(coeff,intra_bgrp_comm) 
           CALL ZGEMM( 'N', 'N', lda, m, nbnd, (1.0_dp,0.0_dp), evcc(:,1:nbnd), lda, coeff, nbnd, (1.0_dp,0.0_dp), hpsi, lda )
           DEALLOCATE(coeff)
           esci = -nelec*sci_vb/ry2eV
        END IF
        !
        ! ... for SIC calculations 
        ! 
        IF (sic) THEN
          IF (sci_vb .ne. 0.d0) THEN
             CALL vb_cb_indexes(ik,0,ibnd_1,ibnd_2)
             ALLOCATE(coeff(ibnd_2-ibnd_1+1,m))
             CALL ZGEMM( 'C', 'N', ibnd_2-ibnd_1+1, m, lda, (1.0_dp,0.0_dp), evcc(:,ibnd_1:ibnd_2), &
                          lda, psi(:,1:m), lda, (0.0_dp,0.0_dp), coeff, ibnd_2-ibnd_1+1 )
             CALL mp_sum(coeff,intra_bgrp_comm)
             CALL mp_sum(coeff,inter_bgrp_comm)
             coeff(:,:) = coeff(:,:) * sci_vb / ry2eV 
             CALL ZGEMM( 'N', 'N', lda, m, ibnd_2-ibnd_1+1, (1.0_dp,0.0_dp), evcc(:,ibnd_1:ibnd_2), &
                         lda, coeff, ibnd_2-ibnd_1+1, (1.0_dp,0.0_dp), hpsi, lda )
             DEALLOCATE(coeff)
             IF(pol_type == 'ep') esci = -(nelec-1)*sci_vb/ry2eV
             IF(pol_type == 'hp') esci = -nelec*sci_vb/ry2eV
          END IF
          IF (sci_cb .ne. 0.d0) THEN
             CALL vb_cb_indexes(ik,1,ibnd_1,ibnd_2)
             ALLOCATE(coeff(ibnd_2-ibnd_1+1,m))
             CALL ZGEMM( 'C', 'N', ibnd_2-ibnd_1+1, m, lda, (1.0_dp,0.0_dp), evcc(:,ibnd_1:ibnd_2), &
                         lda, psi(:,1:m), lda, (0.0_dp,0.0_dp), coeff, ibnd_2-ibnd_1+1 )
             CALL mp_sum(coeff,intra_bgrp_comm)
             CALL mp_sum(coeff,inter_bgrp_comm)
             coeff(:,:) = coeff(:,:) * sci_cb / ry2eV 
             CALL ZGEMM( 'N', 'N', lda, m, ibnd_2-ibnd_1+1, (1.0_dp,0.0_dp), evcc(:,ibnd_1:ibnd_2), &
                         lda, coeff, ibnd_2-ibnd_1+1, (1.0_dp,0.0_dp), hpsi, lda )
             DEALLOCATE(coeff)
          END IF
        END IF
     END IF
     !  
   END SUBROUTINE p_psi

   !-------------------------------------------------
   SUBROUTINE vb_cb_indexes(ik, band, ibnd_1, ibnd_2)
   !-------------------------------------------------
      !
      ! ... return start and end indexes of VB (band = 0), or CB (band=1)
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: ik, band
      INTEGER, INTENT(OUT) :: ibnd_1, ibnd_2 
      INTEGER :: is
      REAL(DP) :: ev2ry = 0.073498810939358
      !
      is =  isk(ik)
      IF(pol_type == 'e') THEN
         IF(band == 0) THEN
            ibnd_1 = 1
            ibnd_2 = nelec/2
         ELSEIF(band == 1) THEN
            IF (is == 1) ibnd_1 = nelec/2+2
            IF (is == 2) ibnd_1 = nelec/2+1
            ibnd_2 = nbnd
         END IF
      ELSEIF(pol_type == 'h') THEN
         IF(band == 0) THEN
            ibnd_1 = 1
            IF (is == 2) ibnd_2 = nelec/2
            IF (is == 1) ibnd_2 = nelec/2+1
         ELSEIF(band == 1) THEN
            ibnd_1 = nelec/2+2
            ibnd_2 = nbnd
         END IF
      END IF
      !
   END SUBROUTINE vb_cb_indexes
   ! 
END MODULE sci_mod
