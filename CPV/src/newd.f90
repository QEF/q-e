!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
   SUBROUTINE newd(vr,irb,eigrb,rhovan,fion)
!-----------------------------------------------------------------------
!
!     this routine calculates array deeq:
!         deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!     and the corresponding term in forces
!         fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!     where
!         rho_lm = \sum_j f_j <psi_j|beta_l><beta_m|psi_j>
!
      USE kinds,            ONLY: dp
      USE uspp_param,       ONLY: nh, nhm, upf
      USE uspp,             ONLY: deeq
      USE ions_base,        ONLY: nat, nsp, na, ityp
      USE constants,        ONLY: pi, fpi
      USE smallbox_gvec,    ONLY: ngb, gxb
      USE smallbox_subs,    ONLY: fft_oned2box, fft_add_oned2box, boxdotgrid
      USE small_box,        ONLY: omegab, tpibab
      USE qgb_mod,          ONLY: qgb
      USE electrons_base,   ONLY: nspin
      USE control_flags,    ONLY: iprint, thdyn, tfor, tprnfor
      USE mp,               ONLY: mp_sum
      USE mp_bands,         ONLY: intra_bgrp_comm, inter_bgrp_comm, &
                                  my_bgrp_id, nbgrp 
      USE fft_interfaces,   ONLY: invfft
      USE fft_base,         ONLY: dfftb, dfftp
!
      IMPLICIT NONE
! input
      INTEGER irb(3,nat)
      REAL(DP) rhovan(nhm*(nhm+1)/2,nat,nspin)
      COMPLEX(DP) eigrb(ngb,nat)
      REAL(DP)  vr(dfftp%nnr,nspin)
! output
      REAL(DP)  fion(3,nat)
! local
      INTEGER isup,isdw,iss, iv,ijv,jv, ik, nfft, ia, is, ig
      REAL(DP)  fac, fac1, fac2, res
      COMPLEX(DP) ci, facg1, facg2
      COMPLEX(DP), ALLOCATABLE :: qv(:), fg1(:), fg2(:)
      REAL(DP), ALLOCATABLE :: fvan(:,:)
      INTEGER :: na_bgrp, ia_bgrp

#if defined(_OPENMP)
      INTEGER :: itid, mytid, ntids, omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif

!
      IF ( dfftb%nr1==0 .OR. dfftb%nr2==0 .OR. dfftb%nr3==0 ) THEN
         RETURN
      END IF

      CALL start_clock( 'newd' )

      ALLOCATE( fvan( 3, nat ) )
      ci=(0.d0,1.d0)
      fac=omegab/DBLE(dfftb%nr1*dfftb%nr2*dfftb%nr3)
      deeq (:,:,:,:) = 0.d0
      fvan (:,:) = 0.d0


!$omp parallel default(none) &
!$omp          shared(ngb, nh, qgb, eigrb, dfftb, irb, vr, ci, deeq, &
!$omp                 fac, nspin, my_bgrp_id, nbgrp, ityp, upf, nat ) &
!$omp          private(mytid, ntids, is, ia, nfft, iv, jv, ijv, ig, qv, fg1, fg2, itid, res, iss )

#if defined(_OPENMP)
      mytid = omp_get_thread_num()  ! take the thread ID
      ntids = omp_get_num_threads() ! take the number of threads
      itid  = 0
#endif


      ALLOCATE( qv( dfftb%nnr ) )
      ALLOCATE( fg1( ngb ) )
      ALLOCATE( fg2( ngb ) )
!
! calculation of deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!
      DO ia = 1, nat

         is = ityp(ia)

         IF( .NOT. upf(is)%tvanp ) &
            CYCLE

#if defined(__MPI)
         IF(  ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2( ia ) <= 0 ) .OR. ( my_bgrp_id /= MOD( ia, nbgrp ) ) ) THEN
            CYCLE
         END IF
#endif

         nfft = 1

#if defined(_OPENMP)
         IF ( mytid /= itid ) THEN
            itid = MOD( itid + 1, ntids )
            CYCLE
         ELSE
            itid = MOD( itid + 1, ntids )
         END IF
#endif

         DO iv=1,nh(is)
             DO jv=iv,nh(is)
                ijv = (jv-1)*jv/2 + iv
                IF (nfft.EQ.2) THEN
                   fg1 = eigrb(1:ngb,ia  )*qgb(1:ngb,ijv,is)
                   fg2 = eigrb(1:ngb,ia+1)*qgb(1:ngb,ijv,is)
                   CALL fft_oned2box( qv, fg1, fg2 )
                ELSE
                   fg1 = eigrb(1:ngb,ia  )*qgb(1:ngb,ijv,is)
                   CALL fft_oned2box( qv, fg1 )
                END IF
!
                CALL invfft( qv, dfftb, ia )
!
                DO iss=1,nspin
                   res = boxdotgrid(irb(:,ia),1,qv,vr(:,iss))
                   deeq(iv,jv,ia,iss) = fac * res  
                   IF (iv.NE.jv)                                      &
     &                  deeq(jv,iv,ia,iss)=deeq(iv,jv,ia,iss)
                   IF (nfft.EQ.2) THEN
                      res = boxdotgrid(irb(:,ia+1),2,qv,vr(:,iss))
                      deeq(iv,jv,ia+1,iss) = fac * res 
                      IF (iv.NE.jv)                                   &
     &                     deeq(jv,iv,ia+1,iss)=deeq(iv,jv,ia+1,iss)
                   END IF
                END DO
             END DO
         END DO
      END DO

      DEALLOCATE( qv )
      DEALLOCATE( fg1 )
      DEALLOCATE( fg2 )

!$omp end parallel

      CALL mp_sum( deeq, intra_bgrp_comm )
      CALL mp_sum( deeq, inter_bgrp_comm )

      IF (.NOT.( tfor .OR. thdyn .OR. tprnfor ) ) go to 10
!
! calculation of fion_i = \int V_eff(r) \sum_lm rho_lm (dq_i,lm(r)/dR_i) dr
!

      IF( nspin == 1 ) THEN

         !     =================================================================
         !     case nspin=1: two ffts at the same time, on two atoms (if possible)
         !     -----------------------------------------------------------------

!$omp parallel default(none) &
!$omp          shared(nat, ngb, nh, qgb, eigrb, dfftb, irb, vr, ci, deeq, &
!$omp                 fac, nspin, rhovan, tpibab, gxb, fvan, my_bgrp_id, nbgrp, ityp, upf ) &
!$omp          private(mytid, ntids, is, ia, ik, nfft, iv, jv, ijv, ig, qv, itid, res, iss, &
!$omp                  fac1, fac2, facg1, facg2, fg1, fg2 )


         ALLOCATE( qv( dfftb%nnr ) )
         ALLOCATE( fg1( ngb ) )
         ALLOCATE( fg2( ngb ) )

         iss=1

#if defined(_OPENMP)
         mytid = omp_get_thread_num()  ! take the thread ID
         ntids = omp_get_num_threads() ! take the number of threads
         itid  = 0
#endif

         DO ia = 1, nat

            is = ityp(ia)

            IF( .NOT. upf(is)%tvanp ) &
               CYCLE

#if defined(__MPI)
            IF ( ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2( ia ) <= 0 ) .OR. ( my_bgrp_id /= MOD( ia, nbgrp ) ) ) THEN
               CYCLE
            END IF
#endif
            nfft=1

#if defined(_OPENMP)
            IF ( mytid /= itid ) THEN
               itid = MOD( itid + 1, ntids )
               CYCLE
            ELSE
               itid = MOD( itid + 1, ntids )
            END IF
#endif
            DO ik=1,3
               qv(:) = (0.d0, 0.d0)
               DO iv=1,nh(is)
                  DO jv=iv,nh(is)
                     ijv = (jv-1)*jv/2 + iv
                     IF(iv.NE.jv) THEN
                        fac1=2.d0*fac*tpibab*rhovan(ijv,ia,iss)
                        IF (nfft.EQ.2) fac2=2.d0*fac*tpibab*         &
     &                                        rhovan(ijv,ia+1,iss)
                     ELSE
                        fac1=     fac*tpibab*rhovan(ijv,ia,iss)
                        IF (nfft.EQ.2) fac2=     fac*tpibab*        &
     &                                        rhovan(ijv,ia+1,iss)
                     ENDIF
                     IF (nfft.EQ.2) THEN
                        DO ig=1,ngb
                           facg1 = CMPLX(0.d0,-gxb(ik,ig),kind=DP) * qgb(ig,ijv,is) * fac1
                           facg2 = CMPLX(0.d0,-gxb(ik,ig),kind=DP) * qgb(ig,ijv,is) * fac2
                           fg1(ig) = eigrb(ig,ia  )*facg1
                           fg2(ig) = eigrb(ig,ia+1)*facg2
                        END DO
                        CALL fft_add_oned2box( qv, fg1, fg2 )
                     ELSE
                        DO ig=1,ngb
                           facg1 = CMPLX(0.d0,-gxb(ik,ig),kind=DP) * qgb(ig,ijv,is)*fac1
                           fg1(ig) = eigrb(ig,ia  )*facg1
                        END DO
                        CALL fft_add_oned2box( qv, fg1 )
                     END IF
                  END DO
               END DO
!
               CALL invfft( qv, dfftb, ia)
!
               res = boxdotgrid(irb(:,ia),1,qv,vr(:,iss))
               fvan(ik,ia) = res
!
               IF (nfft.EQ.2) THEN
                  res = boxdotgrid(irb(:,ia+1),2,qv,vr(:,iss))
                  fvan(ik,ia+1) = res
               END IF
            END DO
         END DO

         DEALLOCATE( qv )
         DEALLOCATE( fg1 )
         DEALLOCATE( fg2 )

!$omp end parallel

      ELSE

         !     =================================================================
         !     case nspin=2: up and down spin fft's combined into a single fft
         !     -----------------------------------------------------------------

         ALLOCATE( qv( dfftb%nnr ) )
         ALLOCATE( fg1( ngb ) )
         ALLOCATE( fg2( ngb ) )
         isup=1
         isdw=2

         DO ia = 1, nat

            is = ityp(ia)

            IF( .NOT. upf(is)%tvanp ) &
               CYCLE

#if defined(__MPI)
            IF ( ( dfftb%np3( ia ) <= 0 ) .OR. ( dfftb%np2( ia ) <= 0 ) .OR. ( my_bgrp_id /= MOD( ia, nbgrp ) ) ) THEN
               CYCLE
            END IF
#endif

            DO ik=1,3
                qv(:) = (0.d0, 0.d0)
!
                DO iv=1,nh(is)
                   DO jv=iv,nh(is)
                      ijv = (jv-1)*jv/2 + iv
                      IF(iv.NE.jv) THEN
                         fac1=2.d0*fac*tpibab*rhovan(ijv,ia,isup)
                         fac2=2.d0*fac*tpibab*rhovan(ijv,ia,isdw)
                      ELSE
                         fac1=     fac*tpibab*rhovan(ijv,ia,isup)
                         fac2=     fac*tpibab*rhovan(ijv,ia,isdw)
                      END IF
                      DO ig=1,ngb
                         fg1(ig) = fac1 * CMPLX(0.d0,-gxb(ik,ig),kind=DP) * &
     &                             qgb(ig,ijv,is) * eigrb(ig,ia)
                         fg2(ig) = fac2 * CMPLX(0.d0,-gxb(ik,ig),kind=DP) * &
     &                             qgb(ig,ijv,is) * eigrb(ig,ia)
                      END DO
                      CALL fft_add_oned2box( qv, fg1, fg2 )
                   END DO
                END DO
!
                CALL invfft( qv, dfftb, ia)
!
                fvan(ik,ia) =                                      &
     &                  boxdotgrid(irb(:,ia),isup,qv,vr(:,isup)) + &
     &                  boxdotgrid(irb(:,ia),isdw,qv,vr(:,isdw))
            END DO
         END DO

         DEALLOCATE( qv )
         DEALLOCATE( fg1 )
         DEALLOCATE( fg2 )

      END IF

      CALL mp_sum( fvan, intra_bgrp_comm )
      CALL mp_sum( fvan, inter_bgrp_comm )

      fion(:,:) = fion(:,:) - fvan(:,:)

  10  CONTINUE
!
      CALL stop_clock( 'newd' )
!
      RETURN
      END SUBROUTINE newd

