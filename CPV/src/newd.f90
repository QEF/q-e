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
      USE uspp_param,       ONLY: nh, nhm, nvb
      USE uspp,             ONLY: deeq
      USE ions_base,        ONLY: nat, nsp, na
      USE constants,        ONLY: pi, fpi
      USE smallbox_gvec,            ONLY: ngb, npb, nmb, gxb
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
      INTEGER isup,isdw,iss, iv,ijv,jv, ik, nfft, isa, ia, is, ig
      REAL(DP)  fvan(3,nat,nvb), fac, fac1, fac2, boxdotgrid, res
      COMPLEX(DP) ci, facg1, facg2
      COMPLEX(DP), ALLOCATABLE :: qv(:)
      INTEGER :: na_bgrp, ia_bgrp
      EXTERNAL boxdotgrid

#if defined(__OPENMP)
      INTEGER :: itid, mytid, ntids, omp_get_thread_num, omp_get_num_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads
#endif

!
      IF ( dfftb%nr1==0 .OR. dfftb%nr2==0 .OR. dfftb%nr3==0 ) RETURN
      CALL start_clock( 'newd' )
      ci=(0.d0,1.d0)
      fac=omegab/DBLE(dfftb%nr1*dfftb%nr2*dfftb%nr3)
      deeq (:,:,:,:) = 0.d0
      fvan (:,:,:) = 0.d0


!$omp parallel default(none) &
!$omp          shared(nvb, na, ngb, nh, qgb, eigrb, dfftb, irb, vr, nmb, npb, ci, deeq, &
!$omp                 fac, nspin, my_bgrp_id, nbgrp ) &
!$omp          private(mytid, ntids, is, ia, nfft, iv, jv, ijv, ig, isa, qv, itid, res, iss )

      isa = 1

#if defined(__OPENMP)
      mytid = omp_get_thread_num()  ! take the thread ID
      ntids = omp_get_num_threads() ! take the number of threads
      itid  = 0
#endif


      ALLOCATE( qv( dfftb%nnr ) )
!
! calculation of deeq_i,lm = \int V_eff(r) q_i,lm(r) dr
!

      DO is = 1, nvb

#if defined(__MPI)

         DO ia=1,na(is)
             nfft = 1
             IF ( ( dfftb%np3( isa ) <= 0 ) .OR. ( my_bgrp_id /= MOD( ia, nbgrp ) ) ) THEN
                isa = isa + nfft
                CYCLE
             END IF
#else
         DO ia=1,na(is),2
            nfft=2
#endif

#if defined(__OPENMP)
            IF ( mytid /= itid ) THEN
               isa = isa + nfft
               itid = MOD( itid + 1, ntids )
               CYCLE
            ELSE
               itid = MOD( itid + 1, ntids )
            END IF
#endif

            IF( ia .EQ. na(is) ) nfft=1
!
! two ffts at the same time, on two atoms (if possible: nfft=2)
!
            DO iv=1,nh(is)
               DO jv=iv,nh(is)
                  ijv = (jv-1)*jv/2 + iv
                  qv(:) = (0.d0, 0.d0)
                  IF (nfft.EQ.2) THEN
                     DO ig=1,ngb
                        qv(npb(ig))= eigrb(ig,isa  )*qgb(ig,ijv,is)   &
     &                          + ci*eigrb(ig,isa+1)*qgb(ig,ijv,is)
                        qv(nmb(ig))= CONJG(                             &
     &                               eigrb(ig,isa  )*qgb(ig,ijv,is))  &
     &                          + ci*CONJG(                             &
     &                               eigrb(ig,isa+1)*qgb(ig,ijv,is))
                     END DO
                  ELSE
                     DO ig=1,ngb
                        qv(npb(ig)) = eigrb(ig,isa)*qgb(ig,ijv,is)
                        qv(nmb(ig)) = CONJG(                            &
     &                                eigrb(ig,isa)*qgb(ig,ijv,is))
                     END DO
                  END IF
!
                  CALL invfft( qv, dfftb, isa )
!
                  DO iss=1,nspin
                     res = boxdotgrid(irb(1,isa),1,qv,vr(1,iss))
                     deeq(iv,jv,isa,iss) = fac * res  
                     IF (iv.NE.jv)                                      &
     &                    deeq(jv,iv,isa,iss)=deeq(iv,jv,isa,iss)
                     IF (nfft.EQ.2) THEN
                        res = boxdotgrid(irb(1,isa+1),2,qv,vr(1,iss))
                        deeq(iv,jv,isa+1,iss) = fac * res 
                        IF (iv.NE.jv)                                   &
     &                       deeq(jv,iv,isa+1,iss)=deeq(iv,jv,isa+1,iss)
                     END IF
                  END DO
               END DO
            END DO
            isa=isa+nfft
         END DO
      END DO

      DEALLOCATE( qv )

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
!$omp          shared(nvb, na, ngb, nh, qgb, eigrb, dfftb, irb, vr, nmb, npb, ci, deeq, &
!$omp                 fac, nspin, rhovan, tpibab, gxb, fvan, my_bgrp_id, nbgrp ) &
!$omp          private(mytid, ntids, is, ia, ik, nfft, iv, jv, ijv, ig, isa, qv, itid, res, iss, &
!$omp                  fac1, fac2, facg1, facg2 )


         ALLOCATE( qv( dfftb%nnr ) )

         iss=1
         isa=1

#if defined(__OPENMP)
         mytid = omp_get_thread_num()  ! take the thread ID
         ntids = omp_get_num_threads() ! take the number of threads
         itid  = 0
#endif

         DO is = 1, nvb

#if defined(__MPI)
            DO ia=1,na(is)
               nfft=1
               IF ( ( dfftb%np3( isa ) <= 0 ) .OR. ( my_bgrp_id /= MOD( ia, nbgrp ) ) ) THEN
                  isa = isa + nfft
                  CYCLE
               END IF
#else
            DO ia=1,na(is),2
               nfft=2
#endif

#if defined(__OPENMP)
               IF ( mytid /= itid ) THEN
                  isa = isa + nfft
                  itid = MOD( itid + 1, ntids )
                  CYCLE
               ELSE
                  itid = MOD( itid + 1, ntids )
               END IF
#endif
               IF( ia.EQ.na(is)) nfft=1
               DO ik=1,3
                  qv(:) = (0.d0, 0.d0)
                  DO iv=1,nh(is)
                     DO jv=iv,nh(is)
                        ijv = (jv-1)*jv/2 + iv
                        IF(iv.NE.jv) THEN
                           fac1=2.d0*fac*tpibab*rhovan(ijv,isa,iss)
                           IF (nfft.EQ.2) fac2=2.d0*fac*tpibab*         &
     &                                           rhovan(ijv,isa+1,iss)
                        ELSE
                           fac1=     fac*tpibab*rhovan(ijv,isa,iss)
                           IF (nfft.EQ.2) fac2=     fac*tpibab*        &
     &                                           rhovan(ijv,isa+1,iss)
                        ENDIF
                        IF (nfft.EQ.2) THEN
                           DO ig=1,ngb
                              facg1 = CMPLX(0.d0,-gxb(ik,ig),kind=DP) * &
     &                                   qgb(ig,ijv,is) * fac1
                              facg2 = CMPLX(0.d0,-gxb(ik,ig),kind=DP) * &
     &                                   qgb(ig,ijv,is) * fac2
                              qv(npb(ig)) = qv(npb(ig))                 &
     &                                    +    eigrb(ig,isa  )*facg1    &
     &                                    + ci*eigrb(ig,isa+1)*facg2
                              qv(nmb(ig)) = qv(nmb(ig))                 &
     &                                +   CONJG(eigrb(ig,isa  )*facg1)&
     &                                +ci*CONJG(eigrb(ig,isa+1)*facg2)
                           END DO
                        ELSE
                           DO ig=1,ngb
                              facg1 = CMPLX(0.d0,-gxb(ik,ig),kind=DP) * &
     &                                   qgb(ig,ijv,is)*fac1
                              qv(npb(ig)) = qv(npb(ig))                 &
     &                                    +    eigrb(ig,isa)*facg1
                              qv(nmb(ig)) = qv(nmb(ig))                 &
     &                               +  CONJG( eigrb(ig,isa)*facg1)
                           END DO
                        END IF
                     END DO
                  END DO
!
                  CALL invfft( qv, dfftb, isa)
!
                  res = boxdotgrid(irb(1,isa),1,qv,vr(1,iss))
                  fvan(ik,ia,is) = res
!
                  IF (nfft.EQ.2) THEN
                     res = boxdotgrid(irb(1,isa+1),2,qv,vr(1,iss))
                     fvan(ik,ia+1,is) = res
                  END IF
               END DO
               isa = isa+nfft
            END DO
         END DO

         DEALLOCATE( qv )

!$omp end parallel

      ELSE

         !     =================================================================
         !     case nspin=2: up and down spin fft's combined into a single fft
         !     -----------------------------------------------------------------

         ALLOCATE( qv( dfftb%nnr ) )
         isup=1
         isdw=2
         isa=1
         DO is=1,nvb
            DO ia=1,na(is)
#if defined(__MPI)
               IF ( dfftb%np3( isa ) <= 0 ) go to 25
#endif
               DO ik=1,3
                  qv(:) = (0.d0, 0.d0)
!
                  DO iv=1,nh(is)
                     DO jv=iv,nh(is)
                        ijv = (jv-1)*jv/2 + iv
                        IF(iv.NE.jv) THEN
                           fac1=2.d0*fac*tpibab*rhovan(ijv,isa,isup)
                           fac2=2.d0*fac*tpibab*rhovan(ijv,isa,isdw)
                        ELSE
                           fac1=     fac*tpibab*rhovan(ijv,isa,isup)
                           fac2=     fac*tpibab*rhovan(ijv,isa,isdw)
                        END IF
                        DO ig=1,ngb
                           facg1 = fac1 * CMPLX(0.d0,-gxb(ik,ig),kind=DP) * &
     &                                qgb(ig,ijv,is) * eigrb(ig,isa)
                           facg2 = fac2 * CMPLX(0.d0,-gxb(ik,ig),kind=DP) * &
     &                                qgb(ig,ijv,is) * eigrb(ig,isa)
                           qv(npb(ig)) = qv(npb(ig))                    &
     &                                    + facg1 + ci*facg2
                           qv(nmb(ig)) = qv(nmb(ig))                    &
     &                                    +CONJG(facg1)+ci*CONJG(facg2)
                        END DO
                     END DO
                  END DO
!
                  CALL invfft( qv, dfftb, isa)
!
                  fvan(ik,ia,is) =                                      &
     &                    boxdotgrid(irb(1,isa),isup,qv,vr(1,isup)) + &
     &                    boxdotgrid(irb(1,isa),isdw,qv,vr(1,isdw))
               END DO
25             isa = isa+1
            END DO
         END DO

         DEALLOCATE( qv )

      END IF

      CALL mp_sum( fvan, intra_bgrp_comm )
      CALL mp_sum( fvan, inter_bgrp_comm )

      isa = 0
      DO is = 1, nvb
        DO ia = 1, na(is)
          isa = isa + 1
          fion(:,isa) = fion(:,isa) - fvan(:,ia,is)
        END DO
      END DO

  10  CONTINUE
!
      CALL stop_clock( 'newd' )
!
      RETURN
      END SUBROUTINE newd

