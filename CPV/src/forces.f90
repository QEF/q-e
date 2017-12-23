!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written and revised by Carlo Cavazzoni
! Task Groups parallelization by C. Bekas (IBM Research Zurich).
!
!-------------------------------------------------------------------------
      SUBROUTINE dforce_x ( i, bec, vkb, c, df, da, v, ldv, ispin, f, n, nspin, v1 )
!-----------------------------------------------------------------------
!computes: the generalized force df=cmplx(dfr,dfi,kind=DP) acting on the i-th
!          electron state at the gamma point of the brillouin zone
!          represented by the vector c=cmplx(cr,ci,kind=DP)
!
!     d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
!              sum_i,ij d^q_i,ij (-i)**l beta_i,i(g) 
!                                 e^-ig.r_i < beta_i,j | c_n >}
!
      USE parallel_include
      USE kinds,                  ONLY: dp
      USE control_flags,          ONLY: iprint
      USE uspp,                   ONLY: nhsa=>nkb, dvan, deeq
      USE uspp_param,             ONLY: nhm, nh, ish
      USE constants,              ONLY: pi, fpi
      USE ions_base,              ONLY: nsp, na, nat
      USE gvecw,                  ONLY: ngw, g2kin
      USE cell_base,              ONLY: tpiba2
      USE ensemble_dft,           ONLY: tens
      USE funct,                  ONLY: dft_is_meta, dft_is_hybrid, exx_is_active
      USE fft_base,               ONLY: dffts
      USE fft_interfaces,         ONLY: fwfft, invfft
      USE mp_global,              ONLY: me_bgrp
      USE control_flags,          ONLY: lwfpbe0nscf
      USE exx_module,             ONLY: exx_potential
      USE fft_helper_subroutines
!
      IMPLICIT NONE
!
      INTEGER,     INTENT(IN)    :: i
      REAL(DP)                   :: bec(:,:)
      COMPLEX(DP)                :: vkb(:,:)
      COMPLEX(DP)                :: c(:,:)
      COMPLEX(DP)                :: df(:), da(:)
      INTEGER,     INTENT(IN)    :: ldv
      REAL(DP)                   :: v( ldv, * )
      INTEGER                    :: ispin( : )
      REAL(DP)                   :: f( : )
      INTEGER,     INTENT(IN)    :: n, nspin
      REAL(DP),    OPTIONAL      :: v1( ldv, * )
      !
      ! local variables
      !
      INTEGER     :: iv, jv, ia, is, isa, ism, ios, iss1, iss2, ir, ig, inl, jnl
      INTEGER     :: ivoff, jvoff, igoff, igno, igrp, ierr
      INTEGER     :: idx, eig_offset, nogrp_ , inc, tg_nr3
      REAL(DP)    :: fi, fip, dd, dv
      COMPLEX(DP) :: fp, fm, ci
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: af, aa, psi, exx_a, exx_b
#endif
#endif
      REAL(DP),    ALLOCATABLE :: af( :, : ), aa( :, : )
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      REAL(DP)    :: tmp1, tmp2                      ! Lingzhu Kong
      REAL(DP),    ALLOCATABLE :: exx_a(:), exx_b(:) ! Lingzhu Kong      
      !
      CALL start_clock( 'dforce' ) 
      !
!=======================================================================
!exx_wf related
      IF(dft_is_hybrid().AND.exx_is_active()) THEN
         allocate( exx_a( dffts%nnr ) ); exx_a=0.0_DP
         allocate( exx_b( dffts%nnr ) ); exx_b=0.0_DP
      END IF
!=======================================================================

      nogrp_ = fftx_ntgrp(dffts)
      ALLOCATE( psi( dffts%nnr_tg ) )
      !
      ci = ( 0.0d0, 1.0d0 )
      !
#if defined(__MPI)


!$omp  parallel
!$omp  single

      DO idx = 1, 2*nogrp_ , 2

!$omp task default(none) &
!$omp          firstprivate( idx, i, n, ngw, ci, nogrp_ ) &
!$omp          private( igoff, ig ) &
!$omp          shared( c, dffts, psi )
         !
         !  This loop is executed only ONCE when NOGRP=1.
         !  Equivalent to the case with no task-groups
         !  dfft%nsw(me) holds the number of z-sticks for the current processor per wave-function
         !  We can either send these in the group with an mpi_allgather...or put the
         !  in the PSIS vector (in special positions) and send them with them.
         !  Otherwise we can do this once at the beginning, before the loop.
         !  we choose to do the latter one.
         !
         !  important: if n is odd => c(*,n+1)=0.
         ! 
         IF ( ( idx + i - 1 ) == n ) c( : , idx + i ) = 0.0d0

         igoff = ( idx - 1 )/2 * dffts%nnr 

         psi( igoff + 1 : igoff + dffts%nnr ) = (0.d0, 0.d0)

         IF( idx + i - 1 <= n ) THEN
            DO ig=1,ngw
               psi(dffts%nlm(ig)+igoff) = conjg( c(ig,idx+i-1) - ci * c(ig,idx+i) )
               psi(dffts%nl(ig)+igoff) =         c(ig,idx+i-1) + ci * c(ig,idx+i)
            END DO
         END IF
!$omp end task

      END DO

!$omp  end single
!$omp  end parallel

      CALL invfft('tgWave', psi, dffts)

#else

      psi = 0.0d0
      DO ig=1,ngw
         psi(dffts%nlm(ig)) = conjg( c(ig,i) - ci * c(ig,i+1) )
         psi(dffts%nl(ig)) =  c(ig,i) + ci * c(ig,i+1)
      END DO
      CALL invfft( 'Wave', psi, dffts )

#endif
      !
      ! the following avoids a potential out-of-bounds error
      !
      IF ( i < n ) THEN
         iss1 = ispin(i)
         iss2 = ispin(i+1)
      ELSE
         iss1 = ispin(i)
         iss2 = iss1
      END IF
      !
      IF( dffts%have_task_groups ) THEN
         !
         CALL tg_get_group_nr3( dffts, tg_nr3 )
         !
!===============================================================================
!exx_wf related
         IF(dft_is_hybrid().AND.exx_is_active()) THEN
            !$omp parallel do private(tmp1,tmp2) 
            DO ir = 1, dffts%nr1x*dffts%nr2x*tg_nr3
               tmp1 = v(ir,iss1) * DBLE( psi(ir) )+exx_potential(ir,i/nogrp_+1)
               tmp2 = v(ir,iss2) * AIMAG(psi(ir) )+exx_potential(ir,i/nogrp_+2)
               psi(ir) = CMPLX( tmp1, tmp2, kind=DP)
            END DO
            !$omp end parallel do 
         ELSE
            !$omp parallel do 
            DO ir = 1, dffts%nr1x*dffts%nr2x*tg_nr3
               psi(ir) = CMPLX ( v(ir,iss1) * DBLE( psi(ir) ), &
                                 v(ir,iss2) *AIMAG( psi(ir) ) ,kind=DP)
            END DO
            !$omp end parallel do 
         ENDIF
!===============================================================================
         !
      ELSE
         !
         IF( PRESENT( v1 ) ) THEN
!===============================================================================
!exx_wf related
            IF(dft_is_hybrid().AND.exx_is_active()) THEN
               !
               IF ( (mod(n,2).ne.0 ) .and. (i.eq.n) ) THEN
                 !
                 !$omp parallel do 
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%my_nr3p
                   exx_a(ir) = exx_potential(ir, i)
                   exx_b(ir) = 0.0_DP
                 END DO
                 !$omp end parallel do 
                 !
               ELSE
                 !
                 !$omp parallel do 
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%my_nr3p
                   exx_a(ir) = exx_potential(ir, i)
                   exx_b(ir) = exx_potential(ir, i+1)
                 END DO
                 !$omp end parallel do 
                 !
               ENDIF
               !$omp parallel do private(tmp1,tmp2) 
               DO ir=1,dffts%nnr
                  tmp1 =  v(ir,iss1)* DBLE(psi(ir))+exx_a(ir)
                  tmp2 = v1(ir,iss2)*AIMAG(psi(ir))+exx_b(ir)
                  psi(ir)=CMPLX( tmp1, tmp2, kind=DP )
               END DO
               !$omp end parallel do 
               !
            ELSE
               !
               !$omp parallel do 
               DO ir=1,dffts%nnr
                  psi(ir)=CMPLX ( v(ir,iss1)* DBLE(psi(ir)), &
                                 v1(ir,iss2)*AIMAG(psi(ir)) ,kind=DP)
               END DO
               !$omp end parallel do 
               !
            ENDIF
         ELSE
!===============================================================================
!exx_wf related
            IF(dft_is_hybrid().AND.exx_is_active()) THEN
               IF ( (mod(n,2).ne.0 ) .and. (i.eq.n) ) THEN
                 !$omp parallel do 
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%my_nr3p
                   exx_a(ir) = exx_potential(ir, i)
                   exx_b(ir) = 0.0_DP
                 END DO
                 !$omp end parallel do 
               ELSE
                 !$omp parallel do 
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%my_nr3p
                   exx_a(ir) = exx_potential(ir, i)
                   exx_b(ir) = exx_potential(ir, i+1)
                 END DO
                 !$omp end parallel do 
               ENDIF
               !$omp parallel do private(tmp1,tmp2) 
               DO ir=1,dffts%nnr
                  tmp1 = v(ir,iss1)* DBLE(psi(ir))+exx_a(ir)
                  tmp2 = v(ir,iss2)*AIMAG(psi(ir))+exx_b(ir)
                  psi(ir)=CMPLX( tmp1, tmp2, kind=DP )
               END DO
               !$omp end parallel do 
!===============================================================================
            ELSE
               !$omp parallel do 
               DO ir=1,dffts%nnr
                  psi(ir)=CMPLX( v(ir,iss1)* DBLE(psi(ir)), &
                                 v(ir,iss2)*AIMAG(psi(ir)) ,kind=DP)
               END DO
               !$omp end parallel do 
            ENDIF
         END IF
         !
      END IF
      !
#if defined(__MPI)
      CALL fwfft( 'tgWave', psi, dffts )
#else
      CALL fwfft( 'Wave', psi, dffts )
#endif
      !
      !   note : the factor 0.5 appears 
      !       in the kinetic energy because it is defined as 0.5*g**2
      !       in the potential part because of the logics
      !
      !   Each processor will treat its own part of the eigenstate
      !   assigned to its ORBITAL group
      !
!$omp  parallel
!$omp  single

      eig_offset = 0
      CALL tg_get_recip_inc( dffts, inc )
      igno = 1

      DO idx = 1, 2*nogrp_ , 2

!$omp task default(none)  &
!$omp          private( fi, fip, fp, fm, ig ) &
!$omp          firstprivate( eig_offset, igno, idx, nogrp_, ngw, tpiba2, me_bgrp, i, n, tens ) &
!$omp          shared( f, psi, df, da, c, dffts, g2kin  )

         IF( idx + i - 1 <= n ) THEN
            if (tens) then
               fi = -0.5d0
               fip = -0.5d0
            else
               fi = -0.5d0*f(i+idx-1)
               fip = -0.5d0*f(i+idx)
            endif
            IF( dffts%have_task_groups ) THEN
               DO ig=1,ngw
                  fp= psi(dffts%nl(ig)+eig_offset) +  psi(dffts%nlm(ig)+eig_offset)
                  fm= psi(dffts%nl(ig)+eig_offset) -  psi(dffts%nlm(ig)+eig_offset)
                  df(ig+igno-1)= fi *(tpiba2 * g2kin(ig) * c(ig,idx+i-1) + &
                                 CMPLX(real (fp), aimag(fm), kind=dp ))
                  da(ig+igno-1)= fip*(tpiba2 * g2kin(ig) * c(ig,idx+i  ) + &
                                 CMPLX(aimag(fp),-real (fm), kind=dp ))
               END DO
            ELSE
               DO ig=1,ngw
                  fp= psi(dffts%nl(ig)) + psi(dffts%nlm(ig))
                  fm= psi(dffts%nl(ig)) - psi(dffts%nlm(ig))
                  df(ig)= fi*(tpiba2*g2kin(ig)* c(ig,idx+i-1)+CMPLX(DBLE(fp), AIMAG(fm),kind=DP))
                  da(ig)=fip*(tpiba2*g2kin(ig)* c(ig,idx+i  )+CMPLX(AIMAG(fp),-DBLE(fm),kind=DP))
               END DO
            END IF
         END IF
!$omp end task

         igno = igno + ngw
         eig_offset = eig_offset + inc

         ! We take into account the number of elements received from other members of the orbital group

      ENDDO

!$omp end single
!$omp end parallel 
      !
      IF(dft_is_meta()) THEN
         ! HK/MCA : warning on task groups
         if (nogrp_.gt.1) call errore('forces','metagga force not supporting taskgroup parallelization',1)
         ! HK/MCA : reset occupation numbers since omp private screws it up... need a better fix FIXME
         if (tens) then
            fi = -0.5d0
            fip = -0.5d0
         else
            fi = -0.5d0*f(i)
            fip = -0.5d0*f(i+1)
         endif
         CALL dforce_meta(c(1,i),c(1,i+1),df,da,psi,iss1,iss2,fi,fip) !METAGGA
      END IF


      IF( nhsa > 0 ) THEN
         !
         !     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
         ! 
         ALLOCATE( af( nhsa, nogrp_ ), aa( nhsa, nogrp_ ) )

         af = 0.0d0
         aa = 0.0d0
         !
!$omp parallel
!$omp single
         !
         igrp = 1

         DO idx = 1, 2*nogrp_ , 2

!$omp task default(none) &
!$omp          firstprivate(igrp,idx, nogrp_, ngw, i, n, nsp, na, nh, ish, iss1, iss2, tens ) &
!$omp          private(iv,jv,ivoff,jvoff,dd,dv,inl,jnl,is,isa,ism,fi,fip) &
!$omp          shared( f, deeq, bec, af, aa, dvan )

            IF( idx + i - 1 <= n ) THEN

               IF (tens) THEN
                  fi = 1.0d0
                  fip= 1.0d0
               ELSE
                  fi = f(i+idx-1)
                  fip= f(i+idx)
               END IF
               !
               DO is = 1, nsp
                  DO iv = 1, nh(is)
                        DO jv = 1, nh(is)
                           isa = 0
                           DO ism = 1, is-1
                              isa = isa + na( ism )
                           END DO
                           dv = dvan(iv,jv,is)
                           ivoff = ish(is)+(iv-1)*na(is)
                           jvoff = ish(is)+(jv-1)*na(is)
                           IF( i + idx - 1 /= n ) THEN
                              DO ia=1,na(is)
                                 inl = ivoff + ia
                                 jnl = jvoff + ia
                                 dd = deeq(iv,jv,isa+ia,iss1) + dv
                                 af(inl,igrp) = af(inl,igrp) - fi  * dd * bec(jnl,i+idx-1)
                                 dd = deeq(iv,jv,isa+ia,iss2) + dv
                                 aa(inl,igrp) = aa(inl,igrp) - fip * dd * bec(jnl,i+idx)
                              END DO
                           ELSE
                              DO ia=1,na(is)
                                 inl = ivoff + ia
                                 jnl = jvoff + ia
                                 dd = deeq(iv,jv,isa+ia,iss1) + dv
                                 af(inl,igrp) = af(inl,igrp) - fi * dd * bec(jnl,i+idx-1)
                              END DO
                           END IF
                        END DO
                  END DO
               END DO

            END IF

!$omp end task

            igrp = igrp + 1

         END DO

!$omp end single
!$omp end parallel

         IF( ngw > 0 ) THEN
           CALL dgemm ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, af, nhsa, 1.0d0, df, 2*ngw)
           CALL dgemm ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, aa, nhsa, 1.0d0, da, 2*ngw)
         END IF
         !
         DEALLOCATE( aa, af )
         !
      ENDIF
!
      IF(dft_is_hybrid().AND.exx_is_active()) DEALLOCATE(exx_a, exx_b)
      DEALLOCATE( psi )
!
      CALL stop_clock( 'dforce' ) 
!
      RETURN
   END SUBROUTINE dforce_x
