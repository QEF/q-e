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
      USE gvecs,                  ONLY: nlsm, nls
      USE uspp,                   ONLY: nhsa=>nkb, dvan, deeq
      USE uspp_param,             ONLY: nhm, nh, ish
      USE constants,              ONLY: pi, fpi
      USE ions_base,              ONLY: nsp, na, nat
      USE gvecw,                  ONLY: ngw, ggp
      USE cell_base,              ONLY: tpiba2
      USE ensemble_dft,           ONLY: tens
      USE funct,                  ONLY: dft_is_meta, dft_is_hybrid, exx_is_active
      USE fft_base,               ONLY: dffts
      USE fft_interfaces,         ONLY: fwfft, invfft
      USE mp_global,              ONLY: me_bgrp
      USE control_flags,          ONLY: lwfpbe0nscf
      USE exx_module,             ONLY: exx_potential,exxalfa
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
      INTEGER     :: idx, eig_offset, eig_index, nogrp_
      REAL(DP)    :: fi, fip, dd, dv
      COMPLEX(DP) :: fp, fm, ci
      REAL(DP),    ALLOCATABLE :: af( :, : ), aa( :, : )
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      REAL(DP)    :: tmp1, tmp2                      
      REAL(DP),    ALLOCATABLE :: exx_a(:), exx_b(:)       
      !
      CALL start_clock( 'dforce' ) 
      !
!=======================================================================
!exx_wf related
      !
      !Note that the mixing parameter exxalfa is multiplied here ...
      !
      IF(dft_is_hybrid().AND.exx_is_active()) THEN
         exx_potential=exx_potential*exxalfa 
      END IF
      !
      IF(dft_is_hybrid().AND.exx_is_active()) THEN
         allocate( exx_a( dffts%nnr ) ); exx_a=0.0_DP
         allocate( exx_b( dffts%nnr ) ); exx_b=0.0_DP
      END IF
!=======================================================================
      IF( dffts%have_task_groups ) THEN
         nogrp_ = dffts%nogrp
         ALLOCATE( psi( dffts%tg_nnr * dffts%nogrp ) )
      ELSE
         nogrp_ = 1
         ALLOCATE( psi( dffts%nnr ) )
      END IF
      !
      ci = ( 0.0d0, 1.0d0 )
      !
      psi( : ) = (0.d0, 0.d0)

      igoff = 0

      DO idx = 1, 2*nogrp_ , 2
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

         IF( idx + i - 1 <= n ) THEN
            !$omp parallel do 
            DO ig=1,ngw
               psi(nlsm(ig)+igoff) = conjg( c(ig,idx+i-1) - ci * c(ig,idx+i) )
               psi(nls(ig)+igoff) =        c(ig,idx+i-1) + ci * c(ig,idx+i)
            END DO
            !$omp end parallel do 
         END IF

         igoff = igoff + dffts%tg_nnr

      END DO

      CALL invfft( 'Wave', psi, dffts )
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
!===============================================================================
!exx_wf related
         IF(dft_is_hybrid().AND.exx_is_active()) THEN
            !$omp parallel do private(tmp1,tmp2) 
            DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
               tmp1 = v(ir,iss1) * DBLE( psi(ir) )+exx_potential(ir,i/nogrp_+1)
               tmp2 = v(ir,iss2) * AIMAG(psi(ir) )+exx_potential(ir,i/nogrp_+2)
               psi(ir) = CMPLX( tmp1, tmp2, kind=DP)
            END DO
            !$omp end parallel do 
         ELSE
            !$omp parallel do 
            DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%tg_npp( me_bgrp + 1 )
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
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%npp( me_bgrp + 1 )
                   exx_a(ir) = exx_potential(ir, i)
                   exx_b(ir) = 0.0_DP
                 END DO
                 !$omp end parallel do 
                 !
               ELSE
                 !
                 !$omp parallel do 
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%npp( me_bgrp + 1 )
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
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%npp( me_bgrp + 1 )
                   exx_a(ir) = exx_potential(ir, i)
                   exx_b(ir) = 0.0_DP
                 END DO
                 !$omp end parallel do 
               ELSE
                 !$omp parallel do 
                 DO ir = 1, dffts%nr1x*dffts%nr2x*dffts%npp( me_bgrp + 1 )
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
      CALL fwfft( 'Wave', psi, dffts ) 
      !
      !   note : the factor 0.5 appears 
      !       in the kinetic energy because it is defined as 0.5*g**2
      !       in the potential part because of the logics
      !
      !   Each processor will treat its own part of the eigenstate
      !   assigned to its ORBITAL group
      !
!$omp parallel default(none) &
!$omp          private( eig_offset, igno, fi, fip, idx, fp, fm, ig ) &
!$omp          shared( nogrp_ , f, ngw, psi, df, da, c, tpiba2, tens, dffts, me_bgrp, &
!$omp                  i, n, ggp, nls, nlsm )

      eig_offset = 0
      igno = 1

      DO idx = 1, 2*nogrp_ , 2

         IF( idx + i - 1 <= n ) THEN
            if (tens) then
               fi = -0.5d0
               fip = -0.5d0
            else
               fi = -0.5d0*f(i+idx-1)
               fip = -0.5d0*f(i+idx)
            endif
            IF( dffts%have_task_groups ) THEN
!$omp do 
               DO ig=1,ngw
                  fp= psi(nls(ig)+eig_offset) +  psi(nlsm(ig)+eig_offset)
                  fm= psi(nls(ig)+eig_offset) -  psi(nlsm(ig)+eig_offset)
                  df(ig+igno-1)= fi *(tpiba2 * ggp(ig) * c(ig,idx+i-1) + &
                                 CMPLX(real (fp), aimag(fm), kind=dp ))
                  da(ig+igno-1)= fip*(tpiba2 * ggp(ig) * c(ig,idx+i  ) + &
                                 CMPLX(aimag(fp),-real (fm), kind=dp ))
               END DO
!$omp end do
               igno = igno + ngw
            ELSE
!$omp do 
               DO ig=1,ngw
                  fp= psi(nls(ig)) + psi(nlsm(ig))
                  fm= psi(nls(ig)) - psi(nlsm(ig))
                  df(ig)= fi*(tpiba2*ggp(ig)* c(ig,idx+i-1)+CMPLX(DBLE(fp), AIMAG(fm),kind=DP))
                  da(ig)=fip*(tpiba2*ggp(ig)* c(ig,idx+i  )+CMPLX(AIMAG(fp),-DBLE(fm),kind=DP))
               END DO
!$omp end do
            END IF
         END IF

         eig_offset = eig_offset + dffts%nr3x * dffts%nsw(me_bgrp+1)

         ! We take into account the number of elements received from other members of the orbital group

      ENDDO

!$omp end parallel 
      !
      IF(dft_is_meta()) THEN
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
!$omp parallel default(none) &
!$omp          private(iv,jv,ivoff,jvoff,dd,dv,inl,jnl,is,isa,ism,igrp,idx,fi,fip) &
!$omp          shared( nogrp_ , f, ngw, deeq, bec, af, aa, i, n, nsp, na, nh, dvan, tens, ish, iss1, iss2 )
         !
         igrp = 1

         DO idx = 1, 2*nogrp_ , 2

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
!$omp do
                              DO ia=1,na(is)
                                 inl = ivoff + ia
                                 jnl = jvoff + ia
                                 dd = deeq(iv,jv,isa+ia,iss1) + dv
                                 af(inl,igrp) = af(inl,igrp) - fi  * dd * bec(jnl,i+idx-1)
                                 dd = deeq(iv,jv,isa+ia,iss2) + dv
                                 aa(inl,igrp) = aa(inl,igrp) - fip * dd * bec(jnl,i+idx)
                              END DO
                           ELSE
!$omp do
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

            igrp = igrp + 1

         END DO

!$omp end parallel
!
         CALL dgemm ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, af, nhsa, 1.0d0, df, 2*ngw)

         CALL dgemm ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, aa, nhsa, 1.0d0, da, 2*ngw)
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

