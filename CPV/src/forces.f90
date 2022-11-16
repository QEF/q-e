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
#if defined(__CUDA)
#define PINMEM 
#else
#define PINMEM
#endif
!
!-------------------------------------------------------------------------
      SUBROUTINE dforce_x ( i, bec, vkb, c, df_out, da, v, ldv, ispin, f, n, nspin, v1 )
      !-----------------------------------------------------------------------
      !! Computes: the generalized force df=cmplx(dfr,dfi,kind=DP) acting on the i-th
      !!           electron state at the gamma point of the Brillouin zone
      !!           represented by the vector c=cmplx(cr,ci,kind=DP).
      !
      !! \[ d_n(g) = f_n { 0.5 g^2 c_n(g) + [\text{vc}_n](g) +
      !!             \sum_{i,ij} d^q_{i,ij} (-i)^l \beta_{i,'i}(g) 
      !!                          e^{-ig} r_i \langle \beta_{i,j} | c_n \rangle}  \]
      !
      ! d_n(g) = f_n { 0.5 g^2 c_n(g) + [vc_n](g) +
      !          sum_i,ij d^q_i,ij (-i)**l beta_i,`i(g) 
      !                         e^-ig.r_i < beta_i,j | c_n >}
      !
      USE parallel_include
      USE kinds,                  ONLY: dp
      USE control_flags,          ONLY: iprint
      USE uspp,                   ONLY: nhsa=>nkb, dvan, deeq, ofsbeta
      USE uspp_param,             ONLY: nhm, nh
      USE constants,              ONLY: pi, fpi
      USE ions_base,              ONLY: nsp, na, nat, ityp
      USE gvecw,                  ONLY: ngw, g2kin
      USE cell_base,              ONLY: tpiba2
      USE ensemble_dft,           ONLY: tens
      USE xc_lib,                 ONLY: xclib_dft_is, exx_is_active
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
      COMPLEX(DP)                :: df_out(:), da(:)
      INTEGER,     INTENT(IN)    :: ldv
      REAL(DP)                   :: v( ldv, * )
      INTEGER                    :: ispin( : )
      REAL(DP)                   :: f( : )
      INTEGER,     INTENT(IN)    :: n, nspin
      REAL(DP),    OPTIONAL      :: v1( ldv, * )
      !
      ! local variables
      !
      INTEGER     :: iv, jv, ia, is, iss1, iss2, ir, ig, inl, jnl
      INTEGER     :: igno, igrp, ierr
      INTEGER     :: idx, eig_offset, nogrp_ , inc, tg_nr3, szdf
      REAL(DP)    :: fi, fip, dd, dv
      COMPLEX(DP) :: fp, fm, ci
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
!dir$ attributes align: 4096 :: af, aa, psi, exx_a, exx_b
#endif
#endif
      REAL(DP),    ALLOCATABLE PINMEM :: af( :, : ), aa( :, : )
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      REAL(DP)    :: tmp1, tmp2                      ! Lingzhu Kong
      REAL(DP),    ALLOCATABLE :: exx_a(:), exx_b(:) ! Lingzhu Kong      
      !
      
      COMPLEX(DP), ALLOCATABLE :: df(:,:)
      
      CALL start_clock( 'dforce' ) 
      !
!=======================================================================
!exx_wf related
      IF(xclib_dft_is('hybrid').AND.exx_is_active()) THEN
         allocate( exx_a( dffts%nnr ) ); exx_a=0.0_DP
         allocate( exx_b( dffts%nnr ) ); exx_b=0.0_DP
      END IF
!=======================================================================

      szdf=SIZE(df_out(:))
      ALLOCATE( df(szdf,1) )

      nogrp_ = fftx_ntgrp(dffts)
      IF( nogrp_ > 1 ) THEN
         CALL errore('dforce','Task group not supported',1)
      END IF
      ALLOCATE( psi( dffts%nnr_tg ) )
      !
      CALL fftx_c2psi_gamma( dffts, psi, c(:,i:i), c(:,i+1) )
      !
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
      IF( dffts%has_task_groups ) THEN
         !
         CALL tg_get_group_nr3( dffts, tg_nr3 )
         !
!===============================================================================
!exx_wf related
         IF(xclib_dft_is('hybrid').AND.exx_is_active()) THEN
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
            IF(xclib_dft_is('hybrid').AND.exx_is_active()) THEN
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
            IF(xclib_dft_is('hybrid').AND.exx_is_active()) THEN
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
      CALL fwfft( 'Wave', psi, dffts )
      !
      !   note : the factor 0.5 appears 
      !       in the kinetic energy because it is defined as 0.5*g**2
      !       in the potential part because of the logics
      !
      !   Each processor will treat its own part of the eigenstate
      !   assigned to its ORBITAL group
      !
      eig_offset = 0
      CALL tg_get_recip_inc( dffts, inc )
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
            CALL fftx_psi2c_gamma( dffts, psi(eig_offset+1:eig_offset+inc), &
                                   df(igno:igno+ngw-1,1:1), da(igno:igno+ngw-1))
            IF( dffts%has_task_groups ) THEN
               DO ig=1,ngw
                  df(ig+igno-1,1)= fi *(tpiba2 * g2kin(ig) * c(ig,idx+i-1) + df(ig+igno-1,1))
                  da(ig+igno-1)= fip*(tpiba2 * g2kin(ig) * c(ig,idx+i  ) + da(ig+igno-1))
               END DO
            ELSE
               DO ig=1,ngw
                  df(ig,1)= fi*(tpiba2*g2kin(ig)* c(ig,idx+i-1)+df(ig,1))
                  da(ig)=fip*(tpiba2*g2kin(ig)* c(ig,idx+i  )+da(ig))
               END DO
            END IF
         END IF

         igno = igno + ngw
         eig_offset = eig_offset + inc

         ! We take into account the number of elements received from other members of the orbital group

      ENDDO

      !
      IF(xclib_dft_is('meta')) THEN
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
         CALL dforce_meta(c(1,i),c(1,i+1),df(:,1),da,psi,iss1,iss2,fi,fip) !METAGGA
      END IF


      IF( nhsa > 0 ) THEN
         !
         !     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
         ! 
         ALLOCATE( af( nhsa, nogrp_ ), aa( nhsa, nogrp_ ) )

         af = 0.0d0
         aa = 0.0d0
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
               DO ia = 1, nat
                  is = ityp(ia)
                  DO iv = 1, nh(is)
                     DO jv = 1, nh(is)
                        dv = dvan(iv,jv,is)
                        inl = ofsbeta(ia) + iv
                        jnl = ofsbeta(ia) + jv
                        IF( i + idx - 1 /= n ) THEN
                           dd = deeq(iv,jv,ia,iss1) + dv
                           af(inl,igrp) = af(inl,igrp) - fi  * dd * bec(jnl,i+idx-1)
                           dd = deeq(iv,jv,ia,iss2) + dv
                           aa(inl,igrp) = aa(inl,igrp) - fip * dd * bec(jnl,i+idx)
                        ELSE
                           dd = deeq(iv,jv,ia,iss1) + dv
                           af(inl,igrp) = af(inl,igrp) - fi * dd * bec(jnl,i+idx-1)
                        END IF
                     END DO
                  END DO
               END DO

            END IF

            igrp = igrp + 1

         END DO

         IF( ngw > 0 ) THEN
           CALL dgemm ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, af, nhsa, 1.0d0, df(:,1), 2*ngw)
           CALL dgemm ( 'N', 'N', 2*ngw, nogrp_ , nhsa, 1.0d0, vkb, 2*ngw, aa, nhsa, 1.0d0, da, 2*ngw)
         END IF
         !
         DEALLOCATE( aa, af )
         !
      ENDIF
!
      df_out=df(:,1)
      DEALLOCATE(df)
      IF(xclib_dft_is('hybrid').AND.exx_is_active()) DEALLOCATE(exx_a, exx_b)
      DEALLOCATE( psi )
!
      CALL stop_clock( 'dforce' ) 
!
      RETURN
   END SUBROUTINE dforce_x

#if defined (__CUDA)
!-------------------------------------------------------------------------
      SUBROUTINE dforce_gpu_x ( i, bec, vkb, c, df, da, v, ldv, ispin, f, n, nspin )
      !-----------------------------------------------------------------------
      !! GPU double of \(\texttt{dforce_x}\).
      !
      USE parallel_include
      USE kinds,                  ONLY: dp
      USE control_flags,          ONLY: iprint
      USE uspp,                   ONLY: nhsa=>nkb, dvan, deeq, ofsbeta
      USE uspp_param,             ONLY: nhm, nh
      USE constants,              ONLY: pi, fpi
      USE ions_base,              ONLY: nsp, na, nat, ityp
      USE gvecw,                  ONLY: ngw, g2kin
      USE cell_base,              ONLY: tpiba2
      USE ensemble_dft,           ONLY: tens
      USE xc_lib,                 ONLY: xclib_dft_is, exx_is_active
      USE fft_base,               ONLY: dffts
      USE fft_interfaces,         ONLY: fwfft, invfft
      USE mp_global,              ONLY: me_bgrp
      USE control_flags,          ONLY: many_fft
      USE fft_helper_subroutines
      USE exx_module,             ONLY: exx_potential
      USE device_memcpy_m,          ONLY : dev_memcpy
      USE cudafor
!
      IMPLICIT NONE
!
      INTEGER,     INTENT(IN)    :: i
      REAL(DP)                   :: bec(:,:)
      COMPLEX(DP), DEVICE        :: vkb(:,:)
      COMPLEX(DP), DEVICE        :: c(:,:)
      COMPLEX(DP)  PINMEM        :: df(:), da(:)
      INTEGER,     INTENT(IN)    :: ldv
      REAL(DP), DEVICE           :: v( :, : )
      INTEGER                    :: ispin( : )
      REAL(DP)                   :: f( : )
      INTEGER,     INTENT(IN)    :: n, nspin
      !
      ! local variables
      !
      INTEGER     :: iv, jv, ia, is, iss1, iss2, ir, ig, inl, jnl
      INTEGER     :: igno, igrp, ierr, ii
      INTEGER     :: idx, ioff
      REAL(DP)    :: fi, fip, dd, dv
      REAL(DP)    :: tmp1, tmp2                      
      COMPLEX(DP) :: fp, fm
      complex(DP), parameter :: ci=(0.0d0,1.0d0)

      REAL(DP),    ALLOCATABLE  PINMEM :: af( :, : ), aa( :, : )
      REAL(DP),    ALLOCATABLE, DEVICE :: af_d( :, : ), aa_d( :, : )
      COMPLEX(DP), ALLOCATABLE, DEVICE :: psi(:)
      COMPLEX(DP), ALLOCATABLE, DEVICE :: df_d(:)
      COMPLEX(DP), ALLOCATABLE, DEVICE :: da_d(:)
      REAL(DP),    ALLOCATABLE, DEVICE :: exx_a(:), exx_b(:)       
      REAL(DP),    ALLOCATABLE, DEVICE :: exx_potential_d(:,:)
      INTEGER,     DEVICE, POINTER     :: nl_d(:), nlm_d(:)
      !
      CALL start_clock( 'dforce' ) 
      !
      IF(xclib_dft_is('hybrid').AND.exx_is_active()) THEN
         allocate( exx_a( dffts%nnr ) ); exx_a=0.0_DP
         allocate( exx_b( dffts%nnr ) ); exx_b=0.0_DP
         allocate( exx_potential_d, source=exx_potential )
      END IF

      ALLOCATE( psi( dffts%nnr * many_fft ) )
      ALLOCATE( df_d( SIZE( df ) ) )
      ALLOCATE( da_d( SIZE( da ) ) )
      !
      psi = 0.0d0
      nl_d => dffts%nl_d
      nlm_d => dffts%nlm_d

      ioff = 0
      DO ii = i, i + 2 * many_fft - 1, 2
         IF( ii < n ) THEN
!$cuf kernel do(1)
            do ig = 1, dffts%ngw
               psi( nlm_d( ig ) + ioff) = CONJG( c( ig, ii ) ) + ci * conjg( c( ig, ii+1))
               psi( nl_d( ig )  + ioff) = c( ig, ii ) + ci * c( ig, ii+1)
            end do
         ELSE IF( ii == n ) THEN
!$cuf kernel do(1)
            do ig = 1, dffts%ngw
               psi( nlm_d( ig ) + ioff) = CONJG( c( ig, ii ) )
               psi( nl_d( ig )  + ioff) = c( ig, ii )
            end do
         END IF
         ioff = ioff + dffts%nnr
      END DO
      !
      CALL invfft( 'Wave', psi, dffts, many_fft )
      !
      ioff = 0
      DO ii = i, i + 2 * many_fft - 1, 2
         IF( ii < n ) THEN
            iss1=ispin( ii )
            iss2=ispin( ii + 1 )
         ELSE IF( ii == n ) THEN
            iss1=ispin( ii )
            iss2=iss1
         END IF

         IF (xclib_dft_is('hybrid').AND.exx_is_active()) THEN
           IF (ii.le.n) THEN
             IF ( (mod(n,2).ne.0 ) .and. (ii.eq.n) ) THEN
               !$cuf kernel do(1) 
               DO ir=1,dffts%nnr
                  tmp1 = v(ir,iss1)* DBLE(psi(ir+ioff))+exx_potential_d(ir, ii)
                  tmp2 = v(ir,iss2)*AIMAG(psi(ir+ioff))
                  psi(ir+ioff)=CMPLX( tmp1, tmp2, kind=DP )
               END DO
             ELSE
               !$cuf kernel do(1) 
               DO ir=1,dffts%nnr
                  tmp1 = v(ir,iss1)* DBLE(psi(ir+ioff))+exx_potential_d(ir, ii)
                  tmp2 = v(ir,iss2)*AIMAG(psi(ir+ioff))+exx_potential_d(ir, ii+1)
                  psi(ir+ioff)=CMPLX( tmp1, tmp2, kind=DP )
               END DO
             END IF
           END IF
         ELSE
           !$cuf kernel do(1)
           DO ir=1,dffts%nnr
              psi(ir+ioff)=CMPLX( v(ir,iss1)* DBLE(psi(ir+ioff)), &
                                  v(ir,iss2)*AIMAG(psi(ir+ioff)) ,kind=DP)
           END DO
         END IF
         ioff = ioff + dffts%nnr
      END DO

      CALL fwfft( 'Wave', psi, dffts, many_fft )

      igno = 0
      ioff = 0
      DO idx = 1, 2 * many_fft, 2

         IF( idx + i - 1 <= n ) THEN
            if (tens) then
               fi = -0.5d0
               fip = -0.5d0
            else
               fi = -0.5d0*f(i+idx-1)
               fip = -0.5d0*f(i+idx)
            endif
            CALL fftx_psi2c_gamma_gpu( dffts, psi( 1+ioff : ioff+dffts%nnr ), df_d(1+igno:igno+ngw), da_d(1+igno:igno+ngw))
!$acc parallel loop present(g2kin, da_d, df_d, c)
            DO ig=1,ngw
               df_d(ig+igno)= fi*(tpiba2*g2kin(ig)* c(ig,idx+i-1)+df_d(ig+igno))
               da_d(ig+igno)=fip*(tpiba2*g2kin(ig)* c(ig,idx+i  )+da_d(ig+igno))
            END DO
         END IF

         igno = igno + ngw
         ioff = ioff + dffts%nnr

      ENDDO

      !

      IF( nhsa > 0 ) THEN
         !
         !     aa_i,i,n = sum_j d_i,ij <beta_i,j|c_n>
         ! 
         ALLOCATE( af( nhsa, many_fft ), aa( nhsa, many_fft ) )
         ALLOCATE( af_d( nhsa, many_fft ), aa_d( nhsa, many_fft ) )
         !
!$omp parallel do default(none), &
!$omp shared(many_fft,i,n,tens,f,nat,ityp,nh,dvan,ofsbeta,deeq,af,aa,bec,ispin), &
!$omp private(idx,igrp,fi,fip,ia,is,iv,jv,inl,jnl,dv,dd,iss1,iss2)
         DO idx = 1, 2*many_fft , 2

            igrp = idx/2+1
            af(:,igrp) = 0.0d0
            aa(:,igrp) = 0.0d0

            IF( idx + i - 1 <= n ) THEN

               IF( i + idx - 1 /= n ) THEN
                  fi = f(i+idx-1)
                  fip= f(i+idx)
                  iss1=ispin( i+idx-1 )
                  iss2=ispin( i+idx   )
               ELSE
                  fi = f(i+idx-1)
                  iss1=ispin( i+idx-1 )
               END IF
               IF (tens) THEN
                  fi = 1.0d0
                  fip= 1.0d0
               ENDIF
               !
               DO ia = 1, nat
                  is = ityp(ia)
                  DO iv = 1, nh(is)
                     DO jv = 1, nh(is)
                        dv = dvan(iv,jv,is)
                        inl = ofsbeta(ia) + iv
                        jnl = ofsbeta(ia) + jv
                        IF( i + idx - 1 /= n ) THEN
                           dd = deeq(iv,jv,ia,iss1) + dv
                           af(inl,igrp) = af(inl,igrp) - fi  * dd * bec(jnl,i+idx-1)
                           dd = deeq(iv,jv,ia,iss2) + dv
                           aa(inl,igrp) = aa(inl,igrp) - fip * dd * bec(jnl,i+idx)
                        ELSE
                           dd = deeq(iv,jv,ia,iss1) + dv
                           af(inl,igrp) = af(inl,igrp) - fi * dd * bec(jnl,i+idx-1)
                        END IF
                     END DO
                  END DO
               END DO
            END IF
         END DO
!$omp end parallel do

         CALL dev_memcpy( af_d, af )
         CALL dev_memcpy( aa_d, aa )

         IF( ngw > 0 ) THEN
           CALL MYDGEMM ( 'N', 'N', 2*ngw, many_fft , nhsa, 1.0d0, vkb, 2*ngw, af_d, nhsa, 1.0d0, df_d, 2*ngw)
           CALL MYDGEMM ( 'N', 'N', 2*ngw, many_fft , nhsa, 1.0d0, vkb, 2*ngw, aa_d, nhsa, 1.0d0, da_d, 2*ngw)
         END IF
         !
         DEALLOCATE( aa, af )
         DEALLOCATE( aa_d, af_d )
         !
      ENDIF

      CALL dev_memcpy( df, df_d )
      CALL dev_memcpy( da, da_d )
!
      DEALLOCATE( df_d )
      DEALLOCATE( da_d )
      DEALLOCATE( psi )
      IF (xclib_dft_is('hybrid').AND.exx_is_active()) DEALLOCATE(exx_potential_d)
      NULLIFY(nl_d) 
      NULLIFY(nlm_d)
!
      CALL stop_clock( 'dforce' ) 
!
      RETURN
   END SUBROUTINE dforce_gpu_x
#endif
