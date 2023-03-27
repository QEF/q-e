!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include <cpv_device_macros.h> 
!-------------------------------------------------------------------------
SUBROUTINE gram_bgrp( betae, bec_bgrp, nkbx, cp_bgrp, ngwx )
      !-----------------------------------------------------------------------
      !! Gram-Schmidt orthogonalization of the set of wavefunctions cp.
      !
      USE gvecw,          ONLY : ngw
      USE electrons_base, ONLY : nbspx_bgrp, ibgrp_g2l, nupdwn, iupdwn, nbspx, iupdwn_bgrp, nspin
      USE kinds,          ONLY : DP
      USE mp,             ONLY : mp_sum
      USE gvect,          ONLY : gstart
      USE mp_global,      ONLY : intra_bgrp_comm, inter_bgrp_comm, me_bgrp, nproc_bgrp
      USE uspp_param,     ONLY: nh, upf
      USE uspp,           ONLY: qq_nt, ofsbeta
      USE ions_base,      ONLY: nsp, nat, ityp
      USE mp_world,       ONLY : mpime
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, ngwx
      REAL(DP)      :: bec_bgrp( nkbx, nbspx_bgrp )
      COMPLEX(DP)   :: cp_bgrp( ngwx, nbspx_bgrp )
      COMPLEX(DP), INTENT(IN) :: betae( ngwx, nkbx )
!
      REAL(DP) :: anorm
      REAL(DP), ALLOCATABLE :: csc( : )
      COMPLEX(DP), ALLOCATABLE :: ctmp( : )
      REAL(DP), ALLOCATABLE :: temp(:) 
      COMPLEX(DP), ALLOCATABLE :: cp_tmp(:) 
      REAL(DP), ALLOCATABLE :: bec_tmp(:) 
      REAL(DP), ALLOCATABLE :: csc2( : )
      INTEGER :: i,k,j, ig, ibgrp_k, ibgrp_i, nbgrp_im1, iss
      REAL(DP), PARAMETER :: one  =  1.d0
      REAL(DP), PARAMETER :: mone = -1.d0
      REAL(DP) :: g0
      INTEGER  :: ia_s, ia_e, mykey 
      INTEGER, ALLOCATABLE :: iqq(:,:), nqq(:) 
      LOGICAL  :: tvanp(nsp), all_tvanp
!
      CALL start_clock( 'gram' )

      g0 = 0.0d0
      IF (gstart == 2) g0 = 1.0d0

      ALLOCATE( csc( nbspx ) )
      ALLOCATE( ctmp( ngwx ) )
      ALLOCATE( cp_tmp( ngwx ) )
      ALLOCATE( bec_tmp( nkbx ) )
      ALLOCATE( csc2( SIZE( csc ) ) )
!
      tvanp(1:nsp) = upf(1:nsp)%tvanp
      CALL block_distribute(nat, me_bgrp, nproc_bgrp, ia_s, ia_e, mykey)  
DEV_ACC data copy(cp_bgrp, bec_bgrp) create(ctmp, cp_tmp,bec_tmp,csc,csc2) &
DEV_ACC & copyin(betae, qq_nt, ofsbeta,ityp, ibgrp_g2l, nh, tvanp ) 
      DO iss = 1, nspin
      DO i = iupdwn(iss), iupdwn(iss) + nupdwn(iss) - 1 
         !
         ibgrp_i = ibgrp_g2l( i )
         !
         CALL gracsc_bgrp( i, csc, iss, nbgrp_im1 )
         !
         ! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
         
         IF( ibgrp_i > 0 ) THEN
DEV_ACC kernels present(ctmp, cp_bgrp) 
            ctmp = cp_bgrp( :, ibgrp_i )
DEV_ACC end kernels     
         ELSE
DEV_ACC kernels present(ctmp) 
            ctmp = 0.0d0
DEV_ACC end kernels
         END IF
         !
         IF( nbgrp_im1 > 0 .AND. ngw > 0 ) THEN 
#if defined (__CUDA) && defined (_OPENACC)
DEV_ACC host_data use_device(cp_bgrp, csc, ctmp) 
           CALL mydgemv( 'N', 2*ngw, nbgrp_im1, mone, cp_bgrp(1,iupdwn_bgrp(iss)), 2*ngwx, csc, 1, one, ctmp, 1 )
DEV_ACC end host_data
#else
           CALL dgemv( 'N', 2*ngw, nbgrp_im1, mone, cp_bgrp(1,iupdwn_bgrp(iss)), 2*ngwx, csc, 1, one, ctmp, 1 )
#endif
         END IF 
DEV_ACC host_data  use_device(ctmp) 
         CALL mp_sum( ctmp, inter_bgrp_comm )
DEV_ACC end host_data
         IF( ibgrp_i > 0 ) THEN
DEV_ACC kernels present(cp_bgrp, ctmp) 
            cp_bgrp( :, ibgrp_i ) = ctmp
DEV_ACC end kernels 
            anorm = cscnorm( bec_bgrp, cp_bgrp, ibgrp_i, nbspx_bgrp )
DEV_ACC kernels 
            cp_bgrp(:,ibgrp_i) = cp_bgrp(:,ibgrp_i) / anorm
            bec_bgrp(:,ibgrp_i) = bec_bgrp(:,ibgrp_i) / anorm
DEV_ACC end kernels 
         END IF
      END DO
      END DO
DEV_ACC end data 
      DEALLOCATE( ctmp )
      DEALLOCATE( csc )
      DEALLOCATE( csc2 )
      DEALLOCATE( bec_tmp )
      DEALLOCATE( cp_tmp )

      CALL stop_clock( 'gram' )
!
      RETURN

CONTAINS   
!-----------------------------------------------------------------------
   FUNCTION cscnorm( bec, cp, i, n )
      !-----------------------------------------------------------------------
      !! Compute the norm of the i-th electronic state:
      !! \[ (\langle c|S|c \rangle)^{1/2} \ .\]
      !! Requires in input the updated \(\text{bec}(i)\).
      !
      USE ions_base,          ONLY: nat, ityp
      USE gvecw,              ONLY: ngw
      USE uspp_param,         ONLY: nh, upf
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm
      USE kinds,              ONLY: DP
!
      IMPLICIT NONE
      !
      INTEGER,     INTENT(IN)  :: i, n
      REAL(DP),    INTENT(IN) :: bec( :, : )
      COMPLEX(DP), INTENT(IN) :: cp( :, : )
      !
      REAL(DP) :: cscnorm, ddot
      !
      INTEGER  :: is, iv, jv, ia, indv
      REAL(DP) rsum, rsum_v
      REAL(DP), EXTERNAL  :: myddot 
!
DEV_ACC data present(bec, cp, tvanp,ofsbeta, nh, ityp, qq_nt)  
#if defined(__CUDA) && defined(_OPENACC) 
DEV_ACC host_data use_device(cp) 
      rsum = 2.d0 * myddot(2*ngw,cp(1,i),1,cp(1,i),1) 
DEV_ACC end host_data
#else
       rsum = 2.d0 * ddot(2*ngw,cp(1,i),1,cp(1,i),1) 
#endif  
DEV_ACC kernels present(cp) 
      rsum = rsum - g0 * REAL(CONJG(cp(1,i))*cp(1,i), DP)
DEV_ACC end kernels 
!
DEV_ACC parallel private(ia, is, iv, jv, rsum_v) reduction (+:rsum) vector_length(32) & 
DEV_ACC & present(ityp, ofsbeta, tvanp,nh,qq_nt,bec)  
DEV_ACC loop gang  
      DO ia= ia_s, ia_e 
         IF ( mykey == 0 ) THEN
            is = ityp(ia)
            IF( tvanp( is )) THEN
               rsum_v = 0._DP
               indv = ofsbeta(ia)
DEV_ACC loop vector collapse(2) reduction(+:rsum_v)
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     IF(ABS(qq_nt(iv,jv,is)).GT.1.e-5) THEN
                        rsum_v = rsum_v + qq_nt(iv,jv,is)*bec(indv+iv,i)*bec(indv+jv,i)
                     ENDIF
                  END DO
               END DO
               rsum = rsum + rsum_v 
            END IF
         END IF
      END DO
DEV_ACC end parallel 
DEV_ACC end data
      CALL mp_sum( rsum, intra_bgrp_comm )
!
      cscnorm=SQRT(rsum)
!
      RETURN
      END FUNCTION cscnorm
!
!
!-------------------------------------------------------------------------
      SUBROUTINE gracsc_bgrp( i, csc, iss, nk )
!-----------------------------------------------------------------------
      !! Requires in input the updated \(\text{bec}(k)\) for \(k<i\).
      !! On output: \(\text{bec}(i)\) is recalculated.
!
      USE ions_base,      ONLY: na, nat, ityp
      USE uspp,           ONLY: qq_nt, ofsbeta
      USE electrons_base, ONLY: ispin, ispin_bgrp, nbspx_bgrp, ibgrp_g2l, iupdwn, nupdwn, nbspx
      USE gvecw,          ONLY: ngw
      USE mp,             ONLY: mp_sum
      USE kinds,          ONLY: DP
      USE gvect, ONLY: gstart
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: i, iss
      INTEGER, INTENT(OUT) :: nk
      REAL(DP)    :: csc( : )
      INTEGER     :: k, kmax_bgrp, kmax,ig, is, iv, jv, ia, inl, jnl, ibgrp_k, ibgrp_i
      REAL(DP)    :: rsum, ddot, rsum_w, rsum_v, bec_tmp_inl
      INTEGER     :: omp_get_thread_num, omp_get_num_threads
      INTEGER     :: iupdwn_iss 
      
      !
      !     calculate csc(k)=<cp(i)|cp(k)>,  k<i
      !
      kmax = i - 1
      !
DEV_ACC data present(csc,csc2, cp_tmp, cp_bgrp, bec_tmp, bec_bgrp) &
DEV_ACC & present(ibgrp_g2l, nh, ofsbeta, qq_nt, ityp, tvanp, betae)  
DEV_ACC kernels present(csc)  
      csc    = 0.0d0
DEV_ACC end kernels 
      ibgrp_i = ibgrp_g2l( i )
DEV_ACC kernels present(cp_tmp) 
      IF( ibgrp_i > 0 ) THEN
         cp_tmp = cp_bgrp( :, ibgrp_i )
      ELSE
         cp_tmp = 0.0d0
      END IF
DEV_ACC end kernels 
!!DEV_ACC host_data use_device(cp_tmp) 

!!DEV_ACC update self(cp_tmp)
DEV_ACC host_data use_device(cp_tmp)  
      CALL mp_sum( cp_tmp, inter_bgrp_comm )
!!DEV_ACC update device(cp_tmp) 
DEV_ACC end host_data
!!DEV_ACC end host_data
      kmax_bgrp = 0
      nk = 0
      DO k = iupdwn( iss ), kmax
         IF( ibgrp_g2l( k ) > 0 ) THEN
            kmax_bgrp = ibgrp_g2l( k )
            nk = nk + 1
         END IF
      END DO
      kmax_bgrp = kmax_bgrp - iupdwn_bgrp(iss) + 1

      IF( kmax_bgrp > 0 .AND. ngw > 0 ) THEN
#if defined(__CUDA) && defined (_OPENACC)
DEV_ACC host_data use_device(cp_bgrp, cp_tmp, csc2)  
        CALL mydgemv( 'T', 2*ngw, kmax_bgrp, 1.0d0, cp_bgrp(1,iupdwn_bgrp(iss)), 2*ngwx, cp_tmp, 1, 0.0d0, csc2, 1 )
DEV_ACC end host_data
#else
        CALL dgemv( 'T', 2*ngw, kmax_bgrp, 1.0d0, cp_bgrp(1,iupdwn_bgrp(iss)), 2*ngwx, cp_tmp, 1, 0.0d0, csc2, 1 )
#endif 
      END IF 
      nk = 0
      iupdwn_iss = iupdwn( iss) 
DEV_ACC serial  present(ibgrp_g2l)      
      DO k = iupdwn_iss, kmax
         ibgrp_k = ibgrp_g2l( k )
         IF( ibgrp_k > 0 ) THEN
            nk = nk + 1
            csc(k) = 2.0d0 * csc2(nk) - g0 * DBLE( cp_bgrp(1,ibgrp_k) * CONJG(cp_tmp(1)) )
         END IF
      END DO
DEV_ACC end serial  

      IF(  ibgrp_i > 0 ) THEN
DEV_ACC parallel private(bec_tmp_inl, inl, is, ia, iv) vector_length(32) 
DEV_ACC loop gang 
         DO ia = 1, nat
            is = ityp(ia)
            DO iv=1,nh(is)
               inl=ofsbeta(ia)+iv
#if defined __USE_DDOT 
               bec_tmp_inl = 2.d0 * DDOT( 2*ngw, cp_bgrp(1,ibgrp_i), 1, betae(1,inl), 1) &
                              - g0 * DBLE(cp_bgrp(1,ibgrp_i) * CONJG(betae(1,inl)))
#else 
               bec_tmp_inl = 0._DP
DEV_ACC loop vector reduction(+:bec_tmp_inl)  
               DO ig =1, ngw
                  bec_tmp_inl = bec_tmp_inl + DBLE(CONJG(cp_bgrp(ig,ibgrp_i))*betae(ig,inl)) 
               END DO 
               bec_tmp_inl = 2._DP * bec_tmp_inl & 
                             - g0 * DBLE(CONJG(cp_bgrp(1,ibgrp_i)) * betae(1,inl))
#endif 
               bec_tmp(inl) = bec_tmp_inl 
            END DO
         END DO
DEV_ACC end parallel 
!!!DEV_ACC host_data use_device(bec_tmp)
DEV_ACC host_data use_device (bec_tmp)  
         CALL mp_sum( bec_tmp, intra_bgrp_comm )  ! parallel sum over G vectors within a band group
DEV_ACC end host_data
!!!DEV_ACC end host_data
DEV_ACC kernels 
         bec_bgrp( : , ibgrp_i ) = bec_tmp( : )
DEV_ACC end kernels 
      ELSE
DEV_ACC kernels 
         bec_tmp = 0.0d0
DEV_ACC end kernels 
      END IF
!!!DEV_ACC host_data use_device(bec_tmp) 
!!DEV_ACC update self(bec_tmp) 
DEV_ACC host_data use_device (bec_tmp) 
      CALL mp_sum( bec_tmp, inter_bgrp_comm )
DEV_ACC end host_data
!!DEV_ACC update device(bec_tmp) 
!!!DEV_ACC end host_data
!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i
!
DEV_ACC kernels 
      csc2    = 0.0d0
DEV_ACC end kernels 
iupdwn_iss = iupdwn(iss) 
DEV_OMP parallel if( (kmax - iupdwn_iss ) > omp_get_num_threads() ) default(none), &
DEV_OMP shared(iupdwn_iss,iss,kmax,nproc_bgrp,me_bgrp,nbspx,i,ibgrp_g2l,nh), &
DEV_OMP shared(ofsbeta,qq_nt,na,bec_tmp,bec_bgrp,csc2,nat,ityp,tvanp, ia_s, ia_e, mykey), &
DEV_OMP private( k, is, iv, jv, ia, inl, jnl, rsum, ibgrp_k, rsum_w, rsum_v )
DEV_OMP do
!civn
DEV_ACC parallel present(ibgrp_g2l, nh, ofsbeta, qq_nt, bec_tmp, bec_bgrp, csc2, ityp, tvanp) &
!!!DEV_ACC & (iupdwn_iss, iss, kmax, nproc_bgrp, me_bgrp, nbspx, i, ia_s, ia_e, mykey) &
!DEV_ACC & private(k, is, iv, jv, ia, inl, jnl, rsum, ibgrp_k) &
!DEV_ACC & num_workers(MIN(32,ia_e-ia_s+1)) vector_length(32)
DEV_ACC & vector_length(32)
DEV_ACC loop gang  private(rsum, ibgrp_k, rsum_w, is, inl, rsum_v) reduction (+:rsum_w, rsum) 
      DO k = iupdwn_iss, kmax
         rsum=0.d0
         ibgrp_k = ibgrp_g2l( k )
         IF( ibgrp_k > 0 ) THEN
            rsum_w = 0._DP 
!DEV_ACC loop worker reduction (+:rsum_w) private(rsum_v)   
            DO ia = ia_s, ia_e
               IF ( mykey  == 0 ) THEN
                  is=ityp(ia)
                  IF( tvanp(is) ) THEN
                     inl = ofsbeta(ia)
                     rsum_v = 0._DP 
DEV_ACC loop vector collapse(2) reduction (+:rsum_v) 
                     DO iv=1,nh(is)
                        DO jv=1,nh(is)
                           IF(ABS(qq_nt(iv,jv,is)).GT.1.e-5) THEN
                              rsum_v = rsum_v + qq_nt(iv,jv,is)*bec_tmp(inl+iv)*bec_bgrp(inl+jv,ibgrp_k)
                           ENDIF
                        END DO
                     END DO
                     rsum_w = rsum_w + rsum_v
                  END IF
               END IF
            END DO
            rsum = rsum + rsum_w
         ENDIF
         csc2(k)=csc2(k)+rsum
      END DO
DEV_OMP end do
DEV_OMP end parallel
DEV_ACC end parallel 
!
!     orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
!
!     corresponing bec:  bec(i)=<cp(i)|beta>-csc(k)<cp(k)|beta>
!
DEV_ACC host_data use_device(csc, csc2) 
      CALL mp_sum( csc, intra_bgrp_comm )
      CALL mp_sum( csc2, intra_bgrp_comm )
      CALL mp_sum( csc, inter_bgrp_comm )
      CALL mp_sum( csc2, inter_bgrp_comm )
DEV_ACC end host_data 
DEV_ACC kernels present(csc, csc2)  
      csc = csc + csc2
DEV_ACC end kernels 


      nk = 0
DEV_ACC serial present(ibgrp_g2l, csc) 
      DO k = iupdwn_iss, kmax
         ibgrp_k = ibgrp_g2l( k )
         IF( ibgrp_k > 0 ) THEN
            nk = nk + 1 
            csc( nk ) = csc( k )
         END IF
      END DO
DEV_ACC end serial 

      IF( nk > 0 .AND. ngw > 0 ) THEN
#if defined (__CUDA) && (_OPENACC)
DEV_ACC data copyin(bec_bgrp, csc) copyout(bec_tmp) 
DEV_ACC host_data use_device(bec_bgrp, csc, bec_tmp) 
        CALL mydgemv( 'N', nkbx, nk, -1.0d0, bec_bgrp(1,iupdwn_bgrp(iss)), nkbx, csc, 1, 0.0d0, bec_tmp, 1 )
DEV_ACC end host_data
DEV_ACC end data 
#else
        CALL dgemv( 'N', nkbx, nk, -1.0d0, bec_bgrp(1,iupdwn_bgrp(iss)), nkbx, csc, 1, 0.0d0, bec_tmp, 1 )
#endif
      ELSE
DEV_ACC kernels present(bec_tmp) 
        bec_tmp = 0.0d0
DEV_ACC end kernels 
      END IF
DEV_ACC host_data use_device(bec_tmp)  
      CALL mp_sum( bec_tmp, inter_bgrp_comm )
DEV_ACC end host_data 
DEV_ACC kernels  present(bec_bgrp, bec_tmp) 
      IF( ibgrp_i > 0 ) bec_bgrp(:,ibgrp_i ) = bec_bgrp(:,ibgrp_i ) + bec_tmp
DEV_ACC end kernels 
DEV_ACC end data 
!
      RETURN
      END SUBROUTINE gracsc_bgrp

END SUBROUTINE gram_bgrp
