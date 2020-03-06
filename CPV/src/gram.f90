!
! Copyright (C) 2002-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-------------------------------------------------------------------------
SUBROUTINE gram_bgrp( betae, bec_bgrp, nkbx, cp_bgrp, ngwx )
!-----------------------------------------------------------------------
!     gram-schmidt orthogonalization of the set of wavefunctions cp
!
      USE uspp,           ONLY : nkb, nkbus
      USE uspp,           ONLY : qq_nt
      USE gvecw,          ONLY : ngw
      USE electrons_base, ONLY : nbspx_bgrp, ibgrp_g2l, nupdwn, iupdwn, nbspx, iupdwn_bgrp, nspin
      USE kinds,          ONLY : DP
      USE mp_global,      ONLY : inter_bgrp_comm
      USE mp,             ONLY : mp_sum
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: nkbx, ngwx
      REAL(DP)      :: bec_bgrp( nkbx, nbspx_bgrp )
      COMPLEX(DP)   :: cp_bgrp( ngwx, nbspx_bgrp ), betae( ngwx, nkbx )
!
      REAL(DP) :: anorm
      REAL(DP), ALLOCATABLE :: csc( : )
      COMPLEX(DP), ALLOCATABLE :: ctmp( : )
      REAL(DP), ALLOCATABLE :: temp(:) 
      COMPLEX(DP), ALLOCATABLE :: cp_tmp(:) 
      COMPLEX(DP), ALLOCATABLE :: cp(:,:) 
      REAL(DP), ALLOCATABLE :: bec_tmp(:) 
      REAL(DP), ALLOCATABLE :: csc2( : )
      INTEGER :: i,k,j, ig, ibgrp_k, ibgrp_i, nbgrp_im1, iss, kk
      REAL(DP), PARAMETER :: one  =  1.d0
      REAL(DP), PARAMETER :: mone = -1.d0
!
      CALL start_clock( 'gram' )

      ALLOCATE( csc( nbspx ) )
      ALLOCATE( ctmp( ngwx ) )
      ALLOCATE( cp_tmp( ngwx ) )
!      ALLOCATE( cp( ngw, MAXVAL(nupdwn)  ) )
      ALLOCATE( bec_tmp( nkbx ) )
      ALLOCATE( csc2( SIZE( csc ) ) )
      ALLOCATE( temp( ngw ) )
!
      DO iss = 1, nspin

!      kk = 0
!      DO k = iupdwn(iss), iupdwn(iss) + nupdwn(iss) - 1
!         ibgrp_k = ibgrp_g2l( k )
!         IF( ibgrp_k > 0 ) THEN
!            kk = kk + 1
!            cp(:,kk) =  cp_bgrp(:,ibgrp_k)
!         END IF
!      END DO

      DO i = iupdwn(iss), iupdwn(iss) + nupdwn(iss) - 1 
         !
         ibgrp_i = ibgrp_g2l( i )
         !
         CALL gracsc_bgrp( bec_bgrp, betae, cp_bgrp, i, csc, iss, nbgrp_im1 )
         !
         ! calculate orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
         !
         IF( ibgrp_i > 0 ) THEN
            ctmp = cp_bgrp( :, ibgrp_i )
         ELSE
            ctmp = 0.0d0
         END IF
         !
         IF( nbgrp_im1 > 0 .AND. ngw > 0 ) &
            CALL dgemv( 'N', 2*ngw, nbgrp_im1, mone, cp_bgrp(1,iupdwn_bgrp(iss)), 2*ngwx, csc, 1, one, ctmp, 1 )

         CALL mp_sum( ctmp, inter_bgrp_comm )

         IF( ibgrp_i > 0 ) THEN
            cp_bgrp( :, ibgrp_i ) = ctmp
            anorm = cscnorm( bec_bgrp, cp_bgrp, ibgrp_i, nbspx_bgrp )
            CALL dscal( 2*ngw, 1.0d0/anorm, cp_bgrp(1,ibgrp_i), 1 )
            CALL dscal( nkbx, 1.0d0/anorm, bec_bgrp(1,ibgrp_i), 1 )
         END IF
      END DO
      END DO
!
      DEALLOCATE( temp )
      DEALLOCATE( csc2 )
      DEALLOCATE( bec_tmp )
      DEALLOCATE( cp_tmp )
!      DEALLOCATE( cp )
      DEALLOCATE( ctmp )
      DEALLOCATE( csc )

      CALL stop_clock( 'gram' )
!
      RETURN

CONTAINS

!-----------------------------------------------------------------------
   FUNCTION cscnorm( bec, cp, i, n )
!-----------------------------------------------------------------------
!
!     Compute the norm of the i-th electronic state = (<c|S|c>)^(1/2) 
!     requires in input the updated bec(i)
!
      USE ions_base,          ONLY: nat, ityp
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE uspp_param,         ONLY: nh, upf
      USE uspp,               ONLY: qq_nt, indv_ijkb0
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm, me_bgrp, nproc_bgrp
      USE kinds,              ONLY: DP
!
      IMPLICIT NONE
      !
      INTEGER,     INTENT(IN)  :: i, n
      REAL(DP),    INTENT(IN) :: bec( :, : )
      COMPLEX(DP), INTENT(IN) :: cp( :, : )
      !
      REAL(DP) :: cscnorm
      !
      INTEGER ig, is, iv, jv, ia, inl, jnl, indv
      REAL(DP) rsum
!
      rsum = 0.0d0
      DO ig=1,ngw
         rsum = rsum + DBLE(CONJG(cp(ig,i))*cp(ig,i))
      END DO
      rsum = 2.d0 * rsum
      IF (gstart == 2) rsum=rsum-DBLE(CONJG(cp(1,i))*cp(1,i))
!
      DO ia=1,nat
         IF ( MOD( ia, nproc_bgrp ) == me_bgrp ) THEN
            is = ityp(ia)
            IF( upf(is)%tvanp ) THEN
               indv = indv_ijkb0(ia)
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     IF(ABS(qq_nt(iv,jv,is)).GT.1.e-5) THEN 
                        rsum = rsum + qq_nt(iv,jv,is)*bec(indv+iv,i)*bec(indv+jv,i)
                     ENDIF
                  END DO
               END DO
            END IF
         END IF
      END DO

      CALL mp_sum( rsum, intra_bgrp_comm )
!
      cscnorm=SQRT(rsum)

      RETURN
      END FUNCTION cscnorm
!
!
!-------------------------------------------------------------------------
      SUBROUTINE gracsc_bgrp( bec_bgrp, betae, cp_bgrp, i, csc, iss, nk )
!-----------------------------------------------------------------------
!     requires in input the updated bec(k) for k<i
!     on output: bec(i) is recalculated
!
      USE ions_base,      ONLY: na, nat, ityp
      USE uspp,           ONLY: nkb, nkbus, qq_nt, indv_ijkb0
      USE uspp_param,     ONLY: nh, upf
      USE electrons_base, ONLY: ispin, ispin_bgrp, nbspx_bgrp, ibgrp_g2l, iupdwn, nupdwn, nbspx
      USE gvecw,          ONLY: ngw
      USE mp,             ONLY: mp_sum
      USE mp_global,      ONLY: intra_bgrp_comm, inter_bgrp_comm, me_bgrp, nproc_bgrp
      USE kinds,          ONLY: DP
      USE gvect, ONLY: gstart
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: i, iss
      INTEGER, INTENT(OUT) :: nk
      COMPLEX(DP) :: betae( :, : )
      REAL(DP)    :: bec_bgrp( :, : )
      COMPLEX(DP) :: cp_bgrp( :, : )
      REAL(DP)    :: csc( : )
      INTEGER     :: k, kmax,ig, is, iv, jv, ia, inl, jnl, ibgrp_k, ibgrp_i
      REAL(DP)    :: rsum
#if defined(_OPENMP)
      INTEGER :: mytid, ntids, omp_get_thread_num, omp_get_num_threads
#endif
      !
      !     calculate csc(k)=<cp(i)|cp(k)>,  k<i
      !
      kmax = i - 1
      !

      csc    = 0.0d0

      ibgrp_i = ibgrp_g2l( i )
      IF( ibgrp_i > 0 ) THEN
         cp_tmp = cp_bgrp( :, ibgrp_i )
      ELSE
         cp_tmp = 0.0d0
      END IF

      CALL mp_sum( cp_tmp, inter_bgrp_comm )

      DO k = iupdwn( iss ), kmax
         ibgrp_k = ibgrp_g2l( k )
         IF( ibgrp_k > 0 ) THEN
            DO ig = 1, ngw
               temp(ig) = DBLE( cp_bgrp(ig,ibgrp_k) * CONJG(cp_tmp(ig)) )
            END DO
            csc(k) = 2.0d0 * SUM(temp)
            IF (gstart == 2) csc(k) = csc(k) - temp(1)
         END IF
      END DO
      !
      !
      !     calculate bec(i)=<cp(i)|beta>
      !
      ibgrp_i = ibgrp_g2l( i )
      !
      IF(  ibgrp_i > 0 ) THEN
         DO ia = 1, nat
            is = ityp(ia)
            IF( upf(is)%tvanp ) THEN
               DO iv=1,nh(is)
                  inl=indv_ijkb0(ia)+iv
                  DO ig=1,ngw
                     temp(ig) = DBLE( cp_bgrp(ig,ibgrp_i) * CONJG(betae(ig,inl)) )
                  END DO
                  bec_bgrp(inl,ibgrp_i)=2.d0*SUM(temp)
                  IF (gstart == 2) bec_bgrp(inl,ibgrp_i)= bec_bgrp(inl,ibgrp_i)-temp(1)
               END DO
            END IF
         END DO
      END IF

      IF(  ibgrp_i > 0 ) THEN
         DO ia = 1, nat
            is = ityp(ia)
            inl= indv_ijkb0(ia)
            IF( upf(is)%tvanp ) THEN
               bec_tmp( inl + 1: inl + nh(is) ) = bec_bgrp( inl + 1: inl + nh(is), ibgrp_i )
            ELSE
               bec_tmp( inl + 1: inl + nh(is) ) = 0.0d0
            END IF
         END DO
         CALL mp_sum( bec_tmp, intra_bgrp_comm )
         DO ia = 1, nat
            is = ityp(ia)
            IF( upf(is)%tvanp ) THEN
               inl=indv_ijkb0(ia)
               bec_bgrp( inl + 1: inl + nh(is), ibgrp_i ) = bec_tmp( inl + 1: inl + nh(is) )
            END IF
         END DO
      ELSE
         bec_tmp = 0.0d0
      END IF

      CALL mp_sum( bec_tmp, inter_bgrp_comm )
!
!     calculate csc(k)=<cp(i)|S|cp(k)>,  k<i
!
      csc2    = 0.0d0

      DO k = iupdwn( iss ), kmax
         IF ( MOD( k, nproc_bgrp ) == me_bgrp ) THEN
            rsum=0.d0
            ibgrp_k = ibgrp_g2l( k )
            IF( ibgrp_k > 0 ) THEN
               DO ia = 1, nat
                  is=ityp(ia)
                  IF( upf(is)%tvanp ) THEN
                     DO iv=1,nh(is)
                        inl=indv_ijkb0(ia)+iv
                        DO jv=1,nh(is)
                           jnl=indv_ijkb0(ia)+jv
                           IF(ABS(qq_nt(iv,jv,is)).GT.1.e-5) THEN 
                              rsum = rsum + qq_nt(iv,jv,is)*bec_tmp(inl)*bec_bgrp(jnl,ibgrp_k)
                           ENDIF
                        END DO
                     END DO
                  END IF
               END DO
            END IF
            csc2(k)=csc2(k)+rsum
         ENDIF
      END DO
!
!     orthogonalized cp(i) : |cp(i)>=|cp(i)>-\sum_k<i csc(k)|cp(k)>
!
!     corresponing bec:  bec(i)=<cp(i)|beta>-csc(k)<cp(k)|beta>
!
      CALL mp_sum( csc, intra_bgrp_comm )
      CALL mp_sum( csc2, intra_bgrp_comm )
      CALL mp_sum( csc, inter_bgrp_comm )
      CALL mp_sum( csc2, inter_bgrp_comm )
      csc = csc + csc2

      bec_tmp = 0.0d0
      DO k = iupdwn(iss), kmax
         ibgrp_k = ibgrp_g2l( k )
         IF( ibgrp_k > 0 ) THEN
            DO inl=1,nkbx
               bec_tmp(inl)=bec_tmp(inl)-csc(k)*bec_bgrp(inl,ibgrp_k)
            END DO
         END IF
      END DO
      nk = 0
      DO k = iupdwn(iss), kmax
         ibgrp_k = ibgrp_g2l( k )
         IF( ibgrp_k > 0 ) THEN
            nk = nk + 1 
            csc( nk ) = csc( k )
         END IF
      END DO
      CALL mp_sum( bec_tmp, inter_bgrp_comm )
      IF( ibgrp_i > 0 ) bec_bgrp(:,ibgrp_i ) = bec_bgrp(:,ibgrp_i ) + bec_tmp

!
      RETURN
      END SUBROUTINE gracsc_bgrp

END SUBROUTINE gram_bgrp
