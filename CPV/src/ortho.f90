!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined(__CUDA)
#define DEVICEATTR ,DEVICE
#else
#define DEVICEATTR
#endif

MODULE local_ortho_memory
#if defined(__CUDA)
   USE cudafor
#endif
   USE kinds,              ONLY: DP
   IMPLICIT NONE
   SAVE

   REAL(DP),   ALLOCATABLE DEVICEATTR :: s(:,:), sig(:,:), tau(:,:), stmp(:,:)
   REAL(DP),   ALLOCATABLE DEVICEATTR :: wrk(:,:), rhoa(:,:), rhos(:,:), rhod(:)
#if defined(__CUDA)
   REAL(DP),   ALLOCATABLE, DEVICE :: qbecp_d(:,:), qbephi_d(:,:)
#endif

   REAL(DP), ALLOCATABLE DEVICEATTR :: xloc(:,:)

CONTAINS

   SUBROUTINE allocate_local_ortho_memory( nss, nx0 )
      INTEGER, INTENT(IN) :: nss, nx0
      INTEGER :: info
      IF( ALLOCATED( rhos ) ) THEN
         IF( nx0 == SIZE( rhos, 1 ) ) THEN
            RETURN
         ELSE
            DEALLOCATE( rhos, rhoa, s, sig, tau, rhod )
            IF(ALLOCATED(wrk)) DEALLOCATE(wrk) 
            IF(ALLOCATED(stmp)) DEALLOCATE(stmp) 
         END IF 
      END IF
      ALLOCATE( rhos( nx0, nx0 ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating rhos ', ABS( info ) )
      rhos = 0
      ALLOCATE( rhoa( nx0, nx0 ), STAT = info )   !   antisymmetric part of rho
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating rhoa ', ABS( info ) )
      rhoa = 0
      ALLOCATE( s( nx0, nx0 ), STAT = info ) 
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating s ', ABS( info ) )
      s = 0
      ALLOCATE( sig( nx0, nx0 ), STAT = info ) 
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating sig ', ABS( info ) )
      sig = 0
      ALLOCATE( tau( nx0, nx0 ), STAT = info ) 
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating tau ', ABS( info ) )
      tau = 0
      ALLOCATE( rhod( nss ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating tau ', ABS( info ) )
      rhod = 0
#if defined(__CUDA)
      ALLOCATE( wrk( nss, nss ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho_gamma ', ' allocating wrk ', 1 )
      ALLOCATE( stmp( nss, nss ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho_gamma ', ' allocating stmp ', 1 )
#endif
   END SUBROUTINE

   SUBROUTINE sync_device_ortho_memory( qbecp, qbephi )
      REAL(DP), INTENT(IN)    :: qbephi( :, : ), qbecp( :, : )
      INTEGER :: info
#if defined(__CUDA)
      IF(ALLOCATED(qbecp_d)) THEN
         IF( SIZE(qbecp,1) == SIZE(qbecp_d,1) .AND. SIZE(qbecp,2) == SIZE(qbecp_d,2) ) THEN
            qbecp_d = qbecp
            qbephi_d = qbephi
            RETURN
         ELSE
            DEALLOCATE( qbecp_d, qbephi_d )
         END IF
      END IF
      ALLOCATE( qbecp_d, source=qbecp, STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating qbecp_d ', ABS( info ) )
      ALLOCATE( qbephi_d, source=qbephi, STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating qbephi_d ', ABS( info ) )
      info = cudaDeviceSynchronize()
#endif

   END SUBROUTINE

   SUBROUTINE deallocate_local_ortho_memory()
     IF(ALLOCATED(s)) DEALLOCATE(s) 
     IF(ALLOCATED(sig)) DEALLOCATE(sig) 
     IF(ALLOCATED(tau)) DEALLOCATE(tau) 
     IF(ALLOCATED(stmp)) DEALLOCATE(stmp) 
     IF(ALLOCATED(wrk)) DEALLOCATE(wrk) 
     IF(ALLOCATED(rhoa)) DEALLOCATE(rhoa)
     IF(ALLOCATED(rhos)) DEALLOCATE(rhos) 
     IF(ALLOCATED(rhod)) DEALLOCATE(rhod) 
     IF(ALLOCATED(xloc)) DEALLOCATE(xloc) 
#if defined(__CUDA)
     IF(ALLOCATED(qbecp_d)) DEALLOCATE(qbecp_d) 
     IF(ALLOCATED(qbephi_d)) DEALLOCATE(qbephi_d) 
#endif
   END SUBROUTINE deallocate_local_ortho_memory

   SUBROUTINE x0_to_xloc( x0, nx0, ccc_, idesc )
      REAL(DP), INTENT(IN) :: x0(:,:)
      REAL(DP), INTENT(IN) :: ccc_
      INTEGER, INTENT(IN) :: nx0
      INTEGER,  INTENT(IN)  :: idesc(:)
      include 'laxlib.fh'
      INTEGER :: i, j, info
      REAL(DP) DEVICEATTR :: ccc
      IF( ALLOCATED(xloc) ) THEN
        IF( nx0 /= SIZE(xloc,1) ) THEN
           DEALLOCATE(xloc)
        END IF
      END IF
      IF( .NOT. ALLOCATED(xloc) ) THEN
         ALLOCATE( xloc( nx0, nx0 ), STAT = info )
         IF( info /= 0 ) &
            CALL errore( ' ortho ', ' allocating xloc ', ABS( info ) )
      END IF
      IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
         RETURN
      ENDIF
      !
      xloc = x0
      ccc = ccc_
!$cuf kernel do(2) <<<*,*>>>
      DO j = 1, SIZE(xloc,2)
         DO i = 1, SIZE(xloc,1)
            xloc(i,j) = xloc(i,j) * ccc
         END DO
      END DO
   END SUBROUTINE x0_to_xloc
   
   SUBROUTINE xloc_to_x0( x0, nx0, ccc_, idesc )
      REAL(DP), INTENT(OUT) :: x0(:,:)
      INTEGER, INTENT(IN) :: nx0
      REAL(DP), INTENT(IN) :: ccc_
      INTEGER,  INTENT(IN)  :: idesc(:)
      include 'laxlib.fh'
      INTEGER :: i, j
      REAL(DP) DEVICEATTR :: byccc
      IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
         RETURN
      ENDIF
      byccc = 1.0d0 / ccc_
!$cuf kernel do(2) <<<*,*>>>
      DO j = 1, SIZE(xloc,2)
         DO i = 1, SIZE(xloc,1)
            xloc(i,j) = xloc(i,j) * byccc
         END DO
      END DO
      x0 = xloc
   END SUBROUTINE xloc_to_x0

   SUBROUTINE distribute_matrix( a, b, ir, nr, ic, nc, comm )
      USE mp, ONLY: mp_bcast
      REAL(DP) DEVICEATTR :: a(:,:), b(:,:)
      INTEGER, INTENT(IN) :: ir, nr, ic, nc, comm
      INTEGER :: i, j, info
      CALL qe_sync()
      CALL mp_bcast( a, 0, comm )
!$cuf kernel do(2) <<<*,*>>>
      DO j = 1, nc
         DO i = 1, nr
            b( i, j ) = a( i + ir - 1, j + ic - 1 )
         END DO
      END DO
      CALL qe_sync()
      RETURN
   END SUBROUTINE

   SUBROUTINE collect_matrix( a, b, ir, nr, ic, nc, comm )
      USE mp, ONLY: mp_sum
      REAL(DP) DEVICEATTR :: a(:,:), b(:,:)
      INTEGER, INTENT(IN) :: ir, nr, ic, nc, comm
      INTEGER :: i, j, info
      CALL qe_sync()
      a = 0.0d0
!$cuf kernel do(2) <<<*,*>>>
      DO j = 1, nc
         DO i = 1, nr
            a( ir + i - 1, ic + j - 1 ) = b( i, j )
         END DO
      END DO
      CALL mp_sum( a, comm )
      CALL qe_sync()
      RETURN
   END SUBROUTINE

   SUBROUTINE consistency_check( a, idesc )
      REAL(DP) DEVICEATTR, INTENT(IN) :: a(:,:)
      INTEGER,  INTENT(IN)  :: idesc(:)
      INTEGER :: i, j
      include 'laxlib.fh'
      !
      ! on some machines (IBM RS/6000 for instance) the following test allows
      ! to distinguish between Numbers and Sodium Nitride (NaN, Not a Number).
      ! If a matrix of Not-Numbers is passed to rs, the most likely outcome is
      ! that the program goes on forever doing nothing and writing nothing.
      !
#if ! defined(__CUDA)
      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         DO j = 1, idesc(LAX_DESC_NC)
            DO i = 1, idesc(LAX_DESC_NR)
               IF (a(i,j) /= a(i,j)) &
                  CALL errore(' ortho ',' ortho went bananas ',1)
            END DO
         END DO
      END IF
#endif
      RETURN
   END SUBROUTINE

END MODULE local_ortho_memory

!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_gamma_x( cp, ngwx, phi, becp_dist, qbecp, nkbx, bephi, qbephi, &
                           nx0, idesc, diff, iter, n, nss, istart )
!=----------------------------------------------------------------------------=!
      !
#if defined(__CUDA)
      USE cudafor
#endif
      USE kinds,              ONLY: DP
      USE orthogonalize_base, ONLY: rhoset, sigset, tauset, ortho_iterate,   &
                                    use_parallel_diag
      USE mp_global,          ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, my_bgrp_id, inter_bgrp_comm, nbgrp
      USE mp,                 ONLY: mp_sum, mp_bcast
      USE mp_world,           ONLY: mpime

      USE local_ortho_memory

      IMPLICIT  NONE

      include 'laxlib.fh'

      ! ... Arguments

      INTEGER,  INTENT(IN)  :: ngwx, nkbx, nx0
      INTEGER,  INTENT(IN)  :: n, nss, istart
      COMPLEX(DP) :: phi( :, : ), cp( :, : )
      REAL(DP)    :: bephi( :, : )
      REAL(DP)    :: becp_dist( :, : )
      REAL(DP)    :: qbephi( :, : ), qbecp( :, : )
      INTEGER,  INTENT(IN)  :: idesc(:)
      INTEGER,  INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff
#if defined (__CUDA)
      ATTRIBUTES( DEVICE ) :: bephi, becp_dist, phi, cp
#endif

      ! ... Locals

      INTEGER  :: i, j, info, nr, nc, ir, ic
      INTEGER, SAVE :: icnt = 1
      !
      ! ...   Subroutine body
      !

      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         !
         IF( nx0 /= idesc(LAX_DESC_NRCX) ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nx0 ' , nx0 )
         !
         nr = idesc(LAX_DESC_NR)
         nc = idesc(LAX_DESC_NC)
         !
         ir = idesc(LAX_DESC_IR)
         ic = idesc(LAX_DESC_IC)
         !
      ELSE
         !
         nr = 1
         nc = 1
         !
         IF( nx0 /= 1 ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nx0, should be 1 ' , nx0 )
         !
      END IF
      !
      CALL allocate_local_ortho_memory(nss, nx0)
      !
      CALL sync_device_ortho_memory( qbecp, qbephi )
      !
      !     rho = <s'c0|s|cp>
      !
      CALL start_clock( 'rhoset' )
      !
#if defined(__CUDA)
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp_d, n, nss, istart, rhos, rhoa, nx0, idesc )
#else
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, istart, rhos, rhoa, nx0, idesc )
#endif
      !
      CALL stop_clock( 'rhoset' )
      !
      !     sig = 1-<cp|s|cp>
      !
      CALL start_clock( 'sigset' )
#if defined(__CUDA)
      CALL sigset( cp, ngwx, becp_dist, nkbx, qbecp_d, n, nss, istart, sig, nx0, idesc )
#else
      CALL sigset( cp, ngwx, becp_dist, nkbx, qbecp, n, nss, istart, sig, nx0, idesc )
#endif
      CALL stop_clock( 'sigset' )
      !
      !     tau = <s'c0|s|s'c0>
      !
      CALL start_clock( 'tauset' )
#if defined(__CUDA)
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi_d, n, nss, istart, tau, nx0, idesc )
#else
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, istart, tau, nx0, idesc )
#endif
      CALL stop_clock( 'tauset' )
      !
      CALL consistency_check(rhos,idesc) 

      CALL start_clock( 'rsg' )
      !
      ! ...   Diagonalize symmetric part of rho (rhos)
      ! ...   "s" is the matrix of eigenvectors, "rhod" is the array of eigenvalues
      !
      IF( use_parallel_diag ) THEN
         !
#if defined(__CUDA)
         IF( idesc(LAX_DESC_NR) == idesc(LAX_DESC_NC) .AND. idesc(LAX_DESC_NR) == idesc(LAX_DESC_N) ) THEN
            CALL laxlib_diagonalize( nss, rhos, rhod, s, info )
         ELSE IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
            CALL collect_matrix( wrk, rhos, ir, nr, ic, nc, idesc(LAX_DESC_COMM) )
            IF( idesc(LAX_DESC_IC) == 1 .AND. idesc(LAX_DESC_IR) == 1 ) THEN
               CALL laxlib_diagonalize( nss, wrk, rhod, stmp, info )
            END IF 
            CALL distribute_matrix( stmp, s, ir, nr, ic, nc, idesc(LAX_DESC_COMM) )
            CALL mp_bcast( rhod, 0, idesc(LAX_DESC_COMM) )
         END IF
#else
         CALL laxlib_diagonalize( nss, rhos, rhod, s, idesc )
#endif
         !
      ELSE
         !
         IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
            !
            IF( idesc(LAX_DESC_NR) == idesc(LAX_DESC_NC) .AND. idesc(LAX_DESC_NR) == idesc(LAX_DESC_N) ) THEN
               !
               !  rhos and s matrixes, are replicated, no need of collect them
#if defined(__CUDA)
               CALL laxlib_diagonalize( nss, rhos, rhod, s, info )
#else
               s = rhos
               CALL laxlib_diagonalize( nss, s, rhod )
#endif
            ELSE
               !
               CALL collect_matrix( wrk, rhos, ir, nr, ic, nc, idesc(LAX_DESC_COMM) )
               !
#if defined(__CUDA)
               CALL laxlib_diagonalize( nss, wrk, rhod, s, info )
#else
               CALL laxlib_diagonalize( nss, wrk, rhod )
#endif
               !
               CALL distribute_matrix( wrk, s, ir, nr, ic, nc, idesc(LAX_DESC_COMM) )
               !
            END IF
            !
         END IF
         !
      END IF
      !
      CALL stop_clock( 'rsg' )
      CALL start_clock( 'ortho_iter' )
      !
      IF( my_bgrp_id == 0 ) THEN
         !
         !  Matrices and orthogonalization are replicated on all band groups,  there is no
         !  need to keep all processors busy with this task. The processors of the first
         !  group are enough. Moreover replicating the computation across groups could leads
         ! to small numerical differences and weird numerical effects.
         !
         CALL ortho_iterate( iter, diff, s, nx0, rhod, xloc, nx0, sig, rhoa, rhos, tau, nss, idesc)
         !
      END IF
      !
      IF( nbgrp > 1 ) THEN
         !
         !  All groups must have the same lambda matrix, in order to avoid weird
         !  numerical side effects.
         !
         CALL mp_bcast( xloc, 0, inter_bgrp_comm )
         CALL mp_bcast( iter, 0, inter_bgrp_comm )
         CALL mp_bcast( diff, 0, inter_bgrp_comm )
         !
      END IF
      !
      CALL consistency_check( xloc,idesc )

      CALL stop_clock( 'ortho_iter' )

      RETURN

   END SUBROUTINE ortho_gamma_x




!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_x( eigr, cp_bgrp, phi_bgrp, x0, idesc, diff, iter, ccc, bephi, becp_bgrp )
!=----------------------------------------------------------------------------=!
      !
      !     input = cp (non-orthonormal), beta
      !     input = phi |phi>=s'|c0>
      !     output= cp (orthonormal with s( r(t+dt) ) )
      !     output= bephi, becp
      !     the method used is similar to the version in les houches 1988
      !     'simple molecular systems at..'  p. 462-463  (18-22)
      !      xcx + b x + b^t x^t + a = 1
      !     where c = <s'c0|s|s'c0>   b = <s'c0|s cp>   a = <cp|s|cp>
      !     where s=s(r(t+dt)) and s'=s(r(t))  
      !     for vanderbilt pseudo pot - kl & ap
      !
      USE kinds,          ONLY: DP
      USE ions_base,      ONLY: na, nat, nsp, ityp
      USE uspp,           ONLY: nkb, qq_nt, indv_ijkb0, nkbus
      USE uspp_param,     ONLY: nh, upf
      USE electrons_base, ONLY: f, nbsp_bgrp, iupdwn_bgrp, nupdwn_bgrp, i2gupdwn_bgrp, nbsp, nspin, nupdwn, iupdwn
      USE gvecw,          ONLY: ngw
      USE control_flags,  ONLY: iprint, iverbosity, ortho_max
      USE control_flags,  ONLY: force_pairing
      USE io_global,      ONLY: stdout, ionode
      USE cp_interfaces,  ONLY: ortho_gamma, c_bgrp_expand, c_bgrp_pack, nlsm1, collect_bec, beta_eigr, nlsm1us
      USE mp_global,          ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, inter_bgrp_comm! DEBUG
      USE mp_world,           ONLY: mpime
      USE orthogonalize_base, ONLY: bec_bgrp2ortho
      USE mp,                 ONLY : mp_sum
      USE local_ortho_memory, ONLY : xloc, x0_to_xloc, xloc_to_x0
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(IN) :: idesc(:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP) :: cp_bgrp(:,:), phi_bgrp(:,:)
      REAL(DP)    :: x0(:,:,:), diff, ccc
      INTEGER     :: iter
      REAL(DP)    :: bephi(:,:)
      REAL(DP)    :: becp_bgrp(:,:)
      !
      REAL(DP), ALLOCATABLE :: becp_dist(:,:)
      REAL(DP), ALLOCATABLE :: qbephi(:,:,:), qbecp(:,:,:), bec_col(:,:)
      COMPLEX(DP), ALLOCATABLE :: beigr(:,:)
#if defined (__CUDA)
      ATTRIBUTES( DEVICE ) :: becp_bgrp, bephi, becp_dist, beigr, cp_bgrp, phi_bgrp
#endif

      INTEGER :: nkbx
      INTEGER :: info, i, j, iss, iv, jv, ia, is, inl, jnl
      INTEGER :: n1, n2, m1, m2
      INTEGER :: nspin_sub, nx0, ngwx, nrcx
      REAL(DP) :: qqf, dum, byccc
      !
      nkbx = nkb
      ngwx = SIZE( cp_bgrp, 1 )
      !
      nx0 = SIZE( x0, 1 )
      !
      !     calculation of becp and bephi
      !
      CALL start_clock( 'ortho' )

      nrcx = MAXVAL( idesc( LAX_DESC_NRCX, : ) )

      ALLOCATE( becp_dist( nkbx, nrcx*nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating becp_dist ', ABS( info ) )

      IF( nkbus > 0 ) THEN
         !
         ALLOCATE( beigr(ngw,nkb), STAT=info )
         IF( info /= 0 ) &
            CALL errore( ' ortho ', ' allocating beigr ', ABS( info ) )
         !
         CALL beta_eigr ( beigr, 1, nsp, eigr, 2 )
         CALL nlsm1us ( nbsp_bgrp, beigr, phi_bgrp, becp_bgrp )
         CALL bec_bgrp2ortho( becp_bgrp, bephi, nrcx, idesc )
         !
         CALL nlsm1us ( nbsp_bgrp, beigr, cp_bgrp, becp_bgrp )
         CALL bec_bgrp2ortho( becp_bgrp, becp_dist, nrcx, idesc )
         DEALLOCATE( beigr )
         !
      END IF
      !
      !     calculation of qbephi and qbecp
      !
      ALLOCATE( qbephi( nkbx, nx0, nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating qbephi ', ABS( info ) )
      !
      IF( nkbus > 0 ) THEN
         ALLOCATE( bec_col ( nkbx, nrcx*nspin ), STAT = info )
         IF( info /= 0 ) &
            CALL errore( ' ortho ', ' allocating bec_col ', ABS( info ) )
         CALL redist_row2col( nupdwn(1), bephi, bec_col, nkbx, nrcx, idesc(:,1) )
         IF( nspin == 2 ) THEN
            CALL redist_row2col( nupdwn(2), bephi(:,nrcx+1:), bec_col(:,nrcx+1:), nkbx, nrcx, idesc(:,2) )
         END IF
      END IF
      !
      qbephi = 0.d0
      !
      DO iss = 1, nspin
         IF( idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) THEN
!$omp parallel do default(none) &
!$omp shared(nat,ityp,upf,nh,indv_ijkb0,qq_nt,idesc,qbephi,bec_col,iss,nrcx) &
!$omp private(ia,is,iv,inl,jv,jnl,qqf,i)
            DO ia = 1, nat
               is = ityp(ia)
               IF( upf(is)%tvanp ) THEN
                  DO iv=1,nh(is)
                     inl = indv_ijkb0(ia) + iv 
                     DO jv=1,nh(is)
                        jnl = indv_ijkb0(ia) + jv
                        qqf = qq_nt(iv,jv,is)
                        IF( ABS( qqf ) > 1.D-5 ) THEN
                           DO i = 1, idesc( LAX_DESC_NC, iss )
                              qbephi(inl,i,iss) = qbephi(inl,i,iss) + qqf * bec_col(jnl,i+(iss-1)*nrcx)
                           END DO
                        END IF
                     END DO
                  END DO
               END IF
            END DO
!$omp end parallel do 
         ENDIF
      END DO
      !
      ALLOCATE( qbecp ( nkbx, nx0, nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating qbecp ', ABS( info ) )

      qbecp  = 0.d0

      IF( nkbus > 0 ) THEN
         CALL redist_row2col( nupdwn(1), becp_dist, bec_col, nkbx, nrcx, idesc(:,1) )
         IF( nspin == 2 ) THEN
            CALL redist_row2col( nupdwn(2), becp_dist(:,nrcx+1:), bec_col(:,nrcx+1:), nkbx, nrcx, idesc(:,2) )
         END IF
         DO iss = 1, nspin
            IF( idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) THEN
!$omp parallel do default(none) &
!$omp shared(nat,ityp,upf,nh,indv_ijkb0,qq_nt,idesc,qbecp,bec_col,iss,nrcx) &
!$omp private(ia,is,iv,inl,jv,jnl,qqf,i)
               DO ia = 1, nat
                  is = ityp(ia) 
                  IF( upf(is)%tvanp ) THEN
                     DO iv=1,nh(is)
                        inl = indv_ijkb0(ia) + iv
                        DO jv=1,nh(is)
                           jnl = indv_ijkb0(ia) + jv
                           qqf = qq_nt(iv,jv,is)
                           IF( ABS( qqf ) > 1.D-5 ) THEN
                              DO i = 1, idesc( LAX_DESC_NC, iss )
                                 qbecp(inl,i,iss) = qbecp(inl,i,iss) + qqf * bec_col(jnl,i+(iss-1)*nrcx)
                              END DO
                           ENDIF
                        END DO
                     END DO
                  END IF
               END DO
!$omp end parallel do 
            END IF
         END DO
         DEALLOCATE( bec_col )
      END IF
      !
      ! Expand cp and phi to contain all electronic band
      !
      CALL c_bgrp_expand( cp_bgrp )
      CALL c_bgrp_expand( phi_bgrp )
      !
      nspin_sub = nspin 
      if( force_pairing ) nspin_sub = 1
      !
      DO iss = 1, nspin_sub

         CALL x0_to_xloc( x0(:,:,iss), nx0, ccc, idesc(:,iss) )

         CALL ortho_gamma( cp_bgrp, ngwx, phi_bgrp, becp_dist(:,(iss-1)*nrcx+1:iss*nrcx), qbecp(:,:,iss), nkbx, &
                           bephi(:,((iss-1)*nrcx+1):iss*nrcx), &
                           qbephi(:,:,iss), nx0, idesc(:,iss), diff, iter, nbsp, nupdwn(iss), iupdwn(iss) )

         IF( iter > ortho_max ) THEN
            WRITE( stdout, 100 ) diff, iter
            CALL errore('ortho','max number of iterations exceeded',iter)
         END IF

         IF( iverbosity > 1 ) THEN
            WRITE( stdout, 100 ) diff, iter
         ENDIF
         !     
         CALL xloc_to_x0( x0(:,:,iss), nx0, ccc, idesc(:,iss) )
         !
      END DO

      IF( force_pairing ) cp_bgrp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp_bgrp(:,1:nupdwn(2))
      !
      DEALLOCATE( qbecp )
      DEALLOCATE( qbephi )
      DEALLOCATE( becp_dist )
      !
      ! pack cp so that it contains only the bands in the band subgroup
      !
      CALL c_bgrp_pack( cp_bgrp )
      !
      CALL stop_clock( 'ortho' )
      !
      RETURN
      !
100   FORMAT(3X,'diff = ',D18.10,' iter = ', I5 )
      !
   END SUBROUTINE ortho_x




SUBROUTINE qe_sync()
#if defined(__CUDA)
      USE cudafor
#endif
   INTEGER :: info
#if defined (__CUDA)
   info = cudaDeviceSynchronize()
   IF( info /= 0 ) CALL errore('qe_sync',' error ',ABS(info))
#endif
   RETURN
END SUBROUTINE
