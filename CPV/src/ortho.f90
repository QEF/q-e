!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined(__CUDA)
#define DGEMMDRV cublasDgemm
#define DEVICEATTR ,DEVICE
#else
#define DGEMMDRV dgemm
#define DEVICEATTR
#endif

#if defined(__CUDA)
#define PINMEM 
#else
#define PINMEM
#endif

MODULE ortho_module
   !
#if defined(__CUDA)
   USE cudafor
#endif
   USE kinds,              ONLY: DP
   IMPLICIT NONE
   SAVE

   REAL(DP), ALLOCATABLE DEVICEATTR :: s(:,:), sig(:,:), tau(:,:), stmp(:,:)
   REAL(DP), ALLOCATABLE DEVICEATTR :: wrk(:,:), rhoa(:,:), rhos(:,:), rhod(:)
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
        IF( SIZE(x0,1) /= SIZE(xloc,1) .OR. SIZE(x0,2) /= SIZE(xloc,2) ) THEN
           DEALLOCATE(xloc)
        END IF
      END IF
      IF( .NOT. ALLOCATED(xloc) ) THEN
         ALLOCATE( xloc, MOLD=x0, STAT = info )
         IF( info /= 0 ) &
            CALL errore( ' x0_to_xloc ', ' allocating xloc ', ABS( info ) )
      END IF
      IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) THEN
         RETURN
      ENDIF
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
      IF( .NOT. ALLOCATED(xloc) ) THEN
         CALL errore( ' xloc_to_x0 ', ' xloc not allocated ', 1 )
      END IF
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
      CALL mp_bcast( a, 0, comm )
!$cuf kernel do(2) <<<*,*>>>
      DO j = 1, nc
         DO i = 1, nr
            b( i, j ) = a( i + ir - 1, j + ic - 1 )
         END DO
      END DO
      RETURN
   END SUBROUTINE

   SUBROUTINE collect_matrix( a, b, ir, nr, ic, nc, comm )
      USE mp, ONLY: mp_sum
      REAL(DP) DEVICEATTR :: a(:,:), b(:,:)
      INTEGER, INTENT(IN) :: ir, nr, ic, nc, comm
      INTEGER :: i, j, info
      a = 0.0d0
!$cuf kernel do(2) <<<*,*>>>
      DO j = 1, nc
         DO i = 1, nr
            a( ir + i - 1, ic + j - 1 ) = b( i, j )
         END DO
      END DO
      CALL mp_sum( a, comm )
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

!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_gamma( cp, ngwx, phi, becp_dist, qbecp, nkbx, bephi, qbephi, &
                           nx0, idesc, diff, iter, n, nss, istart )
!=----------------------------------------------------------------------------=!
      !
      USE kinds,              ONLY: DP
      USE orthogonalize_base, ONLY: rhoset, sigset, tauset, ortho_iterate,   &
                                    use_parallel_diag
      USE control_flags,      ONLY: diagonalize_on_host
      USE mp_global,          ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, my_bgrp_id, inter_bgrp_comm, nbgrp
      USE mp,                 ONLY: mp_sum, mp_bcast
      USE mp_world,           ONLY: mpime
      USE device_memcpy_m

      IMPLICIT  NONE

      include 'laxlib.fh'

      ! ... Arguments

      INTEGER,  INTENT(IN)  :: ngwx, nkbx, nx0
      INTEGER,  INTENT(IN)  :: n, nss, istart
      COMPLEX(DP) DEVICEATTR :: phi( :, : ), cp( :, : )
      REAL(DP)    DEVICEATTR :: bephi( :, : )
      REAL(DP)    DEVICEATTR :: becp_dist( :, : )
      REAL(DP)    DEVICEATTR :: qbephi( :, : ), qbecp( :, : )
      INTEGER,  INTENT(IN)  :: idesc(:)
      INTEGER,  INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff

      ! ... Locals

      INTEGER  :: i, j, info, nr, nc, ir, ic
      INTEGER, SAVE :: icnt = 1
      REAL(DP), ALLOCATABLE PINMEM :: rhos_h(:,:), s_h(:,:), rhod_h(:)
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
#if defined(__GPU_MPI)
      !
      ! Workaround for a bug in the MPI for GPUs: "mp_root_sum", called by
      ! "rhoset", fails if the dimension of the array is not the same on
      ! all processes, included those who are not in the communicator
      ! The bug is present in v.22.7 and previous (?) of the NVIDIA HPC SDK 
      !
      if ( nx0 == 1 ) THEN
         deallocate( rhos, rhoa, sig, tau )
         allocate( rhos(idesc(LAX_DESC_NRCX),idesc(LAX_DESC_NRCX)) )
         allocate( rhoa(idesc(LAX_DESC_NRCX),idesc(LAX_DESC_NRCX)) )
         allocate( sig (idesc(LAX_DESC_NRCX),idesc(LAX_DESC_NRCX)) )
         allocate( tau (idesc(LAX_DESC_NRCX),idesc(LAX_DESC_NRCX)) )
      end if
#endif
      !
      !     rho = <s'c0|s|cp>
      !
      CALL start_clock( 'rhoset' )
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, istart, rhos, rhoa, nx0, idesc )
      CALL stop_clock( 'rhoset' )
      !
      !     sig = 1-<cp|s|cp>
      !
      CALL start_clock( 'sigset' )
      CALL sigset( cp, ngwx, becp_dist, nkbx, qbecp, n, nss, istart, sig, nx0, idesc )
      CALL stop_clock( 'sigset' )
      !
      !     tau = <s'c0|s|s'c0>
      !
      CALL start_clock( 'tauset' )
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, istart, tau, nx0, idesc )
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
            IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
               CALL laxlib_diagonalize( nss, rhos, rhod, s, info )
            END IF
         ELSE IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
            IF( diagonalize_on_host ) THEN  !  tune here
               ALLOCATE( rhos_h, SOURCE = rhos )
               ALLOCATE( rhod_h, MOLD = rhod )
               ALLOCATE( s_h, MOLD = s )
               CALL laxlib_diagonalize( nss, rhos_h, rhod_h, s_h, idesc )
               CALL dev_memcpy( s, s_h )
               CALL dev_memcpy( rhod, rhod_h )
               DEALLOCATE( rhos_h, rhod_h, s_h )
            ELSE
               CALL collect_matrix( wrk, rhos, ir, nr, ic, nc, idesc(LAX_DESC_COMM) )
               IF( idesc(LAX_DESC_IC) == 1 .AND. idesc(LAX_DESC_IR) == 1 ) THEN
                  CALL laxlib_diagonalize( nss, wrk, rhod, stmp, info )
               END IF 
               CALL distribute_matrix( stmp, s, ir, nr, ic, nc, idesc(LAX_DESC_COMM) )
               CALL mp_bcast( rhod, 0, idesc(LAX_DESC_COMM) )
            END IF
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
   END SUBROUTINE ortho_gamma

   !

   SUBROUTINE compute_qs_times_betas( bephi, bec_row, qbephi, qbecp, idesc )
      USE uspp,           ONLY: nkb, qq_nt, qq_nt_d, ofsbeta, nkbus
      USE uspp_param,     ONLY: nh, upf
      USE electrons_base, ONLY: nspin, nbsp_bgrp, iupdwn_bgrp, nupdwn_bgrp, nbsp, nupdwn, iupdwn
      USE ions_base,      ONLY: na, nat, nsp, ityp
#if defined (__CUDA)
      USE cublas
#endif
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      REAL(DP), INTENT(OUT) DEVICEATTR :: qbephi(:,:,:), qbecp(:,:,:)
      REAL(DP), INTENT(IN)  DEVICEATTR :: bephi(:,:)
      REAL(DP), INTENT(IN)  DEVICEATTR :: bec_row(:,:)
      INTEGER, INTENT(IN) :: idesc(:,:)

      REAL(DP), ALLOCATABLE DEVICEATTR :: bec_col(:,:), bephi_col(:,:)

      INTEGER :: nkbx, info, nrcx, iss, is, jv, iv, ia
      INTEGER :: i, j, inl, jnl, indv, nc, nhs

      nkbx = nkb
      nrcx = idesc( LAX_DESC_NRCX, 1 )
      IF( nspin > 1 ) nrcx = MAX( nrcx, idesc( LAX_DESC_NRCX, 2 ) )

      qbephi = 0.d0
      qbecp  = 0.d0
 
      IF( nkbus > 0 ) THEN
         !
         ALLOCATE( bec_col ( nkbx, nrcx*nspin ), STAT = info )
         IF( info /= 0 ) &
            CALL errore( ' compute_qs_times_betas ', ' allocating bec_col ', ABS( info ) )
         ALLOCATE( bephi_col ( nkbx, nrcx*nspin ), STAT = info )
         IF( info /= 0 ) &
            CALL errore( ' compute_qs_times_betas ', ' allocating bephi_col ', ABS( info ) )
         !
         CALL redist_row2col( nupdwn(1), bephi, bephi_col, nkbx, nrcx, idesc(:,1) )
         CALL redist_row2col( nupdwn(1), bec_row, bec_col, nkbx, nrcx, idesc(:,1) )
         IF( nspin == 2 ) THEN
            CALL redist_row2col( nupdwn(2), bephi(:,nrcx+1:), bephi_col(:,nrcx+1:), nkbx, nrcx, idesc(:,2) )
            CALL redist_row2col( nupdwn(2), bec_row(:,nrcx+1:), bec_col(:,nrcx+1:), nkbx, nrcx, idesc(:,2) )
         END IF
         !
         DO iss = 1, nspin
            IF( idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) THEN
               nc = idesc( LAX_DESC_NC, iss )
               DO ia = 1, nat
                  is = ityp(ia)
                  IF( upf(is)%tvanp ) THEN
                     indv = ofsbeta(ia)
                     nhs  = nh(is)
#if defined (__CUDA)
                     CALL DGEMMDRV('N', 'N', nhs, nc, nhs, 1.0d0, qq_nt_d(1,1,is), SIZE(qq_nt_d,1), &
                                   bephi_col(indv+1,(iss-1)*nrcx+1), SIZE(bephi_col,1), 0.0d0, qbephi(indv+1,1,iss), SIZE(qbephi,1))
                     CALL DGEMMDRV('N', 'N', nhs, nc, nhs, 1.0d0, qq_nt_d(1,1,is), SIZE(qq_nt_d,1), &
                                   bec_col(indv+1,(iss-1)*nrcx+1), SIZE(bec_col,1), 0.0d0, qbecp(indv+1,1,iss), SIZE(qbecp,1))
#else
                     CALL DGEMMDRV('N', 'N', nhs, nc, nhs, 1.0d0, qq_nt(1,1,is), SIZE(qq_nt,1), &
                                   bephi_col(indv+1,(iss-1)*nrcx+1), SIZE(bephi_col,1), 0.0d0, qbephi(indv+1,1,iss), SIZE(qbephi,1))
                     CALL DGEMMDRV('N', 'N', nhs, nc, nhs, 1.0d0, qq_nt(1,1,is), SIZE(qq_nt,1), &
                                   bec_col(indv+1,(iss-1)*nrcx+1), SIZE(bec_col,1), 0.0d0, qbecp(indv+1,1,iss), SIZE(qbecp,1))
#endif
! !$cuf kernel do (2)
!                     DO iv=1,nhs
!                        DO i = 1, nc
!                           DO jv = 1, nhs
!                              qbephi(indv+iv,i,iss) = qbephi(indv+iv,i,iss) + &
!                                                         qq_nt_d(iv,jv,is) * bephi_col(indv+jv,i+(iss-1)*nrcx)
!                              qbecp(indv+iv,i,iss) = qbecp(indv+iv,i,iss) + &
!                                                        qq_nt_d(iv,jv,is) * bec_col(indv+jv,i+(iss-1)*nrcx)
!                           END DO
!                        END DO
!                     END DO
                  END IF
               END DO
            ENDIF
         END DO
         DEALLOCATE( bec_col )
         DEALLOCATE( bephi_col )
      END IF
   END SUBROUTINE compute_qs_times_betas

   SUBROUTINE keep_only_us(wrk)
      USE uspp,           ONLY: ofsbeta
      USE uspp_param,     ONLY: nh, upf
      USE ions_base,      ONLY: na, nat, nsp, ityp
#if defined (__CUDA)
      USE cublas
#endif
      IMPLICIT NONE
      COMPLEX(DP) DEVICEATTR, INTENT(OUT) :: wrk(:,:)
      INTEGER :: ia, is, inl, nhs, iv
      DO ia = 1, nat
         is  = ityp(ia)
         inl = ofsbeta(ia)
         nhs = nh(is)
         IF( .NOT. upf(is)%tvanp ) THEN
!$cuf kernel do (1)
            DO iv = 1, nhs
               wrk(:,iv+inl) = 0.0d0
            END DO
         END IF
      END DO
   END SUBROUTINE

END MODULE ortho_module


!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_x( betae, cp_bgrp, phi_bgrp, x0, idesc, diff, iter, ccc, bephi, becp_bgrp )
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
      USE uspp,           ONLY: nkb, nkbus
      USE ions_base,      ONLY: nsp
      USE electrons_base, ONLY: f, nbsp_bgrp, iupdwn_bgrp, nupdwn_bgrp, nbsp, nspin, nupdwn, iupdwn
      USE gvecw,          ONLY: ngw
      USE control_flags,  ONLY: iprint, iverbosity, ortho_max
      USE control_flags,  ONLY: force_pairing
      USE io_global,      ONLY: stdout, ionode
      USE cp_interfaces,  ONLY: c_bgrp_expand, c_bgrp_pack, nlsm1, collect_bec, nlsm1us
      USE mp_global,          ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, inter_bgrp_comm! DEBUG
      USE mp_world,           ONLY: mpime
      USE orthogonalize_base, ONLY: bec_bgrp2ortho
      USE mp,                 ONLY : mp_sum
      USE ortho_module
      USE device_memcpy_m
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(IN) :: idesc(:,:)
      COMPLEX(DP) DEVICEATTR :: betae(:,:)
      COMPLEX(DP) DEVICEATTR :: cp_bgrp(:,:), phi_bgrp(:,:)
      REAL(DP)    :: x0(:,:,:), diff, ccc
      INTEGER     :: iter
      REAL(DP)    DEVICEATTR :: bephi(:,:)
      REAL(DP)    DEVICEATTR :: becp_bgrp(:,:)
      !
      REAL(DP),    ALLOCATABLE DEVICEATTR :: bec_row(:,:)
      REAL(DP),    ALLOCATABLE DEVICEATTR :: qbephi(:,:,:), qbecp(:,:,:)
      COMPLEX(DP), ALLOCATABLE DEVICEATTR :: wrk2(:,:)

      INTEGER :: nkbx, info, iss, nspin_sub, nx0, ngwx, nrcx
      !
      CALL start_clock( 'ortho' )
      !
      nkbx = nkb
      ngwx = SIZE( cp_bgrp, 1 )
      nx0  = SIZE( x0, 1 )
      nrcx = MAXVAL( idesc( LAX_DESC_NRCX, : ) )
      !
      !     calculation of becp and bephi
      !
      ALLOCATE( bec_row( nkbx, nrcx*nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating bec_row ', ABS( info ) )

      IF( nkbus > 0 ) THEN
         !
         ALLOCATE( wrk2, MOLD = betae  )
         CALL dev_memcpy( wrk2, betae )
         CALL keep_only_us( wrk2 ) 
         CALL nlsm1us ( nbsp_bgrp, wrk2, phi_bgrp, becp_bgrp )
         CALL bec_bgrp2ortho( becp_bgrp, bephi, nrcx, idesc )
         !
         CALL nlsm1us ( nbsp_bgrp, wrk2, cp_bgrp, becp_bgrp )
         CALL bec_bgrp2ortho( becp_bgrp, bec_row, nrcx, idesc )
         DEALLOCATE( wrk2 )
         !
      END IF
      !
      !     calculation of qbephi and qbecp
      !
      ALLOCATE( qbephi( nkbx, nx0, nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating qbephi ', ABS( info ) )
      ALLOCATE( qbecp ( nkbx, nx0, nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating qbecp ', ABS( info ) )
      !
      CALL compute_qs_times_betas( bephi, bec_row, qbephi, qbecp, idesc )
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

         CALL ortho_gamma( cp_bgrp, ngwx, phi_bgrp, bec_row(:,(iss-1)*nrcx+1:iss*nrcx), qbecp(:,:,iss), nkbx, &
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

      IF( force_pairing ) THEN
         CALL dev_memcpy(cp_bgrp(:,iupdwn(2):), cp_bgrp(:,1:),  [1, ngw], 1 , [1, nupdwn(2)], 1) 
      END IF
      !
      DEALLOCATE( qbecp )
      DEALLOCATE( qbephi )
      DEALLOCATE( bec_row )
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
