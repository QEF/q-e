!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_gamma_x( iopt, cp, ngwx, phi, becp_dist, qbecp, nkbx, bephi, qbephi, &
                           x0, nx0, descla, diff, iter, n, nss, istart )
!=----------------------------------------------------------------------------=!
      !
      USE kinds,              ONLY: DP
      USE orthogonalize_base, ONLY: rhoset, sigset, tauset, ortho_iterate,   &
                                    ortho_alt_iterate, use_parallel_diag
      USE dspev_module,       ONLY: diagonalize_serial, diagonalize_parallel
      USE descriptors,        ONLY: la_descriptor
      USE mp_global,          ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, my_bgrp_id, inter_bgrp_comm, nbgrp
      USE mp,                 ONLY: mp_sum, mp_bcast

      IMPLICIT  NONE

      ! ... Arguments

      INTEGER,  INTENT(IN)  :: iopt
      INTEGER,  INTENT(IN)  :: ngwx, nkbx, nx0
      INTEGER,  INTENT(IN)  :: n, nss, istart
      COMPLEX(DP) :: phi( ngwx, n ), cp( ngwx, n )
      REAL(DP)    :: bephi( :, : )
      REAL(DP)    :: becp_dist( :, : )
      REAL(DP)    :: qbephi( :, : ), qbecp( :, : )
      REAL(DP)    :: x0( nx0, nx0 )
      TYPE(la_descriptor),  INTENT(IN)  :: descla
      INTEGER,  INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff
      ! ... Locals

      REAL(DP),   ALLOCATABLE :: s(:,:), sig(:,:), tau(:,:), rhot(:,:)
      REAL(DP),   ALLOCATABLE :: wrk(:,:), rhoa(:,:), rhos(:,:), rhod(:)
      INTEGER  :: i, j, info, nr, nc, ir, ic
      !
      ! ...   Subroutine body
      !
      IF( descla%active_node > 0 ) THEN
         !
         IF( nx0 /= descla%nrcx ) &
            CALL errore( ' ortho_gamma ', ' inconsistent dimensions nx0 ' , nx0 )
         !
         nr = descla%nr
         nc = descla%nc
         !
         ir = descla%ir
         ic = descla%ic
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
      ALLOCATE( rhos( nx0, nx0 ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating rhos ', ABS( info ) )
      ALLOCATE( rhoa( nx0, nx0 ), STAT = info )   !   antisymmetric part of rho
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating rhoa ', ABS( info ) )
      ALLOCATE( s( nx0, nx0 ), STAT = info ) 
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating s ', ABS( info ) )
      ALLOCATE( sig( nx0, nx0 ), STAT = info ) 
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating sig ', ABS( info ) )
      ALLOCATE( tau( nx0, nx0 ), STAT = info ) 
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating tau ', ABS( info ) )
      !
      ALLOCATE( rhod( nss ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho_gamma ', ' allocating tau ', ABS( rhod ) )
      !
      !     rho = <s'c0|s|cp>
      !
      CALL start_clock( 'rhoset' )
      !
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, istart, rhos, nx0, descla )
      !
      IF( descla%active_node > 0 ) THEN
         !
         ALLOCATE( rhot( nx0, nx0 ), STAT = info )   !   transpose of rho
         IF( info /= 0 ) &
            CALL errore( ' ortho_gamma ', ' allocating rhot ', ABS( rhod ) )
         !
         !    distributed array rhos contains "rho", 
         !    now transpose rhos and store the result in distributed array rhot
         !
         CALL sqr_tr_cannon( nss, rhos, nx0, rhot, nx0, descla )
         !
         !  Compute the symmetric part of rho
         !
         DO j = 1, nc
            DO i = 1, nr
               rhos( i, j ) = 0.5d0 * ( rhos( i, j ) + rhot( i, j ) )
            END DO
         END DO
         !
         !    distributed array rhos now contains symmetric part of "rho", 
         !
         CALL consistency_check( rhos )
         !
         !  Antisymmetric part of rho, alredy distributed across ortho procs.
         !
         DO j = 1, nc
            DO i = 1, nr
               rhoa( i, j ) = rhos( i, j ) - rhot( i, j )
            END DO
         END DO
         !
         DEALLOCATE( rhot )
         !
      END IF

      CALL stop_clock( 'rhoset' )


      CALL start_clock( 'rsg' )
      !
      ! ...   Diagonalize symmetric part of rho (rhos)
      ! ...   "s" is the matrix of eigenvectors, "rhod" is the array of eigenvalues
      !
      IF( use_parallel_diag ) THEN
         !
         CALL diagonalize_parallel( nss, rhos, rhod, s, descla )
         !
      ELSE
         !
         IF( descla%active_node > 0 ) THEN
            !
            ALLOCATE( wrk( nss, nss ), STAT = info )
            IF( info /= 0 ) CALL errore( ' ortho_gamma ', ' allocating wrk ', 1 )
            !
            CALL collect_matrix( wrk, rhos )
            !
            CALL diagonalize_serial( nss, wrk, rhod )
            !
            CALL distribute_matrix( wrk, s )
            !
            DEALLOCATE( wrk )
            !
         END IF
         !
      END IF
      !
      CALL stop_clock( 'rsg' )
      !
      !     sig = 1-<cp|s|cp>
      !
      CALL start_clock( 'sigset' )
      CALL sigset( cp, ngwx, becp_dist, nkbx, qbecp, n, nss, istart, sig, nx0, descla )
      CALL stop_clock( 'sigset' )
      !
      !     tau = <s'c0|s|s'c0>
      !
      CALL start_clock( 'tauset' )
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, istart, tau, nx0, descla )
      CALL stop_clock( 'tauset' )
      !
      CALL start_clock( 'ortho_iter' )
      !
      IF( my_bgrp_id == 0 ) THEN
         !
         !  Matrices and orthogonalization are replicated on all band groups,  there is no
         !  need to keep all processors busy with this task. The processors of the first
         !  group are enough. Moreover replicating the computation across groups could leads
         ! to small numerical differences and weird numerical effects.
         !
         IF( iopt == 0 ) THEN
            !
            CALL ortho_iterate( iter, diff, s, nx0, rhod, x0, nx0, sig, rhoa, rhos, tau, nss, descla)
            !
         ELSE
            !
            CALL ortho_alt_iterate( iter, diff, s, nx0, rhod, x0, nx0, sig, rhoa, tau, nss, descla)
            !
         END IF
         !
      END IF
      !
      IF( nbgrp > 1 ) THEN
         !
         !  All groups must have the same lambda matrix, in order to avoid weird
         !  numerical side effects.
         !
         CALL mp_bcast( x0, 0, inter_bgrp_comm )
         CALL mp_bcast( iter, 0, inter_bgrp_comm )
         CALL mp_bcast( diff, 0, inter_bgrp_comm )
         !
      END IF
      !
      CALL stop_clock( 'ortho_iter' )
      !
      DEALLOCATE( rhoa, rhos, rhod, s, sig, tau )
      !
      IF( descla%active_node > 0 )  CALL consistency_check( x0 )

      RETURN

   CONTAINS

      SUBROUTINE distribute_matrix( a, b )
         REAL(DP) :: a(:,:), b(:,:)
         INTEGER :: i, j
         IF( descla%active_node > 0 ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  b( i, j ) = a( i + ir - 1, j + ic - 1 )
               END DO
            END DO
         END IF
         RETURN
      END SUBROUTINE

      SUBROUTINE collect_matrix( a, b )
         REAL(DP) :: a(:,:), b(:,:)
         INTEGER :: i, j
         a = 0.0d0
         IF( descla%active_node > 0 ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  a( ir + i - 1, ic + j - 1 ) = b( i, j )
               END DO
            END DO
         END IF
         CALL mp_sum( a, descla%comm )
         RETURN
      END SUBROUTINE

      SUBROUTINE consistency_check( a )
         REAL(DP) :: a(:,:)
         INTEGER :: i, j
         !
         ! on some machines (IBM RS/6000 for instance) the following test allows
         ! to distinguish between Numbers and Sodium Nitride (NaN, Not a Number).
         ! If a matrix of Not-Numbers is passed to rs, the most likely outcome is
         ! that the program goes on forever doing nothing and writing nothing.
         !
         DO j = 1, nc
            DO i = 1, nr
               IF (a(i,j) /= a(i,j)) CALL errore(' ortho ',' ortho went bananas ',1)
            END DO
         END DO
         RETURN
      END SUBROUTINE

   END SUBROUTINE ortho_gamma_x




!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_x( eigr, cp_bgrp, phi_bgrp, x0, descla, diff, iter, ccc, bephi, becp_bgrp )
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
      USE ions_base,      ONLY: na, nat
      USE uspp,           ONLY: nkb, qq
      USE uspp_param,     ONLY: nh, ish, nvb
      USE electrons_base, ONLY: f, nbsp_bgrp, iupdwn_bgrp, nupdwn_bgrp, i2gupdwn_bgrp, nbsp, nspin, nupdwn, iupdwn
      USE gvecw,          ONLY: ngw
      USE control_flags,  ONLY: iprint, iverbosity, ortho_max
      USE control_flags,  ONLY: force_pairing
      USE io_global,      ONLY: stdout, ionode
      USE cp_interfaces,  ONLY: ortho_gamma, c_bgrp_expand, c_bgrp_pack, nlsm1, collect_bec
      USE descriptors,    ONLY: la_descriptor
      USE mp_global,          ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, inter_bgrp_comm ! DEBUG
      USE orthogonalize_base, ONLY: bec_bgrp2ortho
      USE mp,                 ONLY : mp_sum
      !
      IMPLICIT NONE
      !
      TYPE(la_descriptor), INTENT(IN) :: descla(:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP) :: cp_bgrp(:,:), phi_bgrp(:,:)
      REAL(DP)    :: x0(:,:,:), diff, ccc
      INTEGER     :: iter
      REAL(DP)    :: bephi(:,:)
      REAL(DP)    :: becp_bgrp(:,:)
      !
      REAL(DP), ALLOCATABLE :: xloc(:,:), becp_dist(:,:)
      REAL(DP), ALLOCATABLE :: qbephi(:,:,:), qbecp(:,:,:), bec_col(:,:)

      INTEGER :: nkbx
      INTEGER :: info, i, j, iss, iv, jv, ia, is, inl, jnl
      INTEGER :: n1, n2, m1, m2
      INTEGER :: nspin_sub, nx0, ngwx, nrcx
      REAL(DP) :: qqf, dum
      !
      nkbx = nkb
      ngwx = SIZE( cp_bgrp, 1 )
      !
      nx0 = SIZE( x0, 1 )
      !
      !     calculation of becp and bephi
      !
      CALL start_clock( 'ortho' )

      nrcx = MAXVAL( descla( : )%nrcx )

      ALLOCATE( becp_dist( nkbx, nrcx*nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating becp_dist ', ABS( info ) )

      IF( nvb > 0 ) THEN
         !
         becp_bgrp = 0.0d0
         !
         CALL nlsm1 ( nbsp_bgrp, 1, nvb, eigr, phi_bgrp, becp_bgrp )
         CALL bec_bgrp2ortho( becp_bgrp, bephi, nrcx, descla )
         !
         becp_bgrp = 0.0d0
         !
         CALL nlsm1 ( nbsp_bgrp, 1, nvb, eigr, cp_bgrp, becp_bgrp )
         CALL bec_bgrp2ortho( becp_bgrp, becp_dist, nrcx, descla )
         !
      END IF
      !
      !     calculation of qbephi and qbecp
      !
      ALLOCATE( qbephi( nkbx, nx0, nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating qbephi ', ABS( info ) )
      !
      IF( nvb > 0 ) THEN
         ALLOCATE( bec_col ( nkbx, nrcx*nspin ), STAT = info )
         IF( info /= 0 ) &
            CALL errore( ' ortho ', ' allocating bec_col ', ABS( info ) )
         CALL redist_row2col( nupdwn(1), bephi, bec_col, nkbx, nrcx, descla(1) )
         IF( nspin == 2 ) THEN
            CALL redist_row2col( nupdwn(2), bephi(1,nrcx+1), bec_col(1,nrcx+1), nkbx, nrcx, descla(2) )
         END IF
      END IF
      !
      qbephi = 0.d0
      !
      DO is=1,nvb
         DO iv=1,nh(is)
            inl = ish(is)+(iv-1)*na(is)
            DO jv=1,nh(is)
               jnl = ish(is)+(jv-1)*na(is)
               qqf = qq(iv,jv,is)
               IF( ABS( qqf ) > 1.D-5 ) THEN
                  DO iss = 1, nspin
                     IF( descla( iss )%active_node > 0 ) THEN
                        DO i = 1, descla( iss )%nc
                           CALL daxpy( na(is), qqf, bec_col(jnl+1,i+(iss-1)*nrcx),1,qbephi(inl+1,i,iss), 1 )
                        END DO
                     END IF
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      !
      ALLOCATE( qbecp ( nkbx, nx0, nspin ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating qbecp ', ABS( info ) )

      qbecp  = 0.d0

      IF( nvb > 0 ) THEN
         CALL redist_row2col( nupdwn(1), becp_dist, bec_col, nkbx, nrcx, descla(1) )
         IF( nspin == 2 ) THEN
            CALL redist_row2col( nupdwn(2), becp_dist(1,nrcx+1), bec_col(1,nrcx+1), nkbx, nrcx, descla(2) )
         END IF
      END IF

      DO is=1,nvb
         DO iv=1,nh(is)
            inl = ish(is)+(iv-1)*na(is)
            DO jv=1,nh(is)
               jnl = ish(is)+(jv-1)*na(is)
               qqf = qq(iv,jv,is)
               IF( ABS( qqf ) > 1.D-5 ) THEN
                  DO iss = 1, nspin
                     IF( descla( iss )%active_node > 0 ) THEN
                        DO i = 1, descla( iss )%nc
                           CALL daxpy( na(is), qqf, bec_col(jnl+1,i+(iss-1)*nrcx),1, qbecp(inl+1,i,iss), 1 )
                        END DO
                     END IF
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      !
      IF( nvb > 0 ) DEALLOCATE( bec_col )
      !
      ! Expand cp and phi to contain all electronic band
      !
      CALL c_bgrp_expand( cp_bgrp )
      CALL c_bgrp_expand( phi_bgrp )
      !
      ALLOCATE( xloc( nx0, nx0 ), STAT = info )
      IF( info /= 0 ) &
         CALL errore( ' ortho ', ' allocating xloc ', ABS( info ) )
      !
      nspin_sub = nspin 
      if( force_pairing ) nspin_sub = 1
      !
      DO iss = 1, nspin_sub

         IF( descla( iss )%active_node > 0 ) xloc = x0(:,:,iss) * ccc

         CALL ortho_gamma( 0, cp_bgrp, ngwx, phi_bgrp, becp_dist(:,(iss-1)*nrcx+1:iss*nrcx), qbecp(:,:,iss), nkbx, &
                           bephi(:,((iss-1)*nrcx+1):iss*nrcx), &
                           qbephi(:,:,iss), xloc, nx0, descla(iss), diff, iter, nbsp, nupdwn(iss), iupdwn(iss) )

         IF( iter > ortho_max ) THEN
            WRITE( stdout, 100 ) diff, iter
            CALL errore('ortho','max number of iterations exceeded',iter)
         END IF

         IF( iverbosity > 1 ) THEN
            WRITE( stdout, 100 ) diff, iter
         ENDIF
         !     
         IF( descla( iss )%active_node > 0 ) x0( :, :, iss ) = xloc / ccc
         !
      END DO

      IF( force_pairing ) cp_bgrp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp_bgrp(:,1:nupdwn(2))
      !
      DEALLOCATE( xloc )
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
