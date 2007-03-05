!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"



!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_m( c0, cp, lambda, descla, ccc, nupdwn, iupdwn, nspin )
      !
      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: force_pairing
      USE cp_main_variables,  ONLY: ema0bg
      USE descriptors,        ONLY: lambda_node_ , nlar_ , nlac_
      USE control_flags,      ONLY: ortho_eps, ortho_max
      USE orthogonalize_base, ONLY: calphi, updatc
      USE cp_interfaces,      ONLY: ortho_gamma
      !
      IMPLICIT NONE

      INTEGER,     INTENT(IN)    :: descla(:,:)
      INTEGER,     INTENT(IN)    :: nupdwn(:), iupdwn(:), nspin
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:), cp(:,:)
      REAL(DP),    INTENT(INOUT) :: lambda(:,:,:)
      REAL(DP),    INTENT(IN)    :: ccc
      !
      COMPLEX(DP), ALLOCATABLE :: phi(:,:)
      INTEGER                  :: iss, nsc, iwfc, nwfc, info
      INTEGER                  :: iter, i, j
      INTEGER                  :: ngwx, n, nr, nc, nx
      REAL(DP)                 :: diff
      REAL(DP)                 :: dum(2,2)
      COMPLEX(DP)              :: cdum(2,2)
      !
      CALL start_clock( 'ortho' )  

      n    = SIZE( c0, 2 )
      ngwx = SIZE( c0, 1 )
      nx   = SIZE( lambda, 1 )

      ALLOCATE( phi( ngwx, n ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho ', ' allocating phi ', 3 )

      CALL calphi( c0, ngwx, dum, 1, cdum, phi, n, ema0bg )
      !
      nsc = nspin
      IF( force_pairing ) nsc = 1
      !
      DO iss = 1, nsc
          !
          nwfc = nupdwn(iss)
          iwfc = iupdwn(iss)
          !
          CALL ortho_gamma( 1, cp, ngwx, phi, dum, dum, 2, dum, dum, lambda(:,:,iss), nx, &
               descla(:,iss), diff, iter, n, nwfc, iwfc )
          !
          IF ( iter > ortho_max ) THEN
             call errore(' ortho ','  itermax ',iter)
          END IF
          !
          CALL updatc( 1.0d0, n, lambda(:,:,iss), nx, phi, ngwx, dum, 1, dum, dum, cp, &
                       nwfc, iwfc, descla(:,iss) )
          !     
          !     lagrange multipliers
          !
          IF( descla( lambda_node_ , iss ) > 0 ) THEN
             DO j = 1, descla( nlac_ , iss )
                DO i = 1, descla( nlar_ , iss )
                   lambda( i, j, iss ) = lambda( i, j, iss ) / ccc
                END DO
             END DO
          END IF
          !
      END DO
      !
      IF( force_pairing ) cp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp(:,1:nupdwn(2))
      !
      DEALLOCATE( phi )
      !
      CALL stop_clock( 'ortho' )
      !
      RETURN
   END SUBROUTINE ortho_m




!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_gamma_x( iopt, cp, ngwx, phi, becp, qbecp, nkbx, bephi, qbephi, &
                           x0, nx0, descla, diff, iter, n, nss, istart )
!=----------------------------------------------------------------------------=!
      !
      USE kinds,              ONLY: DP
      USE orthogonalize_base, ONLY: rhoset, sigset, tauset, ortho_iterate,   &
                                    ortho_alt_iterate, diagonalize_serial,   &
                                    use_parallel_diag, diagonalize_parallel
      USE descriptors,        ONLY: lambda_node_ , nlar_ , nlac_ , ilar_ , ilac_
      USE mp_global,          ONLY: nproc_image, me_image, intra_image_comm
      USE mp,                 ONLY: mp_sum

      IMPLICIT  NONE

      ! ... Arguments

      INTEGER,  INTENT(IN)  :: iopt
      INTEGER,  INTENT(IN)  :: ngwx, nkbx, nx0
      INTEGER,  INTENT(IN)  :: n, nss, istart
      COMPLEX(DP) :: phi( ngwx, n ), cp( ngwx, n )
      REAL(DP)    :: bephi( nkbx, n ), becp( nkbx, n )
      REAL(DP)    :: qbephi( nkbx, nx0 ), qbecp( nkbx, nx0 )
      REAL(DP)    :: x0( nx0, nx0 )
      INTEGER,  INTENT(IN)  :: descla( : )
      INTEGER,  INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff

      ! ... Locals

      REAL(DP),   ALLOCATABLE :: s(:,:), sig(:,:), tau(:,:), rhot(:,:)
      REAL(DP),   ALLOCATABLE :: wrk(:,:), rhoa(:,:), rhos(:,:), rhod(:)
      INTEGER  :: i, j, info, nr, nc, ir, ic
      LOGICAL  :: iter_node
      !
      ! ...   Subroutine body
      !
      IF( descla( lambda_node_ ) > 0 ) THEN
         !
         iter_node = .TRUE.
         !
         nr = descla( nlar_ )
         nc = descla( nlac_ )
         !
         ir = descla( ilar_ )
         ic = descla( ilac_ )
         !
	 ALLOCATE( rhos( nr, nc ) )
         ALLOCATE( rhoa( nr, nc ) )   !   antisymmetric part of rho
         ALLOCATE( s( nr, nc ) ) 
         ALLOCATE( sig( nr, nc ) ) 
         ALLOCATE( tau( nr, nc ) ) 
         !
      ELSE
         !
         iter_node = .FALSE.
         !
	 ALLOCATE( rhos( 1, 1 ) )
         ALLOCATE( rhoa( 1, 1 ) ) 
         ALLOCATE( s( 1, 1 ) ) 
         ALLOCATE( sig( 1, 1 ) ) 
         ALLOCATE( tau( 1, 1 ) ) 
         !
      END IF
      !
      ALLOCATE( rhod( nss ) )
      !
      !     rho = <s'c0|s|cp>
      !
      CALL start_clock( 'rhoset' )
      !
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, istart, rhos, nr, descla )
      !
      IF( iter_node ) THEN
         !
         ALLOCATE( rhot( nr, nc ) )   !   transpose of rho
         !
         !    distributed array rhos contains "rho", 
         !    now transpose rhos and store the result in distributed array rhot
         !
         CALL sqr_tr_cannon( nss, rhos, nr, rhot, nr, descla )
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
         ALLOCATE( wrk( nss, nss ), STAT = info )
         IF( info /= 0 ) CALL errore( ' ortho ', ' allocating matrixes ', 1 )
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
      CALL stop_clock( 'rsg' )
      !
      !     sig = 1-<cp|s|cp>
      !
      CALL sigset( cp, ngwx, becp, nkbx, qbecp, n, nss, istart, sig, nr, descla )
      !
      !     tau = <s'c0|s|s'c0>
      !
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, istart, tau, nr, descla )
      !
      CALL start_clock( 'ortho_iter' )
      !
      IF( iopt == 0 ) THEN
         !
         CALL ortho_iterate( iter, diff, s, nr, rhod, x0, nx0, sig, rhoa, rhos, tau, nss, descla)
         !
      ELSE
         !
         CALL ortho_alt_iterate( iter, diff, s, nr, rhod, x0, nx0, sig, rhoa, tau, nss, descla)
         !
      END IF
      !
      CALL stop_clock( 'ortho_iter' )
      !
      DEALLOCATE( rhoa, rhos, rhod, s, sig, tau )
      !
      IF( iter_node )  CALL consistency_check( x0 )

      RETURN

   CONTAINS

      SUBROUTINE distribute_matrix( a, b )
         REAL(DP) :: a(:,:), b(:,:)
         INTEGER :: i, j
         IF( iter_node ) THEN
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
         IF( iter_node ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  a( ir + i - 1, ic + j - 1 ) = b( i, j )
               END DO
            END DO
         END IF
         CALL mp_sum( a, intra_image_comm )
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
   SUBROUTINE ortho_cp( eigr, cp, phi, ngwx, x0, descla, diff, iter, ccc, &
                        bephi, becp, nbsp, nspin, nupdwn, iupdwn )
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
      USE cvan,           ONLY: ish, nvb
      USE uspp,           ONLY: nkb, qq
      USE uspp_param,     ONLY: nh
      USE electrons_base, ONLY: f
      USE gvecw,          ONLY: ngw
      USE control_flags,  ONLY: iprint, iprsta, ortho_max
      USE control_flags,  ONLY: force_pairing
      USE io_global,      ONLY: stdout, ionode
      USE cp_interfaces,  ONLY: ortho_gamma
      USE descriptors,    ONLY: nlac_ , ilac_ , lambda_node_
      !
      IMPLICIT NONE
      !
      INTEGER,    INTENT(IN) :: ngwx, nbsp, nspin
      INTEGER,    INTENT(IN) :: nupdwn( nspin ), iupdwn( nspin )
      INTEGER,    INTENT(IN) :: descla(:,:)
      COMPLEX(DP) :: cp(ngwx,nbsp), phi(ngwx,nbsp), eigr(ngwx,nat)
      REAL(DP)    :: x0(:,:,:), diff, ccc
      INTEGER     :: iter
      REAL(DP)    :: bephi(nkb,nbsp), becp(nkb,nbsp)
      !
      REAL(DP), ALLOCATABLE :: xloc(:,:)
      REAL(DP), ALLOCATABLE :: qbephi(:,:,:), qbecp(:,:,:)

      INTEGER :: nkbx
      INTEGER :: istart, nss, ifail, i, j, iss, iv, jv, ia, is, inl, jnl
      INTEGER :: nspin_sub, nx0, nc, ic
      REAL(DP) :: qqf

      nkbx = nkb
      nx0  = SIZE( x0, 1 )
      !
      !
      !     calculation of becp and bephi
      !
      CALL start_clock( 'ortho' )

      CALL nlsm1( nbsp, 1, nvb, eigr,  cp,  becp )
      CALL nlsm1( nbsp, 1, nvb, eigr, phi, bephi )
      !
      !     calculation of qbephi and qbecp
      !
      ALLOCATE( qbephi( nkbx, nx0, nspin ) )
      ALLOCATE( qbecp ( nkbx, nx0, nspin ) )
      !
      qbephi = 0.d0
      qbecp  = 0.d0
      !
      DO is=1,nvb
         DO iv=1,nh(is)
            DO jv=1,nh(is)
               IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                  qqf = qq(iv,jv,is)
                  DO ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     DO iss = 1, nspin
                        istart = iupdwn(iss)
                        nc     = descla( nlac_ , iss )
                        ic     = descla( ilac_ , iss )
                        IF( descla( lambda_node_ , iss ) > 0 ) THEN
                           DO i = 1, nc
                              qbephi(inl,i,iss) = qbephi(inl,i,iss)                    &
                              + qqf * bephi(jnl,i+ic-1+istart-1)
                              qbecp (inl,i,iss) = qbecp (inl,i,iss)                     &
                              + qqf * becp (jnl,i+ic-1+istart-1)
                           END DO
                        END IF
                     END DO
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      !
      ALLOCATE( xloc( nx0, nx0 ) )
      !
      nspin_sub = nspin 
      if( force_pairing ) nspin_sub = 1
      !
      DO iss = 1, nspin_sub

         nss    = nupdwn(iss)
         istart = iupdwn(iss)

         xloc = x0(:,:,iss) * ccc

         CALL ortho_gamma( 0, cp, ngwx, phi, becp, qbecp(:,:,iss), nkbx, bephi, qbephi(:,:,iss), &
                           xloc, nx0, descla(:,iss), diff, iter, nbsp, nss, istart )

         IF( iter > ortho_max ) THEN
            WRITE( stdout, * ) ' diff= ',diff,' iter= ',iter
            CALL errore('ortho','max number of iterations exceeded',iter)
         END IF

         IF( iprsta > 2 ) THEN
            WRITE( stdout,*) ' diff= ',diff,' iter= ',iter
         ENDIF
         !     
         x0( :, :, iss ) = xloc / ccc
         !
      END DO

      IF( force_pairing ) cp(:, iupdwn(2):iupdwn(2)+nupdwn(2)-1 ) = cp(:,1:nupdwn(2))
      !
      DEALLOCATE( xloc )
      DEALLOCATE(qbecp )
      DEALLOCATE(qbephi)
      !
      CALL stop_clock( 'ortho' )
      !
      RETURN
   END SUBROUTINE ortho_cp
