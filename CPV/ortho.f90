!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"



!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_m( c0, cp, lambda, ccc, nupdwn, iupdwn, nspin )
!=----------------------------------------------------------------------------=!
      !
      USE kinds,              ONLY: DP
      USE control_flags,      ONLY: force_pairing
      USE cp_main_variables,  ONLY: ema0bg
      USE control_flags,      ONLY: ortho_eps, ortho_max
      USE orthogonalize_base, ONLY: calphi, updatc
      !
      IMPLICIT NONE

      INTEGER,     INTENT(IN)    :: nupdwn(:), iupdwn(:), nspin
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:), cp(:,:)
      REAL(DP),    INTENT(INOUT) :: lambda(:,:,:)
      REAL(DP),    INTENT(IN)    :: ccc
      !
      COMPLEX(DP), ALLOCATABLE :: phi(:,:)
      INTEGER                  :: iss, nss, iwfc, nwfc, info
      INTEGER                  :: iter, i, j
      INTEGER                  :: ngwx, n, nx
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
      nss = nspin
      IF( force_pairing ) nss = 1
      !
      DO iss = 1, nss
          !
          nwfc = nupdwn(iss)
          iwfc = iupdwn(iss)
          !
          CALL ortho_gamma( 1, cp, ngwx, phi, dum, dum, 2, dum, dum, &
                            lambda(:,:,iss), nx, diff, iter, n, nwfc, iwfc )
          !
          IF ( iter > ortho_max ) THEN
             call errore(' ortho ','  itermax ',iter)
          END IF
          !
          CALL updatc( 1.0d0, n, lambda(:,:,iss), nx, phi, ngwx, dum, 1, dum, dum, cp, nwfc, iwfc )
          !     
          !     lagrange multipliers
          !
          DO i=1,nss
             DO j=1,nss
                lambda( i, j, iss ) = lambda( i, j, iss ) / ccc
             END DO
          END DO
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
   SUBROUTINE ortho_gamma( iopt, cp, ngwx, phi, becp, qbecp, nkbx, bephi, qbephi, &
                           x0, nx, diff, iter, n, nss, istart )
!=----------------------------------------------------------------------------=!
      !
      ! 
      ! 

      USE kinds,              ONLY: DP
      USE orthogonalize_base, ONLY: rhoset, sigset, tauset, ortho_iterate, &
                                    ortho_alt_iterate, updatc, &
                                    use_parallel_diag, diagonalize_parallel, &
                                    diagonalize_serial
      USE mp_global,          ONLY: nproc_image, np_ortho, me_ortho, ortho_comm, &
                                    me_image, intra_image_comm

      IMPLICIT  NONE

      ! ... Arguments

      INTEGER,  INTENT(IN)  :: iopt
      INTEGER,  INTENT(IN)  :: ngwx, nx, nkbx
      INTEGER,  INTENT(IN)  :: n, nss, istart
      COMPLEX(DP) :: phi( ngwx, n ), cp( ngwx, n )
      REAL(DP)    :: bephi( nkbx, n ), becp( nkbx, n )
      REAL(DP)    :: qbephi( nkbx, n ), qbecp( nkbx, n )
      REAL(DP)    :: x0( nx, nx )
      INTEGER,  INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff

      ! ... Locals

      REAL(DP),   ALLOCATABLE :: s(:,:), sig(:,:), tau(:,:)
      REAL(DP),   ALLOCATABLE :: wrk(:,:), rhoa(:,:), rhos(:,:), rhod(:)
      INTEGER  :: i, j, info, nr, nc, ir, ic
      LOGICAL  :: iter_node
      !
      INTEGER  :: ldim_block, gind_block
      EXTERNAL :: ldim_block, gind_block

      ! ...   Subroutine body

      IF( me_ortho(1) >= 0 ) THEN
         !
         iter_node = .TRUE.
         !
         nr = ldim_block( nss, np_ortho(1), me_ortho(1) )
         nc = ldim_block( nss, np_ortho(2), me_ortho(2) )
         !
         ir = gind_block( 1, nss, np_ortho(1), me_ortho(1) )
         ic = gind_block( 1, nss, np_ortho(2), me_ortho(2) )
         !
      ELSE
         !
         iter_node = .FALSE.
         !
      END IF


      ALLOCATE( wrk( nx, nx ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho ', ' allocating matrixes ', 1 )
      ALLOCATE( rhos( nx, nx ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho ', ' allocating matrixes ', 2 )
      ALLOCATE( rhod( nx ), STAT = info )
      IF( info /= 0 ) CALL errore( ' ortho ', ' allocating matrixes ', 3 )
      !
      !     rho = <s'c0|s|cp>
      !
      CALL rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, istart, wrk, nx )
      !
      !  Symmetric part of rho
      !
      DO j = 1, nss
        DO i = 1, nss

          rhos(i,j) = 0.5d0*(wrk(i,j)+wrk(j,i))
          !
          ! on some machines (IBM RS/6000 for instance) the following test allows
          ! to distinguish between Numbers and Sodium Nitride (NaN, Not a Number).
          ! If a matrix of Not-Numbers is passed to rs, the most likely outcome is
          ! that the program goes on forever doing nothing and writing nothing.
          !
          IF (rhos(i,j) /= rhos(i,j)) CALL errore('ortho','ortho went bananas',1)

        ENDDO
      ENDDO
      !
      !  Antisymmetric part of rho, alredy distributed across ortho procs.
      !
      IF( iter_node ) THEN
         ALLOCATE( rhoa( nr, nc ) )
         DO j = 1, nc
            DO i = 1, nr
               rhoa( i, j ) = 0.5d0 * ( wrk( i+ir-1, j+ic-1 ) - wrk( j+ic-1, i+ir-1 ) )
            END DO
         END DO
      ELSE
         ALLOCATE( rhoa( 1, 1 ) )
      END IF


      ! ...   Diagonalize Matrix  symmetric part of rho (rhos)


      CALL start_clock( 'rsg' )

      !
      !       Compute s and rho
      ! ...   "s" is the matrix of eigenvectors, "rhod" is the array of eigenvalues
      !

      IF( iter_node ) THEN
         !
         ALLOCATE( s( nr, nc ) ) 
         ALLOCATE( sig( nr, nc ) ) 
         ALLOCATE( tau( nr, nc ) ) 
         !
      ELSE
         !
         ALLOCATE( s( 1, 1 ) ) 
         ALLOCATE( sig( 1, 1 ) ) 
         ALLOCATE( tau( 1, 1 ) ) 
         !
      END IF

      IF( use_parallel_diag ) THEN
         !
         CALL diagonalize_parallel( nss, rhos, rhod, s )
         !
      ELSE
         !
         CALL diagonalize_serial( nss, rhos, rhod, wrk )
         CALL distribute_matrix( wrk, s )
         !
      END IF
      !

      CALL stop_clock( 'rsg' )

      !
      IF( info /= 0 ) CALL errore( ' ortho ', ' allocating matrixes ', 4 )
      !
      !     sig = 1-<cp|s|cp>
      !
      CALL sigset( cp, ngwx, becp, nkbx, qbecp, n, nss, istart, wrk, nx )
      CALL distribute_matrix( wrk, sig )
      !
      !     tau = <s'c0|s|s'c0>
      !
      CALL tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, istart, wrk, nx )
      CALL distribute_matrix( wrk, tau )

      DEALLOCATE( wrk )

      IF( iopt == 0 ) THEN
         !
         CALL ortho_iterate( iter, diff, s, nr, rhod, x0, sig, rhoa, rhos, tau, nx, nss )
         !
      ELSE
         !
         CALL ortho_alt_iterate( iter, diff, s, nr, rhod, x0, sig, rhoa, tau, nx, nss )
         !
      END IF
      !
      DO i=1,nss
        DO j=1,nss
          IF (x0(i,j) /= x0(i,j)) CALL errore('ortho','ortho went bananas',2)
        END DO
      END DO

      DEALLOCATE( rhoa, rhos, rhod, s, sig, tau )

      RETURN

   CONTAINS

      SUBROUTINE distribute_matrix( a, b )
         REAL(DP) :: a(:,:), b(:,:)
         IF( iter_node ) THEN
            DO j = 1, nc
               DO i = 1, nr
                  b( i, j ) = a( i + ir - 1, j + ic - 1 )
               END DO
            END DO
         END IF
         RETURN
      END SUBROUTINE


   END SUBROUTINE ortho_gamma




!=----------------------------------------------------------------------------=!
   SUBROUTINE ortho_cp( eigr, cp, phi, ngwx, x0, nudx, diff, iter, ccc, &
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
      !
      IMPLICIT NONE
!
      INTEGER     :: ngwx, nudx, nbsp, nspin
      INTEGER     :: nupdwn( nspin ), iupdwn( nspin )
      COMPLEX(DP) :: cp(ngwx,nbsp), phi(ngwx,nbsp), eigr(ngwx,nat)
      REAL(DP)    :: x0( nudx, nudx, nspin ), diff, ccc
      INTEGER     :: iter
      REAL(DP)    :: bephi(nkb,nbsp), becp(nkb,nbsp)
!
      REAL(DP), ALLOCATABLE :: xloc(:,:)
      REAL(DP), ALLOCATABLE :: qbephi(:,:), qbecp(:,:)

      INTEGER :: nkbx
      INTEGER :: istart, nss, ifail, i, j, iss, iv, jv, ia, is, inl, jnl
      INTEGER :: nspin_sub

      nkbx = nkb
      !
      !     calculation of becp and bephi
      !
      CALL start_clock( 'ortho' )

      CALL nlsm1( nbsp, 1, nvb, eigr,  cp,  becp )
      CALL nlsm1( nbsp, 1, nvb, eigr, phi, bephi )
      !
      !     calculation of qbephi and qbecp
      !
      ALLOCATE( qbephi( nkbx, nbsp ) )
      ALLOCATE( qbecp ( nkbx, nbsp ) )
      !
      qbephi = 0.d0
      qbecp  = 0.d0
      !
      DO is=1,nvb
         DO iv=1,nh(is)
            DO jv=1,nh(is)
               IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN
                  DO ia=1,na(is)
                     inl=ish(is)+(iv-1)*na(is)+ia
                     jnl=ish(is)+(jv-1)*na(is)+ia
                     DO i=1,nbsp
                        qbephi(inl,i)= qbephi(inl,i)                    &
     &                       +qq(iv,jv,is)*bephi(jnl,i)
                        qbecp (inl,i)=qbecp (inl,i)                     &
     &                       +qq(iv,jv,is)*becp (jnl,i)
                     END DO
                  END DO
               ENDIF
            END DO
         END DO
      END DO
      !
      ALLOCATE( xloc( nudx, nudx ) )
      !
      !
      nspin_sub = nspin 
      if( force_pairing ) nspin_sub = 1
      !
      DO iss = 1, nspin_sub

         nss    = nupdwn(iss)
         istart = iupdwn(iss)

         DO j=1,nss
            DO i=1,nss
               xloc(i,j) = x0( i, j, iss ) * ccc
            END DO
         END DO

         CALL ortho_gamma( 0, cp, ngwx, phi, becp, qbecp, nkbx, bephi, qbephi, &
                           xloc, nudx, diff, iter, nbsp, nss, istart )

         IF( iter > ortho_max ) THEN
            WRITE( stdout, * ) ' diff= ',diff,' iter= ',iter
            CALL errore('ortho','max number of iterations exceeded',iter)
         END IF

         IF( iprsta > 4 ) THEN
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    lambda '
            DO i=1,nss
               WRITE( stdout,'(7f11.6)') (xloc(i,j)/f(i+istart-1),j=1,nss)
            END DO
         ENDIF
         IF( iprsta > 2 ) THEN
            WRITE( stdout,*) ' diff= ',diff,' iter= ',iter
         ENDIF
         !     
         !     lagrange multipliers
         !
         DO i=1,nss
            DO j=1,nss
               x0( i, j, iss ) = xloc(i,j) / ccc
            END DO
         END DO
!
      END DO

      IF( force_pairing .AND. nspin > 1 ) THEN
         !
         x0(1:nupdwn(2), 1:nupdwn(2), 2) = x0(1:nupdwn(2), 1:nupdwn(2), 1)
         x0(nudx, nudx, 2) = 0.d0
         !
      ENDIF
!
      DEALLOCATE( xloc )
      DEALLOCATE(qbecp )
      DEALLOCATE(qbephi)
!
      CALL stop_clock( 'ortho' )
      RETURN
      END SUBROUTINE ortho_cp

