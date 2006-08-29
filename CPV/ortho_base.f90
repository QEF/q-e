!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"


MODULE orthogonalize_base


      USE kinds
      USE parallel_toolkit, ONLY: pdspev_drv, dspev_drv, &
                                  pzhpev_drv, zhpev_drv, &
                                  rep_matmul_drv

      IMPLICIT NONE

      SAVE

      PRIVATE

      REAL(DP) :: one, zero, two, minus_one, minus_two
      PARAMETER ( one = 1.0d0, zero = 0.0d0, two = 2.0d0, minus_one = -1.0d0 )
      PARAMETER ( minus_two = -2.0d0 )
      COMPLEX(DP) :: cone, czero, mcone
      PARAMETER ( cone = (1.0d0, 0.0d0), czero = (0.0d0, 0.0d0) )
      PARAMETER ( mcone = (-1.0d0, 0.0d0) )
      REAL(DP) :: small = 1.0d-14

      INTERFACE sqr_matmul
         MODULE PROCEDURE sqr_dmatmul
      END INTERFACE

      INTERFACE diagonalize_rho
         MODULE PROCEDURE diagonalize_rrho, diagonalize_crho
      END INTERFACE

      PUBLIC :: sigset, rhoset, tauset, diagonalize_rho
      PUBLIC :: ortho_iterate
      PUBLIC :: ortho_alt_iterate
      PUBLIC :: updatc, calphi

CONTAINS


   SUBROUTINE sqr_dmatmul( transa, transb, n, a, b, c )

      ! ...    Multiply square matrices A, B and return the result in C

      USE control_flags, ONLY: iprsta
      USE mp_global,     ONLY: nproc_image, me_image, root_image, intra_image_comm
      USE io_global,     ONLY: ionode, stdout
      USE mp,            ONLY: mp_bcast

      REAL(DP) :: c(:,:), a(:,:), b(:,:)
      CHARACTER(LEN=1), INTENT(IN) :: transa, transb
      INTEGER, INTENT(IN) :: n

      LOGICAL :: lpdrv
      INTEGER, SAVE :: calls_cnt = 0
      REAL(DP) :: t1
      REAL(DP), SAVE :: tser, tpar
      REAL(DP), EXTERNAL :: cclock
      !
      IF( n < 1 ) RETURN
      !
      calls_cnt = calls_cnt + 1

      IF( nproc_image == 1 ) THEN
         lpdrv = .FALSE.          !  with one proc do not use parallel diag
      ELSE IF ( calls_cnt == 1 ) THEN
         lpdrv = .TRUE.           !  use parallel diag the first call to take the time
      ELSE IF ( calls_cnt == 2 ) THEN
         lpdrv = .FALSE.          !  use seria diag the second call to take the time
      ELSE IF ( tpar < tser ) THEN
         lpdrv = .TRUE.           !  use para diag if it is faster
         IF( calls_cnt == 3 .AND. ionode .AND. iprsta > 1 ) WRITE( stdout, 10 ) tpar, tser
      ELSE
         lpdrv = .FALSE.          !  use scalar otherwise
         IF( calls_cnt == 3 .AND. ionode .AND. iprsta > 1 ) WRITE( stdout, 20 ) tpar, tser
      END IF

10    FORMAT(3X,'ortho matmul, time for parallel and serial driver = ', 2D9.2, /, &
             3X,'using parallel driver' )
20    FORMAT(3X,'ortho matmul, time for parallel and serial driver = ', 2D9.2, /, &
             3X,'using serial driver' )

      IF ( lpdrv ) THEN

         IF( calls_cnt < 3 )  t1 = cclock()

         CALL rep_matmul_drv( transa, transb, n, n, n, one, A, SIZE(a,1), B, SIZE(b,1), zero, C, &
                              SIZE(c,1), intra_image_comm )

         IF( calls_cnt < 3 )  THEN
            tpar = cclock() - t1
            CALL mp_bcast( tpar, root_image, intra_image_comm )
         END IF

      ELSE

         IF( calls_cnt < 3 )  t1 = cclock()

         CALL DGEMM( transa, transb, n, n, n, one, a, SIZE(a,1), b, SIZE(b,1), zero, c, SIZE(c,1) )

         IF( calls_cnt < 3 )  THEN
            tser = cclock() - t1
            CALL mp_bcast( tser, root_image, intra_image_comm )
         END IF

      END IF
      !
      RETURN
   END SUBROUTINE sqr_dmatmul


!  ----------------------------------------------


   SUBROUTINE diagonalize_rrho( n, rhos, rhod, s )

         !   Diagonalization of rhos

      USE control_flags, ONLY: iprsta
      USE mp_global, ONLY: nproc_image, me_image, intra_image_comm, root_image
      USE io_global, ONLY: ionode, stdout
      USE mp,        ONLY: mp_sum, mp_bcast
      !
      REAL(DP), INTENT(IN) :: rhos(:,:) !  input symmetric matrix
      REAL(DP)             :: rhod(:)   !  output eigenvalues
      REAL(DP)             :: s(:,:)    !  output eigenvectors
      INTEGER, INTENT(IN)  :: n         !  matrix dimension

      REAL(DP),   ALLOCATABLE :: aux(:)
      REAL(DP),   ALLOCATABLE :: diag(:,:)
      REAL(DP),   ALLOCATABLE :: vv(:,:)
      !
      INTEGER :: nrl
      LOGICAL :: lpdrv
      INTEGER, SAVE :: calls_cnt = 0
      REAL(DP) :: t1
      REAL(DP), SAVE :: tser, tpar
      REAL(DP), EXTERNAL :: cclock

      IF( n < 1 ) RETURN

      calls_cnt = calls_cnt + 1

      IF( nproc_image == 1 ) THEN 
         lpdrv = .FALSE.          !  with one proc do not use parallel diag
      ELSE IF ( calls_cnt == 1 ) THEN 
         lpdrv = .TRUE.           !  use parallel diag the first call to take the time
      ELSE IF ( calls_cnt == 2 ) THEN 
         lpdrv = .FALSE.          !  use seria diag the second call to take the time
      ELSE IF ( tpar < tser ) THEN
         lpdrv = .TRUE.           !  use para diag if it is faster
         IF( calls_cnt == 3 .AND. ionode .AND. iprsta > 1 ) WRITE( stdout, 10 ) tpar, tser
      ELSE
         lpdrv = .FALSE.          !  use scalar otherwise
         IF( calls_cnt == 3 .AND. ionode .AND. iprsta > 1 ) WRITE( stdout, 20 ) tpar, tser
      END IF

10    FORMAT(3X,'ortho diag, time for parallel and serial driver = ', 2D9.2, /, &
             3X,'using parallel driver' )
20    FORMAT(3X,'ortho diag, time for parallel and serial driver = ', 2D9.2, /, &
             3X,'using serial driver' )
         

      IF( SIZE( rhos, 1 ) /= SIZE( s, 1 ) .OR. SIZE( rhos, 2 ) /= SIZE( s, 2 ) ) &
         CALL errore(" diagonalize_rho ", " input matrixes size do not match ", 1 )


      IF ( lpdrv ) THEN

        IF( calls_cnt < 3 )  t1 = cclock()

        !  distribute matrix rows to processors
        !

        nrl = n / nproc_image
        IF( me_image < MOD( n, nproc_image ) ) THEN
          nrl = nrl + 1
        end if

        ALLOCATE( diag( nrl, n ), vv( nrl, n ) )

        CALL prpack( n, diag, rhos)
        CALL pdspev_drv( 'V', diag, nrl, rhod, vv, nrl, nrl, n, nproc_image, me_image)
        CALL prunpack( n, s, vv)

        DEALLOCATE( diag, vv )

        CALL mp_sum( s, intra_image_comm )

        IF( calls_cnt < 3 )  THEN

           tpar = cclock() - t1

           CALL mp_bcast( tpar, root_image, intra_image_comm )

        END IF

      ELSE

        IF( calls_cnt < 3 )  t1 = cclock()

        ALLOCATE( aux( n * ( n + 1 ) / 2 ) )

        CALL rpack( n, aux, rhos )   !  pack lower triangle of rho into aux

        CALL dspev_drv( 'V', 'L', n, aux, rhod, s, SIZE(s,1) )

        DEALLOCATE( aux )

        IF( calls_cnt < 3 )  THEN

           tser = cclock() - t1

           CALL mp_bcast( tser, root_image, intra_image_comm )

        END IF

      END IF

      RETURN
   END SUBROUTINE diagonalize_rrho


!  ----------------------------------------------


   SUBROUTINE diagonalize_crho( n, a, d, ev, use_pdrv )

      !  this routine calls the appropriate Lapack routine for diagonalizing a
      !  complex Hermitian matrix

         USE mp_global, ONLY: nproc_image, me_image, intra_image_comm
         USE mp, ONLY: mp_sum
         IMPLICIT NONE

         REAL(DP)             :: d(:)
         COMPLEX(DP)          :: a(:,:), ev(:,:)
         LOGICAL, OPTIONAL, &
                  INTENT(IN)  :: use_pdrv  ! if true use parallel driver
         INTEGER, INTENT(IN)  :: n

         INTEGER :: nrl

         COMPLEX(DP), ALLOCATABLE :: aloc(:)
         COMPLEX(DP), ALLOCATABLE :: ap(:,:)
         COMPLEX(DP), ALLOCATABLE :: vp(:,:)

! ...   end of declarations
!  ----------------------------------------------
         LOGICAL :: lpdrv

         lpdrv = .FALSE.
         
         IF( PRESENT( use_pdrv ) ) lpdrv = use_pdrv

         IF ( ( nproc_image > 1 ) .AND. use_pdrv ) THEN

           nrl = n/nproc_image
           IF(me_image.LT.MOD(n,nproc_image)) THEN
             nrl = nrl + 1
           END IF

           ALLOCATE(ap(nrl,n), vp(nrl,n))

           CALL pzpack(ap, a)
           CALL pzhpev_drv( 'V', ap, nrl, d, vp, nrl, nrl, n, nproc_image, me_image)
           CALL pzunpack(a, vp)
           CALL mp_sum(a, ev, intra_image_comm)

           DEALLOCATE(ap, vp)

         ELSE

           ALLOCATE(aloc(n*(n+1)/2))

           ! ...      copy the lower-diagonal part of the matrix according to the
           ! ...      Lapack packed storage scheme for Hermitian matrices

           CALL zpack(aloc, a)

           ! ...      call the Lapack routine

           CALL zhpev_drv('V', 'L', n, aloc, d, ev, n)

           DEALLOCATE(aloc)

         END IF

         RETURN
       END SUBROUTINE diagonalize_crho


!=----------------------------------------------------------------------------=!


      SUBROUTINE prpack( n, ap, a)
        USE mp_global, ONLY: me_image, nproc_image
        REAL(DP), INTENT(IN) :: a(:,:)
        REAL(DP), INTENT(OUT) :: ap(:,:)
        INTEGER, INTENT(IN) :: n
        INTEGER :: i, j, jl
        DO i = 1, SIZE( ap, 2)
           j = me_image + 1
           DO jl = 1, SIZE( ap, 1)
             ap(jl,i) = a(j,i)
             j = j + nproc_image
           END DO
        END DO
        RETURN
      END SUBROUTINE prpack

      SUBROUTINE pzpack( ap, a)
        USE mp_global, ONLY: me_image, nproc_image
        COMPLEX(DP), INTENT(IN) :: a(:,:)
        COMPLEX(DP), INTENT(OUT) :: ap(:,:)
        INTEGER :: i, j, jl
        DO i = 1, SIZE( ap, 2)
           j = me_image + 1
           DO jl = 1, SIZE( ap, 1)
             ap(jl,i) = a(j,i)
             j = j + nproc_image
           END DO
        END DO
        RETURN
      END SUBROUTINE pzpack

      SUBROUTINE prunpack( n, a, ap )
        USE mp_global, ONLY: me_image, nproc_image
        REAL(DP), INTENT(IN) :: ap(:,:)
        REAL(DP), INTENT(OUT) :: a(:,:)
        INTEGER, INTENT(IN) :: n
        INTEGER :: i, j, jl
        DO i = 1, n
          DO j = 1, n
            a(j,i) = zero
          END DO
          j = me_image + 1
          DO jl = 1, SIZE(ap, 1)
            a(j,i) = ap(jl,i)
            j = j + nproc_image
          END DO
        END DO
        RETURN
      END SUBROUTINE prunpack

      SUBROUTINE pzunpack( a, ap)
        USE mp_global, ONLY: me_image, nproc_image
        COMPLEX(DP), INTENT(IN) :: ap(:,:)
        COMPLEX(DP), INTENT(OUT) :: a(:,:)
        INTEGER :: i, j, jl
        DO i = 1, SIZE(a, 2)
          DO j = 1, SIZE(a, 1)
            a(j,i) = zero
          END DO
          j = me_image + 1
          DO jl = 1, SIZE(ap, 1)
            a(j,i) = ap(jl,i)
            j = j + nproc_image
          END DO
        END DO
        RETURN
      END SUBROUTINE pzunpack

      SUBROUTINE rpack( n, ap, a)
        INTEGER, INTENT(IN) :: n
        REAL(DP), INTENT(IN) :: a(:,:)
        REAL(DP), INTENT(OUT) :: ap(:)
        INTEGER :: i, j, k
        K = 0
        DO J = 1, n
          DO I = J, n
            K = K + 1
            ap( k ) = a( i, j )
          END DO
        END DO
        RETURN
      END SUBROUTINE rpack

      SUBROUTINE zpack( ap, a)
        COMPLEX(DP), INTENT(IN) :: a(:,:)
        COMPLEX(DP), INTENT(OUT) :: ap(:)
        INTEGER :: i, j, k
        K=0
        DO J = 1, SIZE(a, 2)
          DO I = J, SIZE(a, 1)
            K = K + 1
            ap(k) = a(i,j)
          END DO
        END DO
        RETURN
      END SUBROUTINE zpack


!=----------------------------------------------------------------------------=!


   SUBROUTINE ortho_iterate( iter, diff, u, diag, xloc, sig, rhor, rhos, tau, nx, nss )

      USE kinds,         ONLY: DP
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: ortho_eps, ortho_max

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx, nss
      REAL(DP) :: u( nx, nx )
      REAL(DP) :: diag( nx )
      REAL(DP) :: xloc( nx, nx )
      REAL(DP) :: rhor( nx, nx )
      REAL(DP) :: rhos( nx, nx )
      REAL(DP) :: tau( nx, nx )
      REAL(DP) :: sig( nx, nx )
      INTEGER, INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff 

      INTEGER :: i, j
      REAL(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:), dd(:,:)
      REAL(DP), ALLOCATABLE :: con(:,:), x1(:,:)

      IF( nss < 1 ) RETURN

      ALLOCATE( tmp1(nx,nx), tmp2(nx,nx), dd(nx,nx), x1(nx,nx), con(nx,nx) )


      ITERATIVE_LOOP: DO iter = 1, ortho_max
         !
         !       the following 4 MXMA-calls do the following matrix
         !       multiplications:
         !                       tmp1 = x0*rhor    (1st call)
         !                       dd   = x0*tau*x0  (2nd and 3rd call)
         !                       tmp2 = x0*rhos    (4th call)
         !

         CALL sqr_matmul( 'N', 'N', nss, xloc, rhor, tmp1 )
         CALL sqr_matmul( 'N', 'N', nss, tau,  xloc, tmp2 )
         CALL sqr_matmul( 'N', 'N', nss, xloc, tmp2, dd )
         CALL sqr_matmul( 'N', 'N', nss, xloc, rhos, tmp2 )

         ! CALL MXMA( xloc,1,nx,rhor,1,nx,tmp1,1,nx,nss,nss,nss)
         ! CALL MXMA( tau ,1,nx,xloc,1,nx,tmp2,1,nx,nss,nss,nss)
         ! CALL MXMA( xloc,1,nx,tmp2,1,nx,  dd,1,nx,nss,nss,nss)
         ! CALL MXMA( xloc,1,nx,rhos,1,nx,tmp2,1,nx,nss,nss,nss)
         !
         DO i=1,nss
            DO j=1,nss
               x1(i,j) = sig(i,j)-tmp1(i,j)-tmp1(j,i)-dd(i,j)
               con(i,j)= x1(i,j)-tmp2(i,j)-tmp2(j,i)
            END DO
         END DO
         !
         !         x1      = sig      -x0*rho    -x0*rho^t  -x0*tau*x0
         !
         diff = 0.d0
         DO i=1,nss
            DO j=1,nss
               IF(ABS(con(i,j)).GT.diff) diff=ABS(con(i,j))
            END DO
         END DO

         IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

         !
         !     the following two MXMA-calls do:
         !                       tmp1 = x1*u
         !                       tmp2 = ut*x1*u
         !
         CALL sqr_matmul( 'N', 'N', nss, x1,    u, tmp1 )
         CALL sqr_matmul( 'T', 'N', nss,  u, tmp1, tmp2 )
         !
         ! CALL MXMA(x1,1,nx,   u,1,nx,tmp1,1,nx,nss,nss,nss)
         ! CALL MXMA(u ,nx,1,tmp1,1,nx,tmp2,1,nx,nss,nss,nss)
         !
         !       g=ut*x1*u/d  (g is stored in tmp1)
         !
         DO i=1,nss
            DO j=1,nss
               tmp1(i,j)=tmp2(i,j)/(diag(i)+diag(j))
            END DO
         END DO
         !
         !       the following two MXMA-calls do:
         !                       tmp2 = g*ut
         !                       x0 = u*g*ut
         !
         CALL sqr_matmul( 'N', 'T', nss, tmp1,    u, tmp2 )
         CALL sqr_matmul( 'N', 'N', nss,    u, tmp2, xloc )
         !
         ! CALL MXMA(tmp1,1,nx,  u,nx,1,tmp2,1,nx,nss,nss,nss)
         ! CALL MXMA(   u,1,nx,tmp2,1,nx,xloc,1,nx,nss,nss,nss)
         !
      END DO ITERATIVE_LOOP

      DEALLOCATE( tmp1, tmp2, dd, x1, con )

      RETURN
   END SUBROUTINE ortho_iterate


!=----------------------------------------------------------------------------=!
!
!  Alternative iterative cycle
!
!=----------------------------------------------------------------------------=!


   SUBROUTINE ortho_alt_iterate( iter, diff, u, diag, xloc, sig, rhor, tau, nx, n )

      USE kinds,         ONLY: DP
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: ortho_eps, ortho_max
      USE mp_global,     ONLY: intra_image_comm

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx, n
      REAL(DP) :: u( nx, nx )
      REAL(DP) :: diag( nx )
      REAL(DP) :: xloc( nx, nx )
      REAL(DP) :: rhor( nx, nx )
      REAL(DP) :: tau( nx, nx )
      REAL(DP) :: sig( nx, nx )
      INTEGER, INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff 

      INTEGER :: i, j
      REAL(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:)
      REAL(DP), ALLOCATABLE :: x1(:,:)
      REAL(DP), ALLOCATABLE :: sigd(:)
      REAL(DP) :: den, dx

      IF( n < 1 ) RETURN

      ALLOCATE( tmp1(nx,nx), tmp2(nx,nx), x1(nx,nx), sigd(nx) )


      !
      ! ...   Transform "sig", "rhoa" and "tau" in the new basis through matrix "s"
      !
      CALL sqr_matmul( 'N', 'N', n, sig, u, tmp1 )
      CALL sqr_matmul( 'T', 'N', n, u, tmp1, sig )
      CALL sqr_matmul( 'N', 'N', n, rhor, u, tmp1 )
      CALL sqr_matmul( 'T', 'N', n, u, tmp1, rhor )
      CALL sqr_matmul( 'N', 'N', n, tau, u, tmp1 )
      CALL sqr_matmul( 'T', 'N', n, u, tmp1, tau )
      !
      ! ...   Initialize x0
      !
      DO J = 1, N
        DO I = 1, N
          den = (diag(i)+diag(j))
          IF( ABS( den ) <= small ) den = SIGN( small, den )
          xloc(i,j) = sig(i,j) / den
        ENDDO
      ENDDO

      !
      ! ...   Starting iteration
      !

      ITERATIVE_LOOP: DO iter = 0, ortho_max

        CALL sqr_matmul( 'N', 'N', n, xloc, rhor, tmp2 )
        call mytranspose( tmp2, NX, tmp1, NX, N, N )
        DO J=1,N
          DO I=1,N
            tmp2(I,J) = tmp2(I,J) + tmp1(I,J)
          ENDDO
        ENDDO
!
        CALL sqr_matmul( 'T', 'N', n, tau, xloc, tmp1 )
        !
        DO I = 1, N
          SIGD(I)   =  tmp1(I,I)
          tmp1(I,I) = -SIGD(I)
        ENDDO

        CALL sqr_matmul( 'T', 'N', n, xloc, tmp1, X1 )
        !
        call mytranspose( X1, NX, tmp1, NX, N, N )

        ! ...     X1   = SIG - tmp2 - 0.5d0 * ( X1 + X1^t )

        diff = 0.0d0
        !
        DO j = 1, n
          DO i = 1, n
            !
            den = ( diag(i) + sigd(i) + diag(j) + sigd(j) )
            IF( ABS( den ) <= small ) den = SIGN( small, den )
            x1(i,j) = sig(i,j) - tmp2(i,j) - 0.5d0 * (x1(i,j)+tmp1(i,j))
            x1(i,j) = x1(i,j) / den
            diff = MAX( ABS( x1(i,j) - xloc(i,j) ), diff )
            xloc(i,j) = x1(i,j)
          END DO
        END DO

        IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

      END DO ITERATIVE_LOOP
      !
      ! ...   Transform x0 back to the original basis

      CALL sqr_matmul( 'N', 'N', n, u, xloc, tmp1 )
      CALL sqr_matmul( 'N', 'T', n, u, tmp1, xloc )

      DEALLOCATE( tmp1, tmp2, x1 )

      RETURN
   END SUBROUTINE ortho_alt_iterate



!-------------------------------------------------------------------------
   SUBROUTINE sigset( cp, ngwx, becp, nkbx, qbecp, n, nss, ist, sig, nx )
!-----------------------------------------------------------------------
!     input: cp (non-orthonormal), becp, qbecp
!     computes the matrix
!       sig = 1 - a ,  a = <cp|s|cp> = <cp|cp> + sum q_ij <cp|i><j|cp>
!     where s=s(r(t+dt))
!     routine makes use of c(-q)=c*(q)
!
      USE kinds,              ONLY: DP
      USE uspp,               ONLY: nkbus
      USE cvan,               ONLY: nvb
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE mp,                 ONLY: mp_sum
      USE control_flags,      ONLY: iprsta
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_image_comm
!
      IMPLICIT NONE
!
      INTEGER nss, ist, ngwx, nkbx, n, nx
      COMPLEX(DP) :: cp( ngwx, n )
      REAL(DP)    :: becp( nkbx, n ), qbecp( nkbx, n ), sig( nx, nx )
!
      INTEGER :: i, j
      REAL(DP), ALLOCATABLE :: tmp1(:,:)
!
      IF( nss < 1 ) RETURN

      CALL DGEMM( 'T', 'N',  nss, nss, 2*ngw, -2.0d0, cp( 1, ist ), 2*ngwx, &
                  cp( 1, ist ), 2*ngwx, 0.0d0, sig, nx)
      !
      !     q = 0  components has weight 1.0
      !
      IF ( gstart == 2 ) THEN
         DO j=1,nss
            DO i=1,nss
               sig(i,j) = sig(i,j) +                                    &
     &              DBLE(cp(1,i+ist-1))*DBLE(cp(1,j+ist-1))
            END DO
         END DO
      END IF
      !
      CALL mp_sum( sig, intra_image_comm )
      !
      DO i = 1, nss
         sig(i,i) = sig(i,i) + 1.0d0
      END DO
!
      IF( nvb > 0 ) THEN

         CALL DGEMM( 'T', 'N', nss, nss, nkbus, -1.0d0, becp( 1, ist ), nkbx, &
                  qbecp( 1, ist ), nkbx, 1.0d0, sig, nx )
!
!         ALLOCATE( tmp1( nx, nx ) )
!
!         CALL MXMA( becp( 1, ist ), nkbx, 1, qbecp( 1, ist ), 1, nkbx,                &
!     &              tmp1, 1, nx, nss, nkbus, nss )
!
!         DO j=1,nss
!            DO i=1,nss
!               sig(i,j)=sig(i,j)-tmp1(i,j)
!            END DO
!         END DO
!
!         DEALLOCATE( tmp1 )

      ENDIF

      IF(iprsta.GT.4) THEN
         WRITE( stdout,*)
         WRITE( stdout,'(26x,a)') '    sig '
         DO i=1,nss
            WRITE( stdout,'(7f11.6)') (sig(i,j),j=1,nss)
         END DO
      ENDIF

!
      RETURN
   END SUBROUTINE sigset

!
!-----------------------------------------------------------------------
   SUBROUTINE rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, ist, rho, nx )
!-----------------------------------------------------------------------
!     input: cp (non-orthonormal), phi, bephi, qbecp
!     computes the matrix
!       rho = <s'c0|s cp> = <phi|s cp>
!     where  |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
!     where s=s(r(t+dt)) and s'=s(r(t))
!     routine makes use of  c(-q)=c*(q)
!
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE uspp,               ONLY: nkbus
      USE cvan,               ONLY: nvb
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_sum
      USE mp_global,          ONLY: intra_image_comm
      USE control_flags,      ONLY: iprsta
      USE io_global,          ONLY: stdout
!
      IMPLICIT NONE
!
      INTEGER     :: nss, ist, ngwx, nkbx, nx, n
      COMPLEX(DP) :: cp( ngwx, n ), phi( ngwx, n )
      REAL(DP)    :: bephi( nkbx, n ), qbecp( nkbx, n ), rho( nx, nx )
      !
      INTEGER     :: i, j
      REAL(DP), ALLOCATABLE :: tmp1(:,:)
      !
      !     <phi|cp>
      !
      !
      IF( nss < 1 ) RETURN

      CALL DGEMM( 'T', 'N', nss, nss, 2*ngw, 2.0d0, phi( 1, ist ), 2*ngwx, &
                  cp( 1, ist ), 2*ngwx, 0.0d0, rho, nx)
      !
      !     q = 0  components has weight 1.0
      !
      IF (gstart == 2) THEN
         DO j=1,nss
            DO i=1,nss
               rho(i,j) = rho(i,j) -                                    &
     &              DBLE(phi(1,i+ist-1))*DBLE(cp(1,j+ist-1))
            END DO
         END DO
      END IF

      CALL mp_sum( rho, intra_image_comm )
!
      IF( nvb > 0 ) THEN
         !
         ! rho(i,j) = rho(i,j) + SUM_b bephi( b, i ) * qbecp( b, j ) 
         !
         CALL DGEMM( 'T', 'N', nss, nss, nkbus, 1.0d0, bephi( 1, ist ), nkbx, &
                  qbecp( 1, ist ), nkbx, 1.0d0, rho, nx )

!         ALLOCATE( tmp1( nx, nx ) )
!
!         CALL MXMA( bephi( 1, ist ), nkbx, 1, qbecp( 1, ist ), 1, nkbx,               &
!     &                                tmp1, 1, nx, nss, nkbus, nss )
!
!         DO j=1,nss
!            DO i=1,nss
!               rho(i,j)=rho(i,j)+tmp1(i,j)
!            END DO
!         END DO
!
!         DEALLOCATE( tmp1 )

      ENDIF

      IF(iprsta.GT.4) THEN
         WRITE( stdout,*)
         WRITE( stdout,'(26x,a)') '    rho '
         DO i=1,nss
            WRITE( stdout,'(7f11.6)') (rho(i,j),j=1,nss)
         END DO
      ENDIF
!
      RETURN
   END SUBROUTINE rhoset

!-------------------------------------------------------------------------
   SUBROUTINE tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, ist, tau, nx )
!-----------------------------------------------------------------------
!     input: phi
!     computes the matrix
!        tau = <s'c0|s|s'c0> = <phi|s|phi>,  where  |phi> = s'|c0>
!     where s=s(r(t+dt)) and s'=s(r(t))
!     routine makes use of c(-q)=c*(q)
!
      USE kinds,              ONLY: DP
      USE cvan,               ONLY: nvb
      USE uspp,               ONLY: nkbus
      USE gvecw,              ONLY: ngw
      USE reciprocal_vectors, ONLY: gstart
      USE mp,                 ONLY: mp_sum
      USE control_flags,      ONLY: iprsta
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_image_comm
!
      IMPLICIT NONE
      INTEGER :: nss, ist, ngwx, nkbx, n, nx
      COMPLEX(DP) :: phi( ngwx, n )
      REAL(DP)    :: bephi( nkbx, n ), qbephi( nkbx, n ), tau( nx, nx )
      !
      INTEGER     :: i, j
      REAL(DP), ALLOCATABLE :: tmp1( :, : )
      !
      IF( nss < 1 ) RETURN
      !
      CALL DGEMM( 'T', 'N', nss, nss, 2*ngw, 2.0d0, phi( 1, ist ), 2*ngwx, &
                  phi( 1, ist ), 2*ngwx, 0.0d0, tau, nx)
      !
      !     q = 0  components has weight 1.0
      !
      IF (gstart == 2) THEN
         DO j=1,nss
            DO i=1,nss
               tau(i,j) = tau(i,j) -                                    &
     &              DBLE(phi(1,i+ist-1))*DBLE(phi(1,j+ist-1))
            END DO
         END DO
      END IF

      CALL mp_sum( tau, intra_image_comm )
!
      IF( nvb > 0 ) THEN
         !
         CALL DGEMM( 'T', 'N', nss, nss, nkbus, 1.0d0, bephi( 1, ist ), nkbx, &
                  qbephi( 1, ist ), nkbx, 1.0d0, tau, nx )

!         ALLOCATE( tmp1( nx, nx ) )
!
!         CALL MXMA( bephi( 1, ist ), nkbx, 1, qbephi( 1, ist ), 1, nkbx,              &
!     &              tmp1, 1, nx, nss, nkbus, nss )
!
!         DO j=1,nss
!            DO i=1,nss
!               tau(i,j)=tau(i,j)+tmp1(i,j)
!            END DO
!         END DO
!
!         DEALLOCATE( tmp1 )

      ENDIF

      IF(iprsta.GT.4) THEN
         WRITE( stdout,*)
         WRITE( stdout,'(26x,a)') '    tau '
         DO i=1,nss
            WRITE( stdout,'(7f11.6)') (tau(i,j),j=1,nss)
         END DO
      ENDIF
!
      RETURN
   END SUBROUTINE tauset

!
!-------------------------------------------------------------------------
   SUBROUTINE updatc( ccc, n, x0, nudx, phi, ngwx, bephi, nkbx, becp, bec, cp, nss, istart )
!-----------------------------------------------------------------------
!
      !     input ccc : dt**2/emass OR 1.0d0 demending on ortho
      !     input x0  : converged lambdas from ortho-loop (unchanged in output)
      !     input cp  : non-orthonormal cp=c0+dh/dc*ccc
      !     input bec : <cp|beta_i>
      !     input phi 
      !     output cp : orthonormal cp=cp+lambda*phi
      !     output bec: bec=becp+lambda*bephi
      !
      USE kinds,         ONLY: DP
      USE ions_base,     ONLY: nsp, na
      USE io_global,     ONLY: stdout
      USE cvan,          ONLY: nvb, ish
      USE uspp,          ONLY: nkb, nkbus
      USE uspp_param,    ONLY: nh
      USE gvecw,         ONLY: ngw
      USE control_flags, ONLY: iprint, iprsta
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: n, nudx, ngwx, nkbx, istart, nss
      COMPLEX(DP) :: cp( ngwx, n ), phi( ngwx, n )
      REAL(DP), INTENT(IN) :: ccc
      REAL(DP)    :: bec( nkbx, n ), x0( nudx, nudx )
      REAL(DP)    :: bephi( nkbx, n ), becp( nkbx, n )

      ! local variables

      INTEGER :: i, j, ig, is, iv, ia, inl
      REAL(DP),    ALLOCATABLE :: wtemp(:,:) 
      !
      !     lagrange multipliers
      !
      IF( nss < 1 ) RETURN
      !
      CALL start_clock( 'updatc' )
      
      IF ( ccc /= 1.0d0 ) THEN
         DO j = 1, nss
            CALL DSCAL( nss, ccc, x0(1,j), 1 )
         END DO
      END IF
      !
      CALL DGEMM( 'N', 'N', 2*ngw, nss, nss, 1.0d0, phi(1,istart), 2*ngwx, &
                  x0, nudx, 1.0d0, cp(1,istart), 2*ngwx )
      !    
      !     updating of the <beta|c(n,g)>
      !
      !     bec of vanderbilt species are updated 
      !
      IF( nvb > 0 )THEN

         ALLOCATE( wtemp( nss, nkb ) )

         CALL DGEMM( 'N', 'T', nss, nkbus, nss, 1.0d0, x0, nudx, &
                  bephi( 1, istart ), nkbx, 0.0d0, wtemp, nss )
         !
         DO i = 1, nss
            DO inl = 1, nkbus
               bec( inl, i + istart - 1 ) = wtemp( i, inl ) + becp( inl, i + istart - 1 )
            END DO
         END DO

         DEALLOCATE( wtemp )

      ENDIF
!
      IF ( iprsta > 2 ) THEN
         WRITE( stdout,*)
         DO is = 1, nvb
            IF( nvb .GT. 1 ) THEN
               WRITE( stdout,'(33x,a,i4)') ' updatc: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &            ((bec(ish(is)+(iv-1)*na(is)+1,i+istart-1),iv=1,nh(is)),i=1,nss)
            ELSE
               DO ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' updatc: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &            ((bec(ish(is)+(iv-1)*na(is)+ia,i+istart-1),iv=1,nh(is)),i=1,nss)
               END DO
            END IF
            WRITE( stdout,*)
         END DO
      ENDIF
!
      IF ( ccc /= 1.0d0 ) THEN
         DO j=1,nss
            CALL DSCAL(nss,1.0d0/ccc,x0(1,j),1)
         END DO
      END IF
!
      CALL stop_clock( 'updatc' )
!
      RETURN
   END SUBROUTINE updatc




!-------------------------------------------------------------------------
      SUBROUTINE calphi( c0, ngwx, bec, nkbx, betae, phi, n, ema0bg )
!-----------------------------------------------------------------------
!     input: c0 (orthonormal with s(r(t)), bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
!     where s'=s(r(t))  
!
      USE kinds,          ONLY: DP
      USE ions_base,      ONLY: na, nsp
      USE io_global,      ONLY: stdout
      USE mp_global,      ONLY: intra_image_comm
      USE cvan,           ONLY: ish, nvb
      USE uspp_param,     ONLY: nh
      USE uspp,           ONLY: nhsavb=>nkbus, qq
      USE gvecw,          ONLY: ngw
      USE constants,      ONLY: pi, fpi
      USE control_flags,  ONLY: iprint, iprsta
      USE mp,             ONLY: mp_sum
!
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ngwx, nkbx, n
      COMPLEX(DP)         :: c0( ngwx, n ), phi( ngwx, n ), betae( ngwx, nkbx )
      REAL(DP)            :: bec( nkbx, n ), emtot
      REAL(DP), OPTIONAL  :: ema0bg( ngwx )

      ! local variables
      !
      INTEGER  :: is, iv, jv, ia, inl, jnl, i, j
      REAL(DP), ALLOCATABLE :: qtemp( : , : )
!
      IF( n < 1 ) RETURN
      !
      CALL start_clock( 'calphi' )

      !
      IF ( nvb > 0 ) THEN

         ALLOCATE( qtemp( nhsavb, n ) )

         qtemp (:,:) = 0.d0
         DO is=1,nvb
            DO iv=1,nh(is)
               DO jv=1,nh(is)
                  IF(ABS(qq(iv,jv,is)) > 1.e-5) THEN
                     DO ia=1,na(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        jnl=ish(is)+(jv-1)*na(is)+ia
                        DO i=1,n
                           qtemp(inl,i) = qtemp(inl,i) +                &
     &                                    qq(iv,jv,is)*bec(jnl,i)
                        END DO
                     END DO
                  ENDIF
               END DO
            END DO
         END DO
!
!         CALL MXMA                                                     &
!     &       ( betae, 1, 2*ngwx, qtemp, 1, nhsavb, phi, 1, 2*ngwx, 2*ngw, nhsavb, n )

         CALL DGEMM &
              ( 'N', 'N', 2*ngw, n, nhsavb, 1.0d0, betae, 2*ngwx, qtemp, nhsavb, 0.0d0, phi, 2*ngwx )

         DEALLOCATE( qtemp )

      ELSE

         phi = (0.d0, 0.d0)

      END IF
!
      IF( PRESENT( ema0bg ) ) THEN
         DO j=1,n
            DO i=1,ngw
               phi(i,j)=(phi(i,j)+c0(i,j))*ema0bg(i)
            END DO
         END DO
      ELSE
         DO j=1,n
            DO i=1,ngw
               phi(i,j)=phi(i,j)+c0(i,j)
            END DO
         END DO
      END IF

      !   

      IF(iprsta > 2) THEN
         emtot=0.0d0
         IF( PRESENT( ema0bg ) ) THEN
            DO j=1,n
               DO i=1,ngw
                  emtot=emtot +2.0d0*DBLE(phi(i,j)*CONJG(c0(i,j)))*ema0bg(i)**(-2.0d0)
               END DO
            END DO
         ELSE
            DO j=1,n
               DO i=1,ngw
                  emtot=emtot +2.0d0*DBLE(phi(i,j)*CONJG(c0(i,j)))
               END DO
            END DO
         END IF
         emtot=emtot/n

         CALL mp_sum( emtot, intra_image_comm )

         WRITE( stdout,*) 'in calphi sqrt(emtot)=',SQRT(emtot)
         WRITE( stdout,*)
         DO is = 1, nvb
            IF( nvb > 1 ) THEN
               WRITE( stdout,'(33x,a,i4)') ' calphi: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &            ((bec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
            ELSE
               DO ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calphi: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &               ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
               END DO
            END IF
         END DO
      ENDIF


      CALL stop_clock( 'calphi' )
!
      RETURN
      END SUBROUTINE calphi

   END MODULE orthogonalize_base
