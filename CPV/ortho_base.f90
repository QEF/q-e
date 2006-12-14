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
      LOGICAL :: use_parallel_diag 

      PUBLIC :: sigset, rhoset, tauset
      PUBLIC :: ortho_iterate
      PUBLIC :: ortho_alt_iterate
      PUBLIC :: updatc, calphi
      PUBLIC :: mesure_diag_perf
      PUBLIC :: mesure_mmul_perf
      PUBLIC :: diagonalize_parallel
      PUBLIC :: diagonalize_serial
      PUBLIC :: use_parallel_diag

CONTAINS


!  ----------------------------------------------


   SUBROUTINE diagonalize_serial( n, rhos, rhod, s )
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: rhos(:,:) !  input symmetric matrix
      REAL(DP)              :: rhod(:)   !  output eigenvalues
      REAL(DP)              :: s(:,:)    !  output eigenvectors
      REAL(DP), ALLOCATABLE :: aux(:)
      INTEGER,  INTENT(IN)  :: n

      ALLOCATE( aux( n * ( n + 1 ) / 2 ) )

      CALL rpack( n, aux, rhos )   !  pack lower triangle of rho into aux

      CALL dspev_drv( 'V', 'L', n, aux, rhod, s, SIZE(s,1) )

      DEALLOCATE( aux )

      RETURN
   END SUBROUTINE diagonalize_serial


!  ----------------------------------------------


   SUBROUTINE diagonalize_parallel( tdist, n, rhos, rhod, s )
      USE mp,        ONLY: mp_sum
      USE mp_global, ONLY: nproc_image, me_image, intra_image_comm, np_ortho, me_ortho, ortho_comm
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: rhos(:,:) !  input symmetric matrix
      REAL(DP)              :: rhod(:)   !  output eigenvalues
      REAL(DP)              :: s(:,:)    !  output eigenvectors
      CHARACTER, INTENT(IN) :: tdist     !  if tdist == 'D' matrix "s" is block distributed
                                         !  across 2D processors grid ( ortho_comm )
      INTEGER,   INTENT(IN) :: n         !  size of the global matrix

      REAL(DP),   ALLOCATABLE :: diag(:,:)
      REAL(DP),   ALLOCATABLE :: vv(:,:)

      INTEGER :: nrl, me, np, comm
      INTEGER :: ldim_cyclic
      EXTERNAL :: ldim_cyclic

      !  distribute matrix rows to processors
      !
      ! IF( me_ortho(1) < 0 ) RETURN
      ! np = np_ortho(2) * np_ortho(1)
      ! me = me_ortho(2) + me_ortho(1) * np_ortho(2)
      ! comm = ortho_comm

      np   = nproc_image
      me   = me_image
      comm = intra_image_comm

      nrl = ldim_cyclic( n, np, me )

      ALLOCATE( diag( nrl, n ), vv( nrl, n ) )

      CALL prpack( n, diag, rhos, np, me )
      !
      CALL pdspev_drv( 'V', diag, nrl, rhod, vv, nrl, nrl, n, np, me, comm )
      !
      IF( tdist == 'D' .OR. tdist == 'd' ) THEN
         CALL cyc2blk_redist( n, vv, nrl, np, me, comm, &
                                 s, SIZE(s,1), np_ortho, me_ortho, ortho_comm )
      ELSE
         CALL prunpack( n, s, vv, np, me )
      END IF

      DEALLOCATE( diag, vv )

      RETURN
   END SUBROUTINE diagonalize_parallel


!  ----------------------------------------------


   SUBROUTINE mesure_diag_perf( n )
      !
      USE mp_global, ONLY: nproc_image, me_image, intra_image_comm, root_image
      USE mp_global, ONLY: np_ortho, me_ortho, ortho_comm
      USE io_global, ONLY: ionode, stdout
      USE mp,        ONLY: mp_sum, mp_bcast, mp_barrier, mp_group, mp_group_free
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      REAL(DP), ALLOCATABLE :: s(:,:), a(:,:), d(:)
      REAL(DP) :: t1, tpar, tser
      INTEGER  :: nr, nc
      REAL(DP) :: cclock
      INTEGER  :: ldim_block
      EXTERNAL :: cclock, ldim_block
      !
      ALLOCATE( a( n, n ), d( n ) )
      !
      IF( me_ortho(1) >= 0 ) THEN
         !
         nr = ldim_block( n, np_ortho(1), me_ortho(1) )
         nc = ldim_block( n, np_ortho(2), me_ortho(2) )
         !
         ALLOCATE( s( nr, nc ) )
         !
      ELSE
         !
         ALLOCATE( s( 1, 1 ) )
         !
      END IF
      !
      !
      CALL set_a()
      !
      CALL mp_barrier( intra_image_comm )
      t1 = cclock()
      !
      CALL diagonalize_parallel( 'D', n, a, d, s )
      !
      CALL mp_barrier( intra_image_comm )
      tpar = cclock() - t1

      DEALLOCATE( s )
      ALLOCATE( s( n, n ) )

      CALL set_a()
      !
      t1 = cclock()

      CALL diagonalize_serial( n, a, d, s )

      tser = cclock() - t1

      IF( ionode ) THEN
         use_parallel_diag = .FALSE.
         WRITE( stdout, 100 ) tpar, tser
100      FORMAT(3X,'ortho diag, time for parallel and serial driver = ', 2F9.5)
         IF( tpar < tser ) use_parallel_diag = .TRUE.
      END IF

      CALL mp_bcast( use_parallel_diag, root_image, intra_image_comm )
      
      DEALLOCATE( a, s, d )

      RETURN

   CONTAINS

      SUBROUTINE set_a()
         INTEGER :: i, j
         DO j = 1, n
            DO i = 1, n
               IF( i == j ) THEN
                  a(i,j) = ( DBLE( n-i+1 ) ) / DBLE( n ) + 1.0d0 / ( DBLE( i+j ) - 1.0d0 )
               ELSE       
                  a(i,j) = 1.0d0 / ( DBLE( i+j ) - 1.0d0 )
               END IF
            END DO        
         END DO
         RETURN
      END SUBROUTINE set_a

   END SUBROUTINE mesure_diag_perf
   

!  ----------------------------------------------


   SUBROUTINE mesure_mmul_perf( n )
      !
      USE mp_global, ONLY: nproc_image, me_image, intra_image_comm, root_image, &
                           ortho_comm, np_ortho, me_ortho, init_ortho_group
      USE io_global, ONLY: ionode, stdout
      USE mp,        ONLY: mp_sum, mp_bcast, mp_barrier, mp_group, mp_group_free
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      !
      REAL(DP), ALLOCATABLE :: c(:,:), a(:,:), b(:,:)
      REAL(DP) :: t1, tser, tcan, tbest
      INTEGER  :: nr, nc, ir, ic, np, npx, npbest
      !
      REAL(DP) :: cclock
      EXTERNAL :: cclock
      INTEGER  :: ldim_block, gind_block
      EXTERNAL :: ldim_block, gind_block
      !
      npx    = INT( SQRT( DBLE( nproc_image ) + 0.1d0 ) ) 
      tbest  = 1000
      npbest = 1

      DO np = 1, npx

         CALL init_ortho_group( np, me_image, intra_image_comm )

         IF( me_ortho(1) >= 0 ) THEN
            !
            nr = ldim_block( n, np_ortho(1), me_ortho(1) )
            nc = ldim_block( n, np_ortho(2), me_ortho(2) )
            !
            ir = gind_block( 1, n, np_ortho(1), me_ortho(1) )
            ic = gind_block( 1, n, np_ortho(2), me_ortho(2) )
            !
         ELSE
            nr = 1
            nc = 1
         END IF
   
         ALLOCATE( a( nr, nc ), c( nr, nc ), b( nr, nc ) )
   
         a = 1.0d0 / DBLE( n )
         b = 1.0d0 / DBLE( n )
   
         CALL mp_barrier( intra_image_comm )
         t1 = cclock()
   
         CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nr, b, n, 0.0d0, c, nr, &
                             np_ortho, me_ortho, ortho_comm )
   
         CALL mp_barrier( intra_image_comm )
         tcan = cclock() - t1

         IF( tcan < tbest .OR. np == 1 ) THEN
            tbest  = tcan
            npbest = np
         END IF
         !
         IF( np == 1 ) tser = tcan
   
         DEALLOCATE( a, c, b )

      END DO

      IF( ionode ) THEN
         !
         WRITE( stdout, 100 ) tser
         WRITE( stdout, 110 ) tbest, npbest*npbest
100      FORMAT(3X,'ortho mmul, time for serial driver = ', 1F9.5)
110      FORMAT(3X,'ortho mmul, best time for paralle driver = ', 1F9.5, ' with ', I4, ' procs')
         !
      END IF

      CALL init_ortho_group( npbest, me_image, intra_image_comm )

      RETURN

   END SUBROUTINE mesure_mmul_perf
   

!  ----------------------------------------------


   SUBROUTINE diagonalize_crho( n, a, d, ev )

      !  this routine calls the appropriate Lapack routine for diagonalizing a
      !  complex Hermitian matrix

         USE mp_global, ONLY: nproc_image, me_image, intra_image_comm
         USE mp, ONLY: mp_sum
         IMPLICIT NONE

         REAL(DP)             :: d(:)
         COMPLEX(DP)          :: a(:,:), ev(:,:)
         INTEGER, INTENT(IN)  :: n

         INTEGER :: nrl

         COMPLEX(DP), ALLOCATABLE :: aloc(:)
         COMPLEX(DP), ALLOCATABLE :: ap(:,:)
         COMPLEX(DP), ALLOCATABLE :: vp(:,:)
         
         IF (  nproc_image > 1  ) THEN

           nrl = n/nproc_image
           IF(me_image.LT.MOD(n,nproc_image)) THEN
             nrl = nrl + 1
           END IF

           ALLOCATE(ap(nrl,n), vp(nrl,n))

           CALL pzpack( ap, a )
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


      SUBROUTINE prpack( n, ap, a, np, me )
        IMPLICIT NONE
        REAL(DP), INTENT(IN)  :: a(:,:)
        REAL(DP), INTENT(OUT) :: ap(:,:)
        INTEGER,  INTENT(IN)  :: n, np, me
        INTEGER :: i, j, jl
        DO i = 1, SIZE( ap, 2)
           j = me + 1
           DO jl = 1, SIZE( ap, 1)
             ap( jl, i ) = a( j, i )
             j = j + np
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


      SUBROUTINE prunpack( n, a, ap, np, me )
        IMPLICIT NONE
        REAL(DP), INTENT(IN)  :: ap(:,:)
        REAL(DP), INTENT(OUT) :: a(:,:)
        INTEGER,  INTENT(IN)  :: n, np, me
        INTEGER :: i, j, jl
        DO i = 1, n
          DO j = 1, n
            a(j,i) = zero
          END DO
          j = me + 1
          DO jl = 1, SIZE(ap, 1)
            a(j,i) = ap(jl,i)
            j = j + np
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



   SUBROUTINE ortho_iterate( iter, diff, u, ldx, diag, xloc, sig, rhor, rhos, tau, nx, nss )

      !  this iterative loop uses Cannon's parallel matrix multiplication
      !  matrix are distributed over a square processor grid: 1x1 2x2 3x3 ...
      !  But the subroutine work with any number of processors, when
      !  nproc is not a square, some procs are left idle

      USE kinds,         ONLY: DP
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: ortho_eps, ortho_max
      USE mp_global,     ONLY: intra_image_comm, me_image, nproc_image, ortho_comm, &
                               np_ortho, me_ortho
      USE mp,            ONLY: mp_sum, mp_max

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx, nss, ldx
      REAL(DP) :: u   ( ldx, * )
      REAL(DP) :: diag( nx )
      REAL(DP) :: xloc( nx, nx )
      REAL(DP) :: rhor( ldx, * )
      REAL(DP) :: rhos( nx, nx )
      REAL(DP) :: tau ( ldx, * )
      REAL(DP) :: sig ( ldx, * )
      INTEGER, INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff 

      INTEGER :: i, j
      INTEGER :: nr, nc, ir, ic
      REAL(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:), dd(:,:), tr1(:,:), tr2(:,:)
      REAL(DP), ALLOCATABLE :: con(:,:), x1(:,:)
      REAL(DP), ALLOCATABLE :: xx(:,:)
      !
      integer :: ldim_block, gind_block
      external :: ldim_block, gind_block


      IF( nss < 1 ) RETURN

      !
      !  all processors not involved in the parallel orthogonalization
      !  jump at the end of the subroutine
      !

      IF( me_ortho(1) < 0 ) then
         xloc = 0.0d0
         iter = 0
         go to 100
      endif
      !
      !  Compute the size of the local block
      !
      nr = ldim_block( nss, np_ortho(1), me_ortho(1) )
      nc = ldim_block( nss, np_ortho(2), me_ortho(2) )
      !
      ir = gind_block( 1, nss, np_ortho(1), me_ortho(1) )
      ic = gind_block( 1, nss, np_ortho(2), me_ortho(2) )

      ALLOCATE( xx(nr,nc), tr1(nr,nc), tr2(nr,nc) )
      ALLOCATE( tmp1(nr,nc), tmp2(nr,nc), dd(nr,nc), x1(nr,nc), con(nr,nc) )

      do j = 1, nc
         do i = 1, nr
            xx( i, j ) = xloc( i + ir - 1, j + ic - 1 )
         end do
      end do

      xloc = 0.0d0

      ITERATIVE_LOOP: DO iter = 1, ortho_max
         !
         !       the following 4 MXMA-calls do the following matrix
         !       multiplications:
         !                       tmp1 = x0*rhor    (1st call)
         !                       dd   = x0*tau*x0  (2nd and 3rd call)
         !                       tmp2 = x0*rhos    (4th call)
         !

         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xx, nr, rhor, ldx, 0.0d0, tmp1, nr, &
                             np_ortho, me_ortho, ortho_comm )
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, tau, ldx, xx, nr, 0.0d0, tmp2, nr, &
                             np_ortho, me_ortho, ortho_comm )
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xx, nr, tmp2, nr, 0.0d0, dd, nr, &
                             np_ortho, me_ortho, ortho_comm )
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xx, nr, rhos(ir,ic), nx, 0.0d0, tmp2, nr, &
                             np_ortho, me_ortho, ortho_comm )
         !
         CALL sqr_tr_cannon( nss, tmp1, nr, tr1, nr, np_ortho, me_ortho, ortho_comm )
         CALL sqr_tr_cannon( nss, tmp2, nr, tr2, nr, np_ortho, me_ortho, ortho_comm )
         !
         DO i=1,nr
            DO j=1,nc
               x1(i,j) = sig(i,j)-tmp1(i,j)-tr1(i,j)-dd(i,j)
               con(i,j)= x1(i,j)-tmp2(i,j)-tr2(i,j)
            END DO
         END DO
         !
         !         x1      = sig      -x0*rho    -x0*rho^t  -x0*tau*x0
         !
         diff = 0.d0
         DO i=1,nr
            DO j=1,nc
               IF(ABS(con(i,j)).GT.diff) diff=ABS(con(i,j))
            END DO
         END DO

         CALL mp_max( diff, ortho_comm )

         IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

         !
         !     the following two MXMA-calls do:
         !                       tmp1 = x1*u
         !                       tmp2 = ut*x1*u
         !
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, x1, nr, u, ldx, 0.0d0, tmp1, nr, &
                             np_ortho, me_ortho, ortho_comm )
         CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, u,  ldx, tmp1,     nr, 0.0d0, tmp2, nr, &
                             np_ortho, me_ortho, ortho_comm )
         !
         !       g=ut*x1*u/d  (g is stored in tmp1)
         !
         DO i=1,nr
            DO j=1,nc
               tmp1(i,j)=tmp2(i,j)/(diag(i+ir-1)+diag(j+ic-1))
            END DO
         END DO
         !
         !       the following two MXMA-calls do:
         !                       tmp2 = g*ut
         !                       x0 = u*g*ut
         !
         CALL sqr_mm_cannon( 'N', 'T', nss, 1.0d0, tmp1,     nr, u, ldx, 0.0d0, tmp2, nr, &
                             np_ortho, me_ortho, ortho_comm )
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, u, ldx, tmp2,     nr, 0.0d0, xx,   nr, &
                             np_ortho, me_ortho, ortho_comm )
         !
      END DO ITERATIVE_LOOP

      do j = 1, nc
         do i = 1, nr
            xloc( i + ir - 1, j + ic - 1) = xx( i, j )
         end do
      end do

      DEALLOCATE( tmp1, tmp2, dd, x1, con, xx, tr1, tr2 )
            
100   continue
            
      CALL mp_sum( xloc, intra_image_comm ) 
      CALL mp_max( iter, intra_image_comm ) 


      RETURN
   END SUBROUTINE ortho_iterate


!=----------------------------------------------------------------------------=!
!
!  Alternative iterative cycle
!
!=----------------------------------------------------------------------------=!
!


   SUBROUTINE ortho_alt_iterate( iter, diff, u, ldx, diag, xloc, sig, rhor, tau, nx, n )

      USE kinds,         ONLY: DP
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: ortho_eps, ortho_max
      USE mp_global,     ONLY: intra_image_comm, me_image, nproc_image, ortho_comm, &
                               np_ortho, me_ortho
      USE mp,            ONLY: mp_sum, mp_max

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx, n, ldx
      REAL(DP) :: u   ( ldx, * )
      REAL(DP) :: diag( nx )
      REAL(DP) :: xloc( nx, nx )
      REAL(DP) :: rhor( ldx, * )
      REAL(DP) :: tau ( ldx, * )
      REAL(DP) :: sig ( ldx, * )
      INTEGER, INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff 

      INTEGER :: i, j
      INTEGER :: nr, nc, ir, ic
      REAL(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:)
      REAL(DP), ALLOCATABLE :: x1(:,:), xx(:,:)
      REAL(DP), ALLOCATABLE :: sigd(:)
      REAL(DP) :: den, dx
      !
      integer  :: ldim_block, gind_block
      external :: ldim_block, gind_block

      IF( n < 1 ) RETURN

      xloc = 0.0d0
      iter = 0

      if( me_ortho(1) < 0 ) then
         go to 100
      endif
      !
      !  Compute the size of the local block
      !
      nr = ldim_block( n, np_ortho(1), me_ortho(1) )
      nc = ldim_block( n, np_ortho(2), me_ortho(2) )
      !
      ir = gind_block( 1, n, np_ortho(1), me_ortho(1) )
      ic = gind_block( 1, n, np_ortho(2), me_ortho(2) )

      ALLOCATE( tmp1(nr,nc), tmp2(nr,nc), x1(nr,nc), sigd(nx) )
      ALLOCATE( xx(nr,nc) )

      !
      ! ...   Transform "sig", "rhoa" and "tau" in the new basis through matrix "s"
      !
      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, sig, ldx, u, ldx, 0.0d0, tmp1, nr, &
                          np_ortho, me_ortho, ortho_comm )
      CALL sqr_mm_cannon( 'T', 'N', n, 1.0d0, u, ldx, tmp1, nr, 0.0d0, sig, ldx, &
                          np_ortho, me_ortho, ortho_comm )
      !
      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, rhor, ldx, u, ldx, 0.0d0, tmp1, nr, &
                          np_ortho, me_ortho, ortho_comm )
      CALL sqr_mm_cannon( 'T', 'N', n, 1.0d0, u, ldx, tmp1, nr, 0.0d0, rhor, ldx, &
                          np_ortho, me_ortho, ortho_comm )
      !
      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, tau, ldx, u, ldx, 0.0d0, tmp1, nr, &
                          np_ortho, me_ortho, ortho_comm )
      CALL sqr_mm_cannon( 'T', 'N', n, 1.0d0, u, ldx, tmp1, nr, 0.0d0, tau, ldx, &
                          np_ortho, me_ortho, ortho_comm )
      !
      ! ...   Initialize x0 with preconditioning
      !
      DO J = 1, nc
        DO I = 1, nr
          den = ( diag( i + ir - 1 ) + diag( j + ic - 1 ) )
          IF( ABS( den ) <= small ) den = SIGN( small, den )
          xx( i, j ) = sig( i, j ) / den
        ENDDO
      ENDDO

      !
      ! ...   Starting iteration
      !

      ITERATIVE_LOOP: DO iter = 0, ortho_max

        CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, xx, nr, rhor, ldx, 0.0d0, tmp2, nr, &
                            np_ortho, me_ortho, ortho_comm )

        CALL sqr_tr_cannon( n, tmp2, nr, tmp1, nr, &
                            np_ortho, me_ortho, ortho_comm )

        DO J=1,nc
          DO I=1,nr
            tmp2(I,J) = tmp2(I,J) + tmp1(I,J)
          ENDDO
        ENDDO
!
        CALL sqr_mm_cannon( 'T', 'N', n, 1.0d0, tau, ldx, xx, nr, 0.0d0, tmp1, nr, &
                            np_ortho, me_ortho, ortho_comm )
        !
        sigd = 0.0d0
        IF( me_ortho(1) == me_ortho(2) ) THEN
           DO i = 1, nr
              SIGD( i + ir - 1 )   =  tmp1(i,i)
              tmp1(i,i) = -SIGD( i + ir - 1 )
           ENDDO
        END IF
        CALL mp_sum( sigd, ortho_comm )

        CALL sqr_mm_cannon( 'T', 'N', n, 1.0d0, xx, nr, tmp1, nr, 0.0d0, x1, nr, &
                            np_ortho, me_ortho, ortho_comm )
        !
        CALL sqr_tr_cannon( n, x1, nr, tmp1, nr, &
                            np_ortho, me_ortho, ortho_comm )

        ! ...     X1   = SIG - tmp2 - 0.5d0 * ( X1 + X1^t )

        diff = 0.0d0
        !
        DO j = 1, nc
          DO i = 1, nr
            !
            den = ( diag(i+ir-1) + sigd(i+ir-1) + diag(j+ic-1) + sigd(j+ic-1) )
            IF( ABS( den ) <= small ) den = SIGN( small, den )
            x1(i,j) = sig(i,j) - tmp2(i,j) - 0.5d0 * (x1(i,j)+tmp1(i,j))
            x1(i,j) = x1(i,j) / den
            diff = MAX( ABS( x1(i,j) - xx(i,j) ), diff )
            xx(i,j) = x1(i,j)
            !
          END DO
        END DO

        CALL mp_max( diff, ortho_comm )

        IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

      END DO ITERATIVE_LOOP
      !
      ! ...   Transform x0 back to the original basis

      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, u, ldx, xx, nr, 0.0d0, tmp1, nr, &
                          np_ortho, me_ortho, ortho_comm )
      CALL sqr_mm_cannon( 'N', 'T', n, 1.0d0, u, ldx, tmp1, nr, 0.0d0, xx, nr, &
                          np_ortho, me_ortho, ortho_comm )

      do j = 1, nc
         do i = 1, nr
            xloc( i + ir - 1, j + ic - 1) = xx( i, j )
         end do
      end do

      DEALLOCATE( tmp1, tmp2, x1, sigd )
      DEALLOCATE( xx )

100   continue

      CALL mp_sum( xloc, intra_image_comm )
      CALL mp_max( iter, intra_image_comm ) 


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
                  IF(ABS(qq(iv,jv,is)) > 1.d-5) THEN
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
