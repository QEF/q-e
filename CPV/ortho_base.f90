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
      USE dspev_module, ONLY: pdspev_drv, dspev_drv

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

      PUBLIC :: sigset
      PUBLIC :: tauset
      PUBLIC :: rhoset
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


   SUBROUTINE diagonalize_serial( n, rhos, rhod )
      IMPLICIT NONE
      INTEGER,  INTENT(IN)  :: n
      REAL(DP)              :: rhos(:,:) 
      REAL(DP)              :: rhod(:)   
      !
      ! inputs:
      ! n     size of the eigenproblem
      ! rhos  the symmetric matrix
      ! outputs:  
      ! rhos  eigenvectors
      ! rhod  eigenvalues
      !
      REAL(DP), ALLOCATABLE :: aux(:)
      INTEGER :: i, j, k

      ALLOCATE( aux( n * ( n + 1 ) / 2 ) )

      !  pack lower triangle of rho into aux
      !
      k = 0
      DO j = 1, n
         DO i = j, n
            k = k + 1
            aux( k ) = rhos( i, j )
         END DO
      END DO

      CALL dspev_drv( 'V', 'L', n, aux, rhod, rhos, SIZE(rhos,1) )

      DEALLOCATE( aux )

      RETURN

   END SUBROUTINE diagonalize_serial


!  ----------------------------------------------


   SUBROUTINE diagonalize_parallel( n, rhos, rhod, s, desc )
      USE mp,          ONLY: mp_sum
      USE mp_global,   ONLY: nproc_image, me_image, intra_image_comm
      USE descriptors, ONLY: la_myr_ , la_myc_ , la_comm_ , la_npr_ , la_npc_ , &
                             lambda_node_ , la_nrl_ , la_me_ , nlax_ , la_nrlx_
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: rhos(:,:) !  input symmetric matrix
      REAL(DP)              :: rhod(:)   !  output eigenvalues
      REAL(DP)              :: s(:,:)    !  output eigenvectors
      INTEGER,   INTENT(IN) :: n         !  size of the global matrix
      INTEGER,   INTENT(IN) :: desc(:)

      REAL(DP),   ALLOCATABLE :: diag(:,:)
      REAL(DP),   ALLOCATABLE :: vv(:,:)

      INTEGER :: nrl, me, np, comm, nrlx

      !  distribute matrix rows to processors 
      !  Matrix is distributed on the same processors group
      !  used for parallel matrix multiplication
      !
      np   = desc( la_npr_ ) * desc( la_npc_ )
      me   = desc( la_me_ )
      nrl  = desc( la_nrl_ )
      nrlx = desc( la_nrlx_ )
      comm = desc( la_comm_ )

      IF( SIZE(s,1) /= SIZE(rhos,1) .OR. SIZE(s,2) /= SIZE(rhos,2) ) &
         CALL errore( " diagonalize_parallel ", " inconsistent dimension for s and rhos ", 1 )

      IF ( desc( lambda_node_ ) > 0 ) THEN
         !
         IF( SIZE(s,1) /= desc( nlax_ ) ) &
            CALL errore( " diagonalize_parallel ", " inconsistent dimension ", 1 )
         !
         !  Compute local dimension of the cyclically distributed matrix
         !
         ALLOCATE( diag( nrlx, n ), vv( nrlx, n ) )
         !
         CALL blk2cyc_redist( n, diag, nrlx, rhos, SIZE(s,1), desc(1) ) 
         !
         CALL pdspev_drv( 'V', diag, nrlx, rhod, vv, nrlx, nrl, n, np, me, comm )
         !
         !  Redistribute matrix "vv" into "s"
         !  matrix "s" is block distributed
         !  across 2D processors grid ( ortho_comm )
         !
         CALL cyc2blk_redist( n, vv, nrlx, s, SIZE(s,1), desc(1) ) 
         !
         DEALLOCATE( diag, vv )
         !
      END IF

      RETURN

   END SUBROUTINE diagonalize_parallel


!  ----------------------------------------------


   SUBROUTINE mesure_diag_perf( n )
      !
      USE mp_global,   ONLY: nproc_image, me_image, intra_image_comm, root_image
      USE mp_global,   ONLY: np_ortho, me_ortho, ortho_comm, ortho_comm_id
      USE io_global,   ONLY: ionode, stdout
      USE mp,          ONLY: mp_sum, mp_bcast, mp_barrier
      USE mp,          ONLY: mp_max
      USE descriptors, ONLY: descla_siz_ , descla_init , nlar_ , nlac_ , &
                             ilar_ , ilac_ , nlax_ , lambda_node_ , la_myc_ , la_myr_
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      REAL(DP), ALLOCATABLE :: s(:,:), a(:,:), d(:)
      REAL(DP) :: t1, tpar, tser
      INTEGER  :: nr, nc, ir, ic, nx
      INTEGER  :: desc( descla_siz_ )
      REAL(DP) :: cclock
      EXTERNAL :: cclock
      !
      ALLOCATE( d( n ) )
      !
      CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id )

      nx = 1
      IF( desc( lambda_node_ ) > 0 ) nx = desc( nlax_ )

      nr = desc( nlar_ )
      nc = desc( nlac_ )
      ir = desc( ilar_ )
      ic = desc( ilac_ )

      ALLOCATE( s( nx, nx ) )
      ALLOCATE( a( nx, nx ) )
      !
      CALL set_a()
      !
      CALL mp_barrier( intra_image_comm )
      t1 = cclock()
      !
      CALL diagonalize_parallel( n, a, d, s, desc )
      !
      tpar = cclock() - t1
      CALL mp_max( tpar, intra_image_comm )

      DEALLOCATE( s, a )
      !
      IF( desc( la_myc_ ) == 0 .AND. desc( la_myr_ ) == 0 .AND. desc( lambda_node_ ) > 0 ) THEN

         ALLOCATE( a( n, n ) )
         nr = n
         nc = n
         ir = 1
         ic = 1

         CALL set_a()

         t1 = cclock()

         CALL diagonalize_serial( n, a, d )

         tser = cclock() - t1

         DEALLOCATE( a )

      ELSE

         tser = 0_DP

      END IF

      CALL mp_max( tser, intra_image_comm )

#if defined __PARA

      IF( ionode ) THEN
         use_parallel_diag = .FALSE.
         WRITE( stdout,  90 ) 
         WRITE( stdout, 100 ) tser
         WRITE( stdout, 110 ) tpar, np_ortho(1) * np_ortho(2)
 90      FORMAT(/,3X,'Diagonalization Performances')
100      FORMAT(3X,'ortho diag, time for serial   driver = ', 1F9.5)
110      FORMAT(3X,'ortho diag, time for parallel driver = ', 1F9.5, ' with ', I4, ' procs' )
         IF( tpar < tser ) use_parallel_diag = .TRUE.
      END IF

#else

      use_parallel_diag = .FALSE.

#endif

      CALL mp_bcast( use_parallel_diag, root_image, intra_image_comm )
      
      DEALLOCATE( d )

      RETURN

   CONTAINS

      SUBROUTINE set_a()
         INTEGER :: i, j, ii, jj
         IF( desc( lambda_node_ ) < 0 ) RETURN
         DO j = 1, nc
            DO i = 1, nr
               ii = i + ir - 1
               jj = j + ic - 1
               IF( ii == jj ) THEN
                  a(i,j) = ( DBLE( n-ii+1 ) ) / DBLE( n ) + 1.0d0 / ( DBLE( ii+jj ) - 1.0d0 )
               ELSE       
                  a(i,j) = 1.0d0 / ( DBLE( ii+jj ) - 1.0d0 )
               END IF
            END DO        
         END DO
         RETURN
      END SUBROUTINE set_a

   END SUBROUTINE mesure_diag_perf
   

!  ----------------------------------------------


   SUBROUTINE mesure_mmul_perf( n )
      !
      USE mp_global,   ONLY: nproc_image, me_image, intra_image_comm, root_image, &
                             ortho_comm, np_ortho, me_ortho, init_ortho_group, ortho_comm_id
      USE io_global,   ONLY: ionode, stdout
      USE mp,          ONLY: mp_sum, mp_bcast, mp_barrier
      USE mp,          ONLY: mp_max
      USE descriptors, ONLY: descla_siz_ , descla_init , nlar_ , nlac_ , la_comm_ , lambda_node_
      USE control_flags, ONLY: ortho_para
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      !
      REAL(DP), ALLOCATABLE :: c(:,:), a(:,:), b(:,:)
      REAL(DP) :: t1, tser, tcan, tbest
      INTEGER  :: nr, nc, ir, ic, np, npx, npbest, nps, lnode
      INTEGER  :: desc( descla_siz_ )
      !
      REAL(DP) :: cclock
      EXTERNAL :: cclock
      !
      IF( ortho_para > 0 ) THEN
         !
         !  Here the number of processors is suggested on input
         !
         np     = ortho_para 
         if( np > MIN( n, nproc_image ) ) np = MIN( n, nproc_image )
         nps    = INT( SQRT( DBLE( np ) + 0.1d0 ) ) 
         npx    = nps
         !
      ELSE
         !
         !  Guess a range or value to be tested
         !
         np     = MIN( n, nproc_image )
         nps    = 1
         npx    = INT( SQRT( DBLE( np ) + 0.1d0 ) ) 
         !
      END IF

      ! 
      !  Now test the allowed processors mesh sizes
      !

      tbest  = 1000
      npbest = 1

      DO np = nps, npx

         CALL init_ortho_group( np * np, intra_image_comm )

         CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id )

         nr = desc( nlar_ )
         nc = desc( nlac_ )
   
         ALLOCATE( a( nr, nc ), c( nr, nc ), b( nr, nc ) )
   
         a = 1.0d0 / DBLE( n )
         b = 1.0d0 / DBLE( n )
   
         CALL mp_barrier( intra_image_comm )
         t1 = cclock()

         CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nr, b, nr, 0.0d0, c, nr, desc)
   
         tcan = cclock() - t1
         CALL mp_max( tcan, intra_image_comm )

         IF( np == nps ) THEN
            tser   = tcan
            tbest  = tcan
            npbest = np
         ELSE IF( tcan < tbest ) THEN
            tbest  = tcan
            npbest = np
         END IF
         !
         DEALLOCATE( a, c, b )

      END DO

#if defined __PARA

      IF( ionode ) THEN
         !
         WRITE( stdout, 90 )
         IF( ortho_para == 0 )  THEN
            WRITE( stdout, 100 ) tser
            WRITE( stdout, 110 ) tbest, npbest*npbest
         ELSE
            WRITE( stdout, 120 ) tbest, npbest*npbest
         END IF
 90      FORMAT(/,3X,'Matrix Multiplication Performances')
100      FORMAT(3X,'ortho mmul, time for serial driver        = ', 1F9.5)
110      FORMAT(3X,'ortho mmul, best time for parallel driver = ', 1F9.5, ' with ', I4, ' procs')
120      FORMAT(3X,'ortho mmul, time for parallel driver      = ', 1F9.5, ' with ', I4, ' procs')
         !
      END IF

#else

      npbest = 1

#endif

      CALL init_ortho_group( npbest*npbest, intra_image_comm )

#if defined __PARA
      IF( ionode ) THEN
         WRITE( stdout, '(/,3X,"Constraints matrixes will be distributed block like on")' )
         WRITE( stdout, '(3X,"ortho sub-group = ", I4, "*", I4, " procs",/)' ) np_ortho(1), np_ortho(2)
      END IF
#endif

      RETURN

   END SUBROUTINE mesure_mmul_perf
   



!=----------------------------------------------------------------------------=!



   SUBROUTINE ortho_iterate( iter, diff, u, ldx, diag, xloc, nx0, sig, rhor, rhos, tau, nss, desc )

      !  this iterative loop uses Cannon's parallel matrix multiplication
      !  matrix are distributed over a square processor grid: 1x1 2x2 3x3 ...
      !  But the subroutine work with any number of processors, when
      !  nproc is not a square, some procs are left idle

      USE kinds,             ONLY: DP
      USE io_global,         ONLY: stdout
      USE control_flags,     ONLY: ortho_eps, ortho_max
      USE mp_global,         ONLY: intra_image_comm, me_image, nproc_image
      USE mp,                ONLY: mp_sum, mp_max
      USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , &
                                   la_myr_ , la_myc_ , la_comm_ , descla_siz_ , nlax_

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nss, ldx, nx0
      INTEGER, INTENT(IN) :: desc( descla_siz_ )
      REAL(DP) :: u   ( ldx, ldx )
      REAL(DP) :: diag( nss )
      REAL(DP) :: xloc( nx0, nx0 )
      REAL(DP) :: rhor( ldx, ldx )
      REAL(DP) :: rhos( ldx, ldx )
      REAL(DP) :: tau ( ldx, ldx )
      REAL(DP) :: sig ( ldx, ldx )
      INTEGER, INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff 

      INTEGER :: i, j
      INTEGER :: nr, nc, ir, ic
      REAL(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:), dd(:,:), tr1(:,:), tr2(:,:)
      REAL(DP), ALLOCATABLE :: con(:,:), x1(:,:)
      !
      IF( nss < 1 ) RETURN

      !
      !  all processors not involved in the parallel orthogonalization
      !  jump at the end of the subroutine
      !

      IF( ldx/= nx0 ) &
         CALL errore( " ortho_iterate ", " inconsistent dimensions ldx, nx0 ", nx0 )

      IF( desc( lambda_node_ ) < 0 ) then
         xloc = 0.0d0
         iter = 0
         go to 100
      endif
      !
      !  Compute the size of the local block
      !
      nr = desc( nlar_ )
      nc = desc( nlac_ )
      ir = desc( ilar_ )
      ic = desc( ilac_ )

      IF( ldx/= desc( nlax_ ) ) &
         CALL errore( " ortho_iterate ", " inconsistent dimensions ldx ", ldx )

      ALLOCATE( tr1(ldx,ldx), tr2(ldx,ldx) )
      ALLOCATE( tmp1(ldx,ldx), tmp2(ldx,ldx), dd(ldx,ldx), x1(ldx,ldx), con(ldx,ldx) )

      !  Clear elements not involved in the orthogonalization
      !
      do j = nc + 1, nx0
         do i = 1, nx0
            xloc( i, j ) = 0.0d0
         end do
      end do
      do j = 1, nx0
         do i = nr + 1, nx0
            xloc( i, j ) = 0.0d0
         end do
      end do


      ITERATIVE_LOOP: DO iter = 1, ortho_max
         !
         !       the following 4 MXMA-calls do the following matrix
         !       multiplications:
         !                       tmp1 = x0*rhor    (1st call)
         !                       dd   = x0*tau*x0  (2nd and 3rd call)
         !                       tmp2 = x0*rhos    (4th call)
         !

         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xloc, nx0, rhor, ldx, 0.0d0, tmp1, ldx, desc)
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, tau, ldx, xloc, nx0, 0.0d0, tmp2, ldx, desc)
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xloc, nx0, tmp2, ldx, 0.0d0, dd, ldx, desc)
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xloc, nx0, rhos, ldx, 0.0d0, tmp2, ldx, desc)
         !
         CALL sqr_tr_cannon( nss, tmp1, ldx, tr1, ldx, desc )
         CALL sqr_tr_cannon( nss, tmp2, ldx, tr2, ldx, desc )
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

         CALL mp_max( diff, desc( la_comm_ ) )

         IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

         !
         !     the following two MXMA-calls do:
         !                       tmp1 = x1*u
         !                       tmp2 = ut*x1*u
         !
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, x1, ldx,  u,    ldx, 0.0d0, tmp1, ldx, desc )
         CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, u,  ldx, tmp1, ldx,  0.0d0, tmp2, ldx, desc )
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
         CALL sqr_mm_cannon( 'N', 'T', nss, 1.0d0, tmp1, ldx,  u,    ldx, 0.0d0, tmp2, ldx, desc )
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, u,    ldx, tmp2, ldx,  0.0d0, xloc, nx0, desc) 
         !
      END DO ITERATIVE_LOOP

      DEALLOCATE( tmp1, tmp2, dd, x1, con, tr1, tr2 )
            
100   CONTINUE
            
      CALL mp_max( iter, intra_image_comm ) 

      RETURN
   END SUBROUTINE ortho_iterate


!=----------------------------------------------------------------------------=!
!
!  Alternative iterative cycle
!
!=----------------------------------------------------------------------------=!
!


   SUBROUTINE ortho_alt_iterate( iter, diff, u, ldx, diag, xloc, nx0, sig, rhor, tau, nss, desc )

      USE kinds,             ONLY: DP
      USE io_global,         ONLY: stdout
      USE control_flags,     ONLY: ortho_eps, ortho_max
      USE mp_global,         ONLY: intra_image_comm, me_image, nproc_image
      USE mp,                ONLY: mp_sum, mp_max
      USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , &
                                   la_myr_ , la_myc_ , la_comm_ , descla_siz_ , nlax_

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nss, ldx, nx0
      INTEGER, INTENT(IN) :: desc( descla_siz_ )
      REAL(DP) :: u   ( ldx, ldx )
      REAL(DP) :: diag( nss )
      REAL(DP) :: xloc( nx0, nx0 )
      REAL(DP) :: rhor( ldx, ldx )
      REAL(DP) :: tau ( ldx, ldx )
      REAL(DP) :: sig ( ldx, ldx )
      INTEGER, INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff 

      INTEGER :: i, j
      INTEGER :: nr, nc, ir, ic
      REAL(DP), ALLOCATABLE :: tmp1(:,:), tmp2(:,:)
      REAL(DP), ALLOCATABLE :: x1(:,:)
      REAL(DP), ALLOCATABLE :: sigd(:)
      REAL(DP) :: den, dx
      !
      IF( nss < 1 ) RETURN

      IF( ldx/= nx0 ) &
         CALL errore( " ortho_alt_iterate ", " inconsistent dimensions ldx, nx0 ", nx0 )

      if( desc( lambda_node_ ) < 0 ) then
         xloc = 0.0d0
         iter = 0
         go to 100
      endif
      !
      !  Compute the size of the local block
      !
      nr = desc( nlar_ )
      nc = desc( nlac_ )
      ir = desc( ilar_ )
      ic = desc( ilac_ )

      IF( ldx/= desc( nlax_ ) ) &
         CALL errore( " ortho_alt_iterate ", " inconsistent dimensions ldx ", ldx )

      ALLOCATE( tmp1(ldx,ldx), tmp2(ldx,ldx), x1(ldx,ldx), sigd(nss) )

      !  Clear elements not involved in the orthogonalization
      !
      do j = nc + 1, nx0
         do i = 1, nx0
            xloc( i, j ) = 0.0d0
         end do
      end do
      do j = 1, nx0
         do i = nr + 1, nx0
            xloc( i, j ) = 0.0d0
         end do
      end do
      !
      ! ...   Transform "sig", "rhoa" and "tau" in the new basis through matrix "s"
      !
      CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, sig, ldx, u, ldx, 0.0d0, tmp1, ldx, desc)
      CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, u, ldx, tmp1, ldx, 0.0d0, sig, ldx, desc)
      !
      CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, rhor, ldx, u, ldx, 0.0d0, tmp1, ldx, desc)
      CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, u, ldx, tmp1, ldx, 0.0d0, rhor, ldx, desc)
      !
      CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, tau, ldx, u, ldx, 0.0d0, tmp1, ldx, desc)
      CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, u, ldx, tmp1, ldx, 0.0d0, tau, ldx, desc)
      !
      ! ...   Initialize x0 with preconditioning
      !
      DO J = 1, nc
        DO I = 1, nr
          den = ( diag( i + ir - 1 ) + diag( j + ic - 1 ) )
          IF( ABS( den ) <= small ) den = SIGN( small, den )
          xloc( i, j ) = sig( i, j ) / den
        ENDDO
      ENDDO

      !
      ! ...   Starting iteration
      !

      ITERATIVE_LOOP: DO iter = 0, ortho_max

        CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xloc, nx0, rhor, ldx, 0.0d0, tmp2, ldx, desc)

        CALL sqr_tr_cannon( nss, tmp2, ldx, tmp1, ldx, desc )

        DO J=1,nc
          DO I=1,nr
            tmp2(I,J) = tmp2(I,J) + tmp1(I,J)
          ENDDO
        ENDDO
!
        CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, tau, ldx, xloc, nx0, 0.0d0, tmp1, ldx, desc)
        !
        sigd = 0.0d0
        IF( desc( la_myr_ ) == desc( la_myc_ ) ) THEN
           DO i = 1, nr
              SIGD( i + ir - 1 )   =  tmp1(i,i)
              tmp1(i,i) = -SIGD( i + ir - 1 )
           ENDDO
        END IF
        CALL mp_sum( sigd, desc( la_comm_ ) )

        CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, xloc, nx0, tmp1, ldx, 0.0d0, x1, ldx, desc)
        !
        CALL sqr_tr_cannon( nss, x1, ldx, tmp1, ldx, desc )

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
            diff = MAX( ABS( x1(i,j) - xloc(i,j) ), diff )
            xloc(i,j) = x1(i,j)
            !
          END DO
        END DO

        CALL mp_max( diff, desc( la_comm_ ) )

        IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

      END DO ITERATIVE_LOOP
      !
      ! ...   Transform x0 back to the original basis

      CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, u, ldx, xloc, nx0, 0.0d0, tmp1, ldx, desc)
      CALL sqr_mm_cannon( 'N', 'T', nss, 1.0d0, u, ldx, tmp1, ldx, 0.0d0, xloc, nx0, desc)

      DEALLOCATE( tmp1, tmp2, x1, sigd )

100   CONTINUE

      CALL mp_max( iter, intra_image_comm ) 

      RETURN
   END SUBROUTINE ortho_alt_iterate


!-------------------------------------------------------------------------
   SUBROUTINE sigset( cp, ngwx, becp, nkbx, qbecp, n, nss, ist, sig, ldx, desc )
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
      USE mp,                 ONLY: mp_root_sum
      USE control_flags,      ONLY: iprsta
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_image_comm, leg_ortho
      USE descriptors,        ONLY: lambda_node_ , la_npc_ , la_npr_ , descla_siz_ , &
                                    descla_init , la_comm_ , ilar_ , ilac_ , nlar_ , &
                                    nlac_ , la_myr_ , la_myc_ , la_nx_ , la_n_ , nlax_
!
      IMPLICIT NONE
!
      INTEGER nss, ist, ngwx, nkbx, n, ldx, nx
      COMPLEX(DP) :: cp( ngwx, n )
      REAL(DP)    :: becp( nkbx, n ), qbecp( nkbx, ldx )
      REAL(DP)    :: sig( ldx, ldx )
      INTEGER     :: desc( descla_siz_ )
!
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root
      INTEGER :: desc_ip( descla_siz_ )
      INTEGER :: np( 2 ), coor_ip( 2 )
      !
      REAL(DP), ALLOCATABLE :: sigp(:,:)
!
      IF( nss < 1 ) RETURN

      np(1) = desc( la_npr_ )
      np(2) = desc( la_npc_ )

      nx = desc( nlax_ )

      ALLOCATE( sigp( nx, nx ) ) 

      IF( desc( lambda_node_ ) > 0 ) THEN
         IF( desc( nlax_ ) /= ldx ) &
            CALL errore( " sigset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " sigset ", " inconsistent dimension nx ", nx )
      END IF

      DO ipc = 1, np(2)
         DO ipr = 1, np(1)

            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL descla_init( desc_ip, desc( la_n_ ), desc( la_nx_ ), np, coor_ip, desc( la_comm_ ), 1 )

            nr = desc_ip( nlar_ )
            nc = desc_ip( nlac_ )
            ir = desc_ip( ilar_ )
            ic = desc_ip( ilac_ )
            !
            CALL GRID2D_RANK( 'R', desc_ip( la_npr_ ), desc_ip( la_npc_ ), &
                                   desc_ip( la_myr_ ), desc_ip( la_myc_ ), root )

            root = root * leg_ortho

            CALL DGEMM( 'T', 'N',  nr, nc, 2*ngw, -2.0d0, cp( 1, ist + ir - 1), 2*ngwx, &
                        cp( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, sigp, nx )
            !
            !     q = 0  components has weight 1.0
            !
            IF ( gstart == 2 ) THEN
               DO jj=1,nc
                  DO ii=1,nr
                     i = ii + ir - 1
                     j = jj + ic - 1
                     sigp(ii,jj) = sigp(ii,jj) + DBLE(cp(1,i+ist-1))*DBLE(cp(1,j+ist-1))
                  END DO
               END DO
            END IF
            !
            CALL mp_root_sum( sigp, sig, root, intra_image_comm )
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( sigp )
      !
      IF( desc( lambda_node_ ) > 0 ) THEN
         !
         nr = desc( nlar_ )
         nc = desc( nlac_ )
         ir = desc( ilar_ )
         ic = desc( ilac_ )
         !
         IF( desc( la_myr_ ) == desc( la_myc_ ) ) THEN
            DO i = 1, nr
               sig(i,i) = sig(i,i) + 1.0d0
            END DO
         END IF
         !
         IF( nvb > 0 ) THEN
            CALL DGEMM( 'T', 'N', nr, nc, nkbus, -1.0d0, becp( 1, ist+ir-1 ), &
                         nkbx, qbecp( 1, 1 ), nkbx, 1.0d0, sig, ldx )
         ENDIF
         !
         IF(iprsta.GT.4) THEN
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    sig '
            DO i=1,nr
               WRITE( stdout,'(7f11.6)') (sig(i,j),j=1,nc)
            END DO
         ENDIF
         !
      END IF
      !
      RETURN
   END SUBROUTINE sigset


!
!-----------------------------------------------------------------------
   SUBROUTINE rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, ist, rho, ldx, desc )
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
      USE mp,                 ONLY: mp_root_sum
      USE mp_global,          ONLY: intra_image_comm, me_image, leg_ortho
      USE control_flags,      ONLY: iprsta
      USE io_global,          ONLY: stdout
      USE descriptors,        ONLY: lambda_node_ , la_npc_ , la_npr_ , descla_siz_ , &
                                    descla_init , la_comm_ , ilar_ , ilac_ , nlar_ , &
                                    nlac_ , la_myr_ , la_myc_ , la_nx_ , la_n_ , nlax_

!
      IMPLICIT NONE
!
      INTEGER     :: nss, ist, ngwx, nkbx, ldx, n
      COMPLEX(DP) :: cp( ngwx, n ), phi( ngwx, n )
      REAL(DP)    :: bephi( nkbx, n ), qbecp( nkbx, ldx )
      REAL(DP)    :: rho( ldx, ldx )
      INTEGER     :: desc( descla_siz_ )
      !
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root, nx
      INTEGER :: desc_ip( descla_siz_ )
      INTEGER :: np( 2 ), coor_ip( 2 )

      REAL(DP), ALLOCATABLE :: rhop(:,:)
      !
      !     <phi|cp>
      !
      !

      IF( nss < 1 ) RETURN

      np(1) = desc( la_npr_ )
      np(2) = desc( la_npc_ )

      nx = desc( nlax_ )

      IF( desc( lambda_node_ ) > 0 ) THEN
         IF( desc( nlax_ ) /= ldx ) &
            CALL errore( " rhoset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " rhoset ", " inconsistent dimension nx ", nx )
      END IF

      ALLOCATE( rhop( nx, nx ) )
     
      rhop = 0.0d0

      DO ipc = 1, np(2)
         DO ipr = 1, np(1)

            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL descla_init( desc_ip, desc( la_n_ ), desc( la_nx_ ), np, coor_ip, desc( la_comm_ ), 1 )

            nr = desc_ip( nlar_ )
            nc = desc_ip( nlac_ )
            ir = desc_ip( ilar_ )
            ic = desc_ip( ilac_ )
            !
            CALL GRID2D_RANK( 'R', desc_ip( la_npr_ ), desc_ip( la_npc_ ), &
                                   desc_ip( la_myr_ ), desc_ip( la_myc_ ), root )
            !
            root = root * leg_ortho

            CALL DGEMM( 'T', 'N', nr, nc, 2*ngw, 2.0d0, phi( 1, ist + ir - 1 ), 2*ngwx, &
                  cp( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, rhop, nx )
            !
            !     q = 0  components has weight 1.0
            !
            IF (gstart == 2) THEN
               DO jj=1,nc
                  DO ii=1,nr
                     i = ii + ir - 1
                     j = jj + ic - 1
                     rhop(ii,jj) = rhop(ii,jj) - DBLE(phi(1,i+ist-1))*DBLE(cp(1,j+ist-1))
                  END DO
               END DO
            END IF

            CALL mp_root_sum( rhop, rho, root, intra_image_comm )

         END DO
      END DO
 
      DEALLOCATE( rhop )

      IF( desc( lambda_node_ ) > 0 ) THEN
         !
         nr = desc( nlar_ )
         nc = desc( nlac_ )
         ir = desc( ilar_ )
         ic = desc( ilac_ )
         !
         IF( nvb > 0 ) THEN
            !
            ! rho(i,j) = rho(i,j) + SUM_b bephi( b, i ) * qbecp( b, j ) 
            !
            CALL DGEMM( 'T', 'N', nr, nc, nkbus, 1.0d0, bephi( 1, ist+ir-1 ), &
                        nkbx, qbecp( 1, 1 ), nkbx, 1.0d0, rho, ldx )

         END IF

         IF (iprsta.GT.4) THEN
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    rho '
            DO i=1,nr
               WRITE( stdout,'(7f11.6)') (rho(i,j),j=1,nc)
            END DO
         END IF

      END IF
      !
      RETURN
   END SUBROUTINE rhoset


!-------------------------------------------------------------------------
   SUBROUTINE tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, ist, tau, ldx, desc )
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
      USE mp,                 ONLY: mp_root_sum
      USE control_flags,      ONLY: iprsta
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_image_comm, leg_ortho
      USE descriptors,        ONLY: lambda_node_ , la_npc_ , la_npr_ , descla_siz_ , &
                                    descla_init , la_comm_ , ilar_ , ilac_ , nlar_ , &
                                    nlac_ , la_myr_ , la_myc_ , la_nx_ , la_n_ , nlax_
!
      IMPLICIT NONE
      !
      INTEGER     :: nss, ist, ngwx, nkbx, n, ldx, nx
      COMPLEX(DP) :: phi( ngwx, n )
      REAL(DP)    :: bephi( nkbx, n ), qbephi( nkbx, ldx )
      REAL(DP)    :: tau( ldx, ldx )
      INTEGER     :: desc( descla_siz_ )
      !
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root
      INTEGER :: desc_ip( descla_siz_ )
      INTEGER :: np( 2 ), coor_ip( 2 )

      REAL(DP), ALLOCATABLE :: taup( :, : )
      !
      IF( nss < 1 ) RETURN
      !
      !  get dimensions of the square processor grid
      !
      np(1) = desc( la_npr_ )
      np(2) = desc( la_npc_ )
      !
      nx = desc( nlax_ )
      !
      IF( desc( lambda_node_ ) > 0 ) THEN
         IF( desc( nlax_ ) /= ldx ) &
            CALL errore( " tauset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " tauset ", " inconsistent dimension nx ", nx )
      END IF
      !
      ALLOCATE( taup( nx, nx ) )
      !
      taup = 0.0d0
      !
      !  loop on processors coordinates
      !
      DO ipc = 1, np(2)
         !
         DO ipr = 1, np(1)

            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL descla_init( desc_ip, desc( la_n_ ), desc( la_nx_ ), np, coor_ip, desc( la_comm_ ), 1 )

            nr = desc_ip( nlar_ )
            nc = desc_ip( nlac_ )
            ir = desc_ip( ilar_ )
            ic = desc_ip( ilac_ )
            !
            CALL GRID2D_RANK( 'R', desc_ip( la_npr_ ), desc_ip( la_npc_ ), &
                                   desc_ip( la_myr_ ), desc_ip( la_myc_ ), root )
            !
            root = root * leg_ortho
            !
            !  All processors contribute to the tau block of processor (ipr,ipc)
            !  with their own part of wavefunctions
            !
            CALL DGEMM( 'T', 'N', nr, nc, 2*ngw, 2.0d0, phi( 1, ist + ir - 1 ), 2*ngwx, &
                        phi( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, taup, nx )
            !
            !           q = 0  components has weight 1.0
            !
            IF (gstart == 2) THEN
               DO jj=1,nc
                  DO ii=1,nr
                     i = ii + ir - 1
                     j = jj + ic - 1
                     taup(ii,jj) = taup(ii,jj) - DBLE(phi(1,i+ist-1))*DBLE(phi(1,j+ist-1))
                  END DO
               END DO
            END IF
            !
            CALL mp_root_sum( taup, tau, root, intra_image_comm )
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( taup )
      !
      IF( desc( lambda_node_ ) > 0 ) THEN
         !
         nr = desc( nlar_ )
         nc = desc( nlac_ )
         ir = desc( ilar_ )
         ic = desc( ilac_ )
         !
         IF( nvb > 0 ) THEN
            !
            CALL DGEMM( 'T', 'N', nr, nc, nkbus, 1.0d0, bephi( 1, ist + ir - 1 ), nkbx, &
                  qbephi( 1, 1 ), nkbx, 1.0d0, tau, ldx )
         END IF

         IF(iprsta.GT.4) THEN
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    tau '
            DO i=1,nr
               WRITE( stdout,'(7f11.6)') (tau(i,j),j=1,nc)
            END DO
         ENDIF
         !
      ENDIF
      !
      RETURN
   END SUBROUTINE tauset

!
!-------------------------------------------------------------------------
   SUBROUTINE updatc( ccc, n, x0, nx0, phi, ngwx, bephi, nkbx, becp, bec, cp, nss, istart, desc )
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
      USE kinds,             ONLY: DP
      USE ions_base,         ONLY: nsp, na
      USE io_global,         ONLY: stdout
      USE cvan,              ONLY: nvb, ish
      USE uspp,              ONLY: nkb, nkbus
      USE uspp_param,        ONLY: nh
      USE gvecw,             ONLY: ngw
      USE control_flags,     ONLY: iprint, iprsta
      USE mp,                ONLY: mp_sum, mp_bcast
      USE mp_global,         ONLY: intra_image_comm, leg_ortho
      USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , descla_siz_ , la_comm_ , &
                                   la_npc_ , la_npr_ , nlax_ , la_n_ , la_nx_ , la_myr_ , la_myc_ , &
                                   descla_init
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: n, nx0, ngwx, nkbx, istart, nss
      INTEGER, INTENT(IN) :: desc( descla_siz_ )
      COMPLEX(DP) :: cp( ngwx, n ), phi( ngwx, n )
      REAL(DP), INTENT(IN) :: ccc
      REAL(DP)    :: bec( nkbx, n ), x0( nx0, nx0 )
      REAL(DP)    :: bephi( nkbx, n ), becp( nkbx, n )

      ! local variables

      INTEGER :: i, j, ig, is, iv, ia, inl, nr, nc, ir, ic
      REAL(DP),    ALLOCATABLE :: wtemp(:,:) 
      REAL(DP),    ALLOCATABLE :: xd(:,:) 
      REAL(DP) :: beta
      INTEGER :: ipr, ipc, nx, root
      INTEGER :: np( 2 ), coor_ip( 2 )
      INTEGER :: desc_ip( descla_siz_ )
      !
      !     lagrange multipliers
      !
      IF( nss < 1 ) RETURN
      !
      IF( desc( lambda_node_ ) > 0 ) THEN
         IF( nx0 /= desc( nlax_ ) ) &
            CALL errore( " updatc ", " inconsistent dimension nx0 ", nx0 )
      END IF
      !
      !  size of the local block
      !
      nx = desc( nlax_ )
      !
      np(1) = desc( la_npr_ )
      np(2) = desc( la_npc_ )
      !
      CALL start_clock( 'updatc' )

      ALLOCATE( xd( nx, nx ) )

      IF( nvb > 0 )THEN
         ALLOCATE( wtemp( nx, nkb ) )
      END IF

      DO ipr = 1, np(1)

         beta = 0_DP
         !
         DO ipc = 1, np(2)
            !
            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL descla_init( desc_ip, desc( la_n_ ), desc( la_nx_ ), np, coor_ip, desc( la_comm_ ), 1 )

            nr = desc_ip( nlar_ )
            nc = desc_ip( nlac_ )
            ir = desc_ip( ilar_ )
            ic = desc_ip( ilac_ )
            !
            CALL GRID2D_RANK( 'R', desc_ip( la_npr_ ), desc_ip( la_npc_ ), &
                                   desc_ip( la_myr_ ), desc_ip( la_myc_ ), root )

            root = root * leg_ortho

            IF( desc( la_myr_ ) == ipr - 1 .AND. desc( la_myc_ ) == ipc - 1 .AND. desc( lambda_node_ ) > 0 ) THEN
               xd = x0 * ccc
            END IF

            CALL mp_bcast( xd, root, intra_image_comm )

            CALL DGEMM( 'N', 'N', 2*ngw, nc, nr, 1.0d0, phi(1,istart+ir-1), 2*ngwx, &
                        xd, nx, 1.0d0, cp(1,istart+ic-1), 2*ngwx )

            IF( nvb > 0 )THEN
               CALL DGEMM( 'N', 'T', nr, nkbus, nc, 1.0d0, xd, nx, &
                    bephi( 1, istart + ic - 1), nkbx, beta, wtemp, nx )
            END IF

            beta = 1_DP

         END DO
         !    
         !     updating of the <beta|c(n,g)>
         !
         !     bec of vanderbilt species are updated 
         !
         IF( nvb > 0 )THEN
            !
            ! here nr and ir are still valid, since they are the same for all procs in the same row
            !
            DO i = 1, nr
               DO inl = 1, nkbus
                  bec( inl, i + istart + ir - 2 ) = wtemp( i, inl ) + becp( inl, i + istart + ir - 2 )
               END DO
            END DO

         ENDIF

      END DO

      IF( nvb > 0 )THEN
         DEALLOCATE( wtemp )
      END IF
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
      DEALLOCATE( xd )
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
      USE uspp,           ONLY: nkbus, qq
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

         ALLOCATE( qtemp( nkbus, n ) )

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
         CALL DGEMM ( 'N', 'N', 2*ngw, n, nkbus, 1.0d0, betae, &
                       2*ngwx, qtemp, nkbus, 0.0d0, phi, 2*ngwx )

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
