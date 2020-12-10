!
! Copyright (C) 2002-2009 Quantum ESPRESSO groups
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

MODULE orthogonalize_base

      USE kinds
      USE mytime, ONLY : f_wall
#if defined(__CUDA)
      USE cudafor
      USE cublas
#endif
      
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
      PUBLIC :: updatc, calphi_bgrp
      PUBLIC :: mesure_diag_perf, mesure_mmul_perf
      PUBLIC :: use_parallel_diag
      PUBLIC :: bec_bgrp2ortho

      REAL(DP), ALLOCATABLE DEVICEATTR :: tmp1(:,:), tmp2(:,:), dd(:,:), tr1(:,:), tr2(:,:)
      REAL(DP), ALLOCATABLE DEVICEATTR :: con(:,:), x1(:,:)

CONTAINS

   SUBROUTINE allocate_local_arrays(ldx)
      INTEGER, INTENT(IN) :: ldx
      IF( ALLOCATED( tr1 ) ) THEN
         IF( SIZE( tr1, 1 ) /= ldx ) THEN
            DEALLOCATE( tmp1, tmp2, dd, x1, con )
            DEALLOCATE( tr1, tr2 )
         END IF
      END IF
      IF( .NOT. ALLOCATED( tr1 ) ) THEN
         ALLOCATE( tr1(ldx,ldx), tr2(ldx,ldx) )
         ALLOCATE( tmp1(ldx,ldx), tmp2(ldx,ldx), dd(ldx,ldx), x1(ldx,ldx), con(ldx,ldx) )
      END IF
   END SUBROUTINE allocate_local_arrays

   SUBROUTINE deallocate_local_arrays()
      IF( ALLOCATED( tr1 ) ) DEALLOCATE( tr1 )
      IF( ALLOCATED( tr2 ) ) DEALLOCATE( tr2 )
      IF( ALLOCATED( tmp1 ) ) DEALLOCATE( tmp1 )
      IF( ALLOCATED( tmp2 ) ) DEALLOCATE( tmp2 )
      IF( ALLOCATED( dd ) ) DEALLOCATE( dd )
      IF( ALLOCATED( x1 ) ) DEALLOCATE( x1 )
      IF( ALLOCATED( con ) ) DEALLOCATE( con )
   END SUBROUTINE deallocate_local_arrays

   SUBROUTINE clear_unused_elements( x, idesc )
      !
      !  Clear elements not involved in the orthogonalization
      !
      IMPLICIT NONE
      REAL(DP) DEVICEATTR :: x(:,:)
      INTEGER, INTENT(IN) :: idesc(:)
      INTEGER :: nr, nc, i, j
      INCLUDE 'laxlib.fh'
      IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) then
         x = 0.0d0
      ELSE
         nr = idesc(LAX_DESC_NR)
         nc = idesc(LAX_DESC_NC)
!$cuf kernel do(2) <<<*,*>>>
         do j = nc + 1, SIZE( x, 2 )
            do i = 1, SIZE( x, 1 )
               x( i, j ) = 0.0d0
            end do
         end do
!$cuf kernel do(2) <<<*,*>>>
         do j = 1, SIZE( x, 2 )
            do i = nr + 1, SIZE( x, 1 )
               x( i, j ) = 0.0d0
            end do
         end do
      END IF
   END SUBROUTINE

!  ----------------------------------------------

   SUBROUTINE mesure_diag_perf( n )
      !
      USE mp_bands,    ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, root_bgrp
      USE io_global,   ONLY: ionode, stdout
      USE mp,          ONLY: mp_sum, mp_bcast, mp_barrier
      USE mp,          ONLY: mp_max
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(IN) :: n
      REAL(DP), ALLOCATABLE :: s(:,:), a(:,:), d(:)
      REAL(DP) :: t1, tpar, tser
      INTEGER  :: nr, nc, ir, ic, nx
      INTEGER :: idesc(LAX_DESC_SIZE)
      INTEGER, PARAMETER :: paradim = 1000
      INTEGER :: nproc_ortho, np_ortho(2), me_ortho(2), ortho_comm, ortho_comm_id, ortho_cntx
      !
      CALL laxlib_getval( nproc_ortho = nproc_ortho, np_ortho = np_ortho, me_ortho = me_ortho, ortho_comm = ortho_comm, &
        ortho_comm_id = ortho_comm_id, ortho_cntx = ortho_cntx )
      !
      ! Check if number of PEs for orthogonalization/diagonalization is given from the input
      !
      IF( nproc_ortho > 0 ) THEN
         use_parallel_diag = .TRUE. 
         RETURN
      END IF

      ALLOCATE( d( n ) )
      !
      CALL laxlib_init_desc( idesc, n, n, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )

      nx = 1
      IF( idesc(LAX_DESC_ACTIVE_NODE)> 0 ) nx = idesc(LAX_DESC_NRCX)

      nr = idesc(LAX_DESC_NR)
      nc = idesc(LAX_DESC_NC)
      ir = idesc(LAX_DESC_IR)
      ic = idesc(LAX_DESC_IC)

      ALLOCATE( s( nx, nx ) )
      ALLOCATE( a( nx, nx ) )
      !
      CALL set_a()
      !
      ! some MPIs (OpenMPI) the first time they call a collective routine take too much
      ! time to perform initializations, then perform a dummy call to get meaningful time
      !
      CALL laxlib_diagonalize( n, a, d, s, idesc )
      !
      CALL set_a()
      !
      CALL mp_barrier( intra_bgrp_comm )
      t1 = f_wall()
      !
      CALL laxlib_diagonalize( n, a, d, s, idesc )
      !
      tpar = f_wall() - t1
      CALL mp_max( tpar, intra_bgrp_comm )

      DEALLOCATE( s, a )
      !
      IF( idesc(LAX_DESC_MYC) == 0 .AND. idesc(LAX_DESC_MYR) == 0 .AND. idesc(LAX_DESC_ACTIVE_NODE)> 0  .AND. n < paradim ) THEN

         ! when n >= paradim do not mesure serial perf, go parallel

         ALLOCATE( a( n, n ) )
         nr = n
         nc = n
         ir = 1
         ic = 1

         CALL set_a()

         t1 = f_wall()

         CALL laxlib_diagonalize( n, a, d )

         tser = f_wall() - t1

         DEALLOCATE( a )

      ELSE

         tser = 0_DP

      END IF

      CALL mp_max( tser, intra_bgrp_comm )

#if defined __MPI

      IF( ionode ) THEN
         use_parallel_diag = .FALSE.
         WRITE( stdout,  90 ) 
         IF( n < paradim ) WRITE( stdout, 100 ) tser
         WRITE( stdout, 110 ) tpar, np_ortho(1) * np_ortho(2)
 90      FORMAT(/,3X,'Diagonalization Performances')
100      FORMAT(3X,'ortho diag, time for serial   driver = ', 1F9.5)
110      FORMAT(3X,'ortho diag, time for parallel driver = ', 1F9.5, ' with ', I4, ' procs' )
         IF( n < paradim ) THEN
            IF( tpar < tser ) use_parallel_diag = .TRUE.
         ELSE
            use_parallel_diag = .TRUE.
         END IF
      END IF

#else

      use_parallel_diag = .FALSE.

#endif

      CALL mp_bcast( use_parallel_diag, root_bgrp, intra_bgrp_comm )
      
      DEALLOCATE( d )

      RETURN

   CONTAINS

      SUBROUTINE set_a()
         INTEGER :: i, j, ii, jj
         IF( idesc(LAX_DESC_ACTIVE_NODE)< 0 ) RETURN
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
      USE mp_world,    ONLY: world_comm
      USE mp_bands,    ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, &
                             root_bgrp, my_bgrp_id, nbgrp
      USE mp_images,   ONLY: nimage, my_image_id
      USE io_global,   ONLY: ionode, stdout
      USE mp,          ONLY: mp_sum, mp_bcast, mp_barrier
      USE mp,          ONLY: mp_max
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(IN) :: n
      !
      REAL(DP), ALLOCATABLE :: c(:,:), a(:,:), b(:,:)
      REAL(DP) :: t1, tcan
      INTEGER  :: nr, nc, ir, ic, lnode
      INTEGER :: idesc(LAX_DESC_SIZE)
      INTEGER :: nproc_ortho, np_ortho(2), me_ortho(2), ortho_comm, ortho_comm_id, ortho_cntx
      !
      CALL laxlib_getval( nproc_ortho = nproc_ortho, np_ortho = np_ortho, me_ortho = me_ortho, ortho_comm = ortho_comm, &
        ortho_comm_id = ortho_comm_id, ortho_cntx = ortho_cntx )

      CALL laxlib_init_desc( idesc, n, n, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )

      nr = idesc(LAX_DESC_NR)
      nc = idesc(LAX_DESC_NC)
   
      ALLOCATE( a( nr, nc ), c( nr, nc ), b( nr, nc ) )
   
      a = 1.0d0 / DBLE( n )
      b = 1.0d0 / DBLE( n )
   
      ! some MPIs (OpenMPI) the first time they call a collective routine take too much
      ! time to perform initializations, then perform a dummy call to get meaningful time
      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nr, b, nr, 0.0d0, c, nr, idesc) 

      CALL mp_barrier( intra_bgrp_comm )
      t1 = f_wall()

      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nr, b, nr, 0.0d0, c, nr, idesc)
   
      tcan = f_wall() - t1
      CALL mp_max( tcan, intra_bgrp_comm )

      DEALLOCATE( a, c, b )

#if defined __MPI

      IF( ionode ) THEN
         !
         WRITE( stdout, 90 )
         WRITE( stdout, 120 ) tcan, nproc_ortho
 90      FORMAT(/,3X,'Matrix Multiplication Performances')
120      FORMAT(3X,'ortho mmul, time for parallel driver      = ', 1F9.5, ' with ', I4, ' procs')
         !
      END IF

      IF( ionode ) THEN
         WRITE( stdout, '(/,3X,"Constraints matrixes will be distributed block like on")' )
         WRITE( stdout, '(3X,"ortho sub-group = ", I4, "*", I4, " procs",/)' ) np_ortho(1), np_ortho(2)
      END IF
#endif

      RETURN

   END SUBROUTINE mesure_mmul_perf
   



!=----------------------------------------------------------------------------=!



   SUBROUTINE ortho_iterate( iter, diff, u, ldx, diag, xloc, nx0, sig, rhor, rhos, tau, nss, idesc )

      !  this iterative loop uses Cannon's parallel matrix multiplication
      !  matrix are distributed over a square processor grid: 1x1 2x2 3x3 ...
      !  But the subroutine work with any number of processors, when
      !  nproc is not a square, some procs are left idle

      USE kinds,             ONLY: DP
      USE io_global,         ONLY: stdout
      USE control_flags,     ONLY: ortho_eps, ortho_max
      USE mp_bands,          ONLY: intra_bgrp_comm, me_bgrp, nproc_bgrp
      USE mp_world,          ONLY: mpime
      USE mp,                ONLY: mp_sum, mp_max

      IMPLICIT NONE

      include 'laxlib.fh'

      INTEGER, INTENT(IN) :: nss, ldx, nx0
      INTEGER, INTENT(IN) :: idesc(:)
      REAL(DP) DEVICEATTR :: u   ( :, : )
      REAL(DP) DEVICEATTR :: diag( : )
      REAL(DP) DEVICEATTR :: xloc( :, : )
      REAL(DP) DEVICEATTR :: rhor( :, : )
      REAL(DP) DEVICEATTR :: rhos( :, : )
      REAL(DP) DEVICEATTR :: tau ( :, : )
      REAL(DP) DEVICEATTR :: sig ( :, : )
      INTEGER, INTENT(OUT) :: iter
      REAL(DP), INTENT(OUT) :: diff 

      INTEGER :: i, j
      INTEGER :: nr, nc, ir, ic, info
      !
      IF( nss < 1 ) RETURN

      !
      !  all processors not involved in the parallel orthogonalization
      !  jump at the end of the subroutine
      !

      IF( ldx/= nx0 ) &
         CALL errore( " ortho_iterate ", " inconsistent dimensions ldx, nx0 ", nx0 )

      IF( idesc(LAX_DESC_ACTIVE_NODE) < 0 ) then
         xloc = 0.0d0
         iter = 0
         go to 100
      END IF
      !
      !  Compute the size of the local block
      !
      nr = idesc(LAX_DESC_NR)
      nc = idesc(LAX_DESC_NC)
      ir = idesc(LAX_DESC_IR)
      ic = idesc(LAX_DESC_IC)

      IF( ldx/= idesc(LAX_DESC_NRCX) ) &
         CALL errore( " ortho_iterate ", " inconsistent dimensions ldx ", ldx )

      CALL allocate_local_arrays(ldx)

      CALL clear_unused_elements( xloc, idesc )
      CALL clear_unused_elements( con, idesc )


      ITERATIVE_LOOP: DO iter = 1, ortho_max
         !
         !       the following calls do the following matrix multiplications:
         !                       tmp1 = x0*rhor    (1st call)
         !                       dd   = x0*tau*x0  (2nd and 3rd call)
         !                       tmp2 = x0*rhos    (4th call)
         !
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xloc, nx0, rhor, ldx, 0.0d0, tmp1, ldx, idesc)
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, tau, ldx, xloc, nx0, 0.0d0, tmp2, ldx, idesc)
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xloc, nx0, tmp2, ldx, 0.0d0, dd, ldx, idesc)
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, xloc, nx0, rhos, ldx, 0.0d0, tmp2, ldx, idesc)
         !
         CALL sqr_tr_cannon( nss, tmp1, ldx, tr1, ldx, idesc )
         CALL sqr_tr_cannon( nss, tmp2, ldx, tr2, ldx, idesc )
         !
#if defined (__CUDA)
!$cuf kernel do(2) <<<*,*>>>
#else
!$omp parallel do default(shared), private(j)
#endif
         DO i=1,nr
            DO j=1,nc
               x1(i,j) = sig(i,j)-tmp1(i,j)-tr1(i,j)-dd(i,j)
               con(i,j)= ABS(x1(i,j)-tmp2(i,j)-tr2(i,j))
            END DO
         END DO
         !
         !         x1      = sig      -x0*rho    -x0*rho^t  -x0*tau*x0
         !
         diff = MAXVAL( con )

         CALL mp_max( diff, idesc(LAX_DESC_COMM) )

         IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

         !
         !     the following calls do:
         !                       tmp1 = x1*u
         !                       tmp2 = ut*x1*u
         !
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, x1, ldx,  u,    ldx, 0.0d0, tmp1, ldx, idesc )
         CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, u,  ldx, tmp1, ldx,  0.0d0, tmp2, ldx, idesc )
         !
         !       g=ut*x1*u/d  (g is stored in tmp1)
         !
#if defined (__CUDA)
!$cuf kernel do(2) <<<*,*>>>
#else
!$omp parallel do default(shared), private(j)
#endif
         DO i=1,nr
            DO j=1,nc
               tmp1(i,j)=tmp2(i,j)/(diag(i+ir-1)+diag(j+ic-1))
            END DO
         END DO
         !
         !       the following calls do:
         !                       tmp2 = g*ut
         !                       x0 = u*g*ut
         !
         CALL sqr_mm_cannon( 'N', 'T', nss, 1.0d0, tmp1, ldx,  u,    ldx, 0.0d0, tmp2, ldx, idesc )
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, u,    ldx, tmp2, ldx,  0.0d0, xloc, nx0, idesc) 
         !
      END DO ITERATIVE_LOOP
            
100   CONTINUE
            
      CALL mp_max( iter, intra_bgrp_comm ) 

      RETURN
   END SUBROUTINE ortho_iterate



!-------------------------------------------------------------------------
   SUBROUTINE sigset( cp, ngwx, becp_dist, nkbx, qbecp, n, nss, ist, sig, ldx, idesc )
!-----------------------------------------------------------------------
!     input: cp (non-orthonormal), becp, qbecp
!     computes the matrix
!       sig = 1 - a ,  a = <cp|s|cp> = <cp|cp> + sum q_ij <cp|i><j|cp>
!     where s=s(r(t+dt))
!     routine makes use of c(-q)=c*(q)
!
      USE kinds,              ONLY: DP
      USE uspp,               ONLY: nkb, nkbus
      USE gvecw,              ONLY: ngw
      USE gvect, ONLY: gstart
      USE mp,                 ONLY: mp_root_sum, mp_sum
      USE control_flags,      ONLY: iverbosity
      USE io_global,          ONLY: stdout
      USE mp_bands,           ONLY: intra_bgrp_comm, inter_bgrp_comm, my_bgrp_id, nbgrp
!
      IMPLICIT NONE

      include 'laxlib.fh'
!
      INTEGER     :: nss, ist, ngwx, nkbx, n, ldx, nx
      COMPLEX(DP), INTENT(IN) DEVICEATTR :: cp(:,:) ! cp( ngwx, n )
      REAL(DP)   , INTENT(IN) DEVICEATTR :: qbecp(:,:) ! qbecp( nkbx, ldx )
      REAL(DP)   , INTENT(IN) DEVICEATTR :: becp_dist(:,:) ! becp_dist( nkbx, ldx )
      REAL(DP)   , INTENT(OUT) DEVICEATTR :: sig(:,:) ! sig( ldx, ldx )

      INTEGER, INTENT(IN) :: idesc(:)
!
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root, info
      INTEGER :: idesc_ip(LAX_DESC_SIZE)
      INTEGER :: np( 2 ), coor_ip( 2 ), leg_ortho
      !
      REAL(DP), ALLOCATABLE DEVICEATTR :: sigp(:,:)
!
      IF( nss < 1 ) RETURN

      CALL laxlib_getval( leg_ortho = leg_ortho )

      np(1) = idesc(LAX_DESC_NPR)
      np(2) = idesc(LAX_DESC_NPC)

      nx = idesc(LAX_DESC_NRCX)

      ALLOCATE( sigp( nx, nx ) ) 

      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         IF( idesc(LAX_DESC_NRCX) /= ldx ) &
            CALL errore( " sigset ", " inconsistent dimension ldx ", 1 )
         IF( nx /= ldx ) &
            CALL errore( " sigset ", " inconsistent dimension nx ", 2 )
         IF( SIZE(sig,1) /= ldx .OR.  SIZE(sig,2) /= ldx ) &
            CALL errore( " sigset ", " inconsistent dimension SIZE of sig ", 3 )
      END IF

      IF( nbgrp > 1 ) THEN
         sig = 0.0d0
      END IF

      DO ipc = 1, np(2)

         DO ipr = 1, ipc ! np(1) use symmetry

            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL laxlib_init_desc( idesc_ip, idesc(LAX_DESC_N), idesc(LAX_DESC_NX), np, coor_ip, &
                                   idesc(LAX_DESC_COMM), idesc(LAX_DESC_CNTX), 1 )

            nr = idesc_ip(LAX_DESC_NR)
            nc = idesc_ip(LAX_DESC_NC)
            ir = idesc_ip(LAX_DESC_IR)
            ic = idesc_ip(LAX_DESC_IC)
            !
            CALL GRID2D_RANK( 'R', idesc_ip(LAX_DESC_NPR), idesc_ip(LAX_DESC_NPC), &
                                   idesc_ip(LAX_DESC_MYR), idesc_ip(LAX_DESC_MYC), root )

            IF( MOD( root , nbgrp ) == my_bgrp_id ) THEN

               root = root * leg_ortho

               IF( ngw > 0 ) THEN 
                  CALL DGEMMDRV &
                       ( 'T', 'N',  nr, nc, 2*ngw, -2.0d0, cp( 1, ist + ir - 1), 2*ngwx, &
                           cp( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, sigp, nx )
               ELSE
                  sigp = 0.0d0
               END IF
               !
               !     q = 0  components has weight 1.0
               !
               IF ( gstart == 2 ) THEN
                  CALL MYDGER( nr, nc, 1.D0, cp(1,ist+ir-1), 2*ngwx, cp(1,ist+ic-1), 2*ngwx, sigp, nx )
               END IF
               !
               CALL mp_root_sum( sigp, sig, root, intra_bgrp_comm )
               !
            ENDIF
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( sigp )
      !
      IF( nbgrp > 1 ) THEN
         CALL mp_sum( sig, inter_bgrp_comm )
      END IF
      !
      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN

         CALL laxlib_dsqmsym( nss, sig, nx, idesc )
         !
         !
         nr = idesc(LAX_DESC_NR)
         nc = idesc(LAX_DESC_NC)
         ir = idesc(LAX_DESC_IR)
         ic = idesc(LAX_DESC_IC)
         !
         IF( idesc(LAX_DESC_MYR) == idesc(LAX_DESC_MYC) ) THEN
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, nr
               sig(i,i) = sig(i,i) + 1.0d0
            END DO
         END IF
         !
         IF( nkbus > 0 ) THEN
            CALL DGEMMDRV &
                 ( 'T', 'N', nr, nc, nkb, -1.0d0, becp_dist( 1, 1 ), &
                         nkbx, qbecp( 1, 1 ), nkbx, 1.0d0, sig, ldx )
         ENDIF
         !
      END IF
      !
      RETURN
   END SUBROUTINE sigset


!
!-----------------------------------------------------------------------
   SUBROUTINE rhoset( cp, ngwx, phi, bephi, nkbx, qbecp, n, nss, ist, rho, rhoa, ldx, idesc )
!-----------------------------------------------------------------------
!     input: cp (non-orthonormal), phi, bephi, qbecp
!     computes the matrix
!       rho = <s'c0|s cp> = <phi|s cp>
!     where  |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
!     where s=s(r(t+dt)) and s'=s(r(t))
!     routine makes use of  c(-q)=c*(q)
!
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE uspp,               ONLY: nkb, nkbus
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_root_sum, mp_sum
      USE mp_bands,           ONLY: intra_bgrp_comm, me_bgrp, inter_bgrp_comm, my_bgrp_id, nbgrp
      USE control_flags,      ONLY: iverbosity
      USE io_global,          ONLY: stdout
!
      IMPLICIT NONE

      include 'laxlib.fh'
!
      INTEGER     :: nss, ist, ngwx, nkbx, ldx, n
      COMPLEX(DP) DEVICEATTR :: cp( :, : ), phi( :, : )
      REAL(DP)    DEVICEATTR :: bephi( :, : ), qbecp( :, : )
      REAL(DP)    DEVICEATTR :: rho( :, : )
      REAL(DP)    DEVICEATTR :: rhoa( :, : )
      INTEGER, INTENT(IN) :: idesc(:)
      !
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root, nx, info
      INTEGER :: idesc_ip(LAX_DESC_SIZE)
      INTEGER :: np( 2 ), coor_ip( 2 ), leg_ortho

      REAL(DP), ALLOCATABLE DEVICEATTR :: rhop(:,:)
      !
      !     <phi|cp>
      !
      IF( nss < 1 ) RETURN

      CALL laxlib_getval( leg_ortho = leg_ortho )

      np(1) = idesc(LAX_DESC_NPR)
      np(2) = idesc(LAX_DESC_NPC)

      nx = idesc(LAX_DESC_NRCX)

      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         IF( idesc(LAX_DESC_NRCX) /= ldx ) &
            CALL errore( " rhoset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " rhoset ", " inconsistent dimension nx ", nx )
         IF( SIZE(rho,1) /= ldx ) &
            CALL errore( " rhoset ", " inconsistent rho size ", ldx )
         IF( SIZE(bephi,1) /= nkbx ) &
            CALL errore( " rhoset ", " inconsistent bephi size ", nkbx )
         IF( SIZE(qbecp,1) /= nkbx ) &
            CALL errore( " rhoset ", " inconsistent qbecp size ", nkbx )
      END IF

      ALLOCATE( rhop( nx, nx ) )
     
      rhop = 0.0d0
      IF( nbgrp > 1 ) THEN
         rho = 0.0d0
      END IF

      DO ipc = 1, np(2)
         DO ipr = 1, np(1)

            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL laxlib_init_desc( idesc_ip, idesc(LAX_DESC_N), idesc(LAX_DESC_NX), np, coor_ip, &
                                   idesc(LAX_DESC_COMM), idesc(LAX_DESC_CNTX), 1 )

            nr = idesc_ip(LAX_DESC_NR)
            nc = idesc_ip(LAX_DESC_NC)
            ir = idesc_ip(LAX_DESC_IR)
            ic = idesc_ip(LAX_DESC_IC)
            !
            CALL GRID2D_RANK( 'R', idesc_ip(LAX_DESC_NPR), idesc_ip(LAX_DESC_NPC), &
                                   idesc_ip(LAX_DESC_MYR), idesc_ip(LAX_DESC_MYC), root )
            !
            IF( MOD( root , nbgrp ) == my_bgrp_id ) THEN

               root = root * leg_ortho

               IF( ngw > 0 ) THEN
                  CALL DGEMMDRV &
                       ('T', 'N', nr, nc, 2*ngw, 2.0d0, phi( 1, ist + ir - 1 ), 2*ngwx, &
                              cp( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, rhop, nx )
               ELSE
                  rhop = 0.0d0
               END IF
               !
               !     q = 0  components has weight 1.0
               !
               IF (gstart == 2) THEN
                  CALL mydger ( nr, nc, -1.D0, phi(1,ist+ir-1), 2*ngwx, cp(1,ist+ic-1), 2*ngwx, rhop, nx )
               END IF

               CALL mp_root_sum( rhop, rho, root, intra_bgrp_comm )

            END IF

         END DO
      END DO
 
      IF( nbgrp > 1 ) THEN
         CALL mp_sum( rho, inter_bgrp_comm )
      END IF
      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         !
         nr = idesc(LAX_DESC_NR)
         nc = idesc(LAX_DESC_NC)
         !
         !  bephi is distributed among processor rows
         !  qbephi is distributed among processor columns
         !  tau is block distributed among the whole processor 2D grid
         !
         !
         IF( nkbus > 0 ) THEN
            !
            ! rho(i,j) = rho(i,j) + SUM_b bephi( b, i ) * qbecp( b, j ) 
            !
            CALL DGEMMDRV &
                 ( 'T', 'N', nr, nc, nkb, 1.0d0, bephi, nkbx, qbecp, nkbx, 1.0d0, rho, ldx )

         END IF

      END IF
      !
      ! NOW get the symmetric and antisymmetric part of rho
      !
      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         !
         nr = idesc(LAX_DESC_NR)
         nc = idesc(LAX_DESC_NC)
         !
         !    distributed array rhos contains "rho",
         !    now transpose rhos and store the result in distributed array rhot
         !
         CALL sqr_tr_cannon( nss, rho, SIZE(rho,1), rhop, SIZE(rhop,1), idesc )
         !
         !  Compute the symmetric part of rho
         !
!$cuf kernel do(2) <<<*,*>>>
         DO j = 1, nc
            DO i = 1, nr
               rho( i, j ) = 0.5d0 * ( rho( i, j ) + rhop( i, j ) )
            END DO
         END DO
         !
         !  distributed array rhos now contains symmetric part of "rho",
         !
         !  Antisymmetric part of rho, alredy distributed across ortho procs.
         !
!$cuf kernel do(2) <<<*,*>>>
         DO j = 1, nc
            DO i = 1, nr
               rhoa( i, j ) = rho( i, j ) - rhop( i, j )
            END DO
         END DO
         !
      END IF
      !
      DEALLOCATE( rhop )
      !
      RETURN
   END SUBROUTINE rhoset


!-------------------------------------------------------------------------
   SUBROUTINE tauset( phi, ngwx, bephi, nkbx, qbephi, n, nss, ist, tau, ldx, idesc )
!-----------------------------------------------------------------------
!     input: phi
!     computes the matrix
!        tau = <s'c0|s|s'c0> = <phi|s|phi>,  where  |phi> = s'|c0>
!     where s=s(r(t+dt)) and s'=s(r(t))
!     routine makes use of c(-q)=c*(q)
!
      USE kinds,              ONLY: DP
      USE uspp,               ONLY: nkb, nkbus
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE mp,                 ONLY: mp_root_sum, mp_sum, mp_barrier
      USE control_flags,      ONLY: iverbosity
      USE io_global,          ONLY: stdout
      USE mp_bands,           ONLY: intra_bgrp_comm, inter_bgrp_comm, my_bgrp_id, nbgrp
      USE mp_world,           ONLY: mpime
!
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(IN)    :: nss, ist, ngwx, nkbx, n, ldx
      COMPLEX(DP) DEVICEATTR :: phi( :, : )
      REAL(DP)    DEVICEATTR :: bephi( :, : ), qbephi( :, : )
      REAL(DP)    DEVICEATTR :: tau( :, : )
      INTEGER, INTENT(IN) :: idesc(:)
      !
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root, info, nx
      INTEGER :: idesc_ip(LAX_DESC_SIZE)
      INTEGER :: np( 2 ), coor_ip( 2 ), leg_ortho

      REAL(DP), ALLOCATABLE DEVICEATTR :: taup( :, : )
      !
      IF( nss < 1 ) RETURN

      CALL laxlib_getval( leg_ortho = leg_ortho )
      !
      !  get dimensions of the square processor grid
      !
      np(1) = idesc(LAX_DESC_NPR)
      np(2) = idesc(LAX_DESC_NPC)
      !
      nx = idesc(LAX_DESC_NRCX)
      !
      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         IF( idesc(LAX_DESC_NRCX) /= ldx ) &
            CALL errore( " tauset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " tauset ", " inconsistent dimension nx ", nx )
         IF( nx /= SIZE(tau,1) ) &
            CALL errore( " tauset ", " inconsistent dimension nx ", nx )
      END IF
      !
      ALLOCATE( taup( nx, nx ) )
      !
      taup = 0.0d0
      IF( nbgrp > 1 ) THEN
         tau = 0.0d0
      END IF
      !
      !  loop on processors coordinates
      !
      DO ipc = 1, np(2)
         !
         DO ipr = 1, ipc ! np(1)  use symmetry

            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL laxlib_init_desc( idesc_ip, idesc(LAX_DESC_N), idesc(LAX_DESC_NX), np, coor_ip, &
                                   idesc(LAX_DESC_COMM), idesc(LAX_DESC_CNTX), 1 )

            nr = idesc_ip(LAX_DESC_NR)
            nc = idesc_ip(LAX_DESC_NC)
            ir = idesc_ip(LAX_DESC_IR)
            ic = idesc_ip(LAX_DESC_IC)
            !
            CALL GRID2D_RANK( 'R', idesc_ip(LAX_DESC_NPR), idesc_ip(LAX_DESC_NPC), &
                                   idesc_ip(LAX_DESC_MYR), idesc_ip(LAX_DESC_MYC), root )
            !
            IF( MOD( root , nbgrp ) == my_bgrp_id ) THEN

               root = root * leg_ortho
               !
               !  All processors contribute to the tau block of processor (ipr,ipc)
               !  with their own part of wavefunctions
               !
               IF( ngw > 0 ) THEN
                  CALL DGEMMDRV &
                       ( 'T', 'N', nr, nc, 2*ngw, 2.0d0, phi( 1, ist + ir - 1 ), 2*ngwx, &
                           phi( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, taup, nx )
               ELSE
                  taup = 0.0d0
               END IF
               !
               !           q = 0  components has weight 1.0
               !
               IF (gstart == 2) THEN
                  CALL MYDGER( nr, nc, -1.D0, phi(1,ist+ir-1), 2*ngwx, phi(1,ist+ic-1), 2*ngwx, taup, nx )
               END IF
               !
               CALL mp_root_sum( taup, tau, root, intra_bgrp_comm )
               !
            END IF
            !
         END DO
         !
      END DO
      !
      DEALLOCATE( taup )
      !
      IF( nbgrp > 1 ) THEN
         CALL mp_sum( tau, inter_bgrp_comm )
      END IF
      !
      IF( idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
         !
         CALL laxlib_dsqmsym( nss, tau, nx, idesc )
         !
         nr = idesc(LAX_DESC_NR)
         nc = idesc(LAX_DESC_NC)
         !
         !  bephi is distributed among processor rows
         !  qbephi is distributed among processor columns
         !  tau is block distributed among the whole processor 2D grid
         !
         IF( nkbus > 0 ) THEN
            !
            CALL DGEMMDRV &
                 ( 'T', 'N', nr, nc, nkb, 1.0d0, bephi, nkbx, qbephi, nkbx, 1.0d0, tau, ldx )
            !
         END IF

      ENDIF
      !
      RETURN
   END SUBROUTINE tauset

!
!-------------------------------------------------------------------------
   SUBROUTINE updatc( ccc, x0, phi, bephi, becp_bgrp, bec_bgrp, cp_bgrp, idesc )
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
      USE uspp,              ONLY: nkb, nkbus
      USE uspp_param,        ONLY: nh
      USE gvecw,             ONLY: ngw
      USE control_flags,     ONLY: iverbosity
      USE mp,                ONLY: mp_sum, mp_bcast
      USE mp_bands,          ONLY: intra_bgrp_comm, me_bgrp, inter_bgrp_comm
      USE electrons_base,    ONLY: nbspx_bgrp, ibgrp_g2l, nbsp, nspin,  nupdwn, iupdwn, nbspx
      USE device_memcpy_m,     ONLY: dev_memcpy
!
      USE device_memcpy_m,     ONLY : dev_memcpy
!
      IMPLICIT NONE

      include 'laxlib.fh'
!
      INTEGER, INTENT(IN) :: idesc( :, : )
      COMPLEX(DP) :: cp_bgrp( :, : ), phi( :, : )
      REAL(DP), INTENT(IN) :: ccc
      REAL(DP)    :: bec_bgrp( :, : ), x0( :, :, : )
      REAL(DP)    :: bephi( :, : )
      REAL(DP)    :: becp_bgrp( :, : )

      ! local variables

      INTEGER :: i, j, ig, is, iv, ia, inl, nr, nc, ir, ic, nx0, ngwx, nkbx, iss, nrcx
      INTEGER :: ipr, ipc, root, i1, i2, nss, istart
      INTEGER :: ibgrp_i, ibgrp_i_first, nbgrp_i, i_first
      REAL(DP),    ALLOCATABLE :: xd(:,:)
      REAL(DP),    ALLOCATABLE :: bephi_tmp(:,:) 
      INTEGER,     ALLOCATABLE :: indi(:)
      INTEGER :: np( 2 ), coor_ip( 2 ), leg_ortho
      INTEGER :: idesc_ip(LAX_DESC_SIZE)
#if defined (__CUDA)
      ATTRIBUTES( DEVICE ) :: xd, becp_bgrp, bephi, cp_bgrp, phi, bephi_tmp, bec_bgrp, indi
#endif

      CALL start_clock( 'updatc' )

#if defined (__CUDA)
      IF( nkbus > 0 )THEN
         ALLOCATE(indi, SOURCE=ibgrp_g2l)
      END IF
#endif

      CALL laxlib_getval( leg_ortho = leg_ortho )

      DO iss = 1, nspin
         !
         !  size of the local block
         !
         nrcx = idesc( LAX_DESC_NRCX, iss )
         !
         nss = nupdwn(iss)
         istart = iupdwn(iss)
         i1 = (iss-1)*nrcx+1
         i2 = iss*nrcx
         nx0 = SIZE( x0, 1 )
         ngwx = SIZE( phi, 1 )
         nkbx = SIZE( bephi, 1 )
         !
         !     lagrange multipliers
         !
         IF( nss < 1 ) CYCLE
         !
         IF( idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) THEN
            IF( nx0 /= idesc( LAX_DESC_NRCX, iss ) ) &
               CALL errore( " updatc ", " inconsistent dimension nx0 ", nx0 )
         END IF
         !
         np(1) = idesc( LAX_DESC_NPR, iss )
         np(2) = idesc( LAX_DESC_NPC, iss )
         !
   
         ALLOCATE( xd( nrcx, nrcx ) )
   
         IF( nkbus > 0 )THEN
#if defined (__CUDA)
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, nss
               IF( indi( i + istart - 1 ) > 0 ) THEN
                  bec_bgrp( :, indi( i + istart - 1 ) ) = becp_bgrp( :, indi( i + istart - 1 ) )
               END IF
            END DO
#else
            DO i = 1, nss
               ibgrp_i = ibgrp_g2l( i + istart - 1 )
               IF( ibgrp_i > 0 ) THEN
                  bec_bgrp( :, ibgrp_i ) = becp_bgrp( :, ibgrp_i )
               END IF
            END DO
#endif
            ALLOCATE( bephi_tmp( nkbx, nrcx ) )
         END IF
   
   
         DO ipc = 1, np(2)
            !
            IF( nkbus > 0 )THEN
               ! 
               ! For the inner loop we need the block of bephi( :, ic : ic + nc - 1 )
               ! this is the same of block bephi( :, ir : ir + nr - 1 ) on processor
               ! with coords ipr == ipc
               !
               ! get the right processor owning the block of bephi
               !
               CALL GRID2D_RANK( 'R', np(1), np(2), ipc-1, ipc-1, root )
               root = root * leg_ortho
               !
               ! broadcast the block to all processors 
               ! 
               IF( me_bgrp == root ) THEN
                  !bephi_tmp(:,:) = bephi(:, i1 : i1+nrcx-1 )
                  CALL dev_memcpy(bephi_tmp, bephi(:, i1:), [1, nkbx], 1 , [1, nrcx])
               END IF
               CALL mp_bcast( bephi_tmp, root, intra_bgrp_comm )
               !
            END IF
   
            DO ipr = 1, np(1)
               !
               ! Compute the descriptor of processor with coord: ( ipr-1, ipc-1 ), in the ortho group
               !
               coor_ip(1) = ipr - 1
               coor_ip(2) = ipc - 1
   
               CALL laxlib_init_desc( idesc_ip, idesc( LAX_DESC_N, iss ), idesc( LAX_DESC_NX, iss ), np, coor_ip, &
                                      idesc( LAX_DESC_COMM, iss ), idesc( LAX_DESC_CNTX, iss ), 1 )
   
               nr = idesc_ip(LAX_DESC_NR)
               nc = idesc_ip(LAX_DESC_NC)
               ir = idesc_ip(LAX_DESC_IR)
               ic = idesc_ip(LAX_DESC_IC)
               !
               CALL GRID2D_RANK( 'R', idesc_ip(LAX_DESC_NPR), idesc_ip(LAX_DESC_NPC), &
                                      idesc_ip(LAX_DESC_MYR), idesc_ip(LAX_DESC_MYC), root )
               !
               ! we need to update only states local to the current band group,
               ! so here we compute the overlap between ortho and band group.
               !
               nbgrp_i = 0
               DO i = 1, nc
                  ibgrp_i = ibgrp_g2l( i + istart + ic - 2 )
                  IF( ibgrp_i > 0 ) THEN
                     IF( nbgrp_i == 0 ) THEN
                        ibgrp_i_first = ibgrp_i
                        i_first = i
                     END IF
                     nbgrp_i = nbgrp_i + 1
                  END IF
               END DO
   
               root = root * leg_ortho
   
               IF( idesc( LAX_DESC_MYR, iss ) == ipr - 1 .AND. &
                   idesc( LAX_DESC_MYC, iss ) == ipc - 1 .AND. &
                   idesc( LAX_DESC_ACTIVE_NODE, iss ) > 0 ) THEN
                   ! xd = x0(:,:,iss) * ccc
                   CALL dev_memcpy( xd(:,:) , x0(:, :, iss ) )
!$cuf kernel do(2) <<<*,*>>>
                   DO j = 1, SIZE( xd, 2 )
                      DO i = 1, SIZE( xd, 1 )
                         xd(i,j) = ccc * xd(i,j)
                      END DO
                   END DO
               END IF
   
               CALL mp_bcast( xd, root, intra_bgrp_comm )
   
               IF( ngw > 0 ) THEN
                  CALL DGEMMDRV &
                      ( 'N', 'N', 2*ngw, nbgrp_i, nr, 1.0d0, phi(1,istart+ir-1), 2*ngwx, &
                        xd(1,i_first), nrcx, 1.0d0, cp_bgrp(1,ibgrp_i_first), 2*ngwx )
               END IF
   
               IF( nkbus > 0 )THEN
                  nbgrp_i = 0
                  DO i = 1, nr
                     ibgrp_i = ibgrp_g2l( i + istart + ir - 2 )
                     IF( ibgrp_i > 0 ) THEN
                        IF( nbgrp_i == 0 ) THEN
                           ibgrp_i_first = ibgrp_i
                           i_first = i
                        END IF
                        nbgrp_i = nbgrp_i + 1
                     END IF
                  END DO
                  IF( nbgrp_i > 0 ) THEN
                     CALL DGEMMDRV &
                          ( 'N', 'T', nkb, nbgrp_i, nc, 1.0d0, bephi_tmp(1,1), nkbx, &
                            xd(i_first,1), nrcx, 1.0d0, bec_bgrp( 1, ibgrp_i_first ), SIZE(bec_bgrp,1) )
                  END IF
                  !
               END IF
   
            END DO
            !    
         END DO
   
         IF( nkbus > 0 )THEN
            DEALLOCATE( bephi_tmp )
         END IF
         !
         DEALLOCATE(xd)
         !
      END DO
      !
      IF(ALLOCATED(indi)) DEALLOCATE(indi)
      !
      CALL stop_clock( 'updatc' )
      !
      RETURN
   END SUBROUTINE updatc



!-------------------------------------------------------------------------
   SUBROUTINE calphi_bgrp( c0_bgrp, ngwx, bec_bgrp, nkbx, betae, phi_bgrp, nbspx_bgrp, ema0bg )
!-----------------------------------------------------------------------
!     input: c0 (orthonormal with s(r(t)), bec=<c0|beta>, betae=|beta>
!     computes the matrix phi (with the old positions)
!       where  |phi> = s'|c0> = |c0> + sum q_ij |i><j|c0>
!     where s'=s(r(t))  
!
      USE kinds,          ONLY: DP
      USE ions_base,      ONLY: nat, ityp
      USE io_global,      ONLY: stdout
      USE mp_bands,       ONLY: intra_bgrp_comm, inter_bgrp_comm
      USE uspp_param,     ONLY: nh, upf
      USE uspp,           ONLY: nkb, nkbus, qq_nt, indv_ijkb0
      USE gvecw,          ONLY: ngw
      USE electrons_base, ONLY: nbsp_bgrp, nbsp
      USE constants,      ONLY: pi, fpi
      USE control_flags,  ONLY: iverbosity
      USE mp,             ONLY: mp_sum
#if defined (__CUDA)
      USE cudafor
#endif
!
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ngwx, nkbx, nbspx_bgrp
      COMPLEX(DP)         :: c0_bgrp( :, : ), phi_bgrp( :, : )
      COMPLEX(DP)         :: betae( :, : )
      REAL(DP)            :: bec_bgrp( :, : ), emtot
      REAL(DP), OPTIONAL  :: ema0bg( : )

      ! local variables
      !
      INTEGER  :: is, iv, jv, ia, inl, jnl, i, j, indv
      REAL(DP), ALLOCATABLE :: qtemp( : , : )
      REAL(DP), ALLOCATABLE :: qtemp_d( : , : )
      REAL(DP), ALLOCATABLE :: ema0bg_d( : )
      REAL(DP) :: qqf
#if defined (__CUDA)
      ATTRIBUTES( DEVICE ) :: c0_bgrp, phi_bgrp, ema0bg_d, qtemp_d, betae
#endif
!
      IF( nbsp_bgrp < 1 ) RETURN
      !
      CALL start_clock( 'calphi' )
      !
      ! Note that phi here is computed only for my band group
      !
      IF ( nkbus > 0 ) THEN

         ALLOCATE( qtemp( nkb, nbspx_bgrp ) )

         qtemp (:,:) = 0.d0
!$omp parallel do default(none) &
!$omp shared(nat,ityp,upf,nh,indv_ijkb0,qq_nt,qtemp,bec_bgrp,nbsp_bgrp) &
!$omp private(ia,is,iv,inl,jv,jnl,qqf,i,indv)
         DO ia = 1, nat
            is = ityp(ia)
            indv = indv_ijkb0(ia)
            IF( upf(is)%tvanp ) THEN
               DO iv=1,nh(is)
                  inl = indv + iv 
                  DO jv=1,nh(is)
                     jnl = indv + jv
                     IF(ABS(qq_nt(iv,jv,is)) > 1.d-5) THEN
                        qqf = qq_nt(iv,jv,is)
                        DO i=1,nbsp_bgrp
                           qtemp(inl,i) = qtemp(inl,i) + qqf * bec_bgrp(jnl,i)
                        END DO
                     ENDIF
                  END DO
               END DO
            END IF
         END DO
!$omp end parallel do
!
         IF( ngw > 0 ) THEN
#if defined (__CUDA)
            ALLOCATE(qtemp_d, SOURCE=qtemp)
            CALL MYDGEMM ( 'N', 'N', 2*ngw, nbsp_bgrp, nkb, 1.0d0, betae, &
                       2*ngwx, qtemp_d, nkb, 0.0d0, phi_bgrp, 2*ngwx )
            DEALLOCATE(qtemp_d)
#else
            CALL dgemm ( 'N', 'N', 2*ngw, nbsp_bgrp, nkb, 1.0d0, betae, &
                       2*ngwx, qtemp, nkb, 0.0d0, phi_bgrp, 2*ngwx )
#endif
         ELSE
            phi_bgrp = 0.0d0
         END IF

         DEALLOCATE( qtemp )

      ELSE

         phi_bgrp = (0.d0, 0.d0)

      END IF
!
      IF( PRESENT( ema0bg ) ) THEN
#if defined (__CUDA)
         ALLOCATE(ema0bg_d, SOURCE=ema0bg )
!$cuf kernel do(2) <<<*,*>>>
         DO j=1,nbsp_bgrp
            DO i=1,ngw
               phi_bgrp(i,j)=(phi_bgrp(i,j)+c0_bgrp(i,j))*ema0bg_d(i)
            END DO
         END DO
         DEALLOCATE(ema0bg_d)
#else
!$omp parallel do default(shared), private(i)
         DO j=1,nbsp_bgrp
            DO i=1,ngw
               phi_bgrp(i,j)=(phi_bgrp(i,j)+c0_bgrp(i,j))*ema0bg(i)
            END DO
         END DO
!$omp end parallel do
#endif
      ELSE
#if defined (__CUDA)
!$cuf kernel do(2) <<<*,*>>>
         DO j=1,nbsp_bgrp
            DO i=1,ngw
               phi_bgrp(i,j)=phi_bgrp(i,j)+c0_bgrp(i,j)
            END DO
         END DO
#else
!$omp parallel do default(shared), private(i)
         DO j=1,nbsp_bgrp
            DO i=1,ngw
               phi_bgrp(i,j)=phi_bgrp(i,j)+c0_bgrp(i,j)
            END DO
         END DO
!$omp end parallel do
#endif
      END IF
      !   
      CALL stop_clock( 'calphi' )
!
      RETURN
   END SUBROUTINE calphi_bgrp


   SUBROUTINE bec_bgrp2ortho( bec_bgrp, bec_ortho, nrcx, idesc )
      USE kinds,             ONLY: DP
      USE uspp,              ONLY: nkb
      USE mp,                ONLY: mp_sum
      USE mp_bands,          ONLY: intra_bgrp_comm, me_bgrp, inter_bgrp_comm, nbgrp
      USE electrons_base,    ONLY: nbspx_bgrp, ibgrp_g2l, nspin
      !
      USE device_memcpy_m,     ONLY : dev_memcpy
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER, INTENT(IN) :: nrcx
      INTEGER, INTENT(IN) :: idesc( :, : )
      REAL(DP), INTENT(IN)  :: bec_bgrp(:,:)
      REAL(DP), INTENT(OUT) :: bec_ortho(:,:)
#if defined (__CUDA)
      ATTRIBUTES( DEVICE ) :: bec_bgrp, bec_ortho
#endif
      !
      INTEGER :: ir, nr, i, ibgrp_i, nup, leg_ortho
      !
      CALL laxlib_getval( leg_ortho = leg_ortho )

      bec_ortho = 0.0d0
      !
      IF( idesc( LAX_DESC_ACTIVE_NODE, 1 ) > 0 ) THEN
         ir = idesc(LAX_DESC_IR, 1)
         nr = idesc(LAX_DESC_NR, 1)
         IF( nbgrp == 1 ) THEN
            !bec_ortho(:,1:nr) = bec_bgrp(:,ir:ir+nr-1)
            CALL dev_memcpy(bec_ortho, bec_bgrp(:, ir:), [1, ubound(bec_bgrp)], 1 , [1, nr])
         ELSE
            DO i = 1, nr
               ibgrp_i = ibgrp_g2l( i + ir - 1 )
               IF( ibgrp_i > 0 ) THEN
                  !bec_ortho( :, i ) = bec_bgrp( :, ibgrp_i )
                  CALL dev_memcpy(bec_ortho(:, i), bec_bgrp(:, ibgrp_i), [1, ubound(bec_ortho)])
               END IF
            END DO
         END IF
      END IF
      !
      IF( nspin == 2 ) THEN
         IF( idesc( LAX_DESC_ACTIVE_NODE, 2 ) > 0 ) THEN
            nup = idesc(LAX_DESC_N, 1 )
            ir = idesc(LAX_DESC_IR, 2 )
            nr = idesc(LAX_DESC_NR, 2 )
            IF( nbgrp == 1 ) THEN
               !bec_ortho( :, nrcx+1 : nrcx+nr ) = bec_bgrp( :, nup+ir:nup+ir+nr-1 )
               CALL dev_memcpy(bec_ortho( :, nrcx+1:), bec_bgrp(:, nup+ir:), [1, ubound(bec_ortho)], 1 , [1, nr])
            ELSE
               do i = 1, nr
                  ibgrp_i = ibgrp_g2l( i + ir - 1 + nup )
                  IF( ibgrp_i > 0 ) THEN
                     !bec_ortho( :, i + nrcx ) = bec_bgrp( :, ibgrp_i )
                     CALL dev_memcpy(bec_ortho(:, i+nrcx ), bec_bgrp(:, ibgrp_i), [1, ubound(bec_ortho)])
                  END IF
               end do
            END IF
         END IF
      END IF
      !
      CALL mp_sum( bec_ortho, inter_bgrp_comm )
      !
      RETURN
   END SUBROUTINE bec_bgrp2ortho

END MODULE orthogonalize_base
