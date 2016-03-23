!
! Copyright (C) 2002-2009 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


MODULE orthogonalize_base


      USE kinds
      USE dspev_module, ONLY: diagonalize_serial, diagonalize_parallel

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
      PUBLIC :: updatc, calphi_bgrp
      PUBLIC :: mesure_diag_perf, mesure_mmul_perf
      PUBLIC :: use_parallel_diag
      PUBLIC :: bec_bgrp2ortho

CONTAINS

   SUBROUTINE mesure_diag_perf( n )
      !
      USE mp_global,   ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, root_bgrp
      USE mp_global,   ONLY: nproc_ortho, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id
      USE io_global,   ONLY: ionode, stdout
      USE mp,          ONLY: mp_sum, mp_bcast, mp_barrier
      USE mp,          ONLY: mp_max
      USE descriptors, ONLY: la_descriptor, descla_init
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      REAL(DP), ALLOCATABLE :: s(:,:), a(:,:), d(:)
      REAL(DP) :: t1, tpar, tser
      INTEGER  :: nr, nc, ir, ic, nx
      TYPE(la_descriptor) :: desc
      REAL(DP) :: cclock
      EXTERNAL :: cclock
      INTEGER, PARAMETER :: paradim = 1000
      !
      ! Check if number of PEs for orthogonalization/diagonalization is given from the input
      !
      IF( nproc_ortho > 0 ) THEN
         use_parallel_diag = .TRUE. 
         RETURN
      END IF

      ALLOCATE( d( n ) )
      !
      CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )

      nx = 1
      IF( desc%active_node > 0 ) nx = desc%nrcx

      nr = desc%nr
      nc = desc%nc
      ir = desc%ir
      ic = desc%ic

      ALLOCATE( s( nx, nx ) )
      ALLOCATE( a( nx, nx ) )
      !
      CALL set_a()
      !
      ! some MPIs (OpenMPI) the first time they call a collective routine take too much
      ! time to perform initializations, then perform a dummy call to get meaningful time
      !
      CALL diagonalize_parallel( n, a, d, s, desc )
      !
      CALL set_a()
      !
      CALL mp_barrier( intra_bgrp_comm )
      t1 = cclock()
      !
      CALL diagonalize_parallel( n, a, d, s, desc )
      !
      tpar = cclock() - t1
      CALL mp_max( tpar, intra_bgrp_comm )

      DEALLOCATE( s, a )
      !
      IF( desc%myc == 0 .AND. desc%myr == 0 .AND. desc%active_node > 0  .AND. n < paradim ) THEN

         ! when n >= paradim do not mesure serial perf, go parallel

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
         IF( desc%active_node < 0 ) RETURN
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
      USE mp_bands,    ONLY: nproc_bgrp, me_bgrp, intra_bgrp_comm, &
                             root_bgrp, my_bgrp_id, nbgrp
      USE mp_images,   ONLY: nimage, my_image_id
      USE mp_diag,     ONLY: ortho_comm, nproc_ortho, np_ortho, &
                             me_ortho, init_ortho_group, ortho_comm_id, ortho_cntx
      USE io_global,   ONLY: ionode, stdout
      USE mp,          ONLY: mp_sum, mp_bcast, mp_barrier
      USE mp,          ONLY: mp_max
      USE descriptors, ONLY: descla_init , la_descriptor
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: n
      !
      REAL(DP), ALLOCATABLE :: c(:,:), a(:,:), b(:,:)
      REAL(DP) :: t1, tcan
      INTEGER  :: nr, nc, ir, ic, np, lnode
      TYPE(la_descriptor) :: desc
      !
      REAL(DP) :: cclock
      EXTERNAL :: cclock
      !
      np    = MAX( INT( SQRT( DBLE( nproc_ortho ) + 0.1d0 ) ), 1 ) 
      !
      !  Make ortho group compatible with the number of electronic states
      !
      np    = MIN( np, n )
      !
      !  Now re-define the ortho group and test the performance
      !
      CALL init_ortho_group( np * np, intra_bgrp_comm, nimage*nbgrp, my_bgrp_id + nbgrp * my_image_id )

      CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )

      nr = desc%nr
      nc = desc%nc
   
      ALLOCATE( a( nr, nc ), c( nr, nc ), b( nr, nc ) )
   
      a = 1.0d0 / DBLE( n )
      b = 1.0d0 / DBLE( n )
   
      ! some MPIs (OpenMPI) the first time they call a collective routine take too much
      ! time to perform initializations, then perform a dummy call to get meaningful time
      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nr, b, nr, 0.0d0, c, nr, desc) 

      CALL mp_barrier( intra_bgrp_comm )
      t1 = cclock()

      CALL sqr_mm_cannon( 'N', 'N', n, 1.0d0, a, nr, b, nr, 0.0d0, c, nr, desc)
   
      tcan = cclock() - t1
      CALL mp_max( tcan, intra_bgrp_comm )

      DEALLOCATE( a, c, b )


#if defined __MPI

      IF( ionode ) THEN
         !
         WRITE( stdout, 90 )
         WRITE( stdout, 120 ) tcan, np*np
 90      FORMAT(/,3X,'Matrix Multiplication Performances')
120      FORMAT(3X,'ortho mmul, time for parallel driver      = ', 1F9.5, ' with ', I4, ' procs')
         !
      END IF

#else

      np = 1

#endif

#if defined __MPI
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
      USE mp_global,         ONLY: intra_bgrp_comm, me_bgrp, nproc_bgrp
      USE mp,                ONLY: mp_sum, mp_max
      USE descriptors,       ONLY: la_descriptor

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nss, ldx, nx0
      TYPE(la_descriptor), INTENT(IN) :: desc
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

      IF( desc%active_node < 0 ) then
         xloc = 0.0d0
         iter = 0
         go to 100
      endif
      !
      !  Compute the size of the local block
      !
      nr = desc%nr
      nc = desc%nc
      ir = desc%ir
      ic = desc%ic

      IF( ldx/= desc%nrcx ) &
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
         !       the following calls do the following matrix multiplications:
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
!$omp parallel do default(shared), private(j)
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

         CALL mp_max( diff, desc%comm )


         IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

         !
         !     the following calls do:
         !                       tmp1 = x1*u
         !                       tmp2 = ut*x1*u
         !
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, x1, ldx,  u,    ldx, 0.0d0, tmp1, ldx, desc )
         CALL sqr_mm_cannon( 'T', 'N', nss, 1.0d0, u,  ldx, tmp1, ldx,  0.0d0, tmp2, ldx, desc )
         !
         !       g=ut*x1*u/d  (g is stored in tmp1)
         !
!$omp parallel do default(shared), private(j)
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
         CALL sqr_mm_cannon( 'N', 'T', nss, 1.0d0, tmp1, ldx,  u,    ldx, 0.0d0, tmp2, ldx, desc )
         CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, u,    ldx, tmp2, ldx,  0.0d0, xloc, nx0, desc) 
         !
      END DO ITERATIVE_LOOP

      DEALLOCATE( tmp1, tmp2, dd, x1, con, tr1, tr2 )
            
100   CONTINUE
            
      CALL mp_max( iter, intra_bgrp_comm ) 

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
      USE mp_global,         ONLY: intra_bgrp_comm, me_bgrp, nproc_bgrp
      USE mp,                ONLY: mp_sum, mp_max
      USE descriptors,       ONLY: la_descriptor

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nss, ldx, nx0
      TYPE(la_descriptor), INTENT(IN) :: desc
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

      if( desc%active_node < 0 ) then
         xloc = 0.0d0
         iter = 0
         go to 100
      endif
      !
      !  Compute the size of the local block
      !
      nr = desc%nr
      nc = desc%nc
      ir = desc%ir
      ic = desc%ic

      IF( ldx/= desc%nrcx ) &
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
        IF( desc%myr == desc%myc ) THEN
           DO i = 1, nr
              SIGD( i + ir - 1 )   =  tmp1(i,i)
              tmp1(i,i) = -SIGD( i + ir - 1 )
           ENDDO
        END IF
        CALL mp_sum( sigd, desc%comm )

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

        CALL mp_max( diff, desc%comm )

        IF( diff < ortho_eps ) EXIT ITERATIVE_LOOP

      END DO ITERATIVE_LOOP
      !
      ! ...   Transform x0 back to the original basis

      CALL sqr_mm_cannon( 'N', 'N', nss, 1.0d0, u, ldx, xloc, nx0, 0.0d0, tmp1, ldx, desc)
      CALL sqr_mm_cannon( 'N', 'T', nss, 1.0d0, u, ldx, tmp1, ldx, 0.0d0, xloc, nx0, desc)

      DEALLOCATE( tmp1, tmp2, x1, sigd )

100   CONTINUE

      CALL mp_max( iter, intra_bgrp_comm ) 

      RETURN
   END SUBROUTINE ortho_alt_iterate


!-------------------------------------------------------------------------
   SUBROUTINE sigset( cp, ngwx, becp_dist, nkbx, qbecp, n, nss, ist, sig, ldx, desc )
!-----------------------------------------------------------------------
!     input: cp (non-orthonormal), becp, qbecp
!     computes the matrix
!       sig = 1 - a ,  a = <cp|s|cp> = <cp|cp> + sum q_ij <cp|i><j|cp>
!     where s=s(r(t+dt))
!     routine makes use of c(-q)=c*(q)
!
      USE kinds,              ONLY: DP
      USE uspp,               ONLY: nkbus
      USE uspp_param,         ONLY: nvb
      USE gvecw,              ONLY: ngw
      USE gvect, ONLY: gstart
      USE mp,                 ONLY: mp_root_sum, mp_sum
      USE control_flags,      ONLY: iverbosity
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_bgrp_comm, leg_ortho, inter_bgrp_comm, my_bgrp_id, nbgrp
      USE descriptors,        ONLY: la_descriptor, descla_init
      USE parallel_toolkit,   ONLY: dsqmsym
!
      IMPLICIT NONE
!
      INTEGER     :: nss, ist, ngwx, nkbx, n, ldx, nx
      COMPLEX(DP) :: cp( ngwx, n )
      REAL(DP)    :: qbecp( nkbx, ldx )
      REAL(DP)    :: becp_dist( nkbx, ldx )
      REAL(DP)    :: sig( ldx, ldx )
      TYPE(la_descriptor), INTENT(IN) :: desc
!
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root
      TYPE(la_descriptor):: desc_ip
      INTEGER :: np( 2 ), coor_ip( 2 )
      !
      REAL(DP), ALLOCATABLE :: sigp(:,:)
!
      IF( nss < 1 ) RETURN

      np(1) = desc%npr
      np(2) = desc%npc

      nx = desc%nrcx

      ALLOCATE( sigp( nx, nx ) ) 

      IF( desc%active_node > 0 ) THEN
         IF( desc%nrcx /= ldx ) &
            CALL errore( " sigset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " sigset ", " inconsistent dimension nx ", nx )
      END IF

      IF( nbgrp > 1 ) THEN
         sig = 0.0d0
      END IF

      DO ipc = 1, np(2)

         DO ipr = 1, ipc ! np(1) use symmetry

            coor_ip(1) = ipr - 1
            coor_ip(2) = ipc - 1

            CALL descla_init( desc_ip, desc%n, desc%nx, np, coor_ip, desc%comm, desc%cntx, 1 )

            nr = desc_ip%nr
            nc = desc_ip%nc
            ir = desc_ip%ir
            ic = desc_ip%ic
            !
            CALL GRID2D_RANK( 'R', desc_ip%npr, desc_ip%npc, &
                                   desc_ip%myr, desc_ip%myc, root )

            IF( MOD( root , nbgrp ) == my_bgrp_id ) THEN

               root = root * leg_ortho

               CALL dgemm( 'T', 'N',  nr, nc, 2*ngw, -2.0d0, cp( 1, ist + ir - 1), 2*ngwx, &
                           cp( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, sigp, nx )
               !
               !     q = 0  components has weight 1.0
               !
               IF ( gstart == 2 ) THEN
                  CALL DGER( nr, nc, 1.D0, cp(1,ist+ir-1), 2*ngwx, cp(1,ist+ic-1), 2*ngwx, sigp, nx )
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
      CALL dsqmsym( nss, sig, nx, desc )
      !
      IF( desc%active_node > 0 ) THEN
         !
         nr = desc%nr
         nc = desc%nc
         ir = desc%ir
         ic = desc%ic
         !
         IF( desc%myr == desc%myc ) THEN
            DO i = 1, nr
               sig(i,i) = sig(i,i) + 1.0d0
            END DO
         END IF
         !
         IF( nvb > 0 ) THEN
            CALL dgemm( 'T', 'N', nr, nc, nkbus, -1.0d0, becp_dist( 1, 1 ), &
                         nkbx, qbecp( 1, 1 ), nkbx, 1.0d0, sig, ldx )
         ENDIF
         !
         IF( iverbosity > 2 ) THEN
            WRITE( stdout,*)
            WRITE( stdout,'(26x,a)') '    sig '
            DO i = 1, nr
               WRITE( stdout,'(7f11.6)' ) ( sig(i,j), j=1, nc )
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
      USE gvect,              ONLY: gstart
      USE uspp,               ONLY: nkbus
      USE uspp_param,         ONLY: nvb
      USE kinds,              ONLY: DP
      USE mp,                 ONLY: mp_root_sum, mp_sum
      USE mp_global,          ONLY: intra_bgrp_comm, me_bgrp, leg_ortho
      USE mp_global,          ONLY: inter_bgrp_comm, my_bgrp_id, nbgrp
      USE control_flags,      ONLY: iverbosity
      USE io_global,          ONLY: stdout
      USE descriptors,        ONLY: la_descriptor, descla_init
!
      IMPLICIT NONE
!
      INTEGER     :: nss, ist, ngwx, nkbx, ldx, n
      COMPLEX(DP) :: cp( ngwx, n ), phi( ngwx, n )
      REAL(DP)    :: bephi( nkbx, ldx ), qbecp( nkbx, ldx )
      REAL(DP)    :: rho( ldx, ldx )
      TYPE(la_descriptor), INTENT(IN) :: desc
      !
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root, nx
      TYPE(la_descriptor) :: desc_ip
      INTEGER :: np( 2 ), coor_ip( 2 )

      REAL(DP), ALLOCATABLE :: rhop(:,:)
      !
      !     <phi|cp>
      !
      !

      IF( nss < 1 ) RETURN

      np(1) = desc%npr
      np(2) = desc%npc

      nx = desc%nrcx

      IF( desc%active_node > 0 ) THEN
         IF( desc%nrcx /= ldx ) &
            CALL errore( " rhoset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " rhoset ", " inconsistent dimension nx ", nx )
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

            CALL descla_init( desc_ip, desc%n, desc%nx, np, coor_ip, desc%comm, desc%cntx, 1 )

            nr = desc_ip%nr
            nc = desc_ip%nc
            ir = desc_ip%ir
            ic = desc_ip%ic
            !
            CALL GRID2D_RANK( 'R', desc_ip%npr, desc_ip%npc, &
                                   desc_ip%myr, desc_ip%myc, root )
            !
            IF( MOD( root , nbgrp ) == my_bgrp_id ) THEN

            root = root * leg_ortho

            CALL dgemm( 'T', 'N', nr, nc, 2*ngw, 2.0d0, phi( 1, ist + ir - 1 ), 2*ngwx, &
                  cp( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, rhop, nx )
            !
            !     q = 0  components has weight 1.0
            !
            IF (gstart == 2) THEN
               CALL DGER( nr, nc, -1.D0, phi(1,ist+ir-1), 2*ngwx, cp(1,ist+ic-1), 2*ngwx, rhop, nx )
            END IF

            CALL mp_root_sum( rhop, rho, root, intra_bgrp_comm )

            END IF

         END DO
      END DO
 
      DEALLOCATE( rhop )

      IF( nbgrp > 1 ) THEN
         CALL mp_sum( rho, inter_bgrp_comm )
      END IF

      IF( desc%active_node > 0 ) THEN
         !
         nr = desc%nr
         nc = desc%nc
         !
         !  bephi is distributed among processor rows
         !  qbephi is distributed among processor columns
         !  tau is block distributed among the whole processor 2D grid
         !
         !
         IF( nvb > 0 ) THEN
            !
            ! rho(i,j) = rho(i,j) + SUM_b bephi( b, i ) * qbecp( b, j ) 
            !
            CALL dgemm( 'T', 'N', nr, nc, nkbus, 1.0d0, bephi, nkbx, qbecp, nkbx, 1.0d0, rho, ldx )

         END IF

         IF ( iverbosity > 2 ) THEN
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
      USE uspp_param,         ONLY: nvb
      USE uspp,               ONLY: nkbus
      USE gvecw,              ONLY: ngw
      USE gvect,              ONLY: gstart
      USE mp,                 ONLY: mp_root_sum, mp_sum
      USE control_flags,      ONLY: iverbosity
      USE io_global,          ONLY: stdout
      USE mp_global,          ONLY: intra_bgrp_comm, leg_ortho
      USE mp_global,          ONLY: inter_bgrp_comm, my_bgrp_id, nbgrp
      USE descriptors,        ONLY: la_descriptor, descla_init
      USE parallel_toolkit,   ONLY: dsqmsym
!
      IMPLICIT NONE
      !
      INTEGER     :: nss, ist, ngwx, nkbx, n, ldx, nx
      COMPLEX(DP) :: phi( ngwx, n )
      REAL(DP)    :: bephi( nkbx, ldx ), qbephi( nkbx, ldx )
      REAL(DP)    :: tau( ldx, ldx )
      TYPE(la_descriptor), INTENT(IN) :: desc
      !
      INTEGER :: i, j, ipr, ipc, nr, nc, ir, ic, npr, npc
      INTEGER :: ii, jj, root
      TYPE(la_descriptor) :: desc_ip
      INTEGER :: np( 2 ), coor_ip( 2 )

      REAL(DP), ALLOCATABLE :: taup( :, : )
      !
      IF( nss < 1 ) RETURN
      !
      !  get dimensions of the square processor grid
      !
      np(1) = desc%npr
      np(2) = desc%npc
      !
      nx = desc%nrcx
      !
      IF( desc%active_node > 0 ) THEN
         IF( desc%nrcx /= ldx ) &
            CALL errore( " tauset ", " inconsistent dimension ldx ", ldx )
         IF( nx /= ldx ) &
            CALL errore( " tauset ", " inconsistent dimension nx ", nx )
      END IF
      !
      ALLOCATE( taup( nx, nx ) )
      !
      taup = 0.0d0
      !
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

            CALL descla_init( desc_ip, desc%n, desc%nx, np, coor_ip, desc%comm, desc%cntx, 1 )

            nr = desc_ip%nr
            nc = desc_ip%nc
            ir = desc_ip%ir
            ic = desc_ip%ic
            !
            CALL GRID2D_RANK( 'R', desc_ip%npr, desc_ip%npc, &
                                   desc_ip%myr, desc_ip%myc, root )
            !
            IF( MOD( root , nbgrp ) == my_bgrp_id ) THEN

            root = root * leg_ortho
            !
            !  All processors contribute to the tau block of processor (ipr,ipc)
            !  with their own part of wavefunctions
            !
            CALL dgemm( 'T', 'N', nr, nc, 2*ngw, 2.0d0, phi( 1, ist + ir - 1 ), 2*ngwx, &
                        phi( 1, ist + ic - 1 ), 2*ngwx, 0.0d0, taup, nx )
            !
            !           q = 0  components has weight 1.0
            !
            IF (gstart == 2) THEN
               CALL DGER( nr, nc, -1.D0, phi(1,ist+ir-1), 2*ngwx, phi(1,ist+ic-1), 2*ngwx, taup, nx )
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
      CALL dsqmsym( nss, tau, nx, desc )
      !
      IF( desc%active_node > 0 ) THEN
         !
         nr = desc%nr
         nc = desc%nc
         !
         !  bephi is distributed among processor rows
         !  qbephi is distributed among processor columns
         !  tau is block distributed among the whole processor 2D grid
         !
         IF( nvb > 0 ) THEN
            !
            CALL dgemm( 'T', 'N', nr, nc, nkbus, 1.0d0, bephi, nkbx, qbephi, nkbx, 1.0d0, tau, ldx )
            !
         END IF

         IF( iverbosity > 2 ) THEN
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
   SUBROUTINE updatc( ccc, x0, phi, bephi, becp_bgrp, bec_bgrp, cp_bgrp, desc )
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
      USE uspp_param,        ONLY: nh, nvb, ish
      USE gvecw,             ONLY: ngw
      USE control_flags,     ONLY: iverbosity
      USE mp,                ONLY: mp_sum, mp_bcast
      USE mp_global,         ONLY: intra_bgrp_comm, leg_ortho, me_bgrp, inter_bgrp_comm
      USE electrons_base,    ONLY: nbspx_bgrp, ibgrp_g2l, nbsp, nspin,  nupdwn, iupdwn, nbspx
      USE descriptors,       ONLY: descla_init, la_descriptor
!
      IMPLICIT NONE
!
      TYPE(la_descriptor), INTENT(IN) :: desc( : )
      COMPLEX(DP) :: cp_bgrp( :, : ), phi( :, : )
      REAL(DP), INTENT(IN) :: ccc
      REAL(DP)    :: bec_bgrp( :, : ), x0( :, :, : )
      REAL(DP)    :: bephi( :, : )
      REAL(DP)    :: becp_bgrp( :, : )

      ! local variables

      INTEGER :: i, j, ig, is, iv, ia, inl, nr, nc, ir, ic, nx0, ngwx, nkbx, iss, nrcx
      INTEGER :: ipr, ipc, root, i1, i2, nss, istart
      INTEGER :: ibgrp_i, ibgrp_i_first, nbgrp_i, i_first
      REAL(DP),    ALLOCATABLE :: wtemp(:,:) 
      REAL(DP),    ALLOCATABLE :: xd(:,:) 
      REAL(DP),    ALLOCATABLE :: bephi_tmp(:,:) 
      INTEGER :: np( 2 ), coor_ip( 2 )
      TYPE(la_descriptor) :: desc_ip

      DO iss = 1, nspin
         !
         !  size of the local block
         !
         nrcx = desc( iss )%nrcx
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
         IF( desc( iss )%active_node > 0 ) THEN
            IF( nx0 /= desc( iss )%nrcx ) &
               CALL errore( " updatc ", " inconsistent dimension nx0 ", nx0 )
         END IF
         !
         np(1) = desc( iss )%npr
         np(2) = desc( iss )%npc
         !
         CALL start_clock( 'updatc' )
   
         ALLOCATE( xd( nrcx, nrcx ) )
   
         IF( nvb > 0 )THEN
            DO i = 1, nss
               ibgrp_i = ibgrp_g2l( i + istart - 1 )
               IF( ibgrp_i > 0 ) THEN
                  DO inl = 1, nkbus
                     bec_bgrp( inl, ibgrp_i ) = becp_bgrp( inl, ibgrp_i )
                  END DO
               END IF
            END DO
            ALLOCATE( wtemp( nrcx, nkb ) )
            ALLOCATE( bephi_tmp( nkbx, nrcx ) )
         END IF
   
   
         DO ipc = 1, np(2)
            !
            IF( nvb > 0 )THEN
               ! 
               ! For the inner loop we need the block of bebhi( :, ic : ic + nc - 1 )
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
               IF( me_bgrp == root ) bephi_tmp = bephi(:,i1:i2)
               CALL mp_bcast( bephi_tmp, root, intra_bgrp_comm )
               !
            END IF
   
            DO ipr = 1, np(1)
               !
               ! Compute the descriptor of processor with coord: ( ipr-1, ipc-1 ), in the ortho group
               !
               coor_ip(1) = ipr - 1
               coor_ip(2) = ipc - 1
   
               CALL descla_init( desc_ip, desc( iss )%n, desc( iss )%nx, np, coor_ip, desc( iss )%comm, desc( iss )%cntx, 1 )
   
               nr = desc_ip%nr
               nc = desc_ip%nc
               ir = desc_ip%ir
               ic = desc_ip%ic
               !
               CALL GRID2D_RANK( 'R', desc_ip%npr, desc_ip%npc, &
                                      desc_ip%myr, desc_ip%myc, root )
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
   
               IF( desc( iss )%myr == ipr - 1 .AND. &
                   desc( iss )%myc == ipc - 1 .AND. &
                   desc( iss )%active_node > 0 ) THEN
                  xd = x0(:,:,iss) * ccc
               END IF
   
               CALL mp_bcast( xd, root, intra_bgrp_comm )
   
               CALL dgemm( 'N', 'N', 2*ngw, nbgrp_i, nr, 1.0d0, phi(1,istart+ir-1), 2*ngwx, &
                           xd(1,i_first), nrcx, 1.0d0, cp_bgrp(1,ibgrp_i_first), 2*ngwx )
   
               IF( nvb > 0 )THEN
   
                  !     updating of the <beta|c(n,g)>
                  !
                  !     bec of vanderbilt species are updated 
                  !
                  CALL dgemm( 'N', 'T', nr, nkbus, nc, 1.0d0, xd, nrcx, bephi_tmp, nkbx, 0.0d0, wtemp, nrcx )
                  !
                  ! here nr and ir are still valid, since they are the same for all procs in the same row
                  !
!$omp parallel do default(none) private(ibgrp_i,inl) shared(nr,ibgrp_g2l,istart,ir,nkbus,bec_bgrp,wtemp)
                  DO i = 1, nr
                     ibgrp_i = ibgrp_g2l( i + istart + ir - 2 )
                     IF( ibgrp_i > 0 ) THEN
                        DO inl = 1, nkbus
                           bec_bgrp( inl, ibgrp_i ) = bec_bgrp( inl, ibgrp_i ) + wtemp( i, inl ) 
                        END DO
                     END IF
                  END DO
!$omp end parallel do
                  !
               END IF
   
            END DO
            !    
         END DO
   
         IF( nvb > 0 )THEN
            DEALLOCATE( wtemp )
            DEALLOCATE( bephi_tmp )
         END IF
         !
         IF ( iverbosity > 1 ) THEN
            WRITE( stdout,*)
            DO is = 1, nvb
               IF( nvb > 1 ) THEN
                  WRITE( stdout,'(33x,a,i4)') ' updatc: bec (is)',is
                  WRITE( stdout,'(8f9.4)')                                       &
        &            ((bec_bgrp(ish(is)+(iv-1)*na(is)+1,i+istart-1),iv=1,nh(is)),i=1,nss)
               ELSE
                  DO ia=1,na(is)
                     WRITE( stdout,'(33x,a,i4)') ' updatc: bec (ia)',ia
                     WRITE( stdout,'(8f9.4)')                                    &
        &            ((bec_bgrp(ish(is)+(iv-1)*na(is)+ia,i+istart-1),iv=1,nh(is)),i=1,nss)
                  END DO
               END IF
               WRITE( stdout,*)
            END DO
         ENDIF
         !
         DEALLOCATE( xd )
         !
      END DO
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
      USE ions_base,      ONLY: na, nsp
      USE io_global,      ONLY: stdout
      USE mp_global,      ONLY: intra_bgrp_comm, inter_bgrp_comm
      USE uspp_param,     ONLY: nh, ish, nvb
      USE uspp,           ONLY: nkbus, qq
      USE gvecw,          ONLY: ngw
      USE electrons_base, ONLY: nbsp_bgrp, nbsp
      USE constants,      ONLY: pi, fpi
      USE control_flags,  ONLY: iverbosity
      USE mp,             ONLY: mp_sum
!
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ngwx, nkbx, nbspx_bgrp
      COMPLEX(DP)         :: c0_bgrp( ngwx, nbspx_bgrp ), phi_bgrp( ngwx, nbspx_bgrp ), betae( ngwx, nkbx )
      REAL(DP)            :: bec_bgrp( nkbx, nbspx_bgrp ), emtot
      REAL(DP), OPTIONAL  :: ema0bg( ngwx )

      ! local variables
      !
      INTEGER  :: is, iv, jv, ia, inl, jnl, i, j
      REAL(DP), ALLOCATABLE :: qtemp( : , : )
      REAL(DP) :: qqf
!
      IF( nbsp_bgrp < 1 ) RETURN
      !
      CALL start_clock( 'calphi' )
      !
      ! Note that phi here is computed only for my band group
      !
      IF ( nvb > 0 ) THEN

         ALLOCATE( qtemp( nkbus, nbspx_bgrp ) )

         qtemp (:,:) = 0.d0
         DO is=1,nvb
            DO iv=1,nh(is)
               inl = ish(is)+(iv-1)*na(is)
               DO jv=1,nh(is)
                  jnl = ish(is)+(jv-1)*na(is)
                  IF(ABS(qq(iv,jv,is)) > 1.d-5) THEN
                     qqf = qq(iv,jv,is)
                     DO i=1,nbsp_bgrp
                        CALL daxpy( na(is), qqf, bec_bgrp(jnl+1,i),1,qtemp(inl+1,i), 1 )
                     END DO
                  ENDIF
               END DO
            END DO
         END DO
!
         CALL dgemm ( 'N', 'N', 2*ngw, nbsp_bgrp, nkbus, 1.0d0, betae, &
                       2*ngwx, qtemp, nkbus, 0.0d0, phi_bgrp, 2*ngwx )

         DEALLOCATE( qtemp )

      ELSE

         phi_bgrp = (0.d0, 0.d0)

      END IF
!
      IF( PRESENT( ema0bg ) ) THEN
!$omp parallel do default(shared), private(i)
         DO j=1,nbsp_bgrp
            DO i=1,ngw
               phi_bgrp(i,j)=(phi_bgrp(i,j)+c0_bgrp(i,j))*ema0bg(i)
            END DO
         END DO
!$omp end parallel do
      ELSE
!$omp parallel do default(shared), private(i)
         DO j=1,nbsp_bgrp
            DO i=1,ngw
               phi_bgrp(i,j)=phi_bgrp(i,j)+c0_bgrp(i,j)
            END DO
         END DO
!$omp end parallel do
      END IF

      !   

      IF(iverbosity > 1) THEN
         emtot=0.0d0
         IF( PRESENT( ema0bg ) ) THEN
            DO j=1,nbsp_bgrp
               DO i=1,ngw
                  emtot=emtot +2.0d0*DBLE(phi_bgrp(i,j)*CONJG(c0_bgrp(i,j)))*ema0bg(i)**(-2.0d0)
               END DO
            END DO
         ELSE
            DO j=1,nbsp_bgrp
               DO i=1,ngw
                  emtot=emtot +2.0d0*DBLE(phi_bgrp(i,j)*CONJG(c0_bgrp(i,j)))
               END DO
            END DO
         END IF
         emtot=emtot/nbsp

         CALL mp_sum( emtot, intra_bgrp_comm )
         CALL mp_sum( emtot, inter_bgrp_comm )

         WRITE( stdout,*) 'in calphi sqrt(emtot)=',SQRT(emtot)
         WRITE( stdout,*)
         DO is = 1, nvb
            IF( nvb > 1 ) THEN
               WRITE( stdout,'(33x,a,i4)') ' calphi: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &            ((bec_bgrp(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,nbsp_bgrp)
            ELSE
               DO ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calphi: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &               ((bec_bgrp(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,nbsp_bgrp)
               END DO
            END IF
         END DO
      ENDIF


      CALL stop_clock( 'calphi' )
!
      RETURN
      END SUBROUTINE calphi_bgrp


   SUBROUTINE bec_bgrp2ortho( bec_bgrp, bec_ortho, nrcx, desc )
      USE kinds,             ONLY: DP
      USE uspp,              ONLY: nkb, nkbus
      USE mp,                ONLY: mp_sum
      USE mp_global,         ONLY: intra_bgrp_comm, leg_ortho, me_bgrp, inter_bgrp_comm
      USE electrons_base,    ONLY: nbspx_bgrp, ibgrp_g2l, nspin
      USE descriptors,       ONLY: la_descriptor
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nrcx
      TYPE(la_descriptor), INTENT(IN) :: desc( : )
      REAL(DP), INTENT(IN)  :: bec_bgrp(:,:)
      REAL(DP), INTENT(OUT) :: bec_ortho(:,:)
      !
      INTEGER :: ir, nr, i, ibgrp_i, nup
      !
      bec_ortho = 0.0d0
      !
      IF( desc( 1 )%active_node > 0 ) THEN
         ir = desc( 1 )%ir
         nr = desc( 1 )%nr
         do i = 1, nr
            ibgrp_i = ibgrp_g2l( i + ir - 1 )
            IF( ibgrp_i > 0 ) THEN
               bec_ortho( :, i ) = bec_bgrp( :, ibgrp_i )
            END IF
         end do
      END IF
      !
      IF( nspin == 2 ) THEN
         IF( desc( 2 )%active_node > 0 ) THEN
            nup = desc( 1 )%n
            ir = desc( 2 )%ir
            nr = desc( 2 )%nr
            do i = 1, nr
               ibgrp_i = ibgrp_g2l( i + ir - 1 + nup )
               IF( ibgrp_i > 0 ) THEN
                  bec_ortho( :, i + nrcx ) = bec_bgrp( :, ibgrp_i )
               END IF
            end do
         END IF
      END IF
      !
      CALL mp_sum( bec_ortho, inter_bgrp_comm )
      !
      RETURN
   END SUBROUTINE bec_bgrp2ortho

END MODULE orthogonalize_base
