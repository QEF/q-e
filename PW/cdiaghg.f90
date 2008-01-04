!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both ZHEGV and ZHEGVX
  !
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast, mp_sum
  USE mp_global,        ONLY : me_pool, root_pool, intra_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculate
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! actually intent(in) but compilers don't know and complain
    ! matrix to be diagonalized
    ! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER                  :: lwork, nb, mm, info, i, j, k
    ! mm = number of calculated eigenvectors
  REAL(DP)                 :: abstol
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space
  LOGICAL                  :: all_eigenvalues
 ! REAL(DP), EXTERNAL       :: DLAMCH
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  !
  CALL start_clock( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_pool == root_pool ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = DBLE( s(i,i) )
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 1 ) nb = MAX( 1, n )
     !
     IF ( nb == 1 .OR. nb >= n ) THEN
        !
        lwork = 2*n
        !
     ELSE
        !
        lwork = ( nb + 1 )*n
        !
     END IF
     !
     ALLOCATE( work( lwork ) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ALLOCATE( rwork( 3*n - 2 ) )
        !
        ! ... calculate all eigenvalues (overwritten to v)
        !
        v(:,:) = h(:,:)
        !
        CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                    s, ldh, e, work, lwork, rwork, info )
        !
     ELSE
        !
        ALLOCATE( rwork( 7*n ) )
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
        ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = DBLE( h(i,i) )
        END DO
        !
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
        !
        ! ... calculate only m lowest eigenvalues
        !
        abstol = 0.D0
       ! abstol = 2.D0*DLAMCH( 'S' )
        !
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        DO i = 1, n
           h(i,i) = CMPLX( hdiag(i), 0.0_DP )
           DO j = i + 1, n
              h(i,j) = CONJG( h(j,i) )
           END DO
           DO j = n + 1, ldh
              h(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
        END DO
        !
        DEALLOCATE( hdiag )
        !
     END IF
     !
     DEALLOCATE( rwork )
     DEALLOCATE( work )
     !
     CALL errore( 'cdiaghg', 'info =/= 0', ABS( info ) )
     !
     ! ... restore input S matrix from saved diagonal and lower triangle
     !
     DO i = 1, n
        s(i,i) = CMPLX( sdiag(i), 0.0_DP )
        DO j = i + 1, n
           s(i,j) = CONJG( s(j,i) )
        END DO
        DO j = n + 1, ldh
           s(j,i) = ( 0.0_DP, 0.0_DP )
        END DO
     END DO
     !
     DEALLOCATE( sdiag )
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
  CALL mp_bcast( e, root_pool, intra_pool_comm )
  CALL mp_bcast( v, root_pool, intra_pool_comm )
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE cdiaghg
!
!----------------------------------------------------------------------------
SUBROUTINE pcdiaghg( n, h, s, ldh, e, v, desc )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... Parallel version, with full data distribution
  !
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : root_pool, intra_pool_comm
  USE zhpev_module,     ONLY : pzhpev_drv, zhpev_drv
  USE descriptors,      ONLY : descla_siz_ , lambda_node_ , nlax_ , la_nrl_ , la_nrlx_ , &
                               la_npc_ , la_npr_ , la_me_ , la_comm_ , la_myc_ , la_myr_ , &
                               nlar_
  USE parallel_toolkit, ONLY : zsqmdst, zsqmcll
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, ldh
    ! dimension of the matrix to be diagonalized
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh), s(ldh,ldh)
    ! actually intent(in) but compilers don't know and complain
    ! matrix to be diagonalized
    ! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
    ! eigenvectors (column-wise)
  INTEGER, INTENT(IN) :: desc( descla_siz_ )
  !
  INTEGER             :: nx, nrl, nrlx
    ! local block size
  COMPLEX(DP), ALLOCATABLE :: ss(:,:), hh(:,:)
  COMPLEX(DP), ALLOCATABLE :: diag(:,:), vv(:,:)
    ! work space used only in parallel diagonalization
  INTEGER :: i, j, k
  !
  ! ... input s and h are copied so that they are not destroyed
  !
  CALL start_clock( 'cdiaghg' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     nx   = desc( nlax_ )
     nrl  = desc( la_nrl_ )
     nrlx = desc( la_nrlx_ )
     !
     IF( nx /= ldh ) &
        CALL errore(" pcdiaghg ", " inconsistent leading dimension ", ldh )
     !
     ALLOCATE( hh( nx, nx ) )
     ALLOCATE( ss( nx, nx ) )
     !
     hh(1:nx,1:nx) = h(1:nx,1:nx)
     ss(1:nx,1:nx) = s(1:nx,1:nx)
     !
  END IF

  CALL start_clock( 'cdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of sl ( L is stored in sl )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL pzpotf( ss, nx, n, desc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:choldc' )
  !
  ! ... L is inverted ( sl = L^-1 )
  !
  CALL start_clock( 'cdiaghg:inversion' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL pztrtri( ss, nx, n, desc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:inversion' )
  !
  ! ... vl = L^-1*H
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'N', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, desc )
     !
     !
  END IF
  !
  ! ... hl = ( L^-1*H )*(L^-1)^T
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'N', 'C', n, ONE, v, nx, ss, nx, ZERO, hh, nx, desc )
     !
     IF( desc( la_myc_ ) == desc( la_myr_ ) ) THEN
        !
        ! ensure that "hh" is really Hermitian, it is sufficient to set the diagonal
        ! properly, because only the lower triangle of hh will be used
        ! 
        DO j = 1, desc( nlar_ )
           hh( j, j ) = CMPLX( REAL( hh(j,j) ), 0.0_DP )
        END DO
     END IF
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     ! 
     !  Compute local dimension of the cyclically distributed matrix
     !
#ifndef TEST_DIAG
     ALLOCATE( diag( nrlx, n ) )
     ALLOCATE( vv( nrlx, n ) )
     !
     CALL blk2cyc_zredist( n, diag, nrlx, hh, nx, desc )
     !
     CALL pzhpev_drv( 'V', diag, nrlx, e, vv, nrlx, nrl, n, &
          desc( la_npc_ ) * desc( la_npr_ ), desc( la_me_ ), desc( la_comm_ ) )
     !
     CALL cyc2blk_zredist( n, vv, nrlx, hh, nx, desc )
     !
     DEALLOCATE( vv )
     DEALLOCATE( diag )
#else
     ALLOCATE( vv( n, n ) )
     CALL zsqmcll( n, hh, nx, vv, n, desc, desc( la_comm_ ) )
     IF( desc( la_myc_ ) == 0 .AND. desc( la_myr_ ) == 0 ) THEN
        ALLOCATE( diag( n*(n+1)/2, 1 ) )
        k = 1
        write( 100, fmt="(I5)" ) n
        DO j = 1, n
           DO i = j, n
              diag( k, 1 ) = vv( i, j )
              write( 100, fmt="(2I5,2D18.10)" ) i, j, vv( i, j )
              k = k + 1
           END DO
        END DO
        call zhpev_drv( 'V', 'L', N, diag(:,1), e, vv, n )
        write( 100, * ) 'eigenvalues and eigenvectors'
        DO j = 1, n
           write( 100, fmt="(1I5,1D18.10,A)" ) j, e( j ), 'eval'
           DO i = 1, n
              write( 100, fmt="(2I5,2D18.10)" ) i, j, vv( i, j )
           END DO
        END DO
        close(100)
        DEALLOCATE( diag )
     END IF
     CALL mp_bcast( vv, 0, desc( la_comm_ ) )
     CALL zsqmdst( n, vv, n, hh, nx, desc )
     DEALLOCATE( vv )
     CALL errore('cdiaghg','stop serial',1)
#endif
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'C', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, desc )
     !
  END IF
  !
  CALL mp_bcast( e, root_pool, intra_pool_comm )
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     DEALLOCATE( ss, hh )
  END IF
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE pcdiaghg
!
!----------------------------------------------------------------------------
SUBROUTINE pcdiaghg_nodist( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... Parallel version, matrices are NOT distributed
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : use_para_diag
  USE mp,               ONLY : mp_bcast, mp_sum
  USE mp_global,        ONLY : npool, nproc_pool, me_pool, root_pool, &
                               intra_pool_comm, init_ortho_group, &
                               ortho_comm, np_ortho, me_ortho, ortho_comm_id
  USE zhpev_module,     ONLY : pzhpev_drv
  USE descriptors,      ONLY : descla_siz_ , descla_init , lambda_node_ , &
                               nlax_ , la_nrl_ , ilac_ , ilar_ , nlar_ ,  &
                               nlac_ , la_npc_ , la_npr_ , la_me_ , la_comm_
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculate
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! actually intent(in) but compilers don't know and complain
    ! matrix to be diagonalized
    ! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER                  :: lwork, nb, mm, info, i, j, k
    ! mm = number of calculated eigenvectors
  INTEGER                  :: nr, nc, ir, ic, nx, nrl
    ! local block size
  REAL(DP)                 :: abstol
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space
  COMPLEX(DP), ALLOCATABLE :: sl(:,:), hl(:,:), vl(:,:)
  COMPLEX(DP), ALLOCATABLE :: diag(:,:), vv(:,:)
    ! work space used only in parallel diagonalization
  LOGICAL                  :: all_eigenvalues
 ! REAL(DP), EXTERNAL       :: DLAMCH
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  INTEGER                  :: desc( descla_siz_ )
  !
  !
  CALL start_clock( 'cdiaghg' )
  !
  CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
  !
  ! ... input s and h are copied so that they are not destroyed
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     ir  = desc( ilar_ )
     ic  = desc( ilac_ )
     nr  = desc( nlar_ )
     nc  = desc( nlac_ )
     nx  = desc( nlax_ )
     nrl = desc( la_nrl_ )
     ALLOCATE( sl( nx , nx ) )
     ALLOCATE( vl( nx , nx ) )
     ALLOCATE( hl( nx , nx ) )
     DO j = 1, nc
        DO i = 1, nr
           sl( i, j ) = s( i + ir - 1, j + ic - 1 )
        END DO
        DO i = nr+1, nx
           sl( i, j ) = 0.0d0
        END DO
     END DO
     DO j = nc + 1, nx
        DO i = 1, nx
           sl( i, j ) = 0.0d0
        END DO
     END DO
     DO j = 1, nc
        DO i = 1, nr
           hl( i, j ) = h( i + ir - 1, j + ic - 1 )
        END DO
        DO i = nr+1, nx
           hl( i, j ) = 0.0d0
        END DO
     END DO
     DO j = nc + 1, nx
        DO i = 1, nx
           hl( i, j ) = 0.0d0
        END DO
     END DO
  END IF

  CALL start_clock( 'cdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of sl ( L is stored in sl )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL pzpotf( sl, nx, n, desc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:choldc' )
  !
  ! ... L is inverted ( sl = L^-1 )
  !
  CALL start_clock( 'cdiaghg:inversion' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL pztrtri( sl, nx, n, desc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:inversion' )
  !
  ! ... vl = L^-1*H
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'N', 'N', n, ONE, sl, nx, hl, nx, ZERO, vl, nx, desc )
     !
  END IF
  !
  ! ... hl = ( L^-1*H )*(L^-1)^T
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'N', 'C', n, ONE, vl, nx, sl, nx, ZERO, hl, nx, desc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     ! 
     !  Compute local dimension of the cyclically distributed matrix
     !
     ALLOCATE( diag( nrl, n ) )
     ALLOCATE( vv( nrl, n ) )
     !
     CALL blk2cyc_zredist( n, diag, nrl, hl, nx, desc )
     !
     CALL pzhpev_drv( 'V', diag, nrl, e, vv, nrl, nrl, n, &
          desc( la_npc_ ) * desc( la_npr_ ), desc( la_me_ ), desc( la_comm_ ) )
     !
     CALL cyc2blk_zredist( n, vv, nrl, vl, nx, desc )
     !
     DEALLOCATE( vv )
     DEALLOCATE( diag )
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  v(1:n,1:n) = ZERO
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'C', 'N', n, ONE, sl, nx, vl, nx, ZERO, hl, nx, desc )
     !
     DO j = 1, nc
        DO i = 1, nr
           v( i + ir - 1, j + ic - 1 ) = hl( i, j )
        END DO
     END DO
     !
  END IF
  !
  CALL mp_bcast( e, root_pool, intra_pool_comm )
  CALL mp_sum( v(1:n,1:n), intra_pool_comm )
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc( lambda_node_ ) > 0 ) THEN
     DEALLOCATE( sl, vl, hl )
  END IF
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE pcdiaghg_nodist
