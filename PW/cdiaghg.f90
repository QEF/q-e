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
  ! ... the parallel version is also implemented
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : use_para_diago, para_diago_dim
  USE parallel_toolkit, ONLY : cdiagonalize
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : npool, nproc_pool, me_pool, &
                               root_pool, intra_pool_comm
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
  INTEGER                  :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  INTEGER, EXTERNAL        :: ILAENV
    ! ILAENV returns optimal block size "nb"
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space
  COMPLEX(DP), ALLOCATABLE :: sdum(:,:), hdum(:,:), vdum(:,:)
    ! work space used only in parallel diagonalization
  LOGICAL                  :: all_eigenvalues
  !
  !
  CALL start_clock( 'diaghg' )
  !
  !
  IF ( use_para_diago .AND. n > para_diago_dim ) THEN
     !
     ALLOCATE( hdum( n, n ), sdum( n, n ), vdum( n, n) )
     !
     ! ... input s and h are copied so that they are not destroyed
     !
     sdum(:,:) = s(1:n,:)
     hdum(:,:) = h(1:n,:)
     !
     CALL start_clock( 'choldc' )
     !
     ! ... Cholesky decomposition of sdum ( L is stored in sdum )
     !
     CALL para_zcholdc( n, sdum, n, intra_pool_comm )
     !
     CALL stop_clock( 'choldc' )
     !
     ! ... L is inverted ( sdum = L^-1 )
     !
     CALL start_clock( 'inversion' )
     !
     CALL para_ztrtri( n, sdum, n, intra_pool_comm )
     !
     CALL stop_clock( 'inversion' )
     !
     ! ... vdum = L^-1*H
     !
     CALL start_clock( 'paragemm' )
     !
     CALL para_zgemm( 'N', 'N', n, n, n, ONE, sdum, &
                      n, hdum, n, ZERO, vdum, n, intra_pool_comm )
     !
     ! ... hdum = ( L^-1*H )*(L^-1)^T
     !
     CALL para_zgemm( 'N', 'C', n, n, n, ONE, vdum, n, &
                      sdum, n, ZERO, hdum, n, intra_pool_comm )
     !
     CALL stop_clock( 'paragemm' )
     !
     CALL cdiagonalize( 1, hdum, n, e, vdum, n, n, &
                        nproc_pool, me_pool, intra_pool_comm )
     !
     ! ... v = (L^T)^-1 v
     !
     CALL start_clock( 'paragemm' )
     !
     CALL para_zgemm( 'C', 'N', n, m, n, ONE, sdum, n, &
                      vdum, n, ZERO, v, ldh, intra_pool_comm )
     !
     CALL stop_clock( 'paragemm' )
     !
     DEALLOCATE( vdum, sdum, hdum )
     !
  ELSE
     !
     ! ... only the first processor diagonalizes the matrix
     !
     IF ( me_pool == root_pool ) THEN
        !
        ALLOCATE( sdiag( n ) )
        !
        ! ... save the diagonal of input S (it will be overwritten)
        !
        DO i = 1, n
           sdiag(i) = DBLE ( s(i,i) )
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
              hdiag(i) = DBLE ( h(i,i) )
           END DO
           !
           ALLOCATE( iwork( 5*n ) )
           ALLOCATE( ifail( n ) )
           !
           ! ... calculate only m lowest eigenvalues
           !
           CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                        0.D0, 0.D0, 1, m, 0.D0, mm, e, v, ldh, &
                        work, lwork, rwork, iwork, ifail, info )
           !
           DEALLOCATE( ifail )
           DEALLOCATE( iwork )
           !
           ! ... restore input H matrix from saved diagonal and lower triangle 
           !
           DO i = 1, n
              h(i,i) = CMPLX ( hdiag (i), 0.0_DP )
              DO j = i + 1, n
                 h(i,j) = CONJG ( h(j,i) )
              END DO
              DO j = n + 1, ldh
                 h(i,j) = ( 0.0_DP, 0.0_DP )
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
           s(i,i) = CMPLX ( sdiag (i), 0.0_DP )
           DO j = i + 1, n
              s(i,j) = CONJG ( s(j,i) )
           END DO
           DO j = n + 1, ldh
              s(i,j) = ( 0.0_DP, 0.0_DP )
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
  END IF
  !
  CALL stop_clock( 'diaghg' )
  !
  RETURN
  !
END SUBROUTINE cdiaghg
