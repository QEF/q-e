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
  COMPLEX(DP), INTENT(IN) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  INTEGER                  :: i,j, lwork, nb, mm, info
    ! mm = number of calculated eigenvectors
  INTEGER, EXTERNAL        :: ILAENV
    ! ILAENV returns optimal block size "nb"
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:)
  COMPLEX(DP), ALLOCATABLE :: sdum(:,:), hdum(:,:), vdum(:,:), work(:)
  LOGICAL                  :: all_eigenvalues
  !
  !
  CALL start_clock( 'diaghg' )
  !
  ALLOCATE( sdum( n, n ) )
  !
  ! ... input s and (see below) h are copied so that they are not destroyed
  !
  sdum(:,:) = s(1:n,:)
  !
  IF ( use_para_diago .AND. n > para_diago_dim ) THEN
     !
     PRINT *, me_pool, "USING PARALLEL cdiaghg"
     !
     ALLOCATE( hdum( n, n ), vdum( n, n ) )
     !
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
     DEALLOCATE( hdum, vdum )
     !
  ELSE
     !
     ! ... only the first processor diagonalize the matrix
     !
     IF ( me_pool == root_pool ) THEN
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
           ! ... calculate all eigenvalues
           !
           v(:,:) = h(:,:)
           !
           CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                       sdum, n, e, work, lwork, rwork, info )
           !
        ELSE
           !
           ALLOCATE( rwork( 7*n ) )  
           ALLOCATE( hdum( n, n ) )
           ALLOCATE( iwork( 5*n ) )
           ALLOCATE( ifail( n ) )
           !
           ! ... calculate only m lowest eigenvalues
           !
           hdum(:,:) = h(1:n,:)
           !
           CALL ZHEGVX( 1, 'V', 'I', 'U', n, hdum, n, sdum, n, &
                        0.D0, 0.D0, 1, m, 0.D0, mm, e, v, ldh, &
                        work, lwork, rwork, iwork, ifail, info )
           !
           DEALLOCATE( ifail )
           DEALLOCATE( iwork )
           DEALLOCATE( hdum )
           !
        END IF
        !
        DEALLOCATE( rwork )
        DEALLOCATE( work )
        !
        CALL errore( 'cdiaghg', 'info =/= 0', ABS( info ) )
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
  DEALLOCATE( sdum )
  !
  CALL stop_clock( 'diaghg' )
  !
  RETURN
  !
END SUBROUTINE cdiaghg
