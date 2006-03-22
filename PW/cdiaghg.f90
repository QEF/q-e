!
! Copyright (C) 2001-2204 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix .
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both ZHEGV and ZHEGVX
  !
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : npool, me_pool, root_pool, intra_pool_comm
  USE mp,        ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculate
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  ! ... on OUTPUT
  !
  REAL(DP) :: e(n)
    ! eigenvalues
  COMPLEX(DP) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  ! ... LOCAL variables
  !
  INTEGER :: lwork, nb, mm, info
    ! mm = number of calculated eigenvectors
  INTEGER, EXTERNAL :: ILAENV
    ! ILAENV returns optimal block size "nb"
  INTEGER,          ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:)
  COMPLEX(DP), ALLOCATABLE :: sdum(:,:), hdum(:,:),  work(:)
  LOGICAL :: all_eigenvalues
  !
  !
  CALL start_clock( 'cdiaghg' )
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
     lwork = 2 * n
     !
  ELSE
     !
     lwork = ( nb + 1 ) * n
     !
  END IF
  !
  ! ... allocate workspace
  !
  ALLOCATE( work( lwork ) )    
  ALLOCATE( sdum( ldh, n ) )    
  !
  IF ( all_eigenvalues ) THEN
     !
     ALLOCATE( rwork( 3 * n - 2  ) )
     !
  ELSE
     !
     ALLOCATE( rwork( 7 * n ) )    
     ALLOCATE( hdum( ldh, n ) )    
     ALLOCATE( iwork( 5 * n ) )    
     ALLOCATE( ifail( n ) )    
     !
  END IF
  !
  ! input s and (see below) h are copied so that they are not destroyed
  !
  sdum = s
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me_pool == root_pool ) THEN
     !
     IF ( all_eigenvalues ) THEN
        !
        ! ... calculate all eigenvalues
        !
        v(:,1:n) = h(:,:)
        !
        CALL ZHEGV( 1, 'V', 'U', n, v, ldh, sdum, ldh, e, work, &
                    lwork, rwork, info )
        !
     ELSE
        !
        ! ... calculate only m lowest eigenvalues
        !
        hdum = h
        !
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, hdum, ldh, sdum, ldh,  &
                     0.0D0, 0.0D0, 1, m, 0.D0, mm, e(1), v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
        !
     END IF
     !
     CALL errore( 'cdiaghg', 'info =/= 0', ABS( info ) )
     !
  END IF
  !
  ! ... broadcast the eigenvectors and the eigenvalues
  !
  CALL mp_bcast( e, root_pool, intra_pool_comm )
  CALL mp_bcast( v, root_pool, intra_pool_comm )
  !
  ! ... deallocate workspace
  !
  IF ( .NOT. all_eigenvalues ) THEN
     !
     DEALLOCATE( ifail )
     DEALLOCATE( iwork )
     DEALLOCATE( hdum )
     !
  END IF
  !
  DEALLOCATE( sdum )
  DEALLOCATE( rwork )
  DEALLOCATE( work )
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE cdiaghg
