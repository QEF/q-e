!
! Copyright (C) 2001-2204 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix .
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - may use both ZHEGV and ZHEGVX
  ! ... ZHEGVX should be faster but it is not available on many machines
  ! ... define HAS_ZHEGVX in the preprocessing options to use ZHEGVX
  !
  USE kinds,     ONLY : DP
#if defined (__PARA)
  USE para,      ONLY : me, npool, MPI_COMM_POOL
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_bcast  
#endif
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculate
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(KIND=DP) :: h(ldh,n), s(ldh,n)
    ! matrix to be diagonalized
    ! overlap matrix
  !
  ! ... on OUTPUT
  !
  REAL(KIND=DP) :: e(n)
    ! eigenvalues
  COMPLEX(KIND=DP) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  ! ... LOCAL variables
  !
  INTEGER :: lwork, nb, ILAENV, mm, info
    ! ILAENV returns optimal block size "nb"
    ! mm = number of calculated eigenvectors
  INTEGER,          ALLOCATABLE :: iwork(:), ifail(:)
  REAL(KIND=DP),    ALLOCATABLE :: rwork(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: sdum(:,:), hdum(:,:),  work(:)
  LOGICAL :: all_eigenvalues
  !
  ! ... external methods
  !
  EXTERNAL ZCOPY, ZHEGV, ILAENV
  !
  !
  CALL start_clock( 'cdiaghg' )
  !
#if defined (__PARA)
#if defined (__T3E)
  !
  ! ... NB: 150 has been determined empirically on the T3E as the point
  ! ...     where it is convenient to use a parallel routines.
  !
  IF ( npool == 1 .AND. n > 150 ) THEN
     !
     CALL scala_cdiaghg( n, h, ldh, s, ldh, e, v, ldh )
     !
     CALL stop_clock( 'cdiaghg' )
     !
     RETURN
     !
  END IF
  !
#endif
#endif
#if defined (HAS_ZHEGVX)
  !
  all_eigenvalues = ( m == n )
  !
#else
  !
  all_eigenvalues = .TRUE.
  !
#endif
  !
  ! ... check for optimal block size
  !
  nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
  !
  IF ( nb < 1 ) nb = MAX( 1, n )
  !
  IF ( nb == 1 .OR. nb >= n ) THEN
     !
     lwork = 2 * n - 1
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
#if defined (__PARA)
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me == 1 ) THEN
     !
#endif
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
#if defined (HAS_ZHEGVX)
        !
        ! ... calculate only m lowest eigenvalues
        !
        hdum = h
        !
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, hdum, ldh, sdum, ldh,  &
                     0.0D0, 0.0D0, 1, m, 0.D0, mm, e(1), v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
        !
#endif
     END IF
     !
     CALL errore( 'cdiaghg', 'info =/= 0', ABS( info ) )
     !
#if defined (__PARA)
     !
  END IF
  !
  ! ... broadcast the eigenvectors and the eigenvalues
  !
  CALL mp_bcast( e, ionode_id, MPI_COMM_POOL )
  CALL mp_bcast( v, ionode_id, MPI_COMM_POOL )
  !
#endif
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
