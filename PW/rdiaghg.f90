!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE rdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix .
  ! ... On output both matrix are unchanged
  ! ... Uses LAPACK routines
  !
  USE kinds,     ONLY : DP
  USE para,      ONLY : me
  USE mp,        ONLY : mp_bcast
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : intra_pool_comm
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(KIND=DP) :: h(ldh,n)
    ! matrix to be diagonalized
  REAL(KIND=DP) :: s(ldh,n)
    ! overlap matrix
  !
  ! ... on OUTPUT
  !
  REAL(KIND=DP) :: e(n)
    ! eigenvalues
  REAL(KIND=DP) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  ! ... LOCAL variables
  !
  INTEGER                    :: lwork, nb, ILAENV, mm, info
    ! ILAENV returns optimal block size "nb"
    ! mm = number of calculated eigenvectors
  EXTERNAL                      ILAENV
  INTEGER, ALLOCATABLE       :: iwork(:), ifail(:)
  REAL(KIND=DP), ALLOCATABLE :: sdum(:,:), hdum(:,:),  work(:)
  LOGICAL                    :: all_eigenvalues
  !
  !
  CALL start_clock( 'cdiaghg' )
  !
  all_eigenvalues = ( m == n )
  !
  ! ... check for optimal block size
  !
  nb = ILAENV( 1, 'DSYTRD', 'U', n, -1, -1, -1 )
  !
  IF ( nb < 1 .OR. nb >= n ) THEN
     !
     lwork = 8 * n
     !
  ELSE
     !
     lwork = ( nb + 3 ) * n
     !
  END IF
  !
  ! ... allocate workspace
  !
  ALLOCATE( work( lwork ) )    
  ALLOCATE( sdum( ldh, n ) )
  !    
  IF ( .NOT. all_eigenvalues ) THEN
     !
     ALLOCATE( hdum( ldh, n ) )    
     ALLOCATE( iwork( 5*n ) )    
     ALLOCATE( ifail( n ) )    
     !
  END IF
  !
  ! ... input s and (see below) h are copied so that they are not destroyed
  !
  sdum = s
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me == 1 ) THEN
     !
     IF ( all_eigenvalues ) THEN
        !
        ! ... calculate all eigenvalues
        !
        v(:,1:n) = h(:,:)
        !
#if defined (__AIX)
        !
        ! ... there is a name conflict between essl and lapack ...
        !
        CALL DSYGV( 1, v, ldh, sdum, ldh, e, v, ldh, n, work, lwork )
        !
        info = 0
        !
#else
        CALL DSYGV( 1, 'V', 'U', n, v, ldh, sdum, ldh, e, work, &
                    lwork, info )
        !
#endif
     ELSE
        !
        ! ... calculate only m lowest eigenvalues
        !
        hdum = h
        !
        CALL DSYGVX( 1, 'V', 'I', 'U', n, hdum, ldh, sdum, ldh, &
                     0.D0, 0.D0, 1, m, 0.D0, mm, e, v, ldh, work, lwork, &
                     iwork, ifail, info )
        !             
     END IF
     !
     CALL errore( 'rdiaghg', 'info =/= 0', ABS( info ) )
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
  CALL mp_bcast( e, ionode_id, intra_pool_comm )
  CALL mp_bcast( v, ionode_id, intra_pool_comm )
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
  DEALLOCATE( work )
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE rdiaghg
