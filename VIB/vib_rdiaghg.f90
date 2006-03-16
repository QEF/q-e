!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE vib_rdiaghg( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H symmetric matrix, S overlap matrix .
  ! ... On output both matrix are unchanged
  ! ... Uses LAPACK routines
  !
  !USE kinds,     ONLY : 8
  !USE mp,        ONLY : mp_bcast
  !USE mp_global, ONLY : npool, me_pool, root_pool, intra_pool_comm, my_image_id
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculated
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(KIND=8) :: h(ldh,n)
    ! matrix to be diagonalized
  REAL(KIND=8) :: s(ldh,n)
    ! overlap matrix
  !
  ! ... on OUTPUT
  !
  REAL(KIND=8) :: e(n)
    ! eigenvalues
  REAL(KIND=8) :: v(ldh,m)
    ! eigenvectors (column-wise)
  !
  ! ... LOCAL variables
  !
  INTEGER                    :: lwork, nb, ILAENV, mm, info
    ! ILAENV returns optimal block size "nb"
    ! mm = number of calculated eigenvectors
  EXTERNAL                      ILAENV
  INTEGER, ALLOCATABLE       :: iwork(:), ifail(:)
  REAL(KIND=8), ALLOCATABLE :: sdum(:,:), hdum(:,:),  work(:)
  LOGICAL                    :: all_eigenvalues
  !
  !
  !CALL start_clock( 'cdiaghg' )
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
  !IF ( me_pool == root_pool ) THEN
     !
     IF ( all_eigenvalues ) THEN
        !
        ! ... calculate all eigenvalues
        !
        v(:,1:n) = h(:,:)
        !
#if defined (__ESSL)
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
     !CALL errore( 'rdiaghg', 'info =/= 0', ABS( info ) )
     if (info.ne.0) then
        write (6,*) 'error: rdiaghg - info =/= 0i  ',ABS( info )
     end if
     !
  !END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
  !CALL mp_bcast( e, root_pool, intra_pool_comm )
  !CALL mp_bcast( v, root_pool, intra_pool_comm )
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
  !CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE vib_rdiaghg
