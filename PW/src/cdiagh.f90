!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE cdiagh( n, h, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates all the eigenvalues and eigenvectors of a complex
  ! ... hermitean matrix H. On output, the matrix is unchanged
  !
  USE kinds,            ONLY : DP
  USE mp_bands,         ONLY : nbgrp, me_bgrp, root_bgrp, intra_bgrp_comm
  USE mp,               ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, ldh
    ! dimension of the matrix to be diagonalized
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP) :: h(ldh,n)
    ! matrix to be diagonalized
  !
  ! ... on OUTPUT
  !
  REAL(DP)    :: e(n)      ! eigenvalues
  COMPLEX(DP) :: v(ldh,n)  ! eigenvectors (column-wise)
  !
  ! ... local variables for LAPACK 
  !
  INTEGER                  :: lwork, nb, info
  REAL(DP),    ALLOCATABLE :: rwork(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
  INTEGER, EXTERNAL :: ILAENV
  ! ILAENV returns optimal block size "nb"
  !
  CALL start_clock( 'diagh' )  
  !
  ! ... check for the block size
  !
  nb = ILAENV( 1, 'ZHETRD', 'U', n, - 1, - 1, - 1 )
  IF ( nb < 1 .OR. nb >= n ) THEN
     lwork = 2*n
  ELSE
     lwork = ( nb + 1 )*n
  END IF
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     ! ... allocate workspace
     !
#if defined(__PGI)
     !     workaround for PGI compiler bug
     !
     v(1:ldh,1:n) = h(1:ldh,1:n)
#else
     v = h
#endif
     !
     ALLOCATE( work( lwork ) )    
     ALLOCATE( rwork( 3 * n - 2 ) )    
     !
     CALL ZHEEV( 'V', 'U', n, v, ldh, e, work, lwork, rwork, info )
     !
     CALL errore( 'cdiagh', 'diagonalization (ZHEEV) failed', ABS( info ) )
     !
     ! ... deallocate workspace
     !
     DEALLOCATE( rwork )
     DEALLOCATE( work )
     !
  END IF
  !
#if defined(__PGI)
  !      workaround for PGI compiler bug
  !
  CALL mp_bcast( e(1:n), root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v(1:ldh,1:n), root_bgrp, intra_bgrp_comm )      
#else
  CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )      
#endif
  !
  CALL stop_clock( 'diagh' )
  !
  RETURN
  !
END SUBROUTINE cdiagh
