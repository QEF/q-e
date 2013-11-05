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
  CALL start_clock( 'diagh' )  
  !
#if defined (__ESSL)
  CALL cdiagh_aix()
#else
  CALL cdiagh_lapack( v, e )
#endif
  !
  CALL stop_clock( 'diagh' )
  !
  RETURN
  !
  CONTAINS  
    !
    ! ... internal procedures
    !
#if defined (__ESSL)
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE cdiagh_aix()
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! ... local variables (ESSL version)
      !
      INTEGER                  :: naux, i, j, ij
      COMPLEX(DP), ALLOCATABLE :: hp(:), aux(:)
      !
      !
      naux = 4 * n
      !
      ALLOCATE( hp(  n * (n + 1) / 2 ) )    
      ALLOCATE( aux( naux ) )    
      !
      ! ... copy to upper triangular packed matrix
      !
      ij = 0
      DO j = 1, n
         DO i = 1, j
            ij = ij + 1
            hp(ij) = h(i,j)
         END DO
      END DO
      !
      ! ... only the first processor diagonalize the matrix
      !
      IF ( me_bgrp == root_bgrp ) THEN
         !
         CALL ZHPEV( 21, hp, e, v, ldh, n, aux, naux )
         !
      END IF
      !
      CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
      CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )
      !
      DEALLOCATE( aux )
      DEALLOCATE( hp )
      !
      RETURN
      ! 
    END SUBROUTINE cdiagh_aix 
    !
#else 
    !
    !-----------------------------------------------------------------------
    SUBROUTINE cdiagh_lapack( v, e )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP)    :: e(n)      ! eigenvalues
      COMPLEX(DP) :: v(ldh,n)
      !
      ! ... local variables (LAPACK version)
      !
      INTEGER                  :: lwork, nb, info
      REAL(DP),    ALLOCATABLE :: rwork(:)
      COMPLEX(DP), ALLOCATABLE :: work(:)
      !
      INTEGER, EXTERNAL :: ILAENV
        ! ILAENV returns optimal block size "nb"
      !
      ! ... check for the block size
      !
      nb = ILAENV( 1, 'ZHETRD', 'U', n, - 1, - 1, - 1 )
      !
      IF ( nb < 1 .OR. nb >= n ) THEN
         !
         lwork = 2*n
         !
      ELSE
         !
         lwork = ( nb + 1 )*n
         !
      END IF
      !
      ! ... only the first processor diagonalize the matrix
      !
      IF ( me_bgrp == root_bgrp ) THEN
         !
         ! ... allocate workspace
         !
#ifdef __PGI
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
#ifdef __PGI
      !      workaround for PGI compiler bug
      !
      CALL mp_bcast( e(1:n), root_bgrp, intra_bgrp_comm )
      CALL mp_bcast( v(1:ldh,1:n), root_bgrp, intra_bgrp_comm )      
#else
      CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
      CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )      
#endif
      !
      RETURN
      !
    END SUBROUTINE cdiagh_lapack
    !
#endif
    !
END SUBROUTINE cdiagh
