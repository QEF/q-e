!
! Copyright (C) 2001-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rdiagh( n, h, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates all the eigenvalues and eigenvectors of a real
  ! ... simmetric matrix H . On output, the matrix is unchanged
  !
  USE kinds,            ONLY : DP
  USE mp_bands,         ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
  USE mp,               ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, ldh
    ! dimension of the matrix to be diagonalized
    ! leading dimension of h, as declared in the calling pgm unit
  REAL(DP) :: h(ldh,n)
    ! matrix to be diagonalized
  !
  ! ... on OUTPUT
  !
  REAL(DP) :: e(n)       ! eigenvalues
  REAL(DP) :: v(ldh,n)   ! eigenvectors (column-wise)
  !
  !
  CALL start_clock( 'diagh' )  
  !
#if defined (__ESSL)
  CALL rdiagh_aix()
#else
  CALL rdiagh_lapack()
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
    SUBROUTINE rdiagh_aix()
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
         CALL DSPEV( 21, hp, e, v, ldh, n, aux, naux )
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
    END SUBROUTINE rdiagh_aix 
    !
#else 
    !
    !-----------------------------------------------------------------------
    SUBROUTINE rdiagh_lapack( )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! ... local variables (LAPACK version)
      !
      INTEGER :: lwork, nb, info
      INTEGER, EXTERNAL :: ILAENV
        ! ILAENV returns optimal block size "nb"
      REAL (KIND=DP), ALLOCATABLE :: work(:)
      !
      !
      ! ... check for the block size
      !
      nb = ILAENV( 1, 'DSYTRD', 'U', n, - 1, - 1, - 1 )
      !
      IF ( nb < 1 .OR. nb >= n ) THEN
         !
         lwork = 3*n
         !
      ELSE
         !
         lwork = ( nb + 2 ) * n
         !
      END IF
      !
      ! ... only the first processor diagonalize the matrix
      !
      IF ( me_bgrp == root_bgrp ) THEN
         !
         ! ... allocate workspace
         !
         v = h
         !
         ALLOCATE( work( lwork ) )    
         !
         CALL DSYEV( 'V', 'U', n, v, ldh, e, work, lwork, info )
         !
         CALL errore( 'rdiagh', 'diagonalization (DSYEV) failed', ABS( info ) )
         !
         ! ... deallocate workspace
         !
         DEALLOCATE( work )
         !
      END IF
      !
      CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
      CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )      
      !
      RETURN
      !
    END SUBROUTINE rdiagh_lapack
    !
#endif
    !
END SUBROUTINE rdiagh
