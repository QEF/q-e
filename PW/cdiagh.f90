!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cdiagh( n, h, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates all the eigenvalues and eigenvectors of a complex
  ! ... hermitean matrix H . On output, the matrix is unchanged
  !
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : npool, me_pool, root_pool, intra_pool_comm, my_image_id
  USE mp,        ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, &               ! dimension of the matrix to be diagonalized
             ldh                ! leading dimension of h, as declared in
                                ! the calling pgm unit
  COMPLEX(DP) :: &
           h(ldh,n)             ! matrix to be diagonalized
  !
  ! ... on OUTPUT
  !
  REAL(DP)    :: e(n)      ! eigenvalues
  COMPLEX(DP) :: v(ldh,n)  ! eigenvectors (column-wise)
  !
  !
  CALL start_clock( 'cdiagh' )  
  !
#if defined (__AIX)
  !
  CALL cdiagh_aix()
  !
#else
  !   
  ! CALL cdiagh_lapack( )
  ! workaround for Intel ifc8 bug:
  CALL cdiagh_lapack( v )
  !
#endif
  !
  CALL stop_clock( 'cdiagh' )
  !
  RETURN
  !
  CONTAINS  
    !
    ! ... internal procedures
    !
#if defined (__AIX)
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE cdiagh_aix()
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! ... local variables (ESSL version)
      !
      INTEGER                       :: naux, i, j, ij
      COMPLEX(DP), ALLOCATABLE :: hp(:), aux(:)
      !
      !
      naux = 4 * n
      !
      ALLOCATE( hp(  n * (n + 1) / 2 ) )    
      ALLOCATE( aux(naux) )    
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
# if defined (__PARA)
      !
      ! ... only the first processor diagonalize the matrix
      !
      IF ( me_pool == root_pool ) THEN
         !
# endif
         !
         CALL ZHPEV( 21, hp, e, v, ldh, n, aux, naux )
         !
# if defined (__PARA)
         !
      END IF
      !
      CALL mp_bcast( e, root_pool, intra_pool_comm )
      CALL mp_bcast( v, root_pool, intra_pool_comm )
      !
# endif
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
    ! workaround for Intel ifc8 bug:
    !    SUBROUTINE cdiagh_lapack( )
    SUBROUTINE cdiagh_lapack( v )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      ! workaround for Intel ifc8 bug:
      COMPLEX(DP) :: v(ldh,n)
      !
      ! ... local variables (LAPACK version)
      !
      INTEGER :: lwork, nb, info
      INTEGER, EXTERNAL :: ILAENV
        ! ILAENV returns optimal block size "nb"
      REAL(DP),    ALLOCATABLE :: rwork(:)
      COMPLEX(DP), ALLOCATABLE :: work(:)
      !
      !
      ! ... check for the block size
      !
      nb = ILAENV( 1, 'ZHETRD', 'U', n, - 1, - 1, - 1 )
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
#  if defined (__PARA)
      !
      ! ... only the first processor diagonalize the matrix
      !
      IF ( me_pool == root_pool ) THEN
#  endif
         !
         ! ... allocate workspace
         !
         v = h
         !
         ALLOCATE( work( lwork ) )    
         ALLOCATE( rwork( 3 * n - 2 ) )    
         !
         CALL ZHEEV( 'V', 'U', n, v, ldh, e, work, lwork, rwork, info )
         !
         CALL errore( 'cdiagh', 'info =/= 0', ABS( info ) )
         !
         ! ... deallocate workspace
         !
         DEALLOCATE( rwork )
         DEALLOCATE( work )
         !
#  if defined (__PARA)
         !
      END IF
      !
      CALL mp_bcast( e, root_pool, intra_pool_comm )
      CALL mp_bcast( v, root_pool, intra_pool_comm )      
      !
#  endif
      !
      RETURN
      !
    END SUBROUTINE cdiagh_lapack
    !
#endif
    !
END SUBROUTINE cdiagh
