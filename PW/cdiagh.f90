!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cdiagh( n, h, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates all the eigenvalues and eigenvectors of a complex
  ! ... hermitean matrix H . On output, the matrix is unchanged
  !
  USE kinds,     ONLY : DP
#if defined (__PARA)
  USE para,      ONLY : me, mypool, npool, MPI_COMM_POOL
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_bcast  
#endif
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, &               ! dimension of the matrix to be diagonalized
             ldh                ! leading dimension of h, as declared in
                                ! the calling pgm unit
  COMPLEX(KIND=DP) :: &
           h(ldh,n)             ! matrix to be diagonalized
  !
  ! ... on OUTPUT
  !
  REAL(KIND=DP)    :: e(n)      ! eigenvalues
  COMPLEX(KIND=DP) :: v(ldh,n)  ! eigenvectors (column-wise)
  !
  !
  CALL start_clock( 'cdiagh' )  
  !
#if defined (__AIX)
  !
  CALL cdiagh_aix()
  !
#else 
# if defined (CRAYY)
  !
  CALL cdiagh_crayy() 
  !
# else
  !
  CALL cdiagh_lapack()
  !
# endif
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
      COMPLEX(KIND=DP), ALLOCATABLE :: hp(:), aux(:)
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
      IF ( me == 1 ) THEN
         !
# endif
         !
         CALL ZHPEV( 21, hp, e, v, ldh, n, aux, naux )
         !
# if defined (__PARA)
         !
      END IF
      !
      CALL mp_bcast( e, ionode_id, MPI_COMM_POOL )
      CALL mp_bcast( v, ionode_id, MPI_COMM_POOL )
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
# if defined (CRAYY)
    !
    !-----------------------------------------------------------------------
    SUBROUTINE cdiagh_crayy()
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      ! ... local variables (Cray Eispack/Scilib version)
      !
      INTEGER :: i, j, k, info
      REAL(KIND=DP) :: ar(ldh,n), ai(ldh,n), zr(ldh,n), zi(ldh,n)
        ! real and imaginary part of  h(ldh,n) and of  v(ldh,n)
        ! (used as auxiliary arrays)
      REAL(KIND=DP) :: rwork(2,ldh), work(ldh)
      !
      !
      ar = REAL( h )
      ai = AIMAG( h )
      !
      CALL ch( ldh, n, ar, ai, e, 1, zr, zi, work, work, rwork, info )
      !
      CALL errore( 'cdiagh', 'info =/= 0', ABS( info ) )
      !
      v = CMPLX( zr, zi )      
      !
      RETURN      
      !      
    END SUBROUTINE cdiagh_crayy 
    !
# else        
    !
    !-----------------------------------------------------------------------
    SUBROUTINE cdiagh_lapack()
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      !
      ! ... local variables (LAPACK version)
      !
      INTEGER :: lwork, nb, info
      INTEGER :: ILAENV
      EXTERNAL   ILAENV
        ! ILAENV returns optimal block size "nb"
      REAL(KIND=DP),    ALLOCATABLE :: rwork(:)
      COMPLEX(KIND=DP), ALLOCATABLE :: work(:)
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
      ! ... if scalapack library is present and we have just one pool
      ! ... and the matrix is larger than 130 we use the scalapack driver
      !
#   if defined (__T3E)
      !
      IF ( npool == 1 .AND. n > 130 ) THEN
         !
         CALL scala_cdiag( n, h, ldh, e, v, ldh )
         !
         RETURN
         !
      END IF
      !
#   endif
      !
      ! ... else only the first processor diagonalize the matrix
      !
      IF ( me == 1 ) THEN
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
      CALL mp_bcast( e, ionode_id, MPI_COMM_POOL )
      CALL mp_bcast( v, ionode_id, MPI_COMM_POOL )      
      !
#  endif
      !
      RETURN
      !
    END SUBROUTINE cdiagh_lapack
    !
# endif
#endif
    !
END SUBROUTINE cdiagh
