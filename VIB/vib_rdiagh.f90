!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE vib_rdiagh( n, h, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates all the eigenvalues and eigenvectors of a complex
  ! ... hermitean matrix H . On output, the matrix is unchanged
  !
  !USE kinds,     ONLY : 8
  !USE mp_global, ONLY : npool, me_pool, root_pool, intra_pool_comm, my_image_id
  !USE mp,        ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  ! ... on INPUT
  !
  INTEGER :: n, &               ! dimension of the matrix to be diagonalized
             ldh                ! leading dimension of h, as declared in
                                ! the calling pgm unit
  REAL (KIND=8) :: h(ldh,n)    ! matrix to be diagonalized
  !
  ! ... on OUTPUT
  !
  REAL (KIND=8) :: e(n)      ! eigenvalues
  REAL (KIND=8) :: v(ldh,n)  ! eigenvectors (column-wise)
  !
  !
  !CALL start_clock( 'rdiagh' )  
  !   
  CALL rdiagh_lapack( )
  !
  !CALL stop_clock( 'rdiagh' )
  !
  RETURN
  !
  CONTAINS  
    !
    ! ... internal procedures
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
      REAL (KIND=8), ALLOCATABLE :: work(:)
      !
      !
      ! ... check for the block size
      !
      nb = ILAENV( 1, 'DSYTRD', 'U', n, - 1, - 1, - 1 )
      !
      IF ( nb < 1 ) nb = MAX( 1, n )
      !
      lwork = ( nb + 3 ) * n
      !
      ! ... only the first processor diagonalize the matrix
      !
      !IF ( me_pool == root_pool ) THEN
         !
         ! ... allocate workspace
         !
         v = h
         !
         ALLOCATE( work( lwork ) )    
         !
         CALL DSYEV( 'V', 'U', n, v, ldh, e, work, lwork, info )
         !
         if (abs(info).ne.0) then 
         !CALL errore( 'rdiagh', 'info =/= 0', ABS( info ) )
             write (6,*) 'error: rdiagh - info =/= 0   ', ABS( info )
         end if
         !
         ! ... deallocate workspace
         !
         DEALLOCATE( work )
         !
      !END IF
      !
      !CALL mp_bcast( e, root_pool, intra_pool_comm )
      !CALL mp_bcast( v, root_pool, intra_pool_comm )      
      !
      RETURN
      !
    END SUBROUTINE rdiagh_lapack
    !
END SUBROUTINE vib_rdiagh

