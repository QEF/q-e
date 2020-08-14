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
  !! Calculates all the eigenvalues and eigenvectors of a real
  !! simmetric matrix H. On output, the matrix is unchanged.
  !
  USE kinds,            ONLY : DP
  USE mp_bands,         ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
  USE mp,               ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER :: n
  !! input: dimension of the matrix to be diagonalized
  INTEGER :: ldh
  !! input: leading dimension of h, as declared in the calling pgm unit
  REAL(DP) :: h(ldh,n)
  !! input: matrix to be diagonalized
  REAL(DP) :: e(n)
  !! output: eigenvalues
  REAL(DP) :: v(ldh,n)
  !! output: eigenvectors (column-wise)
  !
  ! ... local variables (LAPACK version)
  !
  INTEGER :: lwork, nb, info
  INTEGER, EXTERNAL :: ILAENV
  ! ILAENV returns optimal block size "nb"
  REAL(KIND=DP), ALLOCATABLE :: work(:)
  !
  CALL start_clock( 'diagh' )  
  !
  ! ... check for the block size
  !
  nb = ILAENV( 1, 'DSYTRD', 'U', n, - 1, - 1, - 1 )
  IF ( nb < 1 .OR. nb >= n ) THEN
     lwork = 3*n
  ELSE
     lwork = ( nb + 2 ) * n
  ENDIF
  !
  ! ... only the first processor diagonalize the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     ! ... allocate workspace
     !
     v = h
     !
     ALLOCATE( work(lwork) )    
     !
     CALL DSYEV( 'V', 'U', n, v, ldh, e, work, lwork, info )
     !
     CALL errore( 'rdiagh', 'diagonalization (DSYEV) failed', ABS( info ) )
     !
     ! ... deallocate workspace
     !
     DEALLOCATE( work )
     !
  ENDIF
  !
  CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )      
  !
  CALL stop_clock( 'diagh' )
  !
  RETURN
  !
END SUBROUTINE rdiagh
