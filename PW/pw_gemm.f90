!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE pw_gemm( sum_over_nodes, na, nb, n, a, lda, b, ldb, c, ldc )
  !----------------------------------------------------------------------------
  !
  ! ... matrix times matrix with summation index running on G-vectors or PWs
  ! ... c(ij)=real(a(ik)*b(kj)) using half G vectors or half PWs
  !
  ! ccalbec( nkb, npwx, npw, nbnd, vkb, psi, bec ) =>
  !    pw_gemm( 'Y', nkb, nbnd, npw, vkb, npwx, psi, npwx, bec, nkb )  
  !
  !
  USE kinds, ONLY : DP
  USE gvect, ONLY : gstart
  !
  IMPLICIT NONE
  !
  ! ... input
  !
  CHARACTER(LEN=1) :: sum_over_nodes
  INTEGER          :: na, nb, n, lda, ldb, ldc
  COMPLEX(DP)      :: a(lda,na), b(ldb,nb)
  !
  ! ... output
  !
  REAL(DP) :: c(ldc,nb)
  !
  !
  IF ( na == 0 .OR. nb == 0 ) RETURN
  !
  CALL start_clock( 'calbec' )
  !
  IF ( nb == 1 ) THEN
     !
     CALL DGEMV( 'C', 2*n, na, 2.D0, a, 2*lda, b, 1, 0.D0, c, 1 )
     !
     IF ( gstart == 2 ) c(:,1) = c(:,1) - a(1,:) * b(1,1)
     !
  ELSE
     !
     CALL DGEMM( 'C', 'N', na, nb, 2*n, 2.D0, a, 2*lda, b, 2*ldb, 0.D0, c, ldc )
     !
     IF ( gstart == 2 ) &
        CALL DGER( na, nb, -1.D0, a, 2*lda, b, 2*ldb, c, ldc )
     !
  END IF
  !
  IF ( sum_over_nodes == 'y' .OR. &
       sum_over_nodes == 'Y' ) CALL reduce( ldc*nb, c )
  !
  CALL stop_clock( 'calbec' )
  !
  RETURN
  !
END SUBROUTINE pw_gemm
