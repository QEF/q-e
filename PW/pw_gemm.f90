!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!-----------------------------------------------------------------------
SUBROUTINE pw_gemm( sum_over_nodes, na, nb, n, a, lda, b, ldb, c, ldc )
  !-----------------------------------------------------------------------
  !
  ! ... matrix times matrix with summation index running on G-vectors or PWs
  ! ... c(ij)=real(a(ik)*b(kj)) using half G vectors or half PWs
  !
  ! ccalbec( nkb, npwx, npw, nbnd, vkb, psi, bec ) =>
  !    pw_gemm( 'Y', nkb, nbnd, npw, vkb, npwx, psi, npwx, bec, nkb )  
  !
  !
  USE parameters, ONLY : DP
  USE gvect,      ONLY : gstart
  !
  IMPLICIT NONE
  !
  ! ... input
  !
  INTEGER          :: na, nb, n, lda, ldb, ldc
  CHARACTER(LEN=1) :: sum_over_nodes
  COMPLEX(KIND=DP) :: a(lda,na), b(ldb,nb)
  !
  ! ... output
  !
  REAL(KIND=DP)    :: c(ldc,nb)
  !
  !
  IF ( na == 0 .OR. nb == 0 ) RETURN
  !
  CALL start_clock( 'pw_gemm' )
  !
  CALL DGEMM( 'C', 'N', na, nb, 2 * n, 2.D0, a, 2 * lda, b, &
              2 * ldb, 0.D0, c, ldc )
  !
  IF ( gstart == 2 ) &
     CALL DGER(  na, nb, -1.D0, a, ( 2 * lda ), b, ( 2 * ldb ), c, ldc )
  !
#ifdef __PARA
  IF ( sum_over_nodes == 'y' .OR. sum_over_nodes == 'Y' ) &
     CALL reduce( ldc * nb, c )
#endif
  !
  CALL stop_clock( 'pw_gemm' )
  !
  RETURN
  !
END SUBROUTINE pw_gemm
