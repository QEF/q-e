!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE set_kplusq( xk, wk, xq, nks, npk )
  !-----------------------------------------------------------------------
  !! This routine sets the k and k+q points (with zero weight) used in
  !! the preparatory run for a linear response calculation.
  !
  !! * on input: xk and wk contain k-points and corresponding weights;
  !! * on output: the number of points is doubled and xk and wk in the
  !!              odd  positions are the original ones; those in the
  !!              even positions are the corresponding k+q values.
  !
  !! The gamma point is treated in a special way. No change is done
  !! to the k-points.
  !
  USE kinds,   ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER :: npk
  !! inout: maximum allowed number of k
  INTEGER :: nks
  !! inout: starting and ending number of
  REAL(DP) :: xk(3,npk)
  !! inout: coordinates of k points
  REAL(DP) :: wk(npk)
  !! inout: weights of k points
  REAL(DP) :: xq(3)
  !! input: coordinates of a q-point
  ! 
  ! ... local variables
  !
  REAL(DP) :: eps
  ! the smallest xq
  LOGICAL :: lgamma
  ! true if xq is the gamma point
  INTEGER :: ik, j
  ! counter on k
  ! counter
  !
  eps = 1.d-12
  !
  ! ... shift the k points in the odd positions and fill the even ones with k+
  !
  lgamma = ABS(xq(1))<eps .AND. ABS(xq(2))<eps .AND. ABS(xq(3))<eps
  !
  IF (.NOT.lgamma) THEN
     !
     IF (2 * nks > npk) CALL errore( 'set_kplusq', 'too many k points', nks )
     DO ik = nks, 1, - 1
        DO j = 1, 3
           xk(j,2*ik-1) = xk(j,ik)
           xk(j,2*ik) = xk(j,ik) + xq(j)
        ENDDO
        wk(2*ik-1) = wk(ik)
        wk(2*ik) = 0.d0
     ENDDO
     nks = 2 * nks
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE set_kplusq
