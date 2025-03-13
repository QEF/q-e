!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE g2_kin ( ik )
  !----------------------------------------------------------------------------
  !! Calculation of kinetic energy - includes the case of the modified
  !! kinetic energy functional for variable-cell calculations.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba2 
  USE gvecw,                ONLY : ecfixed, qcutz, q2sigma
  USE klist,                ONLY : xk, ngk, igk_k
  USE gvect,                ONLY : g
  USE wvfct,                ONLY : g2kin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (IN) :: ik
  !
  ! ... local variables
  !
  INTEGER :: ig, npw,i
  REAL(DP):: xk1,xk2,xk3
  !
  npw = ngk(ik)
  !
  xk1 = xk(1,ik)
  xk2 = xk(2,ik)
  xk3 = xk(3,ik)
  !
!$acc parallel loop present(g2kin, g, igk_k)
  DO i=1,npw
     g2kin(i) = ( ( xk1 + g(1,igk_k(i,ik)) )*( xk1 + g(1,igk_k(i,ik)) ) + &
                  ( xk2 + g(2,igk_k(i,ik)) )*( xk2 + g(2,igk_k(i,ik)) ) + &
                  ( xk3 + g(3,igk_k(i,ik)) )*( xk3 + g(3,igk_k(i,ik)) ) ) * tpiba2
  !
  END DO
  !
  IF ( qcutz > 0.D0 ) THEN
     !
!$acc parallel loop present(g2kin)
     DO ig = 1, npw
        !
        g2kin(ig) = g2kin(ig) + qcutz * &
             ( 1.D0 + erf( ( g2kin(ig) - ecfixed ) / q2sigma ) )
        !
     END DO
     !
  END IF
  !
  !$acc update self(g2kin)
  !
  RETURN
  !
END SUBROUTINE g2_kin
