!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE g2_kin_gpu ( ik )
  !----------------------------------------------------------------------------
  !
  ! ... Calculation of kinetic energy - includes the case of the modified
  ! ... kinetic energy functional for variable-cell calculations
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba2 
  USE gvecw,                ONLY : ecfixed, qcutz, q2sigma
  USE klist,                ONLY : xk, ngk, igk_k_d
  USE wvfct_gpum,           ONLY : g2kin_d, using_g2kin_d
  USE gvect_gpum,           ONLY : g_d, using_g_d
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (IN) :: ik
#if defined(__CUDA)
  !
  ! ... local variables
  !
  INTEGER :: ig, npw,i
  REAL(DP):: xk1,xk2,xk3
  !
  CALL using_g_d(0); CALL using_g2kin_d(2)
  !
  !
  !
  npw = ngk(ik)
  !
  xk1 = xk(1,ik)
  xk2 = xk(2,ik)
  xk3 = xk(3,ik)
  !
!$cuf kernel do(1) <<<*,*>>>
  DO i=1,npw
     g2kin_d(i) = ( ( xk1 + g_d(1,igk_k_d(i,ik)) )*( xk1 + g_d(1,igk_k_d(i,ik)) ) + &
                  ( xk2 + g_d(2,igk_k_d(i,ik)) )*( xk2 + g_d(2,igk_k_d(i,ik)) ) + &
                  ( xk3 + g_d(3,igk_k_d(i,ik)) )*( xk3 + g_d(3,igk_k_d(i,ik)) ) ) * tpiba2
  !
  END DO
  !
  IF ( qcutz > 0.D0 ) THEN
     !
!$cuf kernel do(1) <<<*,*>>>
     DO ig = 1, npw
        !
        g2kin_d(ig) = g2kin_d(ig) + qcutz * &
             ( 1.D0 + erf( ( g2kin_d(ig) - ecfixed ) / q2sigma ) )
        !
     END DO
     !
  END IF
  !
  RETURN
  !
#else
  CALL errore('g2_kin_gpu', 'Trying to use device subroutine but code was not compiled with device support!', 1)
#endif
END SUBROUTINE g2_kin_gpu
