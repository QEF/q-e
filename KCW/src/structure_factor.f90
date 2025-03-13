!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!-----------------------------------------------------------------------
SUBROUTINE structure_factor (ik, struct_fact)
  !-----------------------------------------------------------------------
  !
  !! Structure factor
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at
  USE constants,            ONLY : tpi
  USE control_kcw,          ONLY : mp1, mp2, mp3, x_q
  !
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(IN) :: ik
  COMPLEX(DP), INTENT (OUT) :: struct_fact
  REAL (DP) :: xk_(3), R(3)
  INTEGER :: i, j ,k, nkstot_eff
  !
  nkstot_eff = SIZE(x_q)/3
  !
  xk_ (:) = x_q(:,ik)
  CALL cryst_to_cart(1, xk_, at, -1)
  !! The qpoint in crystall coordinate
  !
  struct_fact = 0.D0 
  !
  DO i = 1, mp1
    DO j = 1, mp2
      DO k = 1, mp3
         R(1) = i-1; R(2) =j-1; R(3) = k-1
         ! WRITE(*,*) xk_, R, sum(R(:)*xk_(:))
         ! WRITE(*,*) sum(R(:)*xk_(:)), exp(ONE*tpi*sum(R(:)*xk_(:)))
         struct_fact = struct_fact + exp(-ONE*tpi*sum(R(:)*xk_(:)))
      ENDDO
    ENDDO
  ENDDO
  !
END subroutine
