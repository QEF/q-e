!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------------
SUBROUTINE set_asr_c(nat,nasr,dyn)
  !---------------------------------------------------------------------------
  !
  ! Impose Acoustic Sum Rule on the dynamical matrix
  ! We assume that (3*nat-1) columns have been calculated
  ! and that the missing column corresponds to atom nasr
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: nat, nasr
  COMPLEX(DP) :: dyn(3*nat,3*nat)
  !
  INTEGER :: na, nb, i,j
  COMPLEX(DP) :: sum

  IF (nasr.LE.0 .OR. nasr.GT.nat) RETURN
  DO j=1,3
     DO i=1,3
        DO nb=1,nat
           sum=(0.d0,0.d0)
           DO na=1,nat
              IF (na.NE.nasr) sum = sum + dyn(3*(na-1)+i,3*(nb-1)+j)
           END DO
           dyn(3*(nasr-1)+i,3*(nb-1)+j)= -sum
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE set_asr_c
