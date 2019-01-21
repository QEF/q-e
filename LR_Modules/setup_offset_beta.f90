!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE setup_offset_beta
  !----------------------------------------------------------
  !
  ! Calculate the offset of beta functions for each atom na.
  ! Ordering: first all betas for atoms of type 1,
  ! then all betas for atoms of type 2, and so on.
  !
  USE uspp_param,   ONLY : nh
  USE ions_base,    ONLY : nat, ityp, ntyp => nsp
  USE control_lr,   ONLY : ofsbeta   
  !   
  IMPLICIT NONE
  !
  INTEGER ::  na, iat, nt, jkb2, ih
  ! 
  jkb2 = 0
  DO nt = 1, ntyp
     DO na = 1, nat 
        IF ( ityp(na).EQ.nt ) THEN 
           ofsbeta(na) = jkb2
           DO ih = 1, nh(nt)
              jkb2 = jkb2 + 1
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE setup_offset_beta
