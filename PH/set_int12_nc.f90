!
! Copyright (C) 2007 QUANTUM-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!----------------------------------------------------------------------------
SUBROUTINE set_int12_nc(iflag)
!----------------------------------------------------------------------------
!
!  This is a driver to call the routines that rotate and multiply
!  by the Pauli matrices the integrals.
!
USE ions_base, ONLY : nat, ntyp => nsp, ityp
USE spin_orb, ONLY : lspinorb, so
USE uspp_param, only: tvanp
USE phus, ONLY : int1, int2, int1_nc, int2_so
IMPLICIT NONE
INTEGER :: iflag
INTEGER :: np, na

int1_nc=(0.d0,0.d0)
IF (lspinorb) int2_so=(0.d0,0.d0)
DO np = 1, ntyp
   IF ( tvanp(np) ) THEN
      DO na = 1, nat
         IF (ityp(na)==np) THEN
            IF (so(np)) THEN
               CALL transform_int1_so(int1,na,iflag)
               CALL transform_int2_so(int2,na,iflag)
            ELSE
               CALL transform_int1_nc(int1,na,iflag)
               IF (lspinorb) CALL transform_int2_nc(int2,na,iflag)
            END IF
         END IF
      END DO
   END IF
END DO
END SUBROUTINE set_int12_nc

!----------------------------------------------------------------------------
SUBROUTINE set_int3_nc(npe)
!----------------------------------------------------------------------------
USE ions_base, ONLY : nat, ntyp => nsp, ityp
USE spin_orb, ONLY : so
USE uspp_param, only: tvanp
USE phus, ONLY : int3, int3_nc
IMPLICIT NONE
INTEGER :: npe
INTEGER :: np, na

int3_nc=(0.d0,0.d0)
DO np = 1, ntyp
   IF ( tvanp(np) ) THEN
      DO na = 1, nat
         IF (ityp(na)==np) THEN
            IF (so(np)) THEN
               CALL transform_int3_so(int3,na,npe)
            ELSE
               CALL transform_int3_nc(int3,na,npe)
            END IF
         END IF
      END DO
   END IF
END DO
END SUBROUTINE set_int3_nc
!
