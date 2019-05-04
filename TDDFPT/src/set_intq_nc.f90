!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!----------------------------------------------------------------------------
SUBROUTINE set_intq_nc()
!----------------------------------------------------------------------------
USE ions_base,  ONLY : nat, ntyp => nsp, ityp
USE uspp_param, ONLY : upf
USE lr_variables,    ONLY : intq, intq_nc
IMPLICIT NONE
INTEGER :: npe
INTEGER :: np, na

intq_nc=(0.d0,0.d0)
DO np = 1, ntyp
   IF ( upf(np)%tvanp ) THEN
      DO na = 1, nat
         IF (ityp(na)==np) THEN
            IF (upf(np)%has_so) THEN
               CALL transform_intq_so(intq,na)
            ELSE
               CALL transform_intq_nc(intq,na)
            END IF
         END IF
      END DO
   END IF
END DO

RETURN
END SUBROUTINE set_intq_nc
!
