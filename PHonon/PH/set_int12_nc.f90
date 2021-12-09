!
! Copyright (C) 2007-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!----------------------------------------------------------------------------
SUBROUTINE set_int12_nc(iflag)
!----------------------------------------------------------------------------
!! This is a driver to call the routines that rotate and multiply
!! by the Pauli matrices the integrals.
!
USE ions_base, ONLY : nat, ntyp => nsp, ityp
USE noncollin_module, ONLY : noncolin, domag, lspinorb
USE uspp_param, only: upf
USE phus, ONLY : int1, int2, int1_nc, int2_so
USE nc_mag_aux, ONLY : int1_nc_save
IMPLICIT NONE
INTEGER :: iflag
INTEGER :: np, na

IF (noncolin.AND.domag) THEN
   int1_nc=(0.d0,0.d0)
   int1 ( :, :, :, :, 2:4)=-int1 ( :, :, :, :, 2:4)
   DO np = 1, ntyp
      IF ( upf(np)%tvanp ) THEN
         DO na = 1, nat
            IF (ityp(na)==np) THEN
               IF (upf(np)%has_so) THEN
                  CALL transform_int1_so(int1,na,iflag)
               ELSE
                  CALL transform_int1_nc(int1,na,iflag)
               END IF
            END IF
         END DO
      END IF
   END DO
   int1 ( :, :, :, :, 2:4)=-int1 ( :, :, :, :, 2:4)
   int1_nc_save(:,:,:,:,:,2) = int1_nc
END IF

int1_nc=(0.d0,0.d0)
IF (lspinorb) int2_so=(0.d0,0.d0)
DO np = 1, ntyp
   IF ( upf(np)%tvanp ) THEN
      DO na = 1, nat
         IF (ityp(na)==np) THEN
            IF (upf(np)%has_so) THEN
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
IF (noncolin.AND.domag) int1_nc_save(:,:,:,:,:,1) = int1_nc

END SUBROUTINE set_int12_nc

