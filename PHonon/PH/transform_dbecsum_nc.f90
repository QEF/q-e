!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE transform_dbecsum_nc(dbecsum_nc,dbecsum,na,modes)
!----------------------------------------------------------------------------
!
! This routine multiply dbecsum_nc by the identity and the Pauli
! matrices and saves it in dbecsum to use it in the calculation of
! the charge and magnetization.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE uspp_param,           ONLY : nh, nhm
USE lsda_mod,             ONLY : nspin
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : domag
!
IMPLICIT NONE

INTEGER :: na, modes
COMPLEX(DP) :: dbecsum_nc( nhm, nhm, nat , nspin , modes)
COMPLEX(DP) :: dbecsum( nhm * (nhm + 1) /2 , nat , nspin_mag, modes)
!
! ... local variables
!
INTEGER :: ih, jh, ijh, np, mode

np=ityp(na)
DO mode=1, modes
   ijh=1
   DO ih = 1, nh(np)
      dbecsum(ijh,na,1,mode)= dbecsum(ijh,na,1,mode)+  &
               dbecsum_nc(ih,ih,na,1,mode)+dbecsum_nc(ih,ih,na,4,mode)
      IF (domag) THEN
         dbecsum(ijh,na,2,mode)= dbecsum(ijh,na,2,mode)+  &
                  dbecsum_nc(ih,ih,na,2,mode)+ &
                            dbecsum_nc(ih,ih,na,3,mode)
         dbecsum(ijh,na,3,mode)= dbecsum(ijh,na,3,mode)+ &
                  (0.d0,-1.d0)*(dbecsum_nc(ih,ih,na,2,mode)- &
                            dbecsum_nc(ih,ih,na,3,mode) )
         dbecsum(ijh,na,4,mode)= dbecsum(ijh,na,4,mode)+  &
                  dbecsum_nc(ih,ih,na,1,mode)-dbecsum_nc(ih,ih,na,4,mode)
      END IF
      ijh=ijh+1
      DO jh = ih+1, nh(np)
         dbecsum(ijh,na,1,mode)= dbecsum(ijh,na,1,mode) +                   &
                   dbecsum_nc(ih,jh,na,1,mode)+dbecsum_nc(ih,jh,na,4,mode)  &
                  +dbecsum_nc(jh,ih,na,1,mode)+dbecsum_nc(jh,ih,na,4,mode)
         IF (domag) THEN
            dbecsum(ijh,na,2,mode)= dbecsum(ijh,na,2,mode) +     &
                      dbecsum_nc(ih,jh,na,2,mode)+                 &
                             dbecsum_nc(ih,jh,na,3,mode)     &
                   +  dbecsum_nc(jh,ih,na,2,mode)+           &
                             dbecsum_nc(jh,ih,na,3,mode)
            dbecsum(ijh,na,3,mode)= dbecsum(ijh,na,3,mode) +        &
                      (0.d0,-1.d0)*(dbecsum_nc(ih,jh,na,2,mode)-    &
                                    dbecsum_nc(ih,jh,na,3,mode)     &
                   +                dbecsum_nc(jh,ih,na,2,mode)-    &
                                    dbecsum_nc(jh,ih,na,3,mode) )
            dbecsum(ijh,na,4,mode)= dbecsum(ijh,na,4,mode) +     &
                      dbecsum_nc(ih,jh,na,1,mode)-dbecsum_nc(ih,jh,na,4,mode)+&
                      dbecsum_nc(jh,ih,na,1,mode)-dbecsum_nc(jh,ih,na,4,mode)
         END IF
         ijh=ijh+1
      END DO
   END DO
END DO

RETURN
END SUBROUTINE transform_dbecsum_nc
