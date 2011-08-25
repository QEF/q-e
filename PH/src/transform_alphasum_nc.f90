!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE transform_alphasum_nc(alphasum_nc,na)
!----------------------------------------------------------------------------
!
! This routine multiply alphasum_nc by the identity and the Pauli
! matrices and saves it in alphasum to use it in the calculation of
! the change of the charge and of the magnetization.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol
USE spin_orb,             ONLY : domag
USE phus,                 ONLY : alphasum
!
IMPLICIT NONE

COMPLEX(DP) :: alphasum_nc(nhm*(nhm+1)/2,3,nat,npol,npol)
INTEGER :: na
!
! ... local variables
!
INTEGER :: ih, jh, ijh, np, ipol

np=ityp(na)
DO ipol=1,3
   ijh=1
   DO ih = 1, nh(np)
      alphasum(ijh,ipol,na,1)= alphasum(ijh,ipol,na,1)+  &
               alphasum_nc(ijh,ipol,na,1,1)+alphasum_nc(ijh,ipol,na,2,2)
      IF (domag) THEN
         alphasum(ijh,ipol,na,2)= alphasum(ijh,ipol,na,2)+  &
                  alphasum_nc(ijh,ipol,na,1,2)+alphasum_nc(ijh,ipol,na,2,1)
         alphasum(ijh,ipol,na,3)= alphasum(ijh,ipol,na,3)+(0.d0,-1.d0)*  &
                 (alphasum_nc(ijh,ipol,na,1,2)-alphasum_nc(ijh,ipol,na,2,1))
         alphasum(ijh,ipol,na,4)= alphasum(ijh,ipol,na,4)+  &
                  alphasum_nc(ijh,ipol,na,1,1)-alphasum_nc(ijh,ipol,na,2,2)
      END IF
      ijh=ijh+1
      DO jh = ih+1, nh(np)
         alphasum(ijh,ipol,na,1)= alphasum(ijh,ipol,na,1) +     &
              (alphasum_nc(ijh,ipol,na,1,1)+alphasum_nc(ijh,ipol,na,2,2))  &
            + CONJG(alphasum_nc(ijh,ipol,na,1,1)+alphasum_nc(ijh,ipol,na,2,2))
         IF (domag) THEN
            alphasum(ijh,ipol,na,2)= alphasum(ijh,ipol,na,2) +     &
               alphasum_nc(ijh,ipol,na,1,2)+alphasum_nc(ijh,ipol,na,2,1)   &
            + CONJG(alphasum_nc(ijh,ipol,na,2,1)+alphasum_nc(ijh,ipol,na,1,2))
            alphasum(ijh,ipol,na,3)= alphasum(ijh,ipol,na,3) +(0.d0,-1.d0)*  &
              (alphasum_nc(ijh,ipol,na,1,2)-alphasum_nc(ijh,ipol,na,2,1)   &
            + CONJG(alphasum_nc(ijh,ipol,na,2,1)-alphasum_nc(ijh,ipol,na,1,2)))
            alphasum(ijh,ipol,na,4)= alphasum(ijh,ipol,na,4) +             &
              (alphasum_nc(ijh,ipol,na,1,1)-alphasum_nc(ijh,ipol,na,2,2))  &
            + CONJG(alphasum_nc(ijh,ipol,na,1,1)-alphasum_nc(ijh,ipol,na,2,2))
         END IF
         ijh=ijh+1
      END DO
   END DO
END DO

RETURN
END SUBROUTINE transform_alphasum_nc
