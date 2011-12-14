!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE transform_becsum_nc(becsum_nc,becsum,na)
!----------------------------------------------------------------------------
!
! This routine multiply becsum_nc by the identity and the Pauli
! matrices and saves it in becsum for the calculation of 
! augmentation charge and magnetization.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE uspp_param,           ONLY : nh, nhm
USE lsda_mod,             ONLY : nspin
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : domag
!
IMPLICIT NONE

COMPLEX(DP) :: becsum_nc(nhm*(nhm+1)/2,nat,npol,npol)
REAL(DP) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
INTEGER :: na
!
! ... local variables
!
INTEGER :: ih, jh, ijh, np

np=ityp(na)
ijh=1
DO ih = 1, nh(np)
   becsum(ijh,na,1)= becsum(ijh,na,1)+  &
               becsum_nc(ijh,na,1,1)+becsum_nc(ijh,na,2,2)
   IF (domag) THEN
      becsum(ijh,na,2)= becsum(ijh,na,2)+  &
               becsum_nc(ijh,na,1,2)+becsum_nc(ijh,na,2,1)
      becsum(ijh,na,3)= becsum(ijh,na,3)+(0.d0,-1.d0)*  &
              (becsum_nc(ijh,na,1,2)-becsum_nc(ijh,na,2,1))
      becsum(ijh,na,4)= becsum(ijh,na,4)+  &
               becsum_nc(ijh,na,1,1)-becsum_nc(ijh,na,2,2)
   END IF
   ijh=ijh+1
   DO jh = ih+1, nh(np)
      becsum(ijh,na,1)= becsum(ijh,na,1) +     &
                 (becsum_nc(ijh,na,1,1)+becsum_nc(ijh,na,2,2))  &
                + CONJG(becsum_nc(ijh,na,1,1)+becsum_nc(ijh,na,2,2))
      IF (domag) THEN
         becsum(ijh,na,2)= becsum(ijh,na,2) +     &
                  becsum_nc(ijh,na,1,2)+becsum_nc(ijh,na,2,1)   &
                + CONJG(becsum_nc(ijh,na,2,1)+becsum_nc(ijh,na,1,2))
         becsum(ijh,na,3)= becsum(ijh,na,3) +(0.d0,-1.d0)*     &
                 (becsum_nc(ijh,na,1,2)-becsum_nc(ijh,na,2,1)   &
                + CONJG(becsum_nc(ijh,na,2,1)-becsum_nc(ijh,na,1,2)) )
         becsum(ijh,na,4)= becsum(ijh,na,4) +     &
                 (becsum_nc(ijh,na,1,1)-becsum_nc(ijh,na,2,2))  &
                + CONJG(becsum_nc(ijh,na,1,1)-becsum_nc(ijh,na,2,2))
      END IF
      ijh=ijh+1
   END DO
END DO

RETURN
END SUBROUTINE transform_becsum_nc
