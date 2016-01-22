!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE transform_int1_nc(int1,na,iflag)
!----------------------------------------------------------------------------
!
! This routine multiply int1 by the identity and the Pauli
! matrices and saves it in int1_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE spin_orb,             ONLY : domag
USE noncollin_module,     ONLY : nspin_mag
USE phus,                 ONLY : int1_nc
!
IMPLICIT NONE

INTEGER :: na, iflag
COMPLEX(DP) :: int1(nhm,nhm,3,nat,nspin_mag)
!
! ... local variables
!
INTEGER :: ih, jh, ipol, np

np=ityp(na)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      DO ipol=1,3
         IF (iflag==0) THEN
            IF (domag) THEN
               int1_nc(ih,jh,ipol,na,1)= &
                                 int1(ih,jh,ipol,na,1)+int1(ih,jh,ipol,na,4)
               int1_nc(ih,jh,ipol,na,2)=                                       &
                  int1(ih,jh,ipol,na,2) - (0.d0, 1.d0) * int1(ih,jh,ipol,na,3)
               int1_nc(ih,jh,ipol,na,3)=                                       &
                  int1(ih,jh,ipol,na,2) + (0.d0, 1.d0) * int1(ih,jh,ipol,na,3)
               int1_nc(ih,jh,ipol,na,4)=                                       &
                  int1(ih,jh,ipol,na,1) - int1(ih,jh,ipol,na,4)
            ELSE
               int1_nc(ih,jh,ipol,na,1)=int1(ih,jh,ipol,na,1)
               int1_nc(ih,jh,ipol,na,4)=int1(ih,jh,ipol,na,1)
            END IF
         ELSE
            IF (domag) THEN
               int1_nc(ih,jh,ipol,na,1)= &
                             CONJG(int1(ih,jh,ipol,na,1)+int1(ih,jh,ipol,na,4))
               int1_nc(ih,jh,ipol,na,2)=CONJG(int1(ih,jh,ipol,na,2)) - &
                           (0.d0, 1.d0)*CONJG(int1(ih,jh,ipol,na,3))
               int1_nc(ih,jh,ipol,na,3)=CONJG(int1(ih,jh,ipol,na,2)) + &
                           (0.d0, 1.d0)*CONJG(int1(ih,jh,ipol,na,3))
               int1_nc(ih,jh,ipol,na,4)=                               &
                  CONJG(int1(ih,jh,ipol,na,1) - int1(ih,jh,ipol,na,4))
            ELSE
               int1_nc(ih,jh,ipol,na,1)=CONJG(int1(ih,jh,ipol,na,1))
               int1_nc(ih,jh,ipol,na,4)=CONJG(int1(ih,jh,ipol,na,1))
            END IF
         END IF
      END DO
   END DO
END DO

RETURN
END SUBROUTINE transform_int1_nc
!
!----------------------------------------------------------------------------
SUBROUTINE transform_int2_nc(int2, nb, iflag)
!----------------------------------------------------------------------------
!
! This routine sets int2_so for the atomic species which do not
! have a spin-orbit pseudopotential
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE phus,                 ONLY : int2_so
!
IMPLICIT NONE
INTEGER :: nb, iflag
COMPLEX(DP) :: int2(nhm,nhm,3,nat,nat)
!
! ... local variables
!
INTEGER :: ih, jh, np, na, ipol

np=ityp(nb)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      DO na=1,nat
         DO ipol=1,3
            IF (iflag==0) THEN
               int2_so(ih,jh,ipol,na,nb,1)=int2(ih,jh,ipol,na,nb)
               int2_so(ih,jh,ipol,na,nb,4)=int2(ih,jh,ipol,na,nb)
            ELSE
               int2_so(ih,jh,ipol,na,nb,1)=CONJG(int2(ih,jh,ipol,na,nb))
               int2_so(ih,jh,ipol,na,nb,4)=CONJG(int2(ih,jh,ipol,na,nb))
            END IF
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE transform_int2_nc

!----------------------------------------------------------------------------
SUBROUTINE transform_int4_nc(int4,na)
!----------------------------------------------------------------------------
!
! This routine multiply int4 by the identity and the Pauli
! matrices and saves it in int4_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE uspp,                 ONLY : ijtoh
USE noncollin_module,     ONLY : nspin_mag
USE spin_orb,             ONLY : domag
USE phus,                 ONLY : int4_nc
!
IMPLICIT NONE

INTEGER :: na
COMPLEX(DP) :: int4(nhm*(nhm+1)/2,3,3,nat,nspin_mag)
!
! ... local variables
!
INTEGER :: ih, jh, ipol, jpol, np
INTEGER :: ijh

np=ityp(na)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      ijh=ijtoh(ih,jh,np)
      DO ipol=1,3
         DO jpol=1,3
            IF (domag) THEN
               int4_nc(ih,jh,ipol,jpol,na,1)=                              &
                   int4(ijh,ipol,jpol,na,1)+int4(ijh,ipol,jpol,na,4)
               int4_nc(ih,jh,ipol,jpol,na,2)=                           &
                 int4(ijh,ipol,jpol,na,2)-(0.d0,1.d0)*int4(ijh,ipol,jpol,na,3)
               int4_nc(ih,jh,ipol,jpol,na,3)=                           &
                 int4(ijh,ipol,jpol,na,2)+(0.d0,1.d0)*int4(ijh,ipol,jpol,na,3)
               int4_nc(ih,jh,ipol,jpol,na,4)=                           &
                      int4(ijh,ipol,jpol,na,1)-int4(ijh,ipol,jpol,na,4)
            ELSE
               int4_nc(ih,jh,ipol,jpol,na,1)= int4(ijh,ipol,jpol,na,1)
               int4_nc(ih,jh,ipol,jpol,na,4)= int4(ijh,ipol,jpol,na,1)
            END IF
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE transform_int4_nc

!----------------------------------------------------------------------------
SUBROUTINE transform_int5_nc(int5, nb)
!----------------------------------------------------------------------------
!
! This routine sets int5_so for the atomic species which do not
! have a spin-orbit pseudopotential
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp,                 ONLY : ijtoh
USE uspp_param,           ONLY : nh, nhm
USE phus,                 ONLY : int5_so
!
IMPLICIT NONE
INTEGER :: nb
COMPLEX(DP) :: int5(nhm*(nhm+1)/2,3,3,nat,nat)
!
! ... local variables
!
INTEGER :: ih, jh, np, na, ipol, jpol
INTEGER :: ijh

np=ityp(nb)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      ijh=ijtoh(ih,jh,np)
      DO na=1,nat
         DO ipol=1,3
            DO jpol=1,3
               int5_so(ih,jh,ipol,jpol,na,nb,1)=int5(ijh,ipol,jpol,na,nb)
               int5_so(ih,jh,ipol,jpol,na,nb,4)=int5(ijh,ipol,jpol,na,nb)
            END DO
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE transform_int5_nc
