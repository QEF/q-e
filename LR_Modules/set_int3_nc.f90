!
! Copyright (C) 2007-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!----------------------------------------------------------------------------
SUBROUTINE set_int3_nc(npe)
!----------------------------------------------------------------------------
USE ions_base, ONLY : nat, ntyp => nsp, ityp
USE uspp_param, only: upf
USE lrus, ONLY : int3, int3_nc
IMPLICIT NONE
INTEGER :: npe
INTEGER :: np, na

int3_nc=(0.d0,0.d0)
DO np = 1, ntyp
   IF ( upf(np)%tvanp ) THEN
      DO na = 1, nat
         IF (ityp(na)==np) THEN
            IF (upf(np)%has_so) THEN
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
!----------------------------------------------------------------------------
SUBROUTINE transform_int3_so(int3,na,npert)
!----------------------------------------------------------------------------
!
! This routine multiply int3 by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in int3_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
USE lrus,                 ONLY : int3_nc
!
IMPLICIT NONE

INTEGER :: na, npert
COMPLEX(DP) :: int3(nhm,nhm,nat,nspin_mag,npert)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ipol, np, is1, is2, ijs
LOGICAL :: same_lj

np=ityp(na)
DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  DO ipol=1,npert
                     ijs=0
                     DO is1=1,npol
                        DO is2=1,npol
                           ijs=ijs+1
                           int3_nc(ih,jh,na,ijs,ipol)=                       &
                               int3_nc(ih,jh,na,ijs,ipol) +                  &
                               int3 (kh,lh,na,1,ipol)*                       &
                             (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)  + &
                             fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)   )
                           IF (domag) THEN
                              int3_nc(ih,jh,na,ijs,ipol)=                     &
                                 int3_nc(ih,jh,na,ijs,ipol) +                 &
                                 int3(kh,lh,na,2,ipol)*                       &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)+ &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 (0.D0,-1.D0) * int3(kh,lh,na,3,ipol)*        &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 int3 (kh,lh,na,4,ipol)*                      &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np))
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
       !
RETURN
END SUBROUTINE transform_int3_so
!
!----------------------------------------------------------------------------
SUBROUTINE transform_int3_nc(int3,na,npert)
!----------------------------------------------------------------------------
!
! This routine multiply int3 by the identity and the Pauli
! matrices and saves it in int3_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : nspin_mag
USE spin_orb,             ONLY : domag
USE lrus,                 ONLY : int3_nc
!
IMPLICIT NONE

INTEGER :: na, npert
COMPLEX(DP) :: int3(nhm,nhm,nat,nspin_mag,npert)
!
! ... local variables
!
INTEGER :: ih, jh, ipol, np

np=ityp(na)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      DO ipol=1,npert
         IF (domag) THEN
            int3_nc(ih,jh,na,1,ipol)=int3(ih,jh,na,1,ipol)+int3(ih,jh,na,4,ipol)
            int3_nc(ih,jh,na,2,ipol)=                                       &
               int3(ih,jh,na,2,ipol) - (0.d0, 1.d0) * int3(ih,jh,na,3,ipol)
            int3_nc(ih,jh,na,3,ipol)=                                       &
               int3(ih,jh,na,2,ipol) + (0.d0, 1.d0) * int3(ih,jh,na,3,ipol)
            int3_nc(ih,jh,na,4,ipol)=                                       &
               int3(ih,jh,na,1,ipol) - int3(ih,jh,na,4,ipol)
         ELSE
            int3_nc(ih,jh,na,1,ipol)=int3(ih,jh,na,1,ipol)
            int3_nc(ih,jh,na,4,ipol)=int3(ih,jh,na,1,ipol)
         END IF
      END DO
   END DO
END DO

RETURN
END SUBROUTINE transform_int3_nc
!
