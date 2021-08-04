!
! Copyright (C) 2016 Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!----------------------------------------------------------------------------
SUBROUTINE set_intq_nc()
!----------------------------------------------------------------------------
USE ions_base,  ONLY : nat, ntyp => nsp, ityp
USE uspp_param, ONLY : upf
USE lrus,       ONLY : intq, intq_nc
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
!
!----------------------------------------------------------------------------
SUBROUTINE transform_intq_so(intq,na)
!----------------------------------------------------------------------------
!
! This routine multiply intq by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in intq_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
USE lrus,                 ONLY : intq_nc
!
IMPLICIT NONE

COMPLEX(DP) :: intq(nhm,nhm,nat)
INTEGER :: na
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, np, npert, is1, is2, ijs
LOGICAL :: same_lj

np=ityp(na)
DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  ijs=0
                  DO is1=1,npol
                     DO is2=1,npol
                        ijs=ijs+1
                        intq_nc(ih,jh,na,ijs)=                           &
                            intq_nc(ih,jh,na,ijs) + intq (kh,lh,na)*     &
                          (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np) + &
                          fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)  )
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO
       !
RETURN
END SUBROUTINE transform_intq_so
!
!
!----------------------------------------------------------------------------
SUBROUTINE transform_intq_nc(intq,na)
!----------------------------------------------------------------------------
!
! This routine multiply intq by the identity and the Pauli
! matrices and saves it in intq_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE lrus,                 ONLY : intq_nc
!
IMPLICIT NONE

INTEGER :: na
COMPLEX(DP) :: intq(nhm,nhm,nat)
!
! ... local variables
!
INTEGER :: ih, jh, np

np=ityp(na)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      intq_nc(ih,jh,na,1)=intq(ih,jh,na)
      intq_nc(ih,jh,na,4)=intq(ih,jh,na)
   END DO
END DO

RETURN
END SUBROUTINE transform_intq_nc

