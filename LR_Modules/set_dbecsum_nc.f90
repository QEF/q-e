!
! Copyright (C) 2007-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!----------------------------------------------------------------------------
SUBROUTINE set_dbecsum_nc(dbecsum_nc, dbecsum, npe)
!----------------------------------------------------------------------------
USE kinds, ONLY : DP
USE ions_base, ONLY : nat, ntyp => nsp, ityp
USE uspp_param, only: upf, nhm
USE noncollin_module, ONLY : nspin_mag
USE lsda_mod, ONLY : nspin
IMPLICIT NONE
INTEGER :: npe
INTEGER :: np, na
COMPLEX(DP), INTENT(IN) :: dbecsum_nc( nhm, nhm, nat, nspin, npe)
COMPLEX(DP), INTENT(OUT) :: dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, npe)

DO np = 1, ntyp
   IF ( upf(np)%tvanp ) THEN
      DO na = 1, nat
         IF (ityp(na)==np) THEN
            IF (upf(np)%has_so) THEN
               CALL transform_dbecsum_so(dbecsum_nc,dbecsum,na, npe)
            ELSE
               CALL transform_dbecsum_nc(dbecsum_nc,dbecsum,na, npe)
            END IF
         END IF
      END DO
   END IF
END DO

RETURN
END SUBROUTINE set_dbecsum_nc
!
!----------------------------------------------------------------------------
SUBROUTINE transform_dbecsum_so(dbecsum_nc,dbecsum,na,modes)
!----------------------------------------------------------------------------
!
! This routine multiply dbecsum_nc by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in dbecsum to use it in the calculation of
! the charge and magnetization.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE uspp_param,           ONLY : nh, nhm
USE lsda_mod,             ONLY : nspin
USE uspp,                 ONLY : ijtoh
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
!
IMPLICIT NONE

INTEGER :: na, modes
COMPLEX(DP) :: dbecsum_nc( nhm, nhm, nat, nspin, modes)
COMPLEX(DP) :: dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, modes)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ijh, np, is1, is2, ijs, mode
COMPLEX(DP) :: fac
LOGICAL :: same_lj

np=ityp(na)
DO mode=1,modes
   DO ih = 1, nh(np)
      DO kh = 1, nh(np)
         IF (same_lj(kh,ih,np)) THEN
            DO jh = 1, nh(np)
               ijh=ijtoh(ih,jh,np)
               DO lh=1,nh(np)
                  IF (same_lj(lh,jh,np)) THEN
                     ijs=0
                     DO is1=1,npol
                        DO is2=1,npol
                           ijs=ijs+1
                           fac=dbecsum_nc(kh,lh,na,ijs,mode)
                           dbecsum(ijh,na,1,mode)=dbecsum(ijh,na,1,mode)+fac* &
                              (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) + &
                               fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                           IF (domag) THEN
                              dbecsum(ijh,na,2,mode)=dbecsum(ijh,na,2,mode)+ &
                                 fac * &
                                (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np)+&
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                              dbecsum(ijh,na,3,mode)=dbecsum(ijh,na,3,mode)+ &
                                               fac*(0.d0,-1.d0)*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                              dbecsum(ijh,na,4,mode)=dbecsum(ijh,na,4,mode) &
                                + fac *     &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         END IF
      END DO
   END DO
END DO
RETURN
END SUBROUTINE transform_dbecsum_so
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
!
