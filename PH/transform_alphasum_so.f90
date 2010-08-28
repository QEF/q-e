!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE transform_alphasum_so(alphasum_nc,na)
!----------------------------------------------------------------------------
!
! This routine multiply alphasum_nc by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in alphasum to use it in the calculation of
! the change of the charge and of the magnetization.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol
USE spin_orb,             ONLY : fcoef, domag
USE phus,                 ONLY : alphasum
!
IMPLICIT NONE

COMPLEX(DP) :: alphasum_nc(nhm*(nhm+1)/2,3,nat,npol,npol)
INTEGER :: na
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ijh, np, is1, is2, ipol
INTEGER, ALLOCATABLE :: ijh_save(:,:)
COMPLEX(DP) :: fac
INTEGER :: find_ijh, ijh_l
LOGICAL :: same_lj

ALLOCATE(ijh_save(nhm,nhm))
np=ityp(na)
DO ih=1,nh(np)
   DO jh=1,nh(np)
      ijh_save(ih,jh)=find_ijh(ih,jh,nh(np))
   END DO
END DO
DO ipol=1,3
   DO ih = 1, nh(np)
      DO kh = 1, nh(np)
         IF (same_lj(kh,ih,np)) THEN
            DO jh = 1, nh(np)
               ijh=ijh_save(ih,jh)
               DO lh=1,nh(np)
                  IF (same_lj(lh,jh,np)) THEN
                     ijh_l=ijh_save(kh,lh)
                     DO is1=1,npol
                        DO is2=1,npol
                           IF (kh <= lh) THEN
                              fac=alphasum_nc(ijh_l,ipol,na,is1,is2)
                           ELSE
                              fac=CONJG(alphasum_nc(ijh_l,ipol,na,is2,is1))
                           ENDIF
                           alphasum(ijh,ipol,na,1)=alphasum(ijh,ipol,na,1)+fac*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) + &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                           IF (domag) THEN
                              alphasum(ijh,ipol,na,2)=alphasum(ijh,ipol,na,2)+&
                                fac*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) +&
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                              alphasum(ijh,ipol,na,3)=alphasum(ijh,ipol,na,3)+&
                                fac*(0.d0,-1.d0)*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                              alphasum(ijh,ipol,na,4)=alphasum(ijh,ipol,na,4) +&
                                fac * &
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
       !
DEALLOCATE(ijh_save)
RETURN
END SUBROUTINE transform_alphasum_so

