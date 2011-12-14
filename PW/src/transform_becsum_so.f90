!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE transform_becsum_so(becsum_nc,becsum,na)
!----------------------------------------------------------------------------
!
! This routine multiply becsum_nc by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in becsum for the calculation of 
! augmentation charge and magnetization.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ntyp => nsp, ityp
USE uspp_param,           ONLY : nh, nhm
USE lsda_mod,             ONLY : nspin
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
!
IMPLICIT NONE

COMPLEX(DP) :: becsum_nc(nhm*(nhm+1)/2,nat,npol,npol)
REAL(DP) :: becsum(nhm*(nhm+1)/2,nat,nspin_mag)
INTEGER :: na
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ijh, np, is1, is2
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
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      ijh=ijh_save(ih,jh)
      DO kh = 1, nh(np)
         IF (same_lj(kh,ih,np)) THEN
            DO lh=1,nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  ijh_l=ijh_save(kh,lh)
                  DO is1=1,npol
                     DO is2=1,npol
                        IF (kh <= lh) THEN
                           fac=becsum_nc(ijh_l,na,is1,is2)
                        ELSE
                           fac=CONJG(becsum_nc(ijh_l,na,is2,is1))
                        ENDIF
                        becsum(ijh,na,1)=becsum(ijh,na,1) + fac * &
                            (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) + &
                             fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                        IF (domag) THEN
                           becsum(ijh,na,2)=becsum(ijh,na,2)+fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) +&
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                           becsum(ijh,na,3)=becsum(ijh,na,3)+fac*(0.d0,-1.d0)*&
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,2,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,1,is2,np)  )
                           becsum(ijh,na,4)=becsum(ijh,na,4) + fac * &
                               (fcoef(kh,ih,is1,1,np)*fcoef(jh,lh,1,is2,np) - &
                                fcoef(kh,ih,is1,2,np)*fcoef(jh,lh,2,is2,np)  )
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END IF
      END DO
   END DO
END DO
       !
DEALLOCATE(ijh_save)
RETURN
END SUBROUTINE transform_becsum_so

FUNCTION same_lj(ih,jh,np)

USE uspp, ONLY : nhtol, nhtoj, indv

IMPLICIT NONE

LOGICAL :: same_lj
INTEGER :: ih, jh, np

same_lj = ((nhtol(ih,np)==nhtol(jh,np)).AND. &
           (ABS(nhtoj(ih,np)-nhtoj(jh,np))<1.d8).AND. &
           (indv(ih,np)==indv(jh,np)) )

RETURN
END FUNCTION same_lj

FUNCTION find_ijh(ih,jh,nh)
IMPLICIT NONE
INTEGER :: find_ijh, ih, jh, nh
INTEGER :: ih0, jh0, ijh, i, j

if (jh > ih) then
   ih0=ih
   jh0=jh
else
   ih0=jh
   jh0=ih
end if

ijh=0
do i=1, ih0-1
   do j=i, nh
      ijh=ijh+1
   enddo
enddo
do j=ih0, jh0
   ijh=ijh+1
enddo
find_ijh=ijh
end function find_ijh

