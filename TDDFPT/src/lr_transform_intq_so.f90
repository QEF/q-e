!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_transform_intq_so(intq,na)
  !--------------------------------------------------------------------------
  !
  ! This routine multiply intq by the identity and the Pauli
  ! matrices, rotate it as appropriate for the spin-orbit case
  ! and saves it in intq_nc.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE noncollin_module,     ONLY : npol
  USE spin_orb,             ONLY : fcoef
!  USE lr_variables,         ONLY : intq_nc
  USE lrus,                 ONLY : intq_nc
 
  IMPLICIT NONE
  INTEGER :: na
  COMPLEX(DP) :: intq(nhm,nhm,nat)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, lh, kh, np, is1, is2, ijs
  LOGICAL :: same_lj ! from PW
  !
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
                          intq_nc(ih,jh,na,ijs) =                            &
                              intq_nc(ih,jh,na,ijs) + intq (kh,lh,na) *      &
                              (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np) + &
                               fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np))
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
  !
END SUBROUTINE lr_transform_intq_so
