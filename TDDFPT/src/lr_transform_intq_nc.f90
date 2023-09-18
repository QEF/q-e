!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_transform_intq_nc(intq,na)
  !----------------------------------------------------------------------------
  !
  ! This routine multiply intq by the identity and the Pauli
  ! matrices and saves it in intq_nc.
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ityp
  USE uspp_param,       ONLY : nh, nhm
!  USE lr_variables,     ONLY : intq_nc
  USE lrus,             ONLY : intq_nc

  IMPLICIT NONE
  INTEGER :: na
  COMPLEX(DP) :: intq(nhm,nhm,nat)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, np
  !
  np = ityp(na)
  DO ih = 1, nh(np)
     DO jh = 1, nh(np)
        intq_nc(ih,jh,na,1) = intq(ih,jh,na)
        intq_nc(ih,jh,na,4) = intq(ih,jh,na)
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE lr_transform_intq_nc
