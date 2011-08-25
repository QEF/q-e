!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE apply_dpot(nrxxs, aux1, dvscfins, current_spin)
  !
  !  This routine applies the change of the self consistent potential to
  !  one wavefunction
  !
  USE kinds, ONLY : DP
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE spin_orb, ONLY : domag

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: current_spin, nrxxs
  COMPLEX(DP), INTENT(IN) :: dvscfins(nrxxs,nspin_mag)
  COMPLEX(DP), INTENT(INOUT) :: aux1(nrxxs,npol)

  COMPLEX(DP) :: sup, sdwn
  INTEGER :: ir

  IF (noncolin) THEN
     IF (domag) then
        DO ir = 1, nrxxs
           sup=aux1(ir,1)*(dvscfins(ir,1)+dvscfins(ir,4))+ &
                aux1(ir,2)*(dvscfins(ir,2)-(0.d0,1.d0)*dvscfins(ir,3))
           sdwn=aux1(ir,2)*(dvscfins(ir,1)-dvscfins(ir,4)) + &
                aux1(ir,1)*(dvscfins(ir,2)+(0.d0,1.d0)*dvscfins(ir,3))
           aux1(ir,1)=sup
           aux1(ir,2)=sdwn
        ENDDO
     ELSE
        DO ir = 1, nrxxs
           aux1(ir,:)=aux1(ir,:)*dvscfins(ir,1)
        ENDDO
     ENDIF
  ELSE
     DO ir = 1, nrxxs
        aux1(ir,1)=aux1(ir,1)*dvscfins(ir,current_spin)
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE apply_dpot
