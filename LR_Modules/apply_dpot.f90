!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE apply_dpot(nrxxs, aux1, dv, current_spin)
  !
  !  This routine applies the change of the self consistent potential to
  !  one wavefunction
  !
  USE kinds, ONLY : DP
  USE noncollin_module, ONLY : noncolin, npol, nspin_mag
  USE spin_orb, ONLY : domag
  USE mp_bands, ONLY : me_bgrp
  USE fft_base,  ONLY : dffts, dtgs

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: current_spin, nrxxs
  COMPLEX(DP), INTENT(IN) :: dv(nrxxs,nspin_mag)
  COMPLEX(DP), INTENT(INOUT) :: aux1(nrxxs,npol)

  COMPLEX(DP) :: sup, sdwn
  INTEGER :: ir

  IF (noncolin) THEN
!
!    Noncollinear part with task groups
!
     IF( dtgs%have_task_groups ) THEN
        IF (domag) THEN
           DO ir=1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
              sup = aux1(ir,1) * (dv(ir,1)+dv(ir,4)) + &
                    aux1(ir,2) * (dv(ir,2)-(0.d0,1.d0)*dv(ir,3))
              sdwn = aux1(ir,2) * (dv(ir,1)-dv(ir,4)) + &
                     aux1(ir,1) * (dv(ir,2)+(0.d0,1.d0)*dv(ir,3))
              aux1(ir,1)=sup
              aux1(ir,2)=sdwn
           ENDDO
        ELSE
           DO ir=1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
              aux1(ir,:) = aux1(ir,:) * dv(ir,1)
           ENDDO
        ENDIF
     ELSE
!
!   Noncollinear part without TG
!
        IF (domag) then
           DO ir = 1, nrxxs
              sup=aux1(ir,1)*(dv(ir,1)+dv(ir,4))+ &
                  aux1(ir,2)*(dv(ir,2)-(0.d0,1.d0)*dv(ir,3))
              sdwn=aux1(ir,2)*(dv(ir,1)-dv(ir,4)) + &
                   aux1(ir,1)*(dv(ir,2)+(0.d0,1.d0)*dv(ir,3))
              aux1(ir,1)=sup
              aux1(ir,2)=sdwn
           ENDDO
        ELSE
           DO ir = 1, nrxxs
              aux1(ir,:)=aux1(ir,:)*dv(ir,1)
           ENDDO
        ENDIF
     ENDIF
  ELSE
!
!  collinear part with Task Groups
!
     IF( dtgs%have_task_groups ) THEN
        !
        DO ir = 1, dffts%nr1x*dffts%nr2x*dtgs%tg_npp( me_bgrp + 1 )
           aux1 (ir,1) = aux1 (ir,1) * dv(ir,1)
        ENDDO
     ELSE
!
!  collinear part with Task Groups
!
        DO ir = 1, nrxxs
           aux1(ir,1)=aux1(ir,1)*dv(ir,current_spin)
        ENDDO
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE apply_dpot
