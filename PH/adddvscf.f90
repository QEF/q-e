!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine adddvscf (ipert, ik)
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !
#include "f_defs.h"

  USE kinds,      ONLY : DP
  USE uspp_param, ONLY : nh, tvanp
  USE uspp,       ONLY : vkb, okvan
! modules from pwcom
  USE lsda_mod,   ONLY : lsda, current_spin, isk
  USE ions_base,  ONLY : ntyp => nsp, nat, ityp
  USE wvfct,      ONLY : nbnd, npwx
  USE noncollin_module, ONLY : noncolin, npol
! modules from phcom
  USE control_ph, ONLY : lgamma
  USE qpoint,     ONLY : npwq
  USE phus,       ONLY : int3, int3_nc, becp1, becp1_nc
  USE eqv,        ONLY : dvpsi
  implicit none
  !
  !   The dummy variables
  !
  integer :: ik, ipert
  ! input: the k point
  ! input: the perturbation
  !
  !   And the local variables
  !
  integer :: na, nt, ibnd, ih, jh, ijkb0, ikk, ikb, jkb, ip
  ! counter on atoms
  ! counter on atomic types
  ! counter on bands
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for indexing
  ! counter on the k points
  ! counter on vkb
  ! counter on vkb
  complex(DP) :: sum, sum_nc(npol)
  ! auxiliary variable

  if (.not.okvan) return
  call start_clock ('adddvscf')
  if (lgamma) then
     ikk = ik
  else
     ikk = 2 * ik - 1
  endif
  if (lsda) current_spin = isk (ikk)
  ijkb0 = 0
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              !
              !   we multiply the integral for the becp term and the beta_n
              !
              do ibnd = 1, nbnd
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    IF (noncolin) THEN
                       sum_nc = (0.d0, 0.d0)
                    ELSE
                       sum = (0.d0, 0.d0)
                    END IF
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       IF (noncolin) THEN
                          sum_nc(1)=sum_nc(1)+                     &
                                     int3_nc(ih,jh,ipert,na,1)*    &
                                     becp1_nc (jkb, 1, ibnd, ik)+  &
                                     int3_nc(ih,jh,ipert,na,2)*    &
                                     becp1_nc (jkb, 2, ibnd, ik)
                          sum_nc(2)=sum_nc(2)+                     &
                                     int3_nc(ih,jh,ipert,na,3)*    &
                                     becp1_nc (jkb, 1, ibnd, ik)+  &
                                     int3_nc(ih,jh,ipert,na,4)*    &
                                     becp1_nc (jkb, 2, ibnd, ik)
                       ELSE
                          sum = sum + int3 (ih, jh, ipert, na, current_spin)*&
                                   becp1 (jkb, ibnd, ik)
                       END IF
                    enddo
                    IF (noncolin) THEN
                       call ZAXPY(npwq,sum_nc(1),vkb(1,ikb),1,dvpsi(1,ibnd),1)
                       call ZAXPY(npwq,sum_nc(2),vkb(1,ikb),1, &
                                                 dvpsi(1+npwx,ibnd),1)
                    ELSE
                       call ZAXPY(npwq,sum,vkb(1,ikb),1,dvpsi(1,ibnd),1)
                    END IF
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     else
        do na = 1, nat
           if (ityp (na) .eq.nt) ijkb0 = ijkb0 + nh (nt)
        enddo
     endif
  enddo

  call stop_clock ('adddvscf')
  return
end subroutine adddvscf
