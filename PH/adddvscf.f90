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
#include "machine.h"

  use pwcom
  use parameters, only : DP
  use phcom
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
  integer :: na, nt, ibnd, ih, jh, ijkb0, ikk, ikb, jkb
  ! counter on atoms
  ! counter on atomic types
  ! counter on bands
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for indexing
  ! counter on the k points
  ! counter on vkb
  ! counter on vkb
  complex(kind=DP) :: sum
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
                    sum = (0.d0, 0.d0)
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       sum = sum + int3 (ih, jh, ipert, na, current_spin) * &
                                   becp1 (jkb, ibnd, ik)
                    enddo
                    !                        sum=sum*eigqts(na)
                    call ZAXPY (npwq, sum, vkb (1, ikb), 1, dvpsi (1, ibnd),1)
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
