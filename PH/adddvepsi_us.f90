!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine adddvepsi_us(becp2,ipol,kpoint)
  ! This subdoutine adds to dvpsi the terms which depend on the augmentation
  ! charge. It assume that the variable dpqq, has been set.
  !
#include "f_defs.h"

  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE kinds, only : DP
  USE uspp_param, only: nh
  use phcom

  implicit none

  integer, intent(in) :: ipol, kpoint
  complex(DP), intent(in) :: becp2(nkb,nbnd)

  real(DP) :: fact
  complex(DP), allocatable :: ps(:)
  integer:: ijkb0, nt, na, ih, jh, ikb, jkb, ibnd

  allocate (ps(nbnd))    

  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp(na).eq.nt) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              ps = (0.d0,0.d0)
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 fact=at(1,ipol)*dpqq(ih,jh,1,nt)+  &
                      at(2,ipol)*dpqq(ih,jh,2,nt)+  &
                      at(3,ipol)*dpqq(ih,jh,3,nt)
                 do ibnd=1, nbnd_occ(kpoint)
                    ps(ibnd) = ps(ibnd)                             &
                         + becp2(jkb,ibnd)*(0.d0,1.d0)*qq(ih,jh,nt)+  &
                         becp1(jkb,ibnd,kpoint)*fact
                 enddo
              enddo
              do ibnd = 1, nbnd_occ (kpoint)
                 call ZAXPY(npw,ps(ibnd),vkb(1,ikb),1,dvpsi(1,ibnd),1)
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        endif
     enddo
  enddo
  if (jkb.ne.nkb) call errore ('adddvepsi_us', 'unexpected error', 1)

  deallocate(ps)

  return
end subroutine adddvepsi_us
