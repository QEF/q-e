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
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds, only : DP
  USE uspp_param, only: nh
  use phcom

  implicit none

  integer, intent(in) :: ipol, kpoint
  complex(DP), intent(in) :: becp2(nkb,npol,nbnd)

  real(DP) :: fact
  complex(DP), allocatable :: ps(:), ps_nc(:,:), fact_so(:)
  integer:: ijkb0, nt, na, ih, jh, ikb, jkb, ibnd, ip, is

  IF (noncolin) THEN
     allocate (ps_nc(nbnd,npol))    
     allocate (fact_so(nspin))    
  ELSE
     allocate (ps(nbnd))    
  END IF

  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp(na).eq.nt) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              IF (noncolin) THEN
                 ps_nc = (0.d0,0.d0)
              ELSE
                 ps = (0.d0,0.d0)
              END IF
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 IF (lspinorb) THEN
                    do is=1,nspin
                       fact_so(is)=at(1,ipol)*dpqq_so(ih,jh,is,1,nt)+  &
                                   at(2,ipol)*dpqq_so(ih,jh,is,2,nt)+  &
                                   at(3,ipol)*dpqq_so(ih,jh,is,3,nt)
                    enddo
                 ELSE
                    fact=at(1,ipol)*dpqq(ih,jh,1,nt)+  &
                         at(2,ipol)*dpqq(ih,jh,2,nt)+  &
                         at(3,ipol)*dpqq(ih,jh,3,nt)
                 END IF
                 do ibnd=1, nbnd_occ(kpoint)
                    IF (noncolin) THEN
                       DO ip=1,npol
                          IF (lspinorb) THEN
                             ps_nc(ibnd,ip)=ps_nc(ibnd,ip) +          &
                                 (0.d0,1.d0)*(becp2(jkb,1,ibnd)*      &
                                 qq_so(ih,jh,1+(ip-1)*2,nt) +         &
                                 becp2(jkb,2,ibnd) *                  &
                                 qq_so(ih,jh,2+(ip-1)*2,nt) )         &
                               + becp1_nc(jkb,1,ibnd,kpoint)*         &
                                 fact_so(1+(ip-1)*2)                  &
                               + becp1_nc(jkb,2,ibnd,kpoint)*         &
                                 fact_so(2+(ip-1)*2)    
                          ELSE
                             ps_nc(ibnd,ip)=ps_nc(ibnd,ip)+           &
                                 becp2(jkb,ip,ibnd)*(0.d0,1.d0)*   &
                                 qq(ih,jh,nt)+becp1_nc(jkb,ip,ibnd,kpoint)    &
                                                *fact
                          END IF
                       END DO
                    ELSE
                       ps(ibnd) = ps(ibnd)                             &
                           + becp2(jkb,1,ibnd)*(0.d0,1.d0)*qq(ih,jh,nt)+  &
                             becp1(jkb,ibnd,kpoint)*fact

                    END IF
                 enddo
              enddo
              do ibnd = 1, nbnd_occ (kpoint)
                 IF (noncolin) THEN
                    CALL ZAXPY(npw,ps_nc(ibnd,1),vkb(1,ikb),1, &
                                                     dvpsi(1,ibnd),1)
                    CALL ZAXPY(npw,ps_nc(ibnd,2),vkb(1,ikb),1, &
                                                     dvpsi(1+npwx,ibnd),1)
                 ELSE
                    CALL ZAXPY(npw,ps(ibnd),vkb(1,ikb),1,dvpsi(1,ibnd),1)
                 END IF
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        endif
     enddo
  enddo
  if (jkb.ne.nkb) call errore ('adddvepsi_us', 'unexpected error', 1)

  IF (noncolin) THEN
     deallocate(ps_nc)
     deallocate(fact_so)
  ELSE
     deallocate(ps)
  END IF
  

  RETURN
END SUBROUTINE adddvepsi_us
