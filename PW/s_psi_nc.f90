!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine s_psi_nc (lda, n, m, psi, spsi )
  !-----------------------------------------------------------------------
  !
  !    This routine applies the S matrix to m wavefunctions psi
  !    and puts the results in spsi.
  !    Requires the products of psi with all beta functions
  !    in array becp_nc(nkb,m) (calculated by ccalbec)
  ! input:
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     spsi  S*psi
  !
  USE ions_base, ONLY: nat, ityp, ntyp => nsp
  USE uspp_param, ONLY: upf, nh
  USE uspp, ONLY: nkb, vkb, qq, qq_so, okvan
  use wvfct, ONLY: igk, g2kin
  use gsmooth, ONLY: nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  use ldaU, ONLY: lda_plus_u
  use becmod
  use noncollin_module, ONLY: npol
  use spin_orb, ONLY: lspinorb
  implicit none
  !
  !     First the dummy variables
  !
  integer :: lda, n, m
  complex(DP) :: psi (lda, npol, m), spsi (lda,npol, m)
  complex(DP), external ::  ZDOTU
  !
  !    here the local variables
  !
  integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ipol
  ! counters
  complex(DP), allocatable :: ps (:,:,:)
  ! the product vkb and psi
  call start_clock ('s_psi')
  !
  !   initialize  spsi
  !
  call ZCOPY (lda * m *npol, psi(1,1,1), 1, spsi(1,1,1), 1)
  !
  !  The product with the beta functions
  !
  if (nkb.eq.0.or..not.okvan) goto 10
  !
  allocate (ps(nkb,npol,m))    
  ps(:,:,:) = (0.D0,0.D0)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     if ( upf(nt)%tvanp ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do ih = 1,nh(nt)
                 ikb = ijkb0 + ih
                 do ibnd = 1, m
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       if (lspinorb) then
                          ps(ikb,1,ibnd)=ps(ikb,1,ibnd) + &
                                 qq_so(ih,jh,1,nt)*becp_nc(jkb,1,ibnd)+ &
                                 qq_so(ih,jh,2,nt)*becp_nc(jkb,2,ibnd)
                          ps(ikb,2,ibnd)=ps(ikb,2,ibnd) + &
                                 qq_so(ih,jh,3,nt)*becp_nc(jkb,1,ibnd)+ &
                                 qq_so(ih,jh,4,nt)*becp_nc(jkb,2,ibnd)
                       else
                          do ipol=1,npol
                             ps(ikb,ipol,ibnd)=ps(ikb,ipol,ibnd) + &
                                        qq(ih,jh,nt)*becp_nc(jkb,ipol,ibnd)
                          enddo
                       endif
                    enddo
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
  call ZGEMM ('N', 'N', n, m*npol, nkb, (1.d0, 0.d0) , vkb, &
       lda, ps, nkb, (1.d0, 0.d0) , spsi(1,1,1), lda)

  DEALLOCATE(ps)

10 call stop_clock ('s_psi')
  return
end subroutine s_psi_nc
