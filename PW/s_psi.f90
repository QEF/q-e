!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine s_psi (lda, n, m, psi, spsi )
  !-----------------------------------------------------------------------
  !
  !    This routine applies the S matrix to m wavefunctions psi
  !    and puts the results in spsi.
  !    Requires the products of psi with all beta functions
  !    in array becp(nkb,m) (calculated by ccalbec)
  ! input:
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     spsi  S*psi
  !
#include "machine.h"
  use pwcom
  use becmod
  implicit none
  !
  !     First the dummy variables
  !
  integer :: lda, n, m
  complex(kind=DP) :: psi (lda, m), spsi (lda, m)
  !
  !    here the local variables
  !
  integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
  ! counters
  complex(kind=DP), allocatable :: ps (:,:)
  ! the product vkb and psi
  call start_clock ('s_psi')
  !
  !   initialize  spsi
  !
  call ZCOPY (lda * m, psi, 1, spsi, 1)
  !
  !  The product with the beta functions
  !
  if (nkb.eq.0.or..not.okvan) goto 10
  !
  allocate (ps(nkb,m))    
  ps(:,:) = (0.d0,0.d0)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do ibnd = 1, m
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       ps(ikb,ibnd)=ps(ikb,ibnd) + qq(ih,jh,nt)*becp(jkb,ibnd)
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
  call ZGEMM ('N', 'N', n, m, nkb, (1.d0, 0.d0) , vkb, &
       lda, ps, nkb, (1.d0, 0.d0) , spsi, lda)

  deallocate(ps)

10 call stop_clock ('s_psi')
  return
end subroutine s_psi

