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
subroutine add_vuspsi_nc (lda, n, m, psi, hpsi )
  !-----------------------------------------------------------------------
  !
  !    This routine applies the Ultra-Soft Hamiltonian to a
  !    vector psi and puts the result in hpsi.
  !    Requires the products of psi with all beta functions
  !    in array becp_nc(nkb,m) (calculated by calbec)
  ! input:
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     hpsi  V_US*psi is added to hpsi
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE uspp_param, ONLY: nh
  USE uspp,       ONLY : vkb, nkb, deeq_nc
  USE becmod
  USE noncollin_module
  implicit none
  !
  !     First the dummy variables
  !
  integer :: lda, n, m
  complex(DP) :: psi (lda, npol, m), hpsi (lda,npol, m)
  !
  !    here the local variables
  !
  integer :: jkb, ikb, ih, jh, na, nt, ijkb0, ibnd
  ! counters
  complex(DP), allocatable :: ps (:,:,:)
  ! the product vkb and psi
  !

  if (nkb.eq.0) return
  call start_clock ('add_vuspsi')

  allocate (ps(  nkb,npol, m))    
  ps (:,:,:) = (0.d0, 0.d0)
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           do ih = 1, nh (nt)             !!do ibnd = 1, m
              ikb = ijkb0 + ih             !!do jh = 1, nh (nt)
              do ibnd = 1, m                !!   jkb = ijkb0 + jh
                 do jh = 1, nh (nt)  !!do ih = 1, nh (nt)
                    jkb = ijkb0 + jh  !!ikb = ijkb0 + ih
                    ps(ikb,1,ibnd) = ps(ikb,1,ibnd) +    & 
                         deeq_nc(ih,jh,na,1)*becp_nc(jkb,1,ibnd)+ & 
                         deeq_nc(ih,jh,na,2)*becp_nc(jkb,2,ibnd) 
                    ps(ikb,2,ibnd) = ps(ikb,2,ibnd)  +   & 
                         deeq_nc(ih,jh,na,3)*becp_nc(jkb,1,ibnd)+&
                         deeq_nc(ih,jh,na,4)*becp_nc(jkb,2,ibnd) 
                 enddo
              enddo
           enddo
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
  call ZGEMM ('N', 'N', n, m*npol, nkb, (1.d0, 0.d0) , vkb, &
       lda, ps, nkb, (1.d0, 0.d0) , hpsi, lda)

  deallocate (ps)
  call stop_clock ('add_vuspsi')
  return
end subroutine add_vuspsi_nc

