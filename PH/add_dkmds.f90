!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
subroutine add_dkmds(kpoint, uact, jpol, dvkb)
  !--------=========-------------------------------------------------------
  !
  ! This subroutine adds to dvpsi the terms which depend on the augmentation
  ! charge. It assume that the variable dpqq, has been set.
  !
#include "f_defs.h"

  use pwcom
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE kinds, only : DP
  USE wavefunctions_module,    ONLY : evc
  USE uspp_param, only: nh
  use phcom

  implicit none

  integer, intent(in) :: kpoint, jpol
  complex(kind=DP), intent(in) :: uact (3 * nat)
  complex(kind=DP), intent(in) :: dvkb (npwx,nkb,3)

  complex(kind=dp) :: becp2(nkb,nbnd)

  real(kind=DP), parameter :: eps = 1.d-12

  integer :: ipol, ijkb0, nt, na, ih, jh, ikb, jkb, ibnd, ig, igg, mu

  real(kind=dp) :: fact

  logical :: ok

  complex(kind=dp), allocatable :: ps1(:,:), ps2(:,:,:) 
  complex(kind=DP), allocatable :: alphadk(:,:,:)
  complex(kind=DP), allocatable :: aux(:), aux1(:,:)


  integer :: i,j 

#ifdef TIMING_ADD_DKMDS
  call start_clock('add_dkmds')
  call start_clock('add_dkmds2')
#endif 
  allocate(aux(npwx))
  allocate(aux1(npwx,nbnd))
  if (nkb.gt.0) then 
     allocate (ps1(nkb,nbnd))
     allocate (ps2(nkb,3,nbnd))
     allocate (alphadk(nkb,nbnd,3))
  end if

  ps1 = (0.d0, 0.d0)
  ps2 = (0.d0, 0.d0)
  !
  !   First we calculate the alphadk = <d/dk d/du beta|psi>
  !   and becp2 = < d/dk beta | psi>
  !
  if (lsda) current_spin = isk (kpoint)
  call ccalbec (nkb, npwx, npw, nbnd, becp2, dvkb(1,1,jpol), evc)

#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds2')
  call start_clock('add_dkmds3')
#endif

  do ipol = 1, 3
     do ibnd = 1, nbnd
        do ig = 1, npw
           aux1 (ig, ibnd) = evc(ig,ibnd) * tpiba * (0.d0,1.d0) * & 
                ( xk(ipol,kpoint) + g(ipol,igk(ig)) )
        enddo
     enddo
     call ccalbec (nkb, npwx, npw, nbnd, alphadk(1,1,ipol), dvkb(1,1,jpol), aux1)
  enddo
#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds3')
  call start_clock('add_dkmds4')
#endif

  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp(na).eq.nt) then
           mu = 3 * (na - 1)
           if ( abs (uact (mu + 1) ) + &
                abs (uact (mu + 2) ) + &
                abs (uact (mu + 3) ) > eps) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    fact = at(1,jpol)*dpqq(ih,jh,1,nt) +  &
                         at(2,jpol)*dpqq(ih,jh,2,nt) +  &
                         at(3,jpol)*dpqq(ih,jh,3,nt)
                    do ipol = 1, 3 
                       do ibnd=1, nbnd_occ(kpoint)
                          !
                          ! first we calculate the part coming from the 
                          ! overlapp matrix S
                          !
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +           &
                               (0.d0,1.d0) * qq (ih, jh, nt) *          &
                               alphadk(jkb, ibnd, ipol) *               &
                               uact (mu + ipol)
                          ps2 (ikb, ipol, ibnd) = ps2 (ikb, ipol, ibnd) +  &
                               qq (ih, jh, nt) *                           &
                               becp2(jkb, ibnd) *                          &
                               uact (mu + ipol) * tpiba
                          ! 
                          !  and here the part of the matrix K(r)
                          !
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +      &
                               fact *                              &
                               alphap(jkb, ibnd, ipol, kpoint) *   &
                               uact (mu + ipol)
                          ps2 (ikb, ipol, ibnd) = ps2 (ikb, ipol, ibnd) + &
                               fact  *  (0.d0,-1.d0) *                    &
                               becp1(jkb, ibnd, kpoint) *                 &
                               uact (mu + ipol) * tpiba
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0=ijkb0+nh(nt)
        endif
     enddo
  enddo
#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds4')
  call start_clock('add_dkmds5')
#endif
  !
  !      This term is proportional to beta(k+q+G)
  !
  dvpsi = matmul(vkb, ps1) + dvpsi 
#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds5')
  call start_clock('add_dkmds6')
#endif
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        do ibnd = 1, nbnd
           ok = ok.or. (abs (ps2 (ikb, ipol, ibnd)).gt.eps )
        enddo
        if (ok) then
           do ig = 1, npw
              igg = igkq (ig)
              aux (ig) =  vkb(ig, ikb) * (xk(ipol, kpoint) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              dvpsi(1:npw,ibnd) =                       &
                   ps2(ikb,ipol,ibnd) * aux(1:npw) +    &
                   dvpsi(1:npwq,ibnd)
           enddo
        endif
     enddo
  enddo

  deallocate (aux)
  deallocate(aux1)
  if (allocated(ps1))     deallocate(ps1)
  if (allocated(ps2))     deallocate(ps2)
  if (allocated(alphadk)) deallocate (alphadk)

#ifdef TIMING_ADD_DKMDS
  call stop_clock('add_dkmds6')
  call stop_clock('add_dkmds')
#endif
  return

end subroutine add_dkmds
