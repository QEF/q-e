!--------------------------------------------------------------------------
subroutine add_dkmds(kpoint, uact, jpol)
!----------=========-------------------------------------------------------
!
!
! This subdoutine adds to dvpsi the terms which depend on the augmentation
! charge. It assume that the variable dpqq, has been set.
!
#include "machine.h"

  use pwcom
  use parameters, only : DP
  USE wavefunctions,    ONLY : evc
  use becmod
  use phcom
  
  implicit none

  integer, intent(in) :: kpoint, jpol

  complex(kind=DP), intent(in) :: uact (3 * nat)

  complex(kind=dp) :: becp2(nkb,nbnd)

  real(kind=DP), parameter :: eps = 1.d-12

  integer :: ipol, ijkb0, nt, na, ih, jh, ikb, jkb, ibnd, ig, igg, mu

  real(kind=dp) :: fact

  logical :: ok

  real(kind=DP), allocatable  :: gk (:,:)

  complex(kind=dp), allocatable :: ps1(:,:), ps2(:,:,:), ps3(:,:), ps4(:,:,:)
  complex(kind=DP), allocatable :: dvkb (:,:), dvkb1 (:,:), work (:,:)
  complex(kind=DP), allocatable :: alphadk(:,:,:)
  complex(kind=DP), allocatable :: aux(:), aux1(:,:)


  integer :: i,j 
  
  if (nkb.gt.0) then
     allocate (work ( npwx, nkb))    
  else
     allocate (work ( npwx, 1))    
  endif
  allocate(aux(npwx))
  allocate(aux1(npwx,nbnd))
  allocate (gk(3, npwx))
  if (nkb.gt.0) then 
     allocate (dvkb ( npwx , nkb))    
     allocate (dvkb1( npwx , nkb))
     allocate (ps1(nkb,nbnd))
     allocate (ps2(nkb,3,nbnd))
     allocate (ps3(nkb,nbnd))
     allocate (ps4(nkb,3,nbnd))
     allocate (alphadk(nkb,nbnd,3))
  end if

  ps1 = (0.d0, 0.d0)
  ps2 = (0.d0, 0.d0)
  ps3 = (0.d0, 0.d0)
  ps4 = (0.d0, 0.d0)
!
!   First we calculate the alphadk = <d/dk d/du beta|psi>
!   and becp2 = < d/dk beta | psi>
!
  do ig = 1, npw
     gk (1, ig) = (xk (1, kpoint) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, kpoint) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, kpoint) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2
  enddo

  if (lsda) current_spin = isk (kpoint)
  call gen_us_dj (kpoint, dvkb)
  call gen_us_dy (kpoint, at (1, jpol), dvkb1)
  do ig = 1, npw
     if (g2kin (ig) .lt.1.0d-10) then
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     else
        gk (1, ig) = gk (1, ig) / sqrt (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / sqrt (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / sqrt (g2kin (ig) )
     endif
  enddo
     
  jkb = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (nt.eq.ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
              do ig = 1, npw
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                      (at (1, jpol) * gk (1, ig) + &
                      at (2, jpol) * gk (2, ig) + &
                      at (3, jpol) * gk (3, ig) )
              enddo
           enddo
        endif
     enddo
  enddo
  call ccalbec (nkb, npwx, npw, nbnd, becp2, work, evc)

  do ipol = 1, 3
     do ibnd = 1, nbnd
        do ig = 1, npw
              aux1 (ig, ibnd) = evc(ig,ibnd) * tpiba * (0.d0,1.d0) * & 
                   ( xk(ipol,kpoint) + g(ipol,igk(ig)) )
        enddo
     enddo
     call ccalbec (nkb, npwx, npw, nbnd, alphadk(1,1,ipol), work, aux1)
  enddo


  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp(na).eq.nt) then
           mu = 3 * (na - 1)
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ipol = 1, 3 
                    fact=at(1,jpol)*dpqq(ih,jh,1,nt)+  &
                         at(2,jpol)*dpqq(ih,jh,2,nt)+  &
                         at(3,jpol)*dpqq(ih,jh,3,nt)
                    do ibnd=1, nbnd_occ(kpoint)
                       if ( abs (uact (mu + 1) ) + &
                            abs (uact (mu + 2) ) + &
                            abs (uact (mu + 3) ) > eps) then
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
                          ps3 (ikb, ibnd) = ps3 (ikb, ibnd) +      &
                               fact *                              &
                               alphap(jkb, ibnd, ipol, kpoint) *   &
                               uact (mu + ipol)
                          ps4 (ikb, ipol, ibnd) = ps4 (ikb, ipol, ibnd) + &
                               fact  *  (0.d0,-1.d0) *                    &
                               becp1(jkb, ibnd, kpoint) *                 &
                               uact (mu + ipol) * tpiba
                       endif
                    enddo
                 enddo
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        endif
     enddo
  enddo
  !
  !      This term is proportional to beta(k+q+G)
  !
  dvpsi = matmul(vkb, ps1) + dvpsi 
  dvpsi = matmul(vkb, ps3) + dvpsi
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        do ibnd = 1, nbnd
           ok = ok.or. (abs (ps2 (ikb, ipol, ibnd)).gt.eps .or. abs (ps4 (ikb, ipol, ibnd)) .gt.eps)
        enddo
        if (ok) then
           do ig = 1, npw
              igg = igkq (ig)
              aux (ig) =  vkb(ig, ikb) * (xk(ipol, kpoint) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              dvpsi(1:npw,ibnd) =                       &
                   ps2(ikb,ipol,ibnd) * aux(1:npw) +    &
                   ps4(ikb,ipol,ibnd) * aux(1:npw) +    &
                   dvpsi(1:npwq,ibnd)
           enddo
        endif
     enddo
  enddo


  deallocate(gk)
  deallocate(work)
  deallocate (aux)
  deallocate(aux1)
  if (allocated(dvkb1))   deallocate (dvkb1)
  if (allocated(dvkb))    deallocate (dvkb)
  if (allocated(ps1))     deallocate(ps1)
  if (allocated(ps2))     deallocate(ps2)
  if (allocated(ps3))     deallocate(ps3)
  if (allocated(ps4))     deallocate(ps4)
  if (allocated(alphadk)) deallocate (alphadk)

  return

end subroutine add_dkmds
