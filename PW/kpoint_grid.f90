!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine kpoint_grid &
     ( nrot, s, bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
!-----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points
!
#include "machine.h"
  use parameters, only: DP
  implicit none
  ! INPUT:
  integer nrot, s(3,3,48), npk, k1, k2, k3, nk1, nk2, nk3
  real(kind=DP) bg(3,3)
  ! OUTPUT:
  integer nks
  real(kind=DP)  xk(3,npk), wk(npk)
  ! LOCAL:
  real(kind=DP), parameter :: degspin=2.d0, eps=1.0e-5
  ! degspin: spin degeneracy used to normalize k-points
  real(kind=DP) xkr(3), deltap(3), deltam(3), fact
  real(kind=DP), allocatable:: xkg(:,:)
  integer nkr, i,j,k, ns, n, nk
  integer, allocatable :: equiv(:)
  !
  nkr=nk1*nk2*nk3
  allocate (xkg( 3,nkr))    
  allocate (equiv( nkr))    
  !
  do i=1,nk1
     do j=1,nk2
        do k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        end do
     end do
  end do

  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)

  do nk=1,nkr
     equiv(nk)=nk
  end do

  do nk=1,nkr
     !  check if this k-point has already been found equivalent to another
     if (equiv(nk).eq.nk) then
        wk(nk)   = 1.0
        do n=nk+1,nkr
           !  check if there are equivalent k-point to this in the list
           !  (excepted those previously found to be equivalent to another)
           if (equiv(n).eq.n) then
              do ns=1,nrot
                 do i=1,3
                    xkr(i) = s(i,1,ns) * xkg(1,n) &
                           + s(i,2,ns) * xkg(2,n) &
                           + s(i,3,ns) * xkg(3,n)
                 end do
                 ! xkr are the components in crystal axis of the rotated k-point
                 do i=1,3
                    deltap(i) = xkr(i)-xkg(i,nk) - &
                          nint (xkr(i)-xkg(i,nk) )
                    deltam(i) = xkr(i)+xkg(i,nk) - &
                          nint (xkr(i)+xkg(i,nk) )
                 end do
                 !  deltap is the difference vector,
                 !  brought back in the first BZ
                 !  deltam is the same but with k => -k (for time reversal)
                 if ( sqrt ( deltap(1)**2 +  &
                             deltap(2)**2 +  &
                             deltap(3)**2 ) .lt. eps .or. &
                      sqrt ( deltam(1)**2 +  &
                             deltam(2)**2 +  &
                             deltam(3)**2 ) .lt. eps      ) then
                    !  equivalent k-point found:
                    !  add 1 to the weight and go to next k-point
                    equiv(n) = nk
                    wk(nk)=wk(nk)+1.0
                    go to 10
                 end if
              end do
           end if
10         continue
        end do
     end if
  end do

  !  count irreducible points and order them

  nks=0
  fact=0.0
  do nk=1,nkr
     if (equiv(nk).eq.nk) then
        nks=nks+1
        if (nks.gt.npk) call error('kpoint_grid','too many k-points',1)
        wk(nks) = wk(nk)
        fact    = fact+wk(nks)
        !  bring back into to the first BZ
        do i=1,3
           xk(i,nks) = xkg(i,nk)-nint(xkg(i,nk))
        end do
     end if
  end do
  !  go to cartesian axis (in units 2pi/a0)
  call cryst_to_cart(nks,xk,bg,1)
  !  normalize weights to degspin (every band can accomodate 2 electrons)
  do nk=1,nks
     wk(nk) = degspin * wk(nk)/fact
  end do

  deallocate(equiv)
  deallocate(xkg)

  return
end subroutine kpoint_grid
!
!-----------------------------------------------------------------------
subroutine tetrahedra ( nsym, s, minus_q, at, bg, npk, k1,k2,k3, &
     nk1,nk2,nk3, nks, xk, wk, ntetra, tetra )
  !-----------------------------------------------------------------------
  !
  ! Tetrahedron method according to P. E. Bloechl et al, PRB49, 16223 (1994)
  !
#include "machine.h"
  use parameters, only: DP
  implicit none
  ! INPUT:
  integer nks, nsym, s(3,3,48), npk, k1, k2, k3, nk1, nk2, nk3, ntetra
  logical minus_q
  real(kind=DP) :: at(3,3), bg(3,3), xk(3,npk), wk(npk)
  ! OUTPUT:
  integer tetra(4,ntetra)
  ! LOCAL:
  real(kind=DP) :: xkr(3), deltap(3), deltam(3)
  real(kind=DP), parameter:: eps=1.0d-5
  real(kind=DP), allocatable :: xkg(:,:)
  integer :: nkr, i,j,k, ns, n, nk, ip1,jp1,kp1, &
       n1,n2,n3,n4,n5,n6,n7,n8
  integer, allocatable:: equiv(:)
  !
  ! Re-generate a uniform grid of k-points xkg
  !
  nkr=nk1*nk2*nk3
  !      ntetra=6*nkr
  allocate (xkg( 3,nkr))    
  allocate (equiv( nkr))    
!
  do i=1,nk1
     do j=1,nk2
        do k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
           xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
           xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
        end do
     end do
  end do

  !  locate k-points of the uniform grid in the list of irreducible k-points
  !  that was previously calculated

  !  bring irreducible k-points to crystal axis
  call cryst_to_cart (nks,xk,at,-1)
  !
  do nk=1,nkr
     do n=1,nks
        do ns=1,nsym
           do i=1,3
              xkr(i) = s(i,1,ns) * xk(1,n) + &
                       s(i,2,ns) * xk(2,n) + &
                       s(i,3,ns) * xk(3,n)
           end do
           !  xkr is the n-th irreducible k-point rotated wrt the ns-th symmetry
           do i=1,3
              deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk) )
              deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk) )
           end do
           !  deltap is the difference vector, brought back in the first BZ
           !  deltam is the same but with k => -k (for time reversal)
           if ( sqrt ( deltap(1)**2 + &
                       deltap(2)**2 + &
                       deltap(3)**2 ) .lt. eps .or. ( minus_q .and. &
                sqrt ( deltam(1)**2 +  &
                       deltam(2)**2 +  &
                       deltam(3)**2 ) .lt. eps ) ) then
              !  equivalent irreducible k-point found
              equiv(nk) = n
              go to 15
           end if
        end do
     end do
     !  equivalent irreducible k-point found - something wrong
     call error('tetrahedra','cannot locate  k point',nk)
15   continue
  end do

  do n=1,nks
     do nk=1,nkr
        if (equiv(nk).eq.n) go to 20
     end do
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     call error('tetrahedra','cannot remap grid on k-point list',n)
20   continue
  end do

  !  bring irreducible k-points back to cartesian axis
  call cryst_to_cart (nks,xk,bg, 1)

  !  construct tetrahedra

  do i=1,nk1
     do j=1,nk2
        do k=1,nk3
           !  n1-n8 are the indices of k-point 1-8 forming a cube
           ip1 = mod(i,nk1)+1
           jp1 = mod(j,nk2)+1
           kp1 = mod(k,nk3)+1
           n1 = (  k-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n2 = (  k-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n3 = (  k-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n4 = (  k-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n5 = (kp1-1) + (  j-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n6 = (kp1-1) + (  j-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           n7 = (kp1-1) + (jp1-1)*nk3 + (  i-1)*nk2*nk3 + 1
           n8 = (kp1-1) + (jp1-1)*nk3 + (ip1-1)*nk2*nk3 + 1
           !  there are 6 tetrahedra per cube (and nk1*nk2*nk3 cubes)
           n  = 6 * ( (k-1) + (j-1)*nk3 + (i-1)*nk3*nk2 )

           tetra (1,n+1) = equiv(n1)
           tetra (2,n+1) = equiv(n2)
           tetra (3,n+1) = equiv(n3)
           tetra (4,n+1) = equiv(n6)

           tetra (1,n+2) = equiv(n2)
           tetra (2,n+2) = equiv(n3)
           tetra (3,n+2) = equiv(n4)
           tetra (4,n+2) = equiv(n6)

           tetra (1,n+3) = equiv(n1)
           tetra (2,n+3) = equiv(n3)
           tetra (3,n+3) = equiv(n5)
           tetra (4,n+3) = equiv(n6)

           tetra (1,n+4) = equiv(n3)
           tetra (2,n+4) = equiv(n4)
           tetra (3,n+4) = equiv(n6)
           tetra (4,n+4) = equiv(n8)

           tetra (1,n+5) = equiv(n3)
           tetra (2,n+5) = equiv(n6)
           tetra (3,n+5) = equiv(n7)
           tetra (4,n+5) = equiv(n8)

           tetra (1,n+6) = equiv(n3)
           tetra (2,n+6) = equiv(n5)
           tetra (3,n+6) = equiv(n6)
           tetra (4,n+6) = equiv(n7)
        end do
     end do
  end do

  !  check

  do n=1,ntetra
     do i=1,4
        if ( tetra(i,n).lt.1 .or. tetra(i,n).gt.nks ) &
             call error ('tetrahedra','something wrong',n)
     end do
  end do

  deallocate(equiv)
  deallocate(xkg)

  return
end subroutine tetrahedra
