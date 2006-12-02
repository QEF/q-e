!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
subroutine kpoint_grid &
     ( nrot, s, bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
!-----------------------------------------------------------------------
!
!  Automatic generation of a uniform grid of k-points
!
  USE kinds, only: DP
  USE noncollin_module, ONLY: noncolin
  USE spin_orb,         ONLY : domag
  USE symme,            ONLY : t_rev
  implicit none
  ! INPUT:
  integer nrot, s(3,3,48), npk, k1, k2, k3, nk1, nk2, nk3
  real(DP) bg(3,3)
  ! OUTPUT:
  integer nks
  real(DP)  xk(3,npk), wk(npk)
  ! LOCAL:
  real(DP), parameter :: eps=1.0d-5
  real(DP) xkr(3), deltap(3), deltam(3), fact, xx, yy, zz
  real(DP), allocatable:: xkg(:,:), wkk(:)
  integer nkr, i,j,k, ns, n, nk
  integer, allocatable :: equiv(:)
  logical :: in_the_list
  !
  nkr=nk1*nk2*nk3
  allocate (xkg( 3,nkr),wkk(nkr))    
  allocate (equiv( nkr))    
  !
  do i=1,nk1
     do j=1,nk2
        do k=1,nk3
           !  this is nothing but consecutive ordering
           n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1,n) = DBLE(i-1)/nk1 + DBLE(k1)/2/nk1
           xkg(2,n) = DBLE(j-1)/nk2 + DBLE(k2)/2/nk2
           xkg(3,n) = DBLE(k-1)/nk3 + DBLE(k3)/2/nk3
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
        wkk(nk)   = 1.0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        do ns=1,nrot
           do i=1,3
              xkr(i) = s(i,1,ns) * xkg(1,nk) &
                     + s(i,2,ns) * xkg(2,nk) &
                     + s(i,3,ns) * xkg(3,nk)
              xkr(i) = xkr(i) - nint( xkr(i) )
           end do
           if(t_rev(ns).eq.1) xkr = -xkr
           xx = xkr(1)*nk1 - 0.5d0*k1
           yy = xkr(2)*nk2 - 0.5d0*k2
           zz = xkr(3)*nk3 - 0.5d0*k3
           in_the_list = abs(xx-nint(xx)).le.eps .and. abs(yy-nint(yy)).le.eps &
                                                 .and. abs(zz-nint(zz)).le.eps 
           if (in_the_list) then
              i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
              j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
              k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              if (n.gt.nk .and. equiv(n).eq.n) then
                 equiv(n) = nk
                 wkk(nk)=wkk(nk)+1.0
              else
                 if (equiv(n).ne.nk .or. n.lt.nk ) call errore('kpoint_grid', &
                    'something wrong in the checking algorithm',1)
              end if
           end if
           if (.not. noncolin.or. .not.domag) then
              xx =-xkr(1)*nk1 - 0.5d0*k1
              yy =-xkr(2)*nk2 - 0.5d0*k2
              zz =-xkr(3)*nk3 - 0.5d0*k3
              in_the_list=abs(xx-nint(xx)).le.eps.and.abs(yy-nint(yy)).le.eps &
                                                 .and. abs(zz-nint(zz)).le.eps 
              if (in_the_list) then
                 i = mod ( nint (-xkr(1)*nk1 - 0.5 * k1 + 2*nk1), nk1 ) + 1
                 j = mod ( nint (-xkr(2)*nk2 - 0.5 * k2 + 2*nk2), nk2 ) + 1
                 k = mod ( nint (-xkr(3)*nk3 - 0.5 * k3 + 2*nk3), nk3 ) + 1
                 n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                 if (n.gt.nk .and. equiv(n).eq.n) then
                    equiv(n) = nk
                    wkk(nk)=wkk(nk)+1.0
                 else
                    if (equiv(n).ne.nk.or.n.lt.nk) call errore('kpoint_grid', &
                    'something wrong in the checking algorithm',2)
                 end if
              end if
           endif
        end do
     end if
  end do

  !  count irreducible points and order them

  nks=0
  fact=0.0
  do nk=1,nkr
     if (equiv(nk).eq.nk) then
        nks=nks+1
        if (nks.gt.npk) call errore('kpoint_grid','too many k-points',1)
        wk(nks) = wkk(nk)
        fact    = fact+wk(nks)
        !  bring back into to the first BZ
        do i=1,3
           xk(i,nks) = xkg(i,nk)-nint(xkg(i,nk))
        end do
     end if
  end do
  !  go to cartesian axis (in units 2pi/a0)
  call cryst_to_cart(nks,xk,bg,1)
  !  normalize weights to one
  do nk=1,nks
     wk(nk) = wk(nk)/fact
  end do

  deallocate(equiv)
  deallocate(xkg,wkk)

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
  USE kinds, only: DP
  implicit none
  ! INPUT:
  integer nks, nsym, s(3,3,48), npk, k1, k2, k3, nk1, nk2, nk3, ntetra
  logical minus_q
  real(DP) :: at(3,3), bg(3,3), xk(3,npk), wk(npk)
  ! OUTPUT:
  integer tetra(4,ntetra)
  ! LOCAL:
  real(DP) :: xkr(3), deltap(3), deltam(3)
  real(DP), parameter:: eps=1.0d-5
  real(DP), allocatable :: xkg(:,:)
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
           xkg(1,n) = DBLE(i-1)/nk1 + DBLE(k1)/2/nk1
           xkg(2,n) = DBLE(j-1)/nk2 + DBLE(k2)/2/nk2
           xkg(3,n) = DBLE(k-1)/nk3 + DBLE(k3)/2/nk3
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
     call errore('tetrahedra','cannot locate  k point',nk)
15   continue
  end do

  do n=1,nks
     do nk=1,nkr
        if (equiv(nk).eq.n) go to 20
     end do
     !  this failure of the algorithm may indicate that the displaced grid
     !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
     call errore('tetrahedra','cannot remap grid on k-point list',n)
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
             call errore ('tetrahedra','something wrong',n)
     end do
  end do

  deallocate(equiv)
  deallocate(xkg)

  return
end subroutine tetrahedra
