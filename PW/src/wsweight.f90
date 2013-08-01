!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine wsinit(rws,nrwsx,nrws,atw)
!-----------------------------------------------------------------------
!
  USE kinds, only : DP
  implicit none
  integer i, ii, ir, jr, kr, nrws, nrwsx, nx
  real(DP) eps, rws(0:3,nrwsx), atw(3,3)
  parameter (eps=1.0d-6,nx=2)
  ii = 1
  do ir=-nx,nx
     do jr=-nx,nx
        do kr=-nx,nx
           do i=1,3
              rws(i,ii) = atw(i,1)*ir + atw(i,2)*jr + atw(i,3)*kr
           end do
           rws(0,ii)=rws(1,ii)*rws(1,ii)+rws(2,ii)*rws(2,ii)+            &
                               rws(3,ii)*rws(3,ii)
           rws(0,ii)=0.5d0*rws(0,ii)
           if (rws(0,ii).gt.eps) ii = ii + 1
           if (ii.gt.nrwsx) call errore('wsinit', 'ii.gt.nrwsx',1)
        end do
     end do
  end do
  nrws = ii - 1
  return
end subroutine wsinit
!
!-----------------------------------------------------------------------
function wsweight(r,rws,nrws)
!-----------------------------------------------------------------------
!
! wsweights assigns this weight:
! - if a point is inside the Wigner-Seitz cell:    weight=1
! - if a point is outside the WS cell:             weight=0
! - if a point q is on the border of the WS cell, it finds the number N 
!   of translationally equivalent point q+G  (where G is a lattice vector)
!   that are also on the border of the cell. Then: weight = 1/N

! I.e. if a point is on the surface of the WS cell of a cubic lattice 
! it will have weight 1/2; on the vertex of the WS it would be 1/8; 
! the K point of an hexagonal lattice has weight 1/3 and so on.

  USE kinds, only : dp
  implicit none
  integer ir, nreq, nrws
  real(DP) r(3), rrt, ck, eps, rws(0:3,nrws), wsweight
  parameter (eps=1.0d-6)
!
  wsweight = 0.d0
  nreq = 1
  do ir =1,nrws
     rrt = r(1)*rws(1,ir) + r(2)*rws(2,ir) + r(3)*rws(3,ir)
     ck = rrt-rws(0,ir)
     if ( ck .gt. eps ) return
     if ( abs(ck) .lt. eps ) nreq = nreq + 1
  end do
  wsweight = 1.d0/DBLE(nreq)
  return
end function wsweight
