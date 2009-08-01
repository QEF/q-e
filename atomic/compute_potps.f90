!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
      SUBROUTINE compute_potps(ik,v_in,v_out,xc)
!--------------------------------------------------------------------------
!
!     This routine computes the pseudized pseudopotential by pseudizing the
!     all_electron potential. In input it receives, the point
!     ik where the cut is done.
!
use kinds, only: dp
use radial_grids, only: ndmx
use ld1inc, only: grid
IMPLICIT NONE

REAL(DP) :: &
         v_in(ndmx), & ! input: the potential to pseudize
         v_out(ndmx),& ! output: the pseudized potential
         xc(8)        ! output: the coefficients of the fit

INTEGER :: &
         ik        ! input: the point corresponding to rc

REAL(DP) :: &
         p1aep1, p1aem1, &
         fae,    & ! the value of the all-electron function
         f1ae,   & ! its first derivative
         f2ae      ! the second derivative

REAL(DP) :: &
         bm(2),  &              ! the derivative of the bessel
         fact(2), &             ! factor of normalization
         j1(ndmx,8)              ! the bessel functions
     
REAL(DP) :: &
         deriv_7pts, deriv2_7pts

INTEGER :: &
         iok, &  ! if 0 there are no problem
         n,   &  ! counter on mesh points
         nc      ! counter on bessel
!
!    compute first and second derivative
!
fae=v_in(ik)
f1ae=deriv_7pts(v_in,ik,grid%r(ik),grid%dx)
f2ae=deriv2_7pts(v_in,ik,grid%r(ik),grid%dx)
!
!    find the q_i of the bessel functions
!      
CALL find_qi(f1ae/fae,xc(4),ik,0,2,0,iok)
IF (iok.NE.0) &
        CALL errore('compute_potps','problems with find_qi',1)
!
!    compute the functions
!
DO nc=1,2
   call sph_bes(ik+1,grid%r,xc(3+nc),0,j1(1,nc))
   fact(nc)=v_in(ik)/j1(ik,nc)
   DO n=1,ik+1
      j1(n,nc)=j1(n,nc)*fact(nc)
   ENDDO
ENDDO
!
!    compute the second derivative and impose continuity of zero, 
!    first and second derivative
!
DO nc=1,2
   p1aep1=(j1(ik+1,nc)-j1(ik,nc))/(grid%r(ik+1)-grid%r(ik))
   p1aem1=(j1(ik,nc)-j1(ik-1,nc))/(grid%r(ik)-grid%r(ik-1))
   bm(nc)=(p1aep1-p1aem1)*2.0_dp/(grid%r(ik+1)-grid%r(ik-1))
ENDDO

xc(2)=(f2ae-bm(1))/(bm(2)-bm(1))
xc(1)=1.0_dp-xc(2)
!
!    define the v_out function
!
DO n=1,ik
   v_out(n)=xc(1)*j1(n,1)+xc(2)*j1(n,2)
ENDDO

DO n=ik+1,grid%mesh
   v_out(n)=v_in(n)
ENDDO

RETURN
END SUBROUTINE compute_potps
