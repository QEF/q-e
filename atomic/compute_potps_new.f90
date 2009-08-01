!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE compute_potps_new(ik,v_in,v_out,xc)
   !--------------------------------------------------------------------------
   !
   !     This routine computes the local pseudopotential by smoothing the
   !     all_electron potential. In input it receives, the point
   !     ik where the cut is done.
   !
   !     The smooth potential has the following form for r < rc
   !
   !     V(r) = A + B * f( q * r)
   !
   !     where f(x) = [a*a*sen(x)/x - sen(a*x)/(a*x)]/[a*a-1]
   !
   !     V,V' and V'' are matched to the ae value in rc and 
   !     both V'(0) and V''(0) vanish.
   !
   !     According to Troulier and Martins this is an optimal choice for
   !     rapidly converging and transferable potentials.
   !
   use kinds, only: dp
   use radial_grids, only: ndmx
   use ld1inc, only: grid
   IMPLICIT NONE
   REAL(DP), PARAMETER :: a =  .707106781186547 !=sqrt(0.5_dp) ! the parameter defining f(x)

   REAL(DP) :: &
         v_in(ndmx), & ! input: the potential to pseudize
         v_out(ndmx),& ! output: the pseudized potential
         xc(8)        ! output: the coefficients of the fit

   INTEGER :: &
         ik        ! input: the point corresponding to rc

   REAL(DP) :: &
         fae,    & ! the value of the all-electron function
         f1ae,   & ! its first derivative
         f2ae      ! the second derivative

   REAL(DP) :: &
         AA,     & ! the matching value
         deriv,  &  ! auxilairy quantities
         j1(ndmx,2) ! auxiliary functions
     
   REAL(DP) :: &
         deriv_7pts, deriv2_7pts

   INTEGER :: &
         n,   &  ! counter on mesh points
         nc      ! counter on bessel
   !
   !    compute first and second derivative
   !
   fae=v_in(ik)
   f1ae=deriv_7pts(v_in,ik,grid%r(ik),grid%dx)
   f2ae=deriv2_7pts(v_in,ik,grid%r(ik),grid%dx)
   !
   ! set the matching value
   !
   AA = f2ae*grid%r(ik)/f1ae
   ! find the solution
   call find_matching_x(AA,xc(5))
   ! rescale the solution
   xc(5) = xc(5)/grid%r(ik)
   xc(4) = sqrt(0.5_dp)*xc(5)
   ! build the auxiliary functions const and f(x)
   DO nc=1,2
      call sph_bes(ik+1,grid%r,xc(3+nc),0,j1(1,nc))
   end do
   j1(1:ik+1,2) = 0.5_dp * j1(1:ik+1,2) - j1(1:ik+1,1)
   j1(1:ik+1,1) = 1.0_dp
   ! determine A and B
   deriv = 0.5_dp *(j1(ik+1,2)-j1(ik,2))/(grid%r(ik+1)-grid%r(ik)) +  &
           0.5_dp *(j1(ik,2)-j1(ik-1,2))/(grid%r(ik)-grid%r(ik-1))
   xc(2) = f1ae / deriv
   xc(1) = fae - xc(2) * j1(ik,2)
   !
   ! define the v_out function
   !
   DO n=1,ik
      v_out(n)=xc(1)*j1(n,1)+xc(2)*j1(n,2)
   ENDDO
   DO n=ik+1,grid%mesh
      v_out(n)=v_in(n)
   ENDDO

   RETURN

CONTAINS

   SUBROUTINE find_matching_x(AA,x)
      real(dp) :: x, AA
      real(dP) :: aa_min, aa_max, aa_test, xmin, xmax, xtest, dx
      !
      dx = 0.1_dp
      !
      xmin = 0._dp
      aa_min = 3._dp - AA
      IF (aa_min < 0._dp) CALL errore('compute_potps_new', &
                               'unable to find a solution.. try lloc=-1',1)
      xmax = dx
      call fz(xmax, AA, aa_max)
      ! bracket the first solution
      do while (aa_min * aa_max >0)
         xmin   = xmax
         aa_min = aa_max
         xmax   = xmin + dx
         call fz(xmax, AA, aa_max)
      end do
      ! refine by bisection
      do while (abs(aa_max-aa_min) > 1.d-8) 
         xtest= (xmax+xmin)/2.d0
         call fz(xtest, AA, aa_test)
         if (aa_test * aa_max > 0 ) then
            aa_max = aa_test
            xmax = xtest
         else
            aa_min = aa_test
            xmin = xtest
         endif
      end do
      x = (xmax+xmin)/2.d0
      return
   END SUBROUTINE find_matching_x
   !
   ! the auxiliary target function  x*f''(x)/f'(x)-AA
   SUBROUTINE fz(x,AA,aa_f)
      use kinds, only: dp
      real(dp) :: x, AA, aa_f
      real(dp) :: fs, f1s, f2s, fac
      real(dp) :: gs, g1s, g2s, xx
      fac = 1._dp/(a*a-1._dp)
      call ffz(x,fs,f1s,f2s)
      xx = a*x
      call ffz(xx,gs,g1s,g2s)
      fs = fac*(a*a*fs-gs)
      f1s= fac*(a*a*f1s-a*g1s)
      f2s= fac*(a*a*f2s-a*a*g2s)
      aa_f = x * f2s / f1s - AA
      return
   END SUBROUTINE
   !
   ! the first bessel function and its derivatives
   SUBROUTINE ffz(x,fs,f1s,f2s)
      use kinds, only: dp
      implicit none
      real(dp) :: x, fs, f1s, f2s
      fs  = (sin(x))/x
      f1s = (x*cos(x)-sin(x))/x**2
      f2s = (-x*x*sin(x) - 2*x*cos(x) + 2*sin(x))/x**3
      return
   END SUBROUTINE ffz

END SUBROUTINE compute_potps_new
