!
!---------------------------------------------------------------
subroutine green(y,lam,e,dvy,chi,vpot,ze2)
!---------------------------------------------------------------
   !
   use kinds, only: DP
   use ld1inc, only: grid
   use radial_grids, only: ndmx, series
   implicit none
   !
   ! I/O variables
   !
   integer :: lam
   real(DP) :: y(ndmx), chi(ndmx), dvy(ndmx), vpot(ndmx)
   real(DP) :: e, ze2
   !
   ! local variables
   !
   integer :: i, l1, ncross, imatch
   real(DP) :: f(ndmx), g(ndmx)
   real(DP) :: work(ndmx), int_0_inf_dr
   real(DP) :: a(0:3), b(0:3), c0, c1, c2, c3, c4, b0e
   real(DP) :: rr1, rr2, r1, r2
   real(DP) :: ddx12, sqlhf, xl1, x4l6, x6l12, x8l20
   real(DP) :: gi, gim1
   real(DP) :: fac
! set up constants and initialize
   ddx12=grid%dx*grid%dx/12.0
   l1=lam+1
   sqlhf=(DBLE(lam)+.5)**2
   xl1=l1
   x4l6=4*lam+6
   x6l12=6*lam+12
   x8l20=8*lam+20
! series developement of the potential and r.h.s. term near the origin
   do i=1,4
      y(i)=vpot(i)-ze2/grid%r(i)
   end do
   call series(y,grid%r,grid%r2,b)
   do i=1,4
      y(i)= dvy(i)/grid%r(i)**l1
   end do
   call series(y,grid%r,grid%r2,a)
!
! determine the position of the last change of sign of the f-function 
! f < 0 (approximatively) means classically allowed   region
! f > 0         "           "        "      forbidden   "
! set up the f-, g- and y-functions
!
   f(1)=ddx12*(grid%r2(1)*(vpot(1)-e)+sqlhf)
   do i=2,grid%mesh
      f(i)=ddx12*(grid%r2(i)*(vpot(i)-e)+sqlhf)
      if( f(i) .ne. sign(f(i),f(i-1)) ) imatch=i
   end do
   do i=1,grid%mesh
      f(i)=1-f(i)
      g(i)=ddx12*dvy(i)*grid%r(i)*grid%sqr(i)
      y(i)=0.d0
   end do
! wave-function in the first two points by series developement
   b0e=b(0)-e
   c0 = 1.0d0
   c1=0.5*ze2*c0/xl1
   c2=(c1*ze2+b0e*c0+a(0))/x4l6
   c3=(c2*ze2+c1*b0e+b(1)*c0+a(1))/x6l12
   c4=(c3*ze2+c2*b0e+c1*b(1)+b(2)*c0+a(2))/x8l20
   r1=grid%r(1)
   r2=grid%r(2)
   rr1=(c0+r1*(c1+r1*(c2+r1*(c3+r1*c4))))*r1**l1
   rr2=(c0+r2*(c1+r2*(c2+r2*(c3+r2*c4))))*r2**l1
   y(1)=rr1/grid%sqr(1)
   y(2)=rr2/grid%sqr(2)
!
! new set up of the g- function
   gim1 = g(1)
   do i=2,grid%mesh-1
      gi=g(i)
      g(i)=gim1+10.0d0*gi+g(i+1)
      gim1=gi
   end do
! outward integration: numerov's algorithm
   call outward(y,f,g,grid%mesh,imatch,ncross)
! inward integration: froese's algorithm
   call inward(y,f,g,grid%mesh,imatch)
!-orthogonalize to the 0^th order eigenvector
   do i=1,grid%mesh
      y(i) = y(i)*grid%sqr(i)
      work(i)=y(i)*chi(i)
   end do
   fac=int_0_inf_dr(work,grid,grid%mesh,2*lam+2)
   do i=1,grid%mesh
      y(i) = y(i) - fac*chi(i)
   end do
!
   return

end subroutine green
