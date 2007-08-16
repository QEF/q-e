!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine vxcgc(ndm,mesh,nspin,r,r2,rho,rhoc,vgc,egc,iflag)
  !---------------------------------------------------------------
  !
  !
  !     This routine computes the exchange and correlation potential and
  !     energy to be added to the local density, to have the first
  !     gradient correction.
  !     In input the density is rho(r) (multiplied by 4*pi*r2).
  !
  !     The units of the potential are Ryd.
  !
  use kinds, only : DP
  use constants, only : fpi
  use funct, only : gcxc, gcx_spin, gcc_spin
  use ld1inc, only : grid
  implicit none
  integer,  intent(in) :: ndm,mesh,nspin,iflag
  real(DP), intent(in) :: r(mesh), r2(mesh), rho(ndm,2), rhoc(ndm)
  real(DP), intent(out):: vgc(ndm,2), egc(ndm)

  integer :: i, is, ierr
  real(DP) :: sx,sc,v1x,v2x,v1c,v2c
  real(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw
  real(DP) :: segno, arho
  real(DP) :: rh, zeta, grh2, grho2(2)
  real(DP),parameter :: eps=1.e-12_dp

  real(DP), allocatable :: grho(:,:), h(:,:), dh(:), rhoaux(:,:)
  !
  !      First compute the charge and the charge gradient, assumed  
  !      to have spherical symmetry. The gradient is the derivative of
  !      the charge with respect to the modulus of r. 
  !
  allocate(rhoaux(mesh,2),stat=ierr)
  allocate(grho(mesh,2),stat=ierr)
  allocate(h(mesh,2),stat=ierr)
  allocate(dh(mesh),stat=ierr)

  egc=0.0_dp
  vgc=0.0_dp

  do is=1,nspin
     do i=1, mesh
        rhoaux(i,is)=(rho(i,is)+rhoc(i)/nspin)/fpi/r2(i)
     enddo
     call gradient(rhoaux(1,is),grho(1,is),r,mesh,iflag)
  enddo

  if (nspin.eq.1) then
     !
     !     GGA case
     !
     do i=1,mesh
        arho=abs(rhoaux(i,1)) 
        segno=sign(1.0_dp,rhoaux(i,1))
        if (arho.gt.eps.and.abs(grho(i,1)).gt.eps) then
           call gcxc(arho,grho(i,1)**2,sx,sc,v1x,v2x,v1c,v2c)
           egc(i)=(sx+sc)*segno
           vgc(i,1)= v1x+v1c
           h(i,1)  =(v2x+v2c)*grho(i,1)*r2(i)
           !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
           !                          rho(i,1), grho(i,1)**2,  &
           !                          vgc(i,1),h(i,1)
        else
           vgc(i,1)=0.0_dp
           egc(i)=0.0_dp
           h(i,1)=0.0_dp
        endif
     end do
  else
     !
     !   this is the \sigma-GGA case
     !       
     do i=1,mesh
        !
        !  NB: the special or wrong cases where one or two charges 
        !      or gradients are zero or negative must
        !      be detected within the gcxc_spin routine
        !
        !    spin-polarised case
        !
        do is = 1, nspin
           grho2(is)=grho(i,is)**2
        enddo

        call gcx_spin (rhoaux(i, 1), rhoaux(i, 2), grho2(1), grho2(2), &
             sx, v1xup, v1xdw, v2xup, v2xdw)
        rh = rhoaux(i, 1) + rhoaux(i, 2)
        if (rh.gt.eps) then
           zeta = (rhoaux (i, 1) - rhoaux (i, 2) ) / rh
           grh2 = (grho (i, 1) + grho (i, 2) ) **2 
           call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
        else
           sc = 0.0_dp
           v1cup = 0.0_dp
           v1cdw = 0.0_dp
           v2c = 0.0_dp
        endif

        egc(i)=sx+sc
        vgc(i,1)= v1xup+v1cup
        vgc(i,2)= v1xdw+v1cdw
        h(i,1)  =((v2xup+v2c)*grho(i,1)+v2c*grho(i,2))*r2(i)
        h(i,2)  =((v2xdw+v2c)*grho(i,2)+v2c*grho(i,1))*r2(i)
        !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
        !                          rho(i,1)*2.0_dp, grho(i,1)**2*4.0_dp, &
        !                          vgc(i,1),  h(i,2)
     enddo
  endif
  !     
  !     We need the gradient of h to calculate the last part of the exchange
  !     and correlation potential.
  !     
  do is=1,nspin
     call gradient(h(1,is),dh,r,mesh,iflag)
     !
     !     Finally we compute the total exchange and correlation energy and
     !     potential. We put the original values on the charge and multiply
     !     by two to have as output Ry units.

     do i=1, mesh
        vgc(i,is)=vgc(i,is)-dh(i)/r2(i)
        vgc(i,is)=2.0_dp*vgc(i,is)
        if (is.eq.1) egc(i)=2.0_dp*egc(i)
        !            if (is.eq.1.and.i.lt.4) write(6,'(3f20.12)') &
        !                                      vgc(i,1)
     enddo
  enddo

  deallocate(dh)
  deallocate(h)
  deallocate(grho)
  deallocate(rhoaux)

  return
end subroutine vxcgc

subroutine gradient(f,gf,r,mesh,iflag)
!
!  This subroutine calculates the derivative with respect to r of a
!  radial function defined on the mesh r. If iflag=0 it uses all mesh
!  points. If iflag=1 it uses only a coarse grained mesh close to the
!  origin, to avoid large errors in the derivative when the function
!  is too smooth.
!
use kinds, only : DP
use ld1inc, only: zed
use radial_grids, only : series
implicit none
integer, intent(in) :: mesh, iflag
real(DP), intent(in) :: f(mesh), r(mesh)
real(DP), intent(out) :: gf(mesh)

integer :: i,j,k,imin,npoint
real(DP) :: delta, b(5), faux(6), raux(6)
!
!  This formula is used in the all-electron case.
!
if (iflag==0) then
   do i=2, mesh-1
      gf(i)=( (r(i+1)-r(i))**2*(f(i-1)-f(i))   &
         -(r(i-1)-r(i))**2*(f(i+1)-f(i)) ) &
          /((r(i+1)-r(i))*(r(i-1)-r(i))*(r(i+1)-r(i-1)))
   enddo
   gf(mesh)=0.0_dp
!     
!     The gradient in the first point is a linear interpolation of the
!     gradient at point 2 and 3. 
!     
   gf(1) = gf(2) + (gf(3)-gf(2)) * (r(1)-r(2)) / (r(3)-r(2))
   return
endif
!
!  If the input function is slowly changing (as the pseudocharge),
!  the previous formula is affected by numerical errors close to the 
!  origin where the r points are too close one to the other. Therefore 
!  we calculate the gradient on a coarser mesh. This gradient is often 
!  more accurate but still does not remove all instabilities observed 
!  with the GGA. 
!  At larger r the distances between points become larger than delta 
!  and this formula coincides with the previous one.
!  (ADC 08/2007)
!

delta=0.00001_dp

imin=1
points: do i=2, mesh
   do j=i+1,mesh
      if (r(j)>r(i)+delta) then
         do k=i-1,1,-1
            if (r(k)<r(i)-delta) then
               gf(i)=( (r(j)-r(i))**2*(f(k)-f(i))   &
                      -(r(k)-r(i))**2*(f(j)-f(i)) ) &
                      /((r(j)-r(i))*(r(k)-r(i))*(r(j)-r(k)))
               cycle points
            endif
         enddo
!
!  if the code arrives here there are not enough points on the left: 
!  r(i)-delta is smaller than r(1). 
!
         imin=i
         cycle points
      endif
   enddo
!
! If the code arrives here there are not enough points on the right.
! It should happen only at mesh.
! NB: the f function is assumed to be vanishing for large r, so the gradient
!     in the last points is taken as zero.
!
   gf(i)=0.0_DP
enddo points
!
!  In the first imin points the previous formula cannot be
!  used. We interpolate with a polynomial the points already found
!  and extrapolate in the points from 1 to imin.
!  Presently we fit 5 points with a 3th degree polynomial.
!
npoint=5
raux=0.0_DP
faux=0.0_DP
faux(1)=gf(imin+1)
raux(1)=r(imin+1)
j=imin+1
points_fit: do k=2,npoint
   do i=j,mesh-1
      if (r(i)>r(imin+1)+(k-1)*delta) then
         faux(k)=gf(i)
         raux(k)=r(i)
         j=i+1
         cycle points_fit
      endif
   enddo
enddo points_fit  
call fit_pol(raux,faux,npoint,3,b)
do i=1,imin
   gf(i)=b(1)+r(i)*(b(2)+r(i)*(b(3)+r(i)*b(4)))
enddo
return
end subroutine gradient

subroutine fit_pol(xdata,ydata,n,degree,b)
!
! This routine finds the coefficients of the least-square polynomial which 
! interpolates the n input data points.
!
use kinds, ONLY : DP
implicit none

integer, intent(in) :: n, degree
real(DP), intent(in) :: xdata(n), ydata(n)
real(DP), intent(out) :: b(degree+1)

integer :: ipiv(degree+1), info, i, j, k
real(DP) :: bmat(degree+1,degree+1), amat(degree+1,n)

amat(1,:)=1.0_DP
do i=2,degree+1
   do j=1,n
      amat(i,j)=amat(i-1,j)*xdata(j)
   enddo
enddo
do i=1,degree+1
   b(i)=0.0_DP
   do k=1,n
      b(i)=b(i)+ydata(k)*xdata(k)**(i-1)
   enddo
enddo
do i=1,degree+1
   do j=1,degree+1
      bmat(i,j)=0.0_DP
      do k=1,n
         bmat(i,j)=bmat(i,j)+amat(i,k)*amat(j,k)
      enddo
   enddo
enddo
!
!  This lapack routine solves the linear system that gives the
!  coefficients of the interpolating polynomial.
!
call DGESV(degree+1, 1, bmat, degree+1, ipiv, b, degree+1, info)

if (info.ne.0) call errore('pol_fit','problems with the linear system', &
                                                              abs(info))
return
end subroutine fit_pol

