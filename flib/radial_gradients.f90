!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine radial_gradient(f,gf,r,mesh,iflag)
  !
  !  This subroutine calculates the derivative with respect to r of a
  !  radial function defined on the mesh r. If iflag=0 it uses all mesh
  !  points. If iflag=1 it uses only a coarse grained mesh close to the
  !  origin, to avoid large errors in the derivative when the function
  !  is too smooth.
  !
  use kinds, only : DP
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
  !  Presently we fit 5 points with a 3rd degree polynomial.
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
end subroutine radial_gradient

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

