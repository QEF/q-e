!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine radial_gradient(f,gf,r,mesh,iflag)
  !
  !! This subroutine calculates the derivative with respect to r of a
  !! radial function defined on the mesh r.
  !
  !! * If \(\text{iflag}=0\) it uses all mesh points.
  !! * If iflag=1 it uses only a coarse grained mesh close to the
  !!   origin, to avoid large errors in the derivative when the function
  !!   is too smooth.
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




!------------------------------------------------------------------------------
SUBROUTINE radial_gradient2( f, gf, r, mesh, iflag, ix_l )
  !---------------------------------------------------------------------------
  !! This subroutine calculates the derivative with respect to r of a
  !! radial function defined on the mesh r.
  !
  !! * If \(\text{iflag}=0\) it uses all mesh points.
  !! * If iflag=1 it uses only a coarse grained mesh close to the
  !!   origin, to avoid large errors in the derivative when the function
  !!   is too smooth.
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: mesh, iflag
  REAL(DP), INTENT(IN) :: f(mesh*ix_l), r(mesh)
  REAL(DP), INTENT(OUT) :: gf(mesh*ix_l)
  INTEGER, INTENT(IN) :: ix_l
  !
  INTEGER :: i, j, k, imin, npoint, i0, ix, ixk
  REAL(DP) :: delta, b(5), faux(6), raux(6)
  !
  ! This formula is used in the all-electron case.
  !
  IF (iflag==0) THEN
     !
     !$acc data present_or_copyin(f,r) present_or_copyout(gf)
     !
#if defined(_OPENACC)
!$acc parallel loop collapse(2)
#else
!$omp parallel default(private), &
!$omp shared(ix_l,mesh,gf,r,f)
!$omp do
#endif
     DO ix = 1, ix_l
       DO i = 2, mesh-1
         ixk = (ix-1)*mesh + i
         gf(ixk) = ( (r(i+1)-r(i))**2*(f(ixk-1)-f(ixk)) &
                  -(r(i-1)-r(i))**2*(f(ixk+1)-f(ixk)) ) &
                 / ((r(i+1)-r(i))*(r(i-1)-r(i))*(r(i+1)-r(i-1)))
       ENDDO
     ENDDO
#if !defined(_OPENACC)
!$omp end do
!$omp end parallel
#endif
     !
     !$acc parallel loop
     DO ix = 1, ix_l
       i0 = (ix-1)*mesh
       gf(i0+mesh) = 0.0_DP
       gf(i0+1) = gf(i0+2) + (gf(i0+3)-gf(i0+2)) * (r(1)-r(2)) / (r(3)-r(2))
     ENDDO
     !
     !$acc end data
     !
  ELSE
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
     !$acc data present_or_copyin(f,r) present_or_copyout(gf)
     !
     delta = 0.00001_DP
     !
     imin = 1
     !$acc parallel loop collapse(3)
     points: DO i = 2, mesh
       DO j = 2, mesh
         DO k = mesh, 1, -1
           !
           IF (j<=i .OR. k>i) CYCLE
           !
           IF (r(j)>r(i)+delta) THEN
             imin = i
             IF (r(k)<r(i)-delta) THEN
               gf(i) = ( (r(j)-r(i))**2*(f(k)-f(i))   &
                        -(r(k)-r(i))**2*(f(j)-f(i)) ) &
                      /((r(j)-r(i))*(r(k)-r(i))*(r(j)-r(k)))
               CYCLE points
             ENDIF
           ENDIF
           !
           IF (j==mesh) gf(i) = 0.0_DP
           !
         ENDDO 
       ENDDO
     ENDDO points
     !
     !  In the first imin points the previous formula cannot be
     !  used. We interpolate with a polynomial the points already found
     !  and extrapolate in the points from 1 to imin.
     !  Presently we fit 5 points with a 3rd degree polynomial.
     !
     !$acc end data
     !
     !
     npoint = 5
     raux = 0.0_DP
     faux = 0.0_DP
     faux(1) = gf(imin+1)
     raux(1) = r(imin+1)
     j = imin+1
     ! $ acc parallel loop collapse(2) copyout(raux,faux)
     points_fit: DO k = 2, npoint
        DO i = j, mesh-1
           IF (r(i)>r(imin+1)+(k-1)*delta) THEN
              faux(k) = gf(i)
              raux(k) = r(i)
              j = i+1
              CYCLE points_fit
           ENDIF
        ENDDO
     ENDDO points_fit
     !
     !
     CALL fit_pol( raux, faux, npoint, 3, b )
     !
     DO i = 1, imin
        gf(i) = b(1) + r(i) * (b(2)+r(i)*(b(3)+r(i)*b(4)))
     ENDDO
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE radial_gradient2
!
!
!-------------------------------------------------------------------
SUBROUTINE fit_pol( xdata, ydata, n, degree, b )
  !-----------------------------------------------------------------
  !! This routine finds the coefficients of the least-square polynomial which 
  !! interpolates the n input data points.
  !
  USE kinds, ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, degree
  REAL(DP), INTENT(IN) :: xdata(n), ydata(n)
  REAL(DP), INTENT(OUT) :: b(degree+1)
  !
  INTEGER :: ipiv(degree+1), info, i, j, k
  REAL(DP) :: bmat(degree+1,degree+1), amat(degree+1,n)
  !
  amat(1,:) = 1.0_DP
  !
  DO i = 2, degree+1
     DO j = 1, n
        amat(i,j) = amat(i-1,j)*xdata(j)
     ENDDO
  ENDDO
  !
  DO i = 1, degree+1
     b(i) = 0.0_DP
     DO k = 1, n
        b(i) = b(i)+ydata(k)*xdata(k)**(i-1)
     ENDDO
  ENDDO
  !
  DO i = 1, degree+1
     DO j = 1, degree+1
        bmat(i,j) = 0.0_DP
        DO k = 1, n
           bmat(i,j) = bmat(i,j)+amat(i,k)*amat(j,k)
        ENDDO
     ENDDO
  ENDDO
  !
  !  This lapack routine solves the linear system that gives the
  !  coefficients of the interpolating polynomial.
  !
  CALL DGESV( degree+1, 1, bmat, degree+1, ipiv, b, degree+1, info )
  !
  IF (info/=0) CALL errore('pol_fit','problems with the linear system', &
       ABS(info))
  !
  RETURN
  !
END SUBROUTINE fit_pol


