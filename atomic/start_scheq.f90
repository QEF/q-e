!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine start_scheq(lam,e,b,grid,ze2,y)
  !---------------------------------------------------------------
  !
  !  determines the wave-function in the first two points by
  !  series developement. It receives as input:
  !  lam the angular momentum 
  !  e   the energy in Ry 
  !  b(0:3) the coefficients of a polynomial that interpolates the 
  !         potential in the first three points
  !  grid  the mesh
  !  ze2   the zed of the mesh
  !  in output y(1:2) contains the solution in the first two points
  !
USE kinds, ONLY : DP
USE radial_grids, ONLY : radial_grid_type
IMPLICIT NONE
TYPE(radial_grid_type), INTENT(IN) :: grid
INTEGER, INTENT(IN) :: lam
REAL(DP), INTENT(IN) :: b(0:3), e
REAL(DP) :: ze2, xl1, x4l6, x6l12, x8l20, b0e, c1, c2, c3, c4, rr1, rr2
REAL(DP) :: y(1:2)
INTEGER :: l1
!
!  set up constants and initialize
!
l1=lam+1
xl1=lam+1.0_DP
x4l6=4.0_dp*lam+6.0_dp
x6l12=6.0_dp*lam+12.0_dp
x8l20=8.0_dp*lam+20.0_dp
!
!
b0e=b(0)-e
c1=0.5_dp*ze2/xl1
c2=(c1*ze2+b0e)/x4l6
c3=(c2*ze2+c1*b0e+b(1))/x6l12
c4=(c3*ze2+c2*b0e+c1*b(1)+b(2))/x8l20
rr1=(1.0_dp+grid%r(1)*(c1+grid%r(1)*(c2+grid%r(1)*(c3+grid%r(1)*c4))))*grid%r(1)**l1
rr2=(1.0_dp+grid%r(2)*(c1+grid%r(2)*(c2+grid%r(2)*(c3+grid%r(2)*c4))))*grid%r(2)**l1
y(1)=rr1/grid%sqr(1)
y(2)=rr2/grid%sqr(2)

return
end subroutine start_scheq
