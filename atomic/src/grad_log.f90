!-----------
subroutine grad_log(v,dv,rm1,dx,nx,method)
  !----------
  ! Computes derivative of a function v(:) on a radial logarithmic mesh
  ! using d/dr = 1/r d/dx, since r \propto exp(x)
  ! Two methods are available (method=2,4), the second (method=4)
  ! is more precise (but more sensitive to noise)
  !

  use kinds, only : DP
  implicit none
  ! input/output variables
  integer, intent(in) :: nx, method
  ! dimensions, method used
  real(DP), intent(in) :: v(nx), rm1(nx), dx
  ! function, 1/r, dx (step of the log. mesh)
  real(DP), intent(out) :: dv(nx)
  ! derivative

  integer :: i

  select case ( method )
  case ( 2 )
     ! "leap frog" derivative O(dx^2)
     ! use 2 points for each derivative
     dv(2:nx-1) = 0.5_DP/dx * (v(3:nx) - v(1:nx-2)) * rm1(2:nx)

     ! this (forward) formula is also O(dx^2),
     ! but the error coefficient in front of dx^2 is different
     !gv(1) = 0.5_DP/dx * (4._DP*v(2) - 3._DP*v(1) - v(3)) * rm1(1) 
     ! while this has the same error coefficient in front of dx^2
     dv(1) = 0.5_DP/dx * (-4._DP*(v(1)+v(3)) + 7._DP*v(2) + v(4)) * rm1(1)

     ! backward O(dx^2) formula (with different error coeff.)
     dv(nx) = -0.5_DP/dx *(4._DP*v(nx-1) - 3._DP*v(nx) - v(nx-2)) * rm1(nx)
  case ( 4 )  
     ! use 4 points for each derivative
     !
     ! calculate dv/dr (taken from subroutine lschps.f90)
     dv(1)=(-50.0_dp*v(1)+96.0_dp*v(2)-72.0_dp*v(3)+32.0_dp*v(4) &
          -6.0_dp*v(5))/(24.0_dp*dx) * rm1(1)
     dv(2)=(-6.0_dp*v(1)-20.0_dp*v(2)+36.0_dp*v(3)-12.0_dp*v(4) &
          +2.0_dp*v(5))/(24.0_dp*dx) * rm1(2)
     !
     do i=3,nx-2
        dv(i)=(2.0_dp*v(i-2)-16.0_dp*v(i-1)+16.0_dp*v(i+1) &
             -2.0_dp*v(i+2))/(24.0_dp*dx) * rm1(i)
     end do
     !
     dv(nx-1)=( 3.0_dp*v(nx)+10.0_dp*v(nx-1)-18.0_dp*v(nx-2)+ &
          6.0_dp*v(nx-3)-v(nx-4))/(12.0_dp*dx) * rm1(nx-1)
     dv(nx)=( 25.0_dp*v(nx)-48.0_dp*v(nx-1)+36.0_dp*v(nx-2)-&
          16.0_dp*v(nx-3)+3.0_dp*v(nx-4))/(12.0_dp*dx) * rm1(nx)
   case default
      call errore('grad_log','method unknown',1)
   end select

end subroutine grad_log
