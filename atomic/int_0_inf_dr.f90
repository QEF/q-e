!
!---------------------------------------------------------------
function int_0_inf_dr(f,r,r2,dx,mesh,nst)
  !---------------------------------------------------------------
  !
  !      integral of f from 0 to infinity
  !      f is given on a logarithmic mesh. 
  !      f(r) is assumed to be proportional to r**nst for small r
  !
  use kinds, only : DP
  implicit none
  integer :: mesh, nst, i
  real(kind=dp):: int_0_inf_dr, f(mesh), r(mesh),  r2(mesh), dx
  real(kind=dp):: fs(4), b(4), sum1
  !
  !      series development: contribution for small r
  !
  do i=1,4
     fs(i)=f(i)/r(i)**nst
  end do
  call series(fs,r,r2,b)
  int_0_inf_dr = ( b(1)/(nst+1) + r(1)*(b(2)/(nst+2) &
       +r(1)*b(3)/(nst+3)) ) * r(1)**(nst+1)
  !
  !      simpson integration (logarithmic mesh: dr ==> r dx)
  !
  sum1=0.0_DP
  do i=1,mesh-2,2
     sum1 = sum1 + f(i)*r(i) + 4.0_DP*f(i+1)*r(i+1) + f(i+2)*r(i+2)
  end do
  int_0_inf_dr = int_0_inf_dr + sum1*dx/3.0_DP

  return
end function int_0_inf_dr
