!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine seriesbes(fun,r,r2,npt,xc)
  !
  !     assume that the input function has the form xc(1)
  !                              +xc(2)*r(n)+xc(3)*r(n)**2
  !     and finds the two coefficients. works with KKR3 beta functions
  !
  use kinds, only : DP
  implicit none
  integer :: &
       npt,  &        ! the number of points  
       npt2          ! intermediate point

  real(DP) :: &
       fun(npt),   &  ! the function
       r(npt),     &  ! the mesh     
       r2(npt),    &  ! the mesh     
       xc(4)         ! the coefficients

  if (npt.lt.3) call errore('seriesbes','at least 3 points',1)
  npt2=npt/2+1

  !      xc(1)=0.5_dp*(fun(1)-xc(3)*r2(1)+fun(npt)-xc(3)*r2(npt))

  xc(3)=((fun(1)-fun(npt2))/(r(1)-r(npt2)) &
       -(fun(npt)-fun(npt2))/(r(npt)-r(npt2)) )/ ( r(1)-r(npt) )
  xc(1)=fun(1)
  xc(2)=( fun(npt)-fun(npt2) ) / (r(npt)-r(npt2)) -xc(3)*(r(npt)+ &
       r(npt2) )
  xc(4)=0.0_dp
  return
end subroutine seriesbes
