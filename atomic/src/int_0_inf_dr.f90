!
! Copyright (C) 2004-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
function int_0_inf_dr(f,grid,mesh,nst)
  !---------------------------------------------------------------
  !
  !      integral of f from 0 to infinity
  !      f is given on a logarithmic mesh. 
  !      f(r) is assumed to be proportional to r**nst for small r
  !
  use kinds, only : DP
  use radial_grids, only: radial_grid_type, series
  implicit none
  !
  ! I/O variables
  !
  integer, intent(in) :: mesh, nst
  real(DP), intent(in):: f(mesh)
  type(radial_grid_type), intent(in) :: grid
  real(DP) :: int_0_inf_dr
  !
  ! local variables
  !
  real(DP):: fs(4), b(4), sum1
  integer :: i
  !
  ! series development: contribution for small r
  !
  if (mesh > grid%mesh) &
     call errore('int_0_inf_dr','value of mesh is larger than expected',mesh)

  do i=1,4
     fs(i)=f(i)/grid%r(i)**nst
  end do
  call series(fs,grid%r,grid%r2,b)
  int_0_inf_dr = ( b(1)/(nst+1) + grid%r(1)* &
                  (b(2)/(nst+2) + grid%r(1)*b(3)/(nst+3)) ) * grid%r(1)**(nst+1)
  !
  ! simpson integration (logarithmic mesh: dr ==> r dx)
  !
  sum1=0.0_DP
  do i=1,mesh-2,2
     sum1 = sum1 + f(i)*grid%r(i)+4.0_dp*f(i+1)*grid%r(i+1)+f(i+2)*grid%r(i+2)
  end do
  int_0_inf_dr = int_0_inf_dr + sum1*grid%dx/3.0_DP

  return
end function int_0_inf_dr
