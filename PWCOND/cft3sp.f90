!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine cft3sp (f,n1,n2,n3,nx1,nx2,nx3,sign)
!----------------------------------------------------------------------
!
! Just usual FFT (even in MPI running)
!
  implicit none
  integer :: n1, n2, n3, nx1, nx2, nx3, sign
  complex(kind(0.d0)) ::  f(nx1*nx2*nx3)

  if (sign.eq.1) then
     call cft_3(f,n1,n2,n3,nx1,nx2,nx3,2, 1)
  else if (sign.eq.-1) then
     call cft_3(f,n1,n2,n3,nx1,nx2,nx3,2,-1)
  endif
  return
end subroutine cft3sp
