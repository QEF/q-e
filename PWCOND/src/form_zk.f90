!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine form_zk(n2d, nrzp, zkr, zk, e, tpiba)
!
! To construct complex wavevectors zk=sqrt(e-E_n)
! for an energy e from eigenvalues E_n of 2d problem
!
  USE kinds, only : DP
  implicit none
  integer :: nrzp, n2d, n, k
  real(DP) :: zkr(n2d,nrzp), e, ed, tpiba
  complex(DP) :: zk(n2d,nrzp)

  do k=1, nrzp
    do n=1, n2d
      ed = e-zkr(n,k)
      zk(n,k)=SQRT(CMPLX(ed,0.d0,kind=DP))/tpiba
    enddo
  enddo

  return
end subroutine form_zk

