!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine form_zk(n2d, nrzp, zkr, zk, e, tpiba)
!
! Just to construct complex zk=sqrt(e-E_n) for energy e  
! from eigenvalues E_n found for some initial energy
!
  use parameters, only : DP 
  implicit none
  integer :: nrzp, n2d, n, k 
  real(kind=DP) :: zkr(n2d,nrzp), e, ed, tpiba
  real(kind=DP), parameter :: eps=1.d-4
  complex(kind=DP) :: zk(n2d,nrzp)

  do k=1, nrzp
    do n=1, n2d
      ed = e-zkr(n,k)
      zk(n,k)=SQRT(CMPLX(ed))/tpiba
    enddo
  enddo  

  return
end subroutine form_zk            

