! FOR GWW
!
! Author: P. Umari
!
!----------------
subroutine local_wannier(matx,maty,matz, spreadtot)
!----------------
! #ifdef __GWW

! this subroutine calculates the centers of the wannier functiona
! assuning that (matx,maty,matz) are on the basis of the
! wannier functions, the spread of the wannier functions
! and the total spread
!
! uses Resta's spread definition
!   x_i = (L_x/2\pi)*Im ln matx_{i,i}
!
!   sigma_x_i = (-L_x ^2/(4\pi^2)*ln|matx_{ii}|^2
!   sigma_i = sigma_x_i + sigma_y_i + sigma_z_i
!
! only orthorombic cells


  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : nbnd
  USE cell_base, ONLY: at, alat, tpiba
  USE wannier_gw, ONLY : wannier_centers
  USE constants, ONLY : pi


  implicit none


  COMPLEX(kind=DP), INTENT(in)  :: matx(nbnd,nbnd)
  COMPLEX(kind=DP), INTENT(in)  :: maty(nbnd,nbnd)
  COMPLEX(kind=DP), INTENT(in)  :: matz(nbnd,nbnd)
  REAL(kind=DP), INTENT(out) :: spreadtot


 !local variables


  real(kind=dp):: lx,ly,lz!length of simulation cell sizes in a.u. (good for casino)
  real(kind=dp):: spreadx, spready,spreadz
  integer i

  spreadtot=0.d0

  lx=at(1,1)*alat
  ly=at(2,2)*alat
  lz=at(3,3)*alat


  do i=1,nbnd

     wannier_centers(1,i)=-(lx/2.d0/pi)*aimag(log(matx(i,i)))!ATTENZIONE al segno
     wannier_centers(2,i)=-(ly/2.d0/pi)*aimag(log(maty(i,i)))
     wannier_centers(3,i)=-(lz/2.d0/pi)*aimag(log(matz(i,i)))

     spreadx = (-lx**2/4.d0/pi**2.)*log(abs(conjg(matx(i,i))*matx(i,i)))
     spready = (-ly**2/4.d0/pi**2.)*log(abs(conjg(maty(i,i))*maty(i,i)))
     spreadz = (-lz**2/4.d0/pi**2.)*log(abs(conjg(matz(i,i))*matz(i,i)))

     spreadtot = spreadtot + spreadx + spready + spreadz
  enddo
! #endif
  return
end subroutine local_wannier
