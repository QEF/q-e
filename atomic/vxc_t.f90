!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine vxc_t(rho,rhoc,lsd,vxc)
  !---------------------------------------------------------------
  !
  !  this function returns the XC potential in LDA or LSDA approximation
  !

  use io_global, only : stdout
  use kinds, only : DP
  use funct, only : xc, xc_spin
  implicit none
  integer:: lsd
  real(DP):: vxc(2), rho(2),rhoc,arho,zeta
  real(DP):: vx(2), vc(2), ex, ec
  !
  real(DP), parameter :: e2=2.0_dp, eps=1.e-30_dp

  vxc(1)=0.0_dp
  if (lsd.eq.1) vxc(2)=0.0_dp

  if (lsd.eq.0) then
     !
     !     LDA case
     !
     arho=abs(rho(1)+rhoc)
     if (arho.gt.eps) then      
        call xc(arho,ex,ec,vx(1),vc(1))
        vxc(1)=e2*(vx(1)+vc(1))
     endif
  else
     !
     !     LSDA case
     !
     arho = abs(rho(1)+rho(2)+rhoc)
     if (arho.gt.eps) then      
        zeta = (rho(1)-rho(2)) / arho
        if (abs(zeta).gt.1.0_dp) then 
           write(stdout,*) 'zeta, rho_up, rho_dw, rhoc', zeta, &
                           rho(1),rho(2),rhoc
        else
           call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
           vxc(1) = e2*(vx(1)+vc(1))
           vxc(2) = e2*(vx(2)+vc(2))
        endif
     endif
  endif

  return
end subroutine vxc_t
