!---------------------------------------------------------------
subroutine vxc_t(rho,rhoc,lsd,vxc)
  !---------------------------------------------------------------
  !
  !  this fuction returns the XC potential in LDA or LSDA approximation
  !

  use kinds, only : DP
  implicit none
  integer:: lsd
  real(kind=dp):: vxc(2), rho(2),rhoc,arho,zeta
  real(kind=dp):: vx(2), vc(2), ex, ec
  !
  real(kind=dp), parameter :: e2=2.0_dp, eps=1.e-30_dp

  vxc(1)=0.0_dp
  if (lsd.eq.1) vxc(2)=0.0_dp

  if (lsd.eq.0) then
     !
     !     LDA case
     !
     arho=abs(rho(1)+rhoc)
     if (arho.gt.eps) then      
        call xc(arho,ex,ec,vx,vc)
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
           write(6,*) 'zeta= me', zeta, rho(1),rho(2),rhoc
        else
           call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
           vxc(1) = e2*(vx(1)+vc(1))
           vxc(2) = e2*(vx(2)+vc(2))
        endif
     endif
  endif

  return
end subroutine vxc_t
