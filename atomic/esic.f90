!
!---------------------------------------------------------------
subroutine esic
!---------------------------------------------------------------
!
use ld1inc
implicit none
! output
! local
integer:: n, i
real(kind=dp), parameter :: fourpi=4.0d0*3.141592653589793d0
real(kind=dp) :: int_0_inf_dr,vxup,vxdw,vcup,vcdw,ex,ec,deksic  
real(kind=dp) :: work1(ndm),v(ndm),vsic(ndm)
real(kind=dp) :: egc(ndm)
external int_0_inf_dr
!
deksic = 0.0d0
dhrsic = 0.0d0
dxcsic = 0.0d0
do n=1,nwf
   call sic_correction(n,v,vsic,work1)
   do i=1,mesh
      if (rel.eq.2) then
         v(i)=v(i)*(psi_dir(i,1,n)**2+psi_dir(i,2,n)**2)
         vsic(i)=vsic(i)*(psi_dir(i,1,n)**2+psi_dir(i,2,n)**2)
      else
         v(i)=v(i)*psi(i,n)**2
         vsic(i)=vsic(i)*psi(i,n)**2
      endif
   enddo
   deksic = deksic +  &
     &      oc(n)*int_0_inf_dr(vsic,r,r2,dx,mesh,2*(ll(n)+1))
   dhrsic = dhrsic - 0.5d0*oc(n)*int_0_inf_dr(v    ,r,r2,dx,mesh,2)
   dxcsic = dxcsic - oc(n)*int_0_inf_dr(work1,r,r2,dx,mesh,2)
enddo
!
ekin=ekin+deksic
ehrt=ehrt
ecxc=ecxc+dxcsic+dhrsic
etot=etot+dhrsic+dxcsic+deksic
!
return
end
