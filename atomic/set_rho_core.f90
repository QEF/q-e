!
!-----------------------------------------------------------------------
      subroutine set_rho_core
!-----------------------------------------------------------------------
!
!      input : all-electron wavefunctions + valence states
!
use ld1inc
implicit none

real(kind=dp) :: int_0_inf_dr, drho, const, br1, br2, pi, &
       eps1, eps2, br12, a, b, eps12, totrho, rhov(ndm),    &
       rhoco(ndm)
integer :: i, ik, n, ns, ios

write(6,'(/,5x,'' Computing core charge for nlcc: '')')
pi = 4.d0*atan(1.d0)
!
!      calculates core charge density
!
do n=1,mesh
   rhoc(n) = 0.d0
   rhov(n)=0.d0
   do ns=1,nwf
      if (core_state(ns)) then
         rhoc(n)=rhoc(n)+oc(ns)*psi(n,ns)**2
         rhoco(n)=rhoc(n)
      else
         rhov(n) = rhov(n) + oc(ns)*psi(n,ns)**2
      endif
   enddo
enddo
!
!      "equivalent" core charge density for nonlinear core correction:
!      rhoc(r) = core charge        if valence charge < fac*core charge
!      rhoc(r) = r^2 a sin(b r)/r   otherwise      (fac =0.5-1.0)
!
do ik=1,mesh
   if (r(ik).gt.rcore) go to 100
enddo
call errore('set_rho_core','rcore too big',-1)
return
100  rcore=r(ik)
drho = ( rhoc(ik+1)/r2(ik+1) - rhoc(ik)/r2(ik) ) / dx / r(ik)
!
!      true_rho = rhoc(r)/r**2/4 pi !  factor 1/r from logarithmic mesh
!
if (drho.gt.0.d0) then
   call errore('set_rho_core','d rho/ d r .gt. 0',-1)
   return
endif
const= r(ik)*drho / ( rhoc(ik)/r2(ik) ) + 1.d0
if (const.gt.0.d0) then
   br1 = 0.00001d0
   br2 = pi/2.d0-0.00001d0
else
   br1 = pi/2.d0+0.00001d0
   br2 = pi
end if
do n=1, 15
   eps1 = br1 - const*tan(br1)
   eps2 = br2 - const*tan(br2)
   br12 = (br1+br2)/2.d0
   eps12 = br12 - const*tan(br12)
   if(eps1*eps12.lt.0.d0) then
      br2 = br12
   else if(eps12*eps2.lt.0.d0) then
      br1 = br12
   else
      call errore('set_rho_core','error in bisection',n)
   end if
end do
b = br12/r(ik)
a = ( rhoc(ik)/r2(ik) ) * r(ik)/sin(br12)
do n=1,ik
   rhoc(n) = a*sin(b*r(n))/r(n) * r2(n)
end do   
write(6,'(/,5x,''  r > '',f4.2,'' : true rho core'')') r(ik)
write(6,110) r(ik), a, b
110   format (5x, '  r < ',f4.2,' : rho core = a sin(br)/r', &
                    '    a=',f7.2,'  b=',f7.2/)
file_core=' '
write(6,*) '***',file_core,'***'
if (file_core .ne. ' ') then
   open(unit=26,file=file_core, status='unknown', iostat=ios, &
                        err=300 )
300      call errore('set_rho_core','opening file '//file_core,abs(ios))
   do n=1,mesh
      write(26,'(4f20.10)') r(n),rhoc(n),rhov(n),rhoco(n)
   enddo
   close(26)
endif

totrho = int_0_inf_dr(rhoc,r,r2,dx,mesh,2)
write(6,'(13x,''integrated charge : '',f6.2)')  totrho
return
end
