!   
!---------------------------------------------------------------
subroutine drho0ofvx(drho,dchi0)
!---------------------------------------------------------------
   use constants, only: e2
   use kinds,     only: DP
   use ld1inc,    only: ndm, mesh, nwf, nwfx, zed, oc, psi, ll, enl, &
                        vpot, isw, &
                        r, r2, dx
   implicit none
   !
   ! I/O variables
   !
   real(kind=DP) :: drho(ndm,2), dchi0(ndm,nwfx)
   ! local variables
   real(kind=DP) :: dvy(ndm), dchi(ndm), int_0_inf_dr,wrk(ndm)
   real(kind=DP) :: ze2, fac, ocs
   integer :: i, nu, is

   ze2 = - zed * e2
   
   drho = 0.d0

   do nu=1,nwf
      do i=1,mesh
         dvy(i)=dchi0(i,nu)
         wrk(i)=dvy(i)*psi(i,1,nu)
      end do
      fac=int_0_inf_dr(wrk,r,r2,dx,mesh,2*ll(nu)+2)
      do i=1,mesh
         dvy(i) = dvy(i) - fac*psi(i,1,nu)
      end do
      is = isw(nu)
      call green(dchi,ll(nu),enl(nu),dvy,psi(1,1,nu),vpot(1,is),ze2)
      do i=1,mesh
         drho(i,is)=drho(i,is) + 2.0d0*oc(nu)*psi(i,1,nu)*dchi(i)
      end do
   end do
!
   return
end subroutine drho0ofvx
