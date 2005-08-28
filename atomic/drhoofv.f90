!
!---------------------------------------------------------------
subroutine drhoofv(drho,v)
!---------------------------------------------------------------
   !  
   use constants, only: e2
   use kinds,     only: DP
   use ld1inc,    only: ndm, mesh, nwf, zed, oc, psi, ll, enl, &
                        isw, nspin, vpot, &
                        r, r2, dx
   implicit none
   !
   ! I/O variables
   !
   real (DP) :: drho(ndm,2), v(ndm,2)
   !
   ! local variables
   !
   real (DP) :: dvy(ndm), dchi(ndm), int_0_inf_dr, wrk(ndm)
   real (DP) :: ze2, fac
   integer :: i, nu, is

   ze2 = -zed*e2

   drho = 0.d0

   do nu=1,nwf
      is = isw(nu)
      dvy(1:mesh) = v(1:mesh,is)*psi(1:mesh,1,nu)
      wrk(1:mesh) = dvy(1:mesh)*psi(1:mesh,1,nu)
      fac=int_0_inf_dr(wrk,r,r2,dx,mesh,2*ll(nu)+2)
      dvy(1:mesh) = dvy(1:mesh) - fac*psi(1:mesh,1,nu)
      call green(dchi,ll(nu),enl(nu),dvy,psi(1,1,nu),vpot(1,is),ze2)
      drho(1:mesh,is) = drho(1:mesh,is) + &
                        2.0d0*oc(nu)*psi(1:mesh,1,nu)*dchi(1:mesh)
   end do
!
   return
end subroutine drhoofv
