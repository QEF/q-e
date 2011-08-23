!
!---------------------------------------------------------------
subroutine drhoofv(drho,v)
!---------------------------------------------------------------
   !  
   use constants, only: e2
   use kinds,     only: DP
   use ld1inc,    only: nwf, zed, oc, psi, ll, enl, isw, nspin, vpot, grid
   use radial_grids, only: ndmx
   implicit none
   !
   ! I/O variables
   !
   real (DP) :: drho(ndmx,2), v(ndmx,2)
   !
   ! local variables
   !
   real (DP) :: dvy(ndmx), dchi(ndmx), int_0_inf_dr, wrk(ndmx)
   real (DP) :: ze2, fac
   integer ::  nu, is

   ze2 = -zed*e2

   drho = 0.d0

   do nu=1,nwf
      is = isw(nu)
      dvy(1:grid%mesh) = v(1:grid%mesh,is)*psi(1:grid%mesh,1,nu)
      wrk(1:grid%mesh) = dvy(1:grid%mesh)*psi(1:grid%mesh,1,nu)
      fac=int_0_inf_dr(wrk,grid,grid%mesh,2*ll(nu)+2)
      dvy(1:grid%mesh) = dvy(1:grid%mesh) - fac*psi(1:grid%mesh,1,nu)
      call green(dchi,ll(nu),enl(nu),dvy,psi(1,1,nu),vpot(1,is),ze2)
      drho(1:grid%mesh,is) = drho(1:grid%mesh,is) + &
                        2.0d0*oc(nu)*psi(1:grid%mesh,1,nu)*dchi(1:grid%mesh)
   end do
!
   return
end subroutine drhoofv
