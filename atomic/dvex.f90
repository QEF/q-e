!--------------------------------------------------------------------
subroutine dvex(nu,dvy)
!--------------------------------------------------------------------
   !  
   USE kinds,     ONLY: DP
   USE constants, ONLY: e2
   USE ld1inc,    ONLY: nwf, psi, ll, oc, sl3, nspin, isw, grid
   use radial_grids, only: ndmx, hartree

   implicit none
   ! 
   ! I/O variables
   !
   integer :: nu
   real (DP) :: dvy(ndmx)
   !
   ! local variables
   !
   integer :: i, mu, l0, l1, l2, l3
   real (DP) :: wrk(ndmx), wrk1(ndmx)
   real (DP) :: fac, sss, ocs, doc, half
   !
   do i=1,grid%mesh
      dvy(i)=0.0d0
   end do
   l1 = ll(nu)
   half = 2.d0 * l1 + 1.d0
   do mu=1,nwf
!
! only wfc with the same spin contribute to exchange term
!
      if (isw(mu) /= isw(nu) ) cycle
      ocs = oc(mu) * (0.5d0 * nspin)
!      write (*,*) mu, oc(mu), ocs
      if ( mu == nu ) then
         doc = 0.d0
         if( (l1 /= 0) .and. (ocs > 0.d0) ) then
           i = int(ocs)
           doc = (i*(2.d0*ocs-i-1.d0)/(half-1.d0) - ocs*ocs/half) * half/ocs
         end if
         ocs = ocs + doc
!         if (doc /= 0.d0) write (*,*) "DOC ",nu, doc
      end if
!
      l2 = ll(mu)
      l0=abs(l1-l2)
      do i=1,grid%mesh
         wrk(i) = psi(i,1,mu)*psi(i,1,nu)
      end do
      do l3=l0,l1+l2
         sss = sl3(l1,l2,l3)
!         write (*,*) l1,l2,l3,sss
         if (abs(sss).gt.1.0d-10) then
            call hartree(l3,l1+l2+2,grid%mesh,grid,wrk,wrk1)
            fac =-e2*ocs*sss/2.0d0
            do i=1,grid%mesh
               dvy(i)= dvy(i) + fac*wrk1(i)*psi(i,1,mu)
            end do
         end if
      end do
!- spurious hartree part 
      if (mu == nu ) then
          call hartree(0,2,grid%mesh,grid,wrk,wrk1)
          fac = doc*e2
          do i=1,grid%mesh
             dvy(i)= dvy(i) + fac*wrk1(i)*psi(i,1,mu)
          end do
        end if
   end do 

  return
end subroutine dvex
