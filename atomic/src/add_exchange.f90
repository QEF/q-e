!
!--------------------------------------------------------------------
subroutine add_exchange ( energy )
!--------------------------------------------------------------------
#undef DEBUG
   !  
   use io_global, only: stdout
   use kinds,  only: DP
   use constants, only: e2
   use ld1_parameters, only : nwfx
   use ld1inc, only: nwf, oc, psi, vx, sl3, ll, nn, enl, el, &
                     isw, nspin, enzero, grid
   use radial_grids, only: ndmx, hartree
   implicit none
   !
   ! I/O variables
   !
   real (DP) :: energy
   !
   ! local variables
   !
   integer :: i, l0, l1, l2, l3, nu, mu, nst, is, half
   real (DP) :: ex_hf, ocs, doc, sss, fac, sxc, sxc1
   real (DP) :: wrk(ndmx), wrk1(ndmx), wrk2(ndmx), int_0_inf_dr, enzhf(nwfx)
!
   ex_hf = 0.0
   do nu=1,nwf
      is = isw(nu)
      l1 = ll(nu)
      half = 2.d0 * l1 + 1.d0
      sxc = 0.0d0
      do mu=1,nwf
!
! only wfc with the same spin contribute to exchange term
!
         if (isw(mu) /= is) cycle
         ocs = oc(mu) * (0.5d0 * nspin )
         if ( mu == nu ) then
            doc = 0.d0
            if( (l1 /= 0) .and. (ocs > 0.d0) ) then
              i = int(ocs)
              doc = (i*(2.d0*ocs-i-1.d0)/(half-1.d0) - ocs*ocs/half) * half/ocs
            end if
            ocs = ocs + doc
!            if (doc /= 0.d0) write (*,*) "DOC ",nu, doc
         end if

         l2 = ll(mu)
         l0=abs(l1-l2)
         do i=1,grid%mesh
            wrk(i) = psi(i,1,mu)*psi(i,1,nu)
            wrk1(i)= 0.0d0
         end do
         do l3=l0,l1+l2
            sss = sl3(l1,l2,l3)
            if (abs(sss).gt.1.0d-10) then
               call hartree(l3,l1+l2+2,grid%mesh,grid,wrk,wrk2)
               fac = -e2*ocs*sss/2.0d0
               do i=1,grid%mesh
                  wrk1(i)= wrk1(i) + fac*wrk2(i)*wrk(i)
               end do
            end if
         end do
!- spurious hartree part 
          if (mu.eq.nu) then
            call hartree(0,2,grid%mesh,grid,wrk,wrk2)
            fac = doc*e2
            do i=1,grid%mesh
               wrk1(i) = wrk1(i) + fac*wrk2(i)*wrk(i)
            end do
          end if
!
         nst = 2 * min(l1,l2) + 2
         sxc = sxc + int_0_inf_dr(wrk1,grid,grid%mesh,nst)
      end do 
!
      do i=1,grid%mesh
         wrk1(i) = vx(i,is)*psi(i,1,nu)*psi(i,1,nu)
      end do
      sxc1 = int_0_inf_dr(wrk1,grid,grid%mesh,2*ll(nu)+2)
      ex_hf = ex_hf + 0.5d0*oc(nu)*sxc
      enzhf(nu)=sxc1-sxc
      if(oc(nu)>0) enzero(is) = enzhf(nu)
   end do 

   energy = energy + ex_hf

#ifdef DEBUG
   write (*,*) enzero
   write(stdout,1100) (nn(nu),ll(nu),el(nu),oc(nu),enl(nu)-enzero(isw(nu)), &
                  enl(nu)-enzhf(nu), nu=1,nwf)

1100 format(4x,2i2,5x,a2,'(',f5.2,')',2f11.4)
#endif

   return

end subroutine add_exchange
