!--------------------------------------------------------------------
subroutine dvex(nu,dvy)
   !--------------------------------------------------------------------
   !
   USE kinds,     ONLY: DP
   USE constants, ONLY: electron_charge_squared => e2
   USE ld1inc,    ONLY: nwf, &
                        psi, &
                        get_angular_momentum => ll, &
                        occupation => oc,&
                        sl3, &
                        nspin, &
                        get_spin => isw, &
                        grid
   use radial_grids, only: ndmx, hartree

   implicit none
   !
   ! I/O variables
   !
   integer,intent(in) :: nu
   real (DP),intent(out) :: dvy(ndmx)



   !
   ! local variables
   !
   integer :: i, mu, l0, l1, l2, l3
   real (DP) :: wrk(ndmx), wrk1(ndmx)
   real (DP) :: fac, sss, &
              ocs, &
              doc, &
              half
   !
   do i=1,grid%mesh
      dvy(i)=0.0d0
   end do
   l1 = get_angular_momentum(nu)
   half = 2.d0 * l1 + 1.d0


   ! iterate over all the KS wavefunctions
   ! get_spin(i) gives the spin of the i wavefunction
   do mu = 1,nwf
      !
      ! only wfc with the same spin contribute to exchange term
      !
      ! check for the spin of the wavefunction
      if (get_spin(mu) /= get_spin(nu) ) cycle  ! skip to next
      ocs = occupation(mu) * (0.5d0 * nspin)
!      write (*,*) mu, oc(mu), ocs
      if ( mu == nu ) then
         doc = 0.d0
         !if(AND((l1 /= 0), (ocs > 0.d0))) then
         if((l1 /= 0).AND.(ocs > 0.d0)) then
           i = int(ocs)
           doc = (i*(2.d0*ocs-i-1.d0)/(half-1.d0) - ocs*ocs/half) * half/ocs
         end if
         ocs = ocs + doc
      !   if (doc /= 0.d0) write (*,*) "DOC ",nu, doc
      end if
!
      l2 = get_angular_momentum(mu)
      l0 = abs(l1 - l2)
      wrk = psi(:,1,mu)*psi(:,1,nu)
      ! do i = 1,grid%mesh

      !    wrk(i) = psi(i,1,mu)*psi(i,1,nu) ! \varphi^{*}_{mu \sigma}\varphi^{*}_{nu \sigma}
      ! end do
      do l3 = l0, l1 + l2
         sss = sl3(l1,l2,l3) ! some kind of symbol
!         write (*,*) l1,l2,l3,sss
         if (abs(sss).gt.1.0d-10) then
            ! call the hartree utility function
            ! from radial_grids
            ! the output goes to vh
            call hartree(l3,l1+l2+2,grid%mesh, grid, wrk, vh = wrk1)

            fac =-electron_charge_squared*ocs*sss/2.0d0

            do i=1,grid%mesh
               dvy(i)= dvy(i) + fac*wrk1(i)*psi(i,1,mu)

            end do
         end if
      end do




!- spurious hartree part

      if (mu == nu ) then
          call hartree(0,2,grid%mesh,grid,wrk,wrk1)
          fac = doc  * electron_charge_squared
          do i=1,grid%mesh
             dvy(i)= dvy(i) + fac*wrk1(i)*psi(i,1,mu)
          end do
      end if

   end do

end subroutine dvex
