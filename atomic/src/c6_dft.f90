!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#undef DEBUG
!---------------------------------------------------------------
subroutine c6_dft (mesh, zed, grid)
   !--------------------------------------------------------------------
   !
   use kinds,      only : DP
   use constants,  only : e2, pi, fpi, BOHR_RADIUS_ANGS
   use ld1inc,     only : lsd, nwf, oc, nn, ll, isw, psi, enl, vpot,vxt,vh, &
                          enne, latt, rho
   use radial_grids, only: radial_grid_type, ndmx
   !
   implicit none
   !
   ! I/O variables
   !
   type(radial_grid_type), intent(in) :: grid
   integer mesh , mesh_save
   real(DP) :: zed
   !
   ! local variables
   !
   logical :: csi, l_add_tf_term
   real(DP) :: vnew(ndmx,2), rhoc1(ndmx), ze2, fac, vme(ndmx)
   real(DP) :: rho_save(ndmx,2)
   real(DP) :: error, error2, e, charge, beta, u, alpha, dalpha, c6, du1, &
               du2, factor, thresh
   real(DP), allocatable :: y(:), yy(:), sqr(:)
   real(DP), allocatable :: dvpot(:), dvscf(:), drho(:), dvhx(:), dvxc(:), pp(:)
   complex(DP), allocatable :: dy(:), drho_old(:)
   integer i, is, n, l, iu, Nu, Nc, counter, nstop, nerr

   allocate ( y(mesh),yy(mesh),sqr(mesh) )
   allocate ( dvpot(mesh),dvscf(mesh),drho(mesh),dvhx(mesh),dvxc(mesh),pp(mesh) )
   allocate ( dy(mesh), drho_old(mesh) )
   !
   write(6,'(/,/,/,5x,20(''-''),'' Compute C6 from polarizability.'',10(''-''),/)')
   !
   if (mesh.ne.grid%mesh) call errore('c6_dft',' mesh dimension is not as expected',1)
   counter = 1
   do i = 1, mesh
      if (rho(i,1) .gt. 1.0d-30) counter = counter + 1
   enddo
   mesh_save = mesh
   mesh = counter

   if (lsd .ne. 0) call errore ('c6_dft', 'implemented only for non-magnetic ions', lsd) 
   csi = .true.
   do i = 1, nwf
      csi = csi .and. ( ((ll(i).eq.0) .and. (oc(i).eq.2 )) .or. &
                        ((ll(i).eq.1) .and. (oc(i).eq.6 )) .or. & 
                        ((ll(i).eq.2) .and. (oc(i).eq.10)) .or. &
                        ((ll(i).eq.3) .and. (oc(i).eq.14)) )
   enddo
   if (.not. csi) call errore ('c6_dft', 'implemented only for closed-shell ions', 1)
! 
 
   n = 1
   l = 0
   e = -1.d-7
   charge=0.d0
   ze2 = - zed * e2
   thresh = 1.d-10
!
   rho_save =  rho
   rho=0.0_dp
   do n=1,nwf
      do i=1,mesh
         rho(i,isw(n))=rho(i,isw(n))+oc(n)*(psi(i,1,n)**2+psi(i,2,n)**2)
      enddo
   enddo

   error = 0.d0
   do i=1,mesh
      error = error + abs( rho(i,1)-rho_save(i,1) ) * grid%r2(i) * grid%dx
      error = error + abs( rho(i,2)-rho_save(i,2) ) * grid%r2(i) * grid%dx
   end do

   if (error > 1.d-8) then
      write (*,*) error
      call errore('c6_dft','charge density rho from last vnew is inaccurate',1)
   end if

   rhoc1=0.d0
   call new_potential(ndmx,mesh,grid,zed,vxt,lsd,.false.,latt,enne,rhoc1,rho,vh,vnew,0)
   error = 0.d0
   do i=1,mesh
      error = error + abs( vpot(i,1)-vnew(i,1) ) * grid%r2(i) * grid%dx
      error = error + abs( vpot(i,2)-vnew(i,2) ) * grid%r2(i) * grid%dx
   end do
   write (*,*) "Vpot-Vnew", error

   nerr = 0
   do n=1,nwf
      if (oc(n) >= 0.0_dp) then
         is=isw(n)
         call ascheq (nn(n),ll(n),enl(n),mesh,grid,vnew(1,is),ze2,&
                      thresh,psi(1,1,n),nstop)
         nerr = nerr + nstop
         write (*,'(4i3,2f10.5,i5)') n, nn(n),ll(n),isw(n),oc(n),enl(n),nstop
      else
         enl(n)=0.0_dp
         psi(:,:,n)=0.0_dp
      end if
   end do

!  from now on rho is the REAL rho w/o the volume element
   do i = 1, mesh
      rho(i,1) = rho(i,1) / (fpi*grid%r(i)**2)
   end do
!
! initialize external perturbation (electric field)
!
   call init_dpot(grid%r,mesh,dvpot)
!
! derivative of xc-potential
!
   call dvxc_dn(mesh, rho, dvxc)
!
!write(*,'(1PE20.12)')sum(abs(dvxc))
!stop
   write(6,'(5x,''Frequency dependent polarizability is written into freq-pol.dat'',/)')

   c6    = 0.0d0
   alpha = 0.0d0

   open(1, file = 'freq-pol-dft.dat')
   write (1,'(15x,"    u          alpha(angstrom)       alpha(a.u.)  ",/)')
   !
   Nu  = 230
   Nc  = 50
   du1 = 0.1d0
   du2 = 0.25d0
   u   = -du1
   !
   do iu=0,Nu
      !
      if (iu .le. 50) then
         u = u + du1
      else
         u = u + du2
      endif
      !
      if (iu.eq.0) then
         do i=1,mesh
            dvscf(i) = dvpot(i)
            drho_old(i) = 0.d0
         end do 
      end if
      beta = 0.05
      dalpha = 1.0d+99
      alpha = 0.d0
      counter = 0
      do while (dalpha > 1.d-9)
         counter =  counter + 1
         !
         ! solve Sternheimer equation for the auxiliary wavefunction
         !
         drho = 0.d0
         do n=1,nwf
            do l = 1 + ll(n), max( 1 - ll(n), 0 ), - 2
!               write (*,*) l, ll(n)
               y(1:mesh) = psi(1:mesh,1,n)/grid%sqr(1:mesh)
               vme(:) = vnew(:,isw(n)) - enl(n)
               call sternheimer(u,l,ll(n),mesh,grid%dx,grid%r,grid%sqr,grid%r2,vme,zed,y,dvscf,dy)
               fac = 2.0d0 * (2.d0 * ll(n) + 1.d0 )
               if (ll(n)==1 .and. l==2) fac = fac * 2.d0/3.d0
               if (ll(n)==1 .and. l==0) fac = fac * 1.d0/3.d0
               if (ll(n)==2 .and. l==3) fac = fac * 3.d0/5.d0
               if (ll(n)==2 .and. l==1) fac = fac * 2.d0/5.d0
               call inc_drho_of_r(mesh, grid%dx, grid%r, grid%r2, y, dy, fac, drho)
#ifdef DEBUG
         write (*,*) "========================", n, l
         write (*,*) "y(1:3)"
         write (*,*) y(1:3)
         write (*,*) "dy(1:3)"
         write (*,*) dy(1:3)
         write (*,*) "drho(1:3)"
         write (*,*) drho(1:3)
#endif
             
            end do
         end do
#ifdef DEBUG
         write (*,*) "========================"
         write (*,*) "drho(1:3)"
         write (*,*) drho(1:3)
         write (*,*) "drho(20:22)"
         write (*,*) drho(20:22)
         write (*,*) "drho(40:42)"
         write (*,*) drho(40:42)
         write (*,*) "drho(mesh-2:mesh)"
         write (*,*) drho(mesh-2:mesh)
#endif 
         !
         ! compute dv of drho (w/o the TF term)
         !
         l_add_tf_term = .false.
         call dv_of_drho(mesh, grid%dx, grid%r,grid%r2,rho,drho,dvhx,dvxc,pp,l_add_tf_term)

#ifdef DEBUG
         write (*,*) "========================"
         write (*,*) "dvhx(1:3)"
         write (*,*) dvhx(1:3)
         write (*,*) "dvhx(20:22)"
         write (*,*) dvhx(20:22)
         write (*,*) "dvhx(40:42)"
         write (*,*) dvhx(40:42)
         write (*,*) "dvhx(mesh-2:mesh)"
         write (*,*) dvhx(mesh-2:mesh)
         write (*,*) "========================"
         write (*,*) "pp(1:3)"
         write (*,*) pp(1:3)
         write (*,*) "pp(20:22)"
         write (*,*) pp(20:22)
         write (*,*) "pp(40:42)"
         write (*,*) pp(40:42)
         write (*,*) "pp(mesh-2:mesh)"
         write (*,*) pp(mesh-2:mesh)
#endif 
         !
         ! mix
         !
         error = 0.d0
         error2 = 0.d0
         do i=1,mesh
            dvscf(i) = dvscf(i) + beta * (dvpot(i)+dvhx(i) -dvscf(i))
            error = error + abs (drho(i) -drho_old(i))
            error2 = error2 + abs (drho(i) -drho_old(i))* grid%r(i) * grid%dx
            drho_old(i) = drho(i)
         end do 
         dalpha = abs(alpha + pp(mesh)) 
         alpha = -pp(mesh)
!         write (*,'(4e16.6)') alpha, dalpha, error, error2
        
      end do

      write (1,'(17x, f8.4, 3x, 1PE14.6, 9x, 1PE14.6)') u, pp(mesh)*BOHR_RADIUS_ANGS**3, pp(mesh)
      if (iu .eq. 0) & 
      write (6,'(5x, "Static polarizability: ", f10.5, " (in angstrom^3)   --->", f10.5,&
                & "  (in e^2a0^3)")') pp(mesh)*BOHR_RADIUS_ANGS**3, pp(mesh)

      if (iu .eq. 0)                  factor = 0.5d0 * du1
      if (iu .gt. 0 .and. iu .lt. Nc) factor = du1
      if (iu .eq. Nc)                 factor = 0.5d0 * ( du1 + du2)
      if (iu .gt. Nc .and. iu .lt. Nu) factor = du2
      if (iu .eq. Nu)                 factor = 0.5d0 * du2
      c6 = c6 + factor*alpha*alpha

   end do

   c6 = c6 * 3.d0 / pi 

   write (*,'(/, 5x, a, f12.6)') "C6 coefficient in units [e2*a0**5]", c6/e2
   !
   write(6,'(/,5x,20(''-''),'' End of C6 calculation '',20(''-''),/)')

   deallocate ( dy )
   deallocate ( y, yy, sqr )
   deallocate ( dvpot, dvscf, drho, dvhx, pp )
   
   return
end subroutine c6_dft
  
!--------------------------------------------------------------------
subroutine inc_drho_of_r(mesh, dx, r, r2, y, dy, fac, drho)
   !--------------------------------------------------------------------
   ! compute the first order variation of the density from
   ! the zeroth and first order auxiliary wavefunctions y and dy
   !
   use constants, only : e2, pi, fpi
   implicit none
   !
   ! I/O vaiables
   !
   integer mesh
   real (kind=8) :: dx, fac, r(mesh), r2(mesh), y(mesh), drho(mesh)
   complex (kind=8) :: dy(mesh)
   ! local variables
   integer i

   do i=1,mesh
      drho(i) = drho(i) + fac * 2.d0 * y(i) * real(dy(i)) * r(i) / (fpi*r2(i))
   end do

   return
end subroutine

