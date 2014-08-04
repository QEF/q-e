!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#undef DEBUG
!---------------------------------------------------------------
subroutine c6_tfvw (mesh, zed, grid, rho_input)
   !--------------------------------------------------------------------
   !
   use kinds,      only : DP
   use constants,  only : e2, pi, fpi, BOHR_RADIUS_ANGS
   use ld1inc,     only : lsd, nwf, oc, ll
   use radial_grids, only: radial_grid_type
   !
   implicit none
   !
   ! I/O variables
   !
   type(radial_grid_type), intent(in):: grid
   integer mesh
   real (kind=8) :: rho_input(mesh)
   real (kind=8) :: zed, rho(mesh)
   !
   !real (kind=8) :: vw_lambda=9.0_dp ! describes exactly the low-q limit
   !real (kind=8) :: vw_lambda=3.0_dp ! interpolate between high and low q's
   real (kind=8) :: vw_lambda=1.0_dp ! describes exactly the high-q limit
   !
   ! local variables
   !
   logical :: csi, l_add_tf_term
   real (kind=8) :: error, error2, e, charge, beta, u, alpha, dalpha, c6, du1, &
                    du2, factor, ze2, thresh
   real (kind=8), allocatable :: veff(:), y(:), yy(:)
   real (kind=8), allocatable :: dvpot(:), dvscf(:), drho(:), dvhx(:), dvxc(:), pp(:)
   complex (kind=8), allocatable :: dy(:), drho_old(:)
   integer i, iter, n, l, ly, iu, Nu, Nc, counter, nstop, mesh_save

   allocate ( veff(mesh),y(mesh),yy(mesh))
   allocate ( dvpot(mesh),dvscf(mesh),drho(mesh),dvhx(mesh),dvxc(mesh),pp(mesh) )
   allocate ( dy(mesh), drho_old(mesh) )
   !
   write(6,'(/,/,/,5x,10(''-''),'' Compute C6 from polarizability with TFvW approx.'',10(''-''),/)')
   if (vw_lambda.ne.1.d0) write(6,*) " value of vw_lambda ", vw_lambda
   !
   if (mesh.ne.grid%mesh) call errore('c6_tfwv',' mesh dimension is not as expected',1)
   do i = 1, mesh
      rho(i) = rho_input(i) / (fpi*grid%r(i)**2)
   end do
   !
   counter = 1
   do i = 1, mesh
      if (rho(i) .gt. 1.0d-30) counter = counter + 1 
   enddo
   mesh_save = mesh
   mesh = counter
#ifdef DEBUG
   write (*,*) "mesh ", mesh
#endif
   !
   if (lsd .ne. 0) call errore ('c6_tfvw', 'implemented only for non-magnetic ions', lsd) 
   csi = .true.
   do i = 1, nwf
      csi = csi .and. ( ((ll(i).eq.0) .and. (oc(i).eq.2 )) .or. &
                        ((ll(i).eq.1) .and. (oc(i).eq.6 )) .or. & 
                        ((ll(i).eq.2) .and. (oc(i).eq.10)) .or. &
                        ((ll(i).eq.3) .and. (oc(i).eq.14)) )
   enddo
   if (.not. csi) call errore ('c6_tfvw', 'implemented only for closed-shell ions', 1)
!   rho = 0.d0
!   open (7,file='rho.out',status='unknown',form='formatted')
!   do i=1,mesh
!      read (7,'(P5E20.12)') r(i), rho(i), y(i), y(i), y(i)
!      write (6,'(P5E15.6)') r(i), rho(i)
!   end do
!   close (7)
!
! compute unperturbed effective potential
!
   call veff_of_rho(mesh,grid%dx,grid%r,grid%r2,rho,y,veff)
#ifdef DEBUG
   write (*,*) "veff(1:3)"
   write (*,*) veff(1:3)
   write (*,*) "veff(mesh-5:mesh)"
   write (*,*) veff(mesh-5:mesh)
#endif
!
! check that veff and y are what we think
!
   n = 1
   l = 0
   e = -1.d-7
   charge=0.d0
   ze2 = - zed * e2
   thresh = 1.d-14

   do i=1,mesh
      charge = charge + rho(i) * fpi * grid%r2(i) * grid%r(i) * grid%dx
   end do
!   call solve_scheq(n,l,e,mesh,dx,r,sqr,r2,veff,zed,yy)
   call ascheq (n, l, e, mesh, grid, veff, ze2, thresh, yy, nstop)

   error = 0.d0
   do i=1,mesh
      error = error + (y(i)-yy(i)/grid%sqr(i)*sqrt(charge))**2 * grid%r2(i) * grid%dx
   end do

#ifdef DEBUG
      write (*,*) "ascheq called with mesh"
      write (*,*) "nstop", nstop, e
      write (*,*) "error ",error
      write (*,*) "y(1:3)"
      write (*,*) y(1:3)
      write (*,*) "y(mesh-2:mesh)"
      write (*,*) y(mesh-2:mesh)
      write (*,*) "yy(1:3)"
      write (*,*) yy(1:3)
      write (*,*) "yy(mesh-2:mesh)"
      write (*,*) yy(mesh-2:mesh)
      write (*,*) grid%sqr(1:3)
      write (*,*) "sqrt(charge)", sqrt(charge)
#endif
   if (error > 1.d-8) then
      call errore('c6_tfvw','auxiliary funtions veff(r) and y(r) are inaccurate',1)
   end if
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

   open(1, file = 'freq-pol.dat')
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
         l = 1
         ly = 0
         call sternheimer(u*vw_lambda,l,ly,mesh,grid%dx,grid%r,grid%sqr,grid%r2,veff,zed,y,dvscf,dy)
         dy = dy*vw_lambda
         ! compute drho of r
         !
         call drho_of_r(mesh, grid%dx, grid%r, grid%r2, y, dy, drho)
#ifdef DEGUG
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
         ! compute dv of drho (including the TF term)
         !
         l_add_tf_term = .true.
         call dv_of_drho(mesh, grid%dx, grid%r,grid%r2,rho,drho,dvhx,dvxc,pp, l_add_tf_term)

#ifdef DEGUG
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
      write (6,'(5x, "Static polarizability: ", f10.5, " (in A^3) --->", f10.5,&
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
   deallocate ( veff, y, yy )
   deallocate ( dvpot, dvscf, drho, dvhx, pp )
   
   mesh = mesh_save
   return
end subroutine
  
!--------------------------------------------------------------------
subroutine veff_of_rho(mesh,dx,r,r2,rho,y,veff)
   !--------------------------------------------------------------------
   ! compute unperturbed auxiliary wavefunction y and
   ! the corresponding effective potential veff
   !
   use constants, only : e2, pi, fpi
   implicit none
   ! 
   ! I/O variables
   !
   integer mesh
   real (kind=8) :: dx, r(mesh), r2(mesh), rho(mesh), veff(mesh), y(mesh)
   ! 
   ! local variables
   !
   real (kind=8), allocatable :: vold(:)
   real (kind=8) :: dx2, error
   integer i, iter, k

!
! compute auxiliary wavefunction y
!
   do i=1,mesh
      y(i) = sqrt(rho(i)*r(i)*fpi)
   end do
!
! compute effective potential veff
!
   allocate (vold(mesh))

   do i=1,mesh
      vold(i) = 0.d0
   end do
   dx2= dx*dx
   error = 1.d0
   k = 0
   do while (error > 1.d-9) 
      !
      k=k+1
      !
      do i=2,mesh-1
         veff(i) = ( y(i+1)/y(i) + y(i-1)/y(i) -2.d0 )/dx2  &
                 - ( vold(i+1)*y(i+1)/y(i) + vold(i-1)*y(i-1)/y(i) -2.d0*vold(i) )/12.d0
      end do
      veff(1) = veff(2) + (veff(3)-veff(2))*(r(1)-r(2))/(r(3)-r(2)) 
      veff(mesh) = (y(mesh-1)/y(mesh) -2.d0 )/dx2 &
                 - (vold(mesh-1)*y(mesh-1)/y(mesh) -2.d0*vold(mesh) )/12.d0
!
! the routine that integrates the Sh.Eq. requires that v(mesh) is an upper bound
!
      veff(mesh) = max(veff(mesh),veff(mesh-1))
!
      error = 0.d0
      do i=1,mesh
         error = error + abs( veff(i) - vold(i) )
         vold(i) = veff(i)
      end do
      error = error / mesh
!      write (*,*) 'iteration # ', k, error
   end do

   deallocate (vold)
   !
   do i=1,mesh
      veff(i) = (veff(i) -0.25d0)/r2(i)
   end do
   !
!   open (7,file='veff.out',status='unknown',form='formatted')
!   write (7,*) "#  r(i),        rho(i),        y(i),          veff(i),       veff(i)*r(i) "
!   do i=0,mesh
!      write (7,'(P6E15.6)') r(i), rho(i), y(i), veff(i), veff(i)*r(i)
!   end do
!   close (7)
   !
   return
end subroutine
!
!--------------------------------------------------------------------
subroutine dv_of_drho(mesh,dx,r,r2,rho,drho,dvhx,dvxc,pp, l_add_tf_term)
   !--------------------------------------------------------------------
   use constants, only : e2, pi, fpi
!   use flags,      only : HartreeFock, rpa
   implicit none
   !
   ! I/O variables
   !
   logical l_add_tf_term
   integer mesh
   real (kind=8) :: dx, r(mesh), r2(mesh)
   real (kind=8) :: rho(mesh), drho(mesh), dvhx(mesh), pp(mesh), dvxc(mesh)
   !
   ! local variables
   !
   real (kind=8) :: dr3, kf2, charge
   real (kind=8), allocatable :: qq(:)
   integer i

   allocate (qq(mesh))

   do i=1,mesh 
      dr3   = fpi * r2(i) * r(i) * dx
      pp(i) = drho(i) * r(i)  * dr3 /3.d0
      qq(i) = drho(i) / r2(i) * dr3 /3.d0
   end do
   do i=2,mesh
      pp(i) = pp(i) + pp(i-1) 
   end do
!   write (*,*) "pp in dv_of_drho"
!   write (*,*) pp(1:6)
!   write (*,*) pp(mesh-5:mesh)
!   write (*,*) "qq -prima in dv_of_drho"
!   write (*,*) qq(1:6)
!   write (*,*) qq(mesh-5:mesh)
   do i=mesh-1,1,-1
      qq(i) = qq(i) + qq(i+1)
   end do
!   write (*,*) "qq -dopo in dv_of_drho"
!   write (*,*) qq(1:6)
!   write (*,*) "r2 in dv_of_drho"
!   write (*,*) r2(1:6)
!   write (*,*) "r in dv_of_drho"
!   write (*,*) r(1:6)
   do i=1,mesh
      dvhx(i) = e2 * ( pp(i) / r2(i) + qq(i) * r(i) ) ! Hartree term
   end do

!   write (*,*) "Hartree in dv_of_drho"
!   write (*,*) dvhx(1:6)
!   write (*,*) dvhx(mesh-5:mesh)
! add TF term
   if (l_add_tf_term) then
      do i=1,mesh
         kf2 = ( 3.d0*pi*pi*rho(i) )**(2.d0/3.d0)
         dvhx(i) = dvhx(i) + e2/3.d0* kf2 / rho(i) * drho(i)
      end do
   end if
!
!add xc term
!         
!   write (*,*) "dvxc in dv_of_drho"
!   write (*,*) dvxc(1:6)
   do i=1,mesh
      dvhx(i) = dvhx(i) + dvxc(i) * drho(i)
   end do
!   write (*,*) "Hartree + dvxc in dv_of_drho"
!   write (*,*) dvhx(1:3)
!   write (*,*) dvhx(mesh-2:mesh)
!
   deallocate (qq)

   return
end subroutine
!----------------------------------------------------------------------
subroutine dvxc_dn(mesh, rho, dvxc)
   !-------------------------------------------------------------------
   ! compute the derivative of xc-potential w.r.t local density.
   ! some routine in PH and flibs will be called
   !
   use funct,  only : dft_is_gradient, dmxc
   !
   implicit none
   !
   ! I/O variables
   !
   integer :: mesh
   real(kind=8) :: rho(mesh), dvxc(mesh)
   !
   ! local variables
   !
   integer :: i
   !
   !
   !
   if ( dft_is_gradient() ) &
      call errore ('dvxc_dn', 'gradient correction to dvxc not yet implemented', 1)
   do i = 1, mesh
      ! LDA only
      dvxc(i) = dmxc (rho(i))
      !
   end do
   !
   return
   !
end subroutine dvxc_dn
!--------------------------------------------------------------------
subroutine drho_of_r(mesh, dx, r, r2, y, dy, drho)
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
   real (kind=8) :: dx, r(mesh), r2(mesh), y(mesh), drho(mesh)
   complex (kind=8) :: dy(mesh)
   ! local variables
   integer i
   do i=1,mesh
      drho(i) = 2.d0 * y(i) * real(dy(i)) * r(i) / (fpi*r2(i))
   end do

   return
end subroutine
!--------------------------------------------------------------------
subroutine init_dpot(r,mesh,dvpot)
   !--------------------------------------------------------------------
   !
   ! initialize external potential
   !
   use constants, only : e2
   implicit none
   !
   ! I/O variables
   !
   integer mesh
   real (kind=8) :: r(mesh), dvpot(mesh)
   !
   ! local variables
   !
   integer i

   do i =1,mesh
      dvpot (i) = - e2*r(i)
   end do

   return
end subroutine

!--------------------------------------------------------------------
subroutine sternheimer(u, l, ll, mesh, dx, r, sqr, r2, vpot, zed, y, dvpot, dy)
   !--------------------------------------------------------------------
   !
   ! solve the sternheimer equation for imaginary frequency 
   ! in radial coordinates by Numerov method
   !
   implicit none
   !
   ! I/O variables
   ! 
   integer mesh, l, ll
   real (kind=8) :: u,  dx, zed
   real (kind=8) :: r(mesh), sqr(mesh), r2(mesh)
   real (kind=8) :: vpot(mesh), y(mesh), dvpot(mesh)
   complex (kind=8) :: dy(mesh)
   !
   ! local variables
   !
   integer i, icl, j
   real (kind=8) :: ddx12, sqlhf, x2l2
   complex (kind=8) :: gg, aa, bb, fac, e
   complex (kind=8), allocatable :: f(:), g(:), yy(:)

   allocate ( f(mesh), g(mesh), yy(mesh) )

   ddx12=dx*dx/12.d0
   sqlhf = (l+0.5d0)**2
   x2l2 = 2*l+2
   e = cmplx (0.d0, u, kind=8)
   !
   ! set up the f-function and determine the position of its last
   ! change of sign
   ! f < 0 (approximatively) means classically allowed   region
   ! f > 0         "           "        "      forbidden   "
   !
   icl = 2
!   f(0) = ddx12 *( sqlhf + r2(0) * (vpot(0)-e) )
   f(1) = ddx12 *( r2(1) * (vpot(1)-e) )
   do i=2,mesh
!      f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
      f(i) = ddx12 * ( r2(i) *(vpot(i)-e) )
      if( real(f(i)) .ne. sign(real(f(i)),real(f(i-1))) &
          .and. real(f(i)).gt.0d0 ) icl=i
   end do
   ! write (*,*) icl

   do i=1,mesh
      f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
   end do

   do i=1,mesh
      f(i)=1.0d0-f(i)
      g(i)= ddx12 * r2(i) * dvpot(i) * y(i) 
   end do 
   ! step 1) determine a solution to the homogeneous equation that
   !         satisfies proper boundary conditions at 0 and infty
   !         NB it will NOT satisfy the homogeneous equation at icl 
   !
   ! determination of the wave-function in the first two points 
   !
   yy(1) = r(1)**(l+1) *(1.d0 - 2.d0*zed*r(1)/x2l2) / sqr(1)
   yy(2) = r(2)**(l+1) *(1.d0 - 2.d0*zed*r(2)/x2l2) / sqr(2)
   !
   ! outward integration 
   !
   do i =2, icl-1
      yy(i+1)=( (12.d0-10.d0*f(i))*yy(i)-f(i-1)*yy(i-1) )/f(i+1)
   end do
   !
   ! rescale to 1 at icl
   !
   fac = 1.d0/yy(icl)
   do i =1, icl
      yy(i)= yy(i) * fac 
   end do
   !
   ! determination of the wave-function in the last two points 
   ! assuming y(mesh+1) = 0 and y(mesh) = dx
   !
   yy(mesh) = dx
   yy(mesh-1) = (12.d0-10.d0*f(mesh))*yy(mesh)/f(mesh-1) 
   !
   ! inward integration 
   !
   do i = mesh-1,icl+1,-1
      yy(i-1)=( (12.d0-10.d0*f(i))*yy(i)-f(i+1)*yy(i+1) )/f(i-1)
      if (abs(yy(i-1)).gt.1.d6) then
         fac = 1.d0/yy(i-1)
         do j=i-1,mesh
            yy(j) = yy(j) * fac
         end do
      end if
   end do
   !
   ! rescale to 1 at icl
   !
   fac = 1.d0 / yy(icl)
   do i = icl, mesh
      yy(i)= yy(i) * fac
   end do
   !
   ! step 2) determine a solution to the inhomogeneous equation that
   !         satisfies proper boundary conditions at 0 and infty
   !         and matches in value at icl
   !         NB it will NOT satisfy the inhomogeneous equation at icl 
   !
   ! determination of the wave-function in the first two points 
   !
   dy(1) = r(1)**(l+1) *(1.d0 - 2.d0*zed*r(1)/x2l2) / sqr(1)
   dy(2) = r(2)**(l+1) *(1.d0 - 2.d0*zed*r(2)/x2l2) / sqr(2)
   !
   ! outward integration 
   !
   do i =2, icl-1
      gg = g(i+1) + 10.d0*g(i) + g(i-1)
      dy(i+1)=( (12.d0-10.d0*f(i))*dy(i)-f(i-1)*dy(i-1) + gg )/f(i+1)
   end do
   !
   ! choose the solution that goes to 1 in icl
   !
   fac =  dy(icl) - 1.d0
   do i = 1,icl
      dy(i) = dy(i) - fac * yy(i)
   end do
   !
   ! determination of the wave-function in the last two points 
   ! assuming y(mesh+1) = 0 and y(mesh) = dx
   !
   dy(mesh) = dx
   dy(mesh-1) = (12.d0-10.d0*f(mesh))*yy(mesh)/f(mesh-1) 
   !
   ! inward integration 
   !
   do i = mesh-1,icl+1,-1
      gg = g(i+1) + 10.d0*g(i) + g(i-1)
      dy(i-1)=( (12.d0-10.d0*f(i))*dy(i)-f(i+1)*dy(i+1) + gg )/f(i-1)
      if (abs(dy(i-1)).gt.1.d6) then
         if (abs(yy(i-1)).gt. 1.d-12) then
            fac = ( dy(i-1) - 1.d0 ) / yy(i-1)
            do j = i-1, mesh
               dy(j) = dy(j) - fac * yy(j)
            end do
         else
            do j = i-1, mesh
               dy(j) = 0.d0
            end do
         end if
      end if
   end do
   !
   ! choose the solution that goes to 1 in icl
   !
   fac =  dy(icl) - 1.d0
   do i = icl, mesh
      dy(i) = dy(i) - fac * yy(i)
   end do
   !
   ! step 3) chose the proper combination of dy and yy so that the
   !         inhomogeneous equation is satisfied also in icl
   i = icl
   gg = g(i+1) + 10.d0*g(i) + g(i-1)
   aa= (12.d0-10.d0*f(i))*dy(i)-f(i+1)*dy(i+1)-f(i-1)*dy(i-1) + gg 
   bb= (12.d0-10.d0*f(i))*yy(i)-f(i+1)*yy(i+1)-f(i-1)*yy(i-1) 
   !   write (*,*) aa, bb
   fac = aa/bb
   do i=1,mesh
      dy(i) = dy(i) - fac * yy(i)
   end do
   deallocate ( f, g, yy )

   return

end subroutine

