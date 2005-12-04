!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine c6_tfvw (mesh, zed, dx, r, r2, rho)
   !--------------------------------------------------------------------
   !
   use kinds,      only : DP
   use constants,  only : e2, pi, fpi
   use ld1inc,     only : lsd, nwf, oc, ll
   !
   implicit none
   !
   ! I/O variables
   !
   integer mesh
   real (kind=8) :: dx, zed, r(0:mesh), r2(0:mesh), rho(0:mesh)
   !
   ! local variables
   !
   logical :: csi
   real (kind=8) :: error, error2, e, charge, beta, u, alpha, dalpha, c6, du1, &
                    du2, factor, ze2, thresh
   real (kind=8), allocatable :: veff(:), y(:), yy(:), sqr(:)
   real (kind=8), allocatable :: dvpot(:), dvscf(:), drho(:), dvhx(:), dvxc(:), pp(:)
   complex (kind=8), allocatable :: dy(:), drho_old(:)
   integer i, iter, n, l, iu, Nu, Nc, counter, nstop

   allocate ( veff(0:mesh),y(0:mesh),yy(0:mesh),sqr(0:mesh) )
   allocate ( dvpot(0:mesh),dvscf(0:mesh),drho(0:mesh),dvhx(0:mesh),dvxc(0:mesh),pp(0:mesh) )
   allocate ( dy(0:mesh), drho_old(0:mesh) )
   !
   write(6,'(/,/,/,5x,20(''-''),'' Compute C6 from polarizability with TFvW approx.'',10(''-''),/)')
   !
   do i = 0, mesh
      rho(i) = rho(i) / (fpi*r(i)**2)
   end do
   !
   counter = 0
   do i = 0, mesh
      if (rho(i) .gt. 1.0d-30) counter = counter + 1 
   enddo
   mesh = counter
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
!   do i=0,mesh
!      read (7,'(P5E20.12)') r(i), rho(i), y(i), y(i), y(i)
!      write (6,'(P5E15.6)') r(i), rho(i)
!   end do
!   close (7)
!
! compute unperturbed effective potential
!
   call veff_of_rho(mesh,dx,r,r2,rho,y,veff)
!
! check that veff and y are what we think
!
   n = 1
   l = 0
   e = -1.d-7
   charge=0.d0
   ze2 = - zed * e2
   thresh = 1.d-14

   do i=0,mesh
      sqr(i) = sqrt(r(i))
      charge = charge + rho(i) * fpi * r2(i) * r(i) * dx
   end do
!   call solve_scheq(n,l,e,mesh,dx,r,sqr,r2,veff,zed,yy)
   call ascheq (n, l, e, mesh, dx, r, r2, sqr, veff(0), ze2, thresh, yy(0), nstop)

   error = 0.d0
   do i=0,mesh
      error = error + (y(i)-yy(i)/sqr(i)*sqrt(charge))**2 * r2(i) * dx
   end do

   if (error > 1.d-8) & 
      call errore('c6_tfvw','auxiliary funtions veff(r) and y(r) are inaccurate',1)
!
! initialize external perturbation (electric field)
!
   call init_dpot(r,mesh,dvpot)
!
! derivative of xc-potential
!
   call dvxc_dn(mesh, rho, dvxc)
!
!write(*,'(PE20.12)')sum(abs(dvxc))
!stop
   write(6,'(5x,''Frequency dependent polarizability is written into freq-pol.dat'',/)')

   c6    = 0.0d0
   alpha = 0.0d0

   open(1, file = 'freq-pol.dat')
   write (1,'(15x,"    u          alpha(angstrong)       alpha(a.u.)  ",/)')
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
         do i=0,mesh
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
         call sternheimer(u, l, mesh, dx, r, sqr, r2, veff ,zed, y, dvscf, dy)
         ! compute drho of r
         !
         call drho_of_r(mesh, dx, r, r2, y, dy, drho)
         !
         ! compute dv of drho
         !
         call dv_of_drho(mesh, dx, r,r2,rho,drho,dvhx,dvxc,pp)
         !
         ! mix
         !
         error = 0.d0
         error2 = 0.d0
         do i=0,mesh
            dvscf(i) = dvscf(i) + beta * (dvpot(i)+dvhx(i) -dvscf(i))
            error = error + abs (drho(i) -drho_old(i))
            error2 = error2 + abs (drho(i) -drho_old(i))* r(i) * dx
            drho_old(i) = drho(i)
         end do 
         dalpha = abs(alpha + pp(mesh)) 
         alpha = -pp(mesh)
!         write (*,'(4e16.6)') alpha, dalpha, error, error2
      end do

      write (1,'(17x, f8.4, 3x, PE14.6, 9x, PE14.6)') u, pp(mesh)*0.529177**3, pp(mesh)
      if (iu .eq. 0) & 
      write (6,'(5x, "Static polarizability: ", f10.5, " (in angstrom^3)   --->", f10.5,&
                & "  (in e^2a0^3)")') pp(mesh)*0.529177**3, pp(mesh)

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
   deallocate ( veff, y, yy, sqr )
   deallocate ( dvpot, dvscf, drho, dvhx, pp )
   
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
   real (kind=8) :: dx, r(0:mesh), r2(0:mesh), rho(0:mesh), &
                    veff(0:mesh), y(0:mesh)
   ! 
   ! local variables
   !
   real (kind=8), allocatable :: vold(:)
	real (kind=8) :: dx2, error
   integer i, iter, k

!
! compute auxiliary wavefunction y
!
   do i=0,mesh
      y(i) = sqrt(rho(i)*r(i)*fpi)
   end do
!
! compute effective potential veff
!
   allocate (vold(0:mesh))

   do i=0,mesh
      vold(i) = 0.d0
   end do
   dx2= dx*dx
   error = 1.d0
   k = 0
   do while (error > 1.d-9) 
      !
      k=k+1
      !
      do i=1,mesh-1
         veff(i) = ( y(i+1)/y(i) + y(i-1)/y(i) -2.d0 )/dx2  &
                 - ( vold(i+1)*y(i+1)/y(i) + vold(i-1)*y(i-1)/y(i) -2.d0*vold(i) )/12.d0
      end do
      veff(0) = veff(1) + (veff(2)-veff(1))*(r(0)-r(1))/(r(2)-r(1)) 
      veff(mesh) = (y(mesh-1)/y(mesh) -2.d0 )/dx2 &
                 - (vold(mesh-1)*y(mesh-1)/y(mesh) -2.d0*vold(mesh) )/12.d0
!
      error = 0.d0
      do i=0,mesh
         error = error + abs( veff(i) - vold(i) )
         vold(i) = veff(i)
      end do
      error = error / mesh
!      write (*,*) 'iteration # ', k, error
   end do

   deallocate (vold)
   !
   do i=0,mesh
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
subroutine dv_of_drho(mesh,dx,r,r2,rho,drho,dvhx,dvxc,pp)
   !--------------------------------------------------------------------
   use constants, only : e2, pi, fpi
!   use flags,      only : HartreeFock, rpa
   implicit none
   !
   ! I/O variables
   !
   integer mesh
   real (kind=8) :: dx, r(0:mesh), r2(0:mesh)
   real (kind=8) :: rho(0:mesh), drho(0:mesh), dvhx(0:mesh), pp(0:mesh), dvxc(0:mesh)
   !
   ! local variables
   !
   real (kind=8) :: dr3, kf2, charge
   real (kind=8), allocatable :: qq(:)
   integer i

   allocate (qq(0:mesh))

   do i=0,mesh 
      dr3   = fpi * r2(i) * r(i) * dx
      pp(i) = drho(i) * r(i)  * dr3 /3.d0
      qq(i) = drho(i) / r2(i) * dr3 /3.d0
   end do
   do i=1,mesh
      pp(i) = pp(i) + pp(i-1) 
   end do
   do i=mesh-1,0,-1
      qq(i) = qq(i) + qq(i+1)
   end do
   do i=0,mesh
      dvhx(i) = e2 * ( pp(i) / r2(i) + qq(i) * r(i) ) ! Hartree term
   end do
! add TF term
   do i=0,mesh
      kf2 = ( 3.d0*pi*pi*rho(i) )**(2.d0/3.d0)
      dvhx(i) = dvhx(i) + e2/3.d0* kf2 / rho(i) * drho(i)
   end do
!
!add xc term
!         
   do i=0,mesh
      dvhx(i) = dvhx(i) + dvxc(i) * drho(i)
   end do
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
   real(kind=8) :: rho(0:mesh), dvxc(0:mesh)
   !
   ! local variables
   !
   integer :: i
   !
   !
   !
   if ( dft_is_gradient() ) &
      call errore ('dvxc_dn', 'gradient correction to dvxc not yet implemented', 1)
   do i = 0, mesh
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
   real (kind=8) :: dx, r(0:mesh), r2(0:mesh), y(0:mesh), drho(0:mesh)
   complex (kind=8) :: dy(0:mesh)
   ! local variables
   integer i
   do i=0,mesh
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
   real (kind=8) :: r(0:mesh), dvpot(0:mesh)
   !
   ! local variables
   !
   integer i

   do i =0,mesh
      dvpot (i) = - e2*r(i)
   end do

   return
end subroutine

!--------------------------------------------------------------------
subroutine sternheimer(u, l, mesh, dx, r, sqr, r2, vpot, zed, y, dvpot, dy)
   !--------------------------------------------------------------------
   !
   ! solve the sternheimer equation for imaginary frequency 
   ! in radial coordinates by Numerov method
   !
   implicit none
   !
   ! I/O variables
   ! 
   integer mesh, l
   real (kind=8) :: u,  dx, zed
   real (kind=8) :: r(0:mesh), sqr(0:mesh), r2(0:mesh)
   real (kind=8) :: vpot(0:mesh), y(0:mesh), dvpot(0:mesh)
   complex (kind=8) :: dy(0:mesh)
   !
   ! local variables
   !
   integer i, icl, j
   real (kind=8) :: ddx12, sqlhf, x2l2
   complex (kind=8) :: gg, aa, bb, fac, e
   complex (kind=8), allocatable :: f(:), g(:), yy(:)

   allocate ( f(0:mesh), g(0:mesh), yy(0:mesh) )

   ddx12=dx*dx/12.d0
   sqlhf = (l+0.5d0)**2
   x2l2 = 2*l+2
   e = cmplx (0.d0, u)
   !
   ! set up the f-function and determine the position of its last
   ! change of sign
   ! f < 0 (approximatively) means classically allowed   region
   ! f > 0         "           "        "      forbidden   "
   !
   icl = 2
!   f(0) = ddx12 *( sqlhf + r2(0) * (vpot(0)-e) )
   f(0) = ddx12 *( r2(0) * (vpot(0)-e) )
   do i=1,mesh
!      f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
      f(i) = ddx12 * ( r2(i) *(vpot(i)-e) )
      if( real(f(i)) .ne. sign(real(f(i)),real(f(i-1))) &
          .and. real(f(i)).gt.0d0 ) icl=i
   end do
!   write (*,*) icl

   do i=0,mesh
      f(i) = ddx12 * ( sqlhf + r2(i) *(vpot(i)-e) )
   end do

   do i=0,mesh
      f(i)=1.0d0-f(i)
      g(i)= ddx12 * r2(i) * dvpot(i) * y(i) 
   end do 
   ! step 1) determine a solution to the homogeneous equation that
   !         satisfies proper boundary conditions at 0 and infty
   !         NB it will NOT satisfy the homogeneous equation at icl 
   !
   ! determination of the wave-function in the first two points 
   !
   yy(0) = r(0)**(l+1) *(1.d0 - 2.d0*zed*r(0)/x2l2) / sqr(0)
   yy(1) = r(1)**(l+1) *(1.d0 - 2.d0*zed*r(1)/x2l2) / sqr(1)
   !
   ! outward integration 
   !
   do i =1, icl-1
      yy(i+1)=( (12.d0-10.d0*f(i))*yy(i)-f(i-1)*yy(i-1) )/f(i+1)
   end do
   !
   ! rescale to 1 at icl
   !
   fac = 1.d0/yy(icl)
   do i =0, icl
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
   dy(0) = r(0)**(l+1) *(1.d0 - 2.d0*zed*r(0)/x2l2) / sqr(0)
   dy(1) = r(1)**(l+1) *(1.d0 - 2.d0*zed*r(1)/x2l2) / sqr(1)
   !
   ! outward integration 
   !
   do i =1, icl-1
      gg = g(i+1) + 10.d0*g(i) + g(i-1)
      dy(i+1)=( (12.d0-10.d0*f(i))*dy(i)-f(i-1)*dy(i-1) + gg )/f(i+1)
   end do
   !
   ! choose the solution that goes to 1 in icl
   !
   fac =  dy(icl) - 1.d0
   do i = 0,icl
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
         fac = ( dy(i-1) - 1.d0 ) / yy(i-1)
         do j = i-1, mesh
            dy(j) = dy(j) - fac * yy(j)
         end do
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
   do i=0,mesh
      dy(i) = dy(i) - fac * yy(i)
   end do
   deallocate ( f, g, yy )

   return

end subroutine

