!
      subroutine lschps(mode,z,amesh,al,mmax,nin,mch, &
                        n,l,e,u,r,v)
!
! integrates radial pauli-type scalar-relativistic equation
! on a logarithmic mesh
! modified routine to be used in finding norm-conserving
! pseudopotential
!
! mode = 1 is for full potential bound state
! mode = 2 is for pseudopotential bound state
! mode = 3 is for full potential to find log derivative
! mode = 4 is for pseudopotential to find energy which produces
!            specified log derivative (fake bound state)
! mode = 5 is for pseudopotential to produce wavefunction beyond
!            radius used for pseudopotential construction
!
      implicit none
      integer,parameter :: dp = kind(1.d0)
      real(kind=dp),parameter:: e2=2.d0
      real(kind=dp):: aei, aeo, aii, aio, al, als, amesh,  cn
      real(kind=dp):: dabs, de, dmax1, dmin1, e, emax, emin
      real(kind=dp):: eps, fss, gamma, ro, sc
      real(kind=dp):: sls, sn, tfapot, uld, uout,  upin, upout
      real(kind=dp):: xkap, z
      integer:: i, it, l, mch, mmax, mode, n, nin, nint, node, ndm, ierr

      real(kind=dp):: r(mmax),v(mmax),u(mmax)
! these arrays are used as work space
      real(kind=dp),allocatable :: up(:),upp(:),cf(:),dv(:),fr(:),frp(:)

      allocate(up(mmax), stat=ierr)
      allocate(upp(mmax), stat=ierr)
      allocate(cf(mmax), stat=ierr)
      allocate(dv(mmax), stat=ierr)
      allocate(fr(mmax), stat=ierr)
      allocate(frp(mmax), stat=ierr)

      uld=0.d0
      v=v/e2
      e=e/e2
!
! convergence factor for solution of schroedinger eq.  if calculated
! correction to eigenvalue is smaller in magnitude than eps times
! the magnitude of the current guess, the current guess is not changed.
      eps=1.0d-8
!
! relativistic - non-relativistic switch
!
      if(mode .eq. 1 .or. mode .eq. 3) then
         fss=(1.0d0/137.036d0)**2
         if(l .eq. 0) gamma=dsqrt(1.0d0-fss*z**2)
         if(l .gt. 0) gamma=(l*dsqrt(l**2-fss*z**2) + &
              (l+1)*dsqrt((l+1)**2-fss*z**2))/(2*l+1)
      else
         fss=1.0d-20
         gamma=l+1
      end if
!
      sls=l*(l+1)
!
      if(mode .eq. 1 .or. mode .eq. 2) then
         emax=v(mmax)+0.5d0*sls/r(mmax)**2
         emin=0.0d0
         do i=1,mmax
            emin=min(emin,v(i)+0.5d0*sls/r(i)**2)
!           if (l.eq.0)  write(6,*) r(i),v(i)  
         end do
!         if (l.eq.0) stop  
         if(e .gt. emax) e=1.25d0*emax
         if(e .lt. emin) e=0.75d0*emin
         if(e .gt. emax) e=0.5d0*(emax+emin)
      else if(mode .eq. 4) then
         emax=e + 10.0d0
         emin=e - 10.0d0
      end if
!
      do i=1,4
         u(i)=0.0d0
         up(i)=0.0d0
         upp(i)=0.0d0
      end do
      nint=0
      als=al**2
!
! return point for bound state convergence
 10   nint=nint+1
      if(nint .gt. 60) then
         print '('' warning: wfc '',2i2,'' not converged'')', n, l
         u=0.d0
         go to 999
      end if
!
! coefficient array for u in differential eq.
      do i=1,mmax
         cf(i)=als*sls + 2.0d0*als*(v(i)-e)*r(i)**2
      end do
!
! calculate dv/dr for darwin correction
      dv(1)=(-50.d0*v(1)+96.d0*v(2)-72.d0*v(3)+32.d0*v(4) &
             -6.d0*v(5))/(24.d0*al*r(1))
      dv(2)=(-6.d0*v(1)-20.d0*v(2)+36.d0*v(3)-12.d0*v(4) &
             +2.d0*v(5))/(24.d0*al*r(2))
!
      do i=3,mmax
         dv(i)=(2.d0*v(i-2)-16.d0*v(i-1)+16.d0*v(i+1) &
               -2.d0*v(i+2))/(24.d0*al*r(i))
      end do
!
!  relativistic coefficient arrays for u (fr) and up (frp).
      do i=1,mmax
        fr(i)=als*(r(i)**2)*(-fss*(v(i)-e)**2 + 0.5d0*fss*dv(i)/ &
        (r(i)*(1.0d0+0.5d0*fss*(e-v(i)))))
        frp(i)=-al*r(i)*0.5d0*fss*dv(i)/(1.0d0+0.5d0*fss*(e-v(i)))
      end do
!
! find classical turning point for matching
      if(mode .eq. 1 .or. mode .eq. 2) then
        do i=mmax,2,-1
           if(cf(i-1) .le. 0.d0 .and. cf(i) .gt. 0.d0) then
              mch=i
              go to 40
           end if
        end do
        print '('' warning: wfc '',2i2,'' no turning point'')', n, l
        e=0.0
        do i=1,mmax
           u (i)=0.0
        end do
        go to 999
      else
         nin=mch
      end if
 40   continue
!
! start wavefunction with series
!
      do i=1,4
         u(i)=r(i)**gamma
         up(i)=al*gamma*r(i)**gamma
         upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
      end do
!
! outward integration using predictor once, corrector
! twice
      node=0
!               
      do i=4,mch-1
         u(i+1)=u(i)+aeo(up,i)
         up(i+1)=up(i)+aeo(upp,i)
         do it=1,2
            upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
            up(i+1)=up(i)+aio(upp,i)
            u(i+1)=u(i)+aio(up,i)
         end do
         if(u(i+1)*u(i) .le. 0.0d0) node=node+1
      end do
!
      uout=u(mch)
      upout=up(mch)
!
!
      if(node-n+l+1 .eq. 0 .or. mode .eq. 3 .or. mode .eq. 5) then
!
         if(mode .eq. 1 .or. mode .eq. 2) then
!
! start inward integration at 10*classical turning
! point with simple exponential
            nin=mch+2.3d0/al
            if(nin+4 .gt. mmax) nin=mmax-4
            xkap=dsqrt(sls/r(nin)**2 + 2.0d0*(v(nin)-e))
!
            do i=nin,nin+4
               u(i)=exp(-xkap*(r(i)-r(nin)))
               up(i)=-r(i)*al*xkap*u(i)
               upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
            end do
!
! integrate inward
!
          do i=nin,mch+1,-1
             u(i-1)=u(i)+aei(up,i)
             up(i-1)=up(i)+aei(upp,i)
             do it=1,2
                upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
                up(i-1)=up(i)+aii(upp,i)
                u(i-1)=u(i)+aii(up,i)
             end do
          end do
!
! scale outside wf for continuity
          sc=uout/u(mch)
!
          do i=mch,nin
             up(i)=sc*up(i)
             u (i)=sc*u (i)
          end do
!
          upin=up(mch)
!
       else
!
          upin=uld*uout
!
       end if
!
! perform normalization sum
!
       ro=r(1)/dsqrt(amesh)
       sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)
!
       do i=1,nin-3
          sn=sn+al*r(i)*u(i)**2
       end do
!
       sn=sn + al*(23.0d0*r(nin-2)*u(nin-2)**2 &
                 + 28.0d0*r(nin-1)*u(nin-1)**2 &
                 +  9.0d0*r(nin  )*u(nin  )**2)/24.0d0
!
! normalize u
       cn=1.0d0/dsqrt(sn)
       uout=cn*uout
       upout=cn*upout
       upin=cn*upin
!
       do i=1,nin
          up(i)=cn*up(i)
          u(i)=cn*u(i)
       end do
       do i=nin+1,mmax
          u(i)=0.0d0
       end do
!
! exit for fixed-energy calculation
!
       if(mode .eq. 3 .or. mode .eq. 5) go to 999

! perturbation theory for energy shift
       de=0.5d0*uout*(upout-upin)/(al*r(mch))
!
! convergence test and possible exit
!
       if ( abs(de) .lt. max(abs(e),0.2d0)*eps) go to 999
!
       if(de .gt. 0.0d0) then 
          emin=e
       else
          emax=e
       end if
       e=e+de
       if(e .gt. emax .or. e .lt. emin) e=0.5d0*(emax+emin)
!
! loop back to converge e
!
       go to 10
!
      else if(node-n+l+1 .lt. 0) then
! too few nodes
         emin=e
         e=0.5d0*(emin+emax)
         go to 10
!
      else
! too many nodes
         emax=e
         e=0.5d0*(emin+emax)
         go to 10
      end if      
!
! deallocate arrays and exit
!
 999  continue
      deallocate(frp)
      deallocate(fr)
      deallocate(dv)
      deallocate(cf)
      deallocate(upp)
      deallocate(up)
      e=e*e2
      v=v*e2
      return
      
      end
!
      function aei(y,j)
!
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp):: y, aei
      integer j
!
      dimension y(600)
      aei=-(4.16666666667d-2)*(55.0d0*y(j)-59.0d0*y(j+1) &
       +37.0d0*y(j+2)-9.0d0*y(j+3))
      return
      end
!
! adams extrapolation and interpolation formulas for
! outward and inward integration, abramowitz and
! stegun, p. 896
      function aeo(y,j)
!
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp):: y, aeo
      integer:: j   
!
      dimension y(600)
      aeo=(4.16666666667d-2)*(55.0d0*y(j)-59.0d0*y(j-1) &
       +37.0d0*y(j-2)-9.0d0*y(j-3))
      return
      end
!
      function aii(y,j)
!
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp) :: y, aii
      integer:: j
!
      dimension y(600)
      aii=-(4.16666666667d-2)*(9.0d0*y(j-1)+19.0d0*y(j) &
       -5.0d0*y(j+1)+y(j+2))
      return
      end
!
      function aio(y,j)
!
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp):: y,aio
      integer :: j
!
      dimension y(600)
      aio=(4.16666666667d-2)*(9.0d0*y(j+1)+19.0d0*y(j) &
       -5.0d0*y(j-1)+y(j-2))
      return
      end
