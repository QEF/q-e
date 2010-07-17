!-----------------------------------------------------------------------
      subroutine kin_e_density &
         (ndm, mesh, nwf, ll, oc, chi, r, r2, dx, tau )
!-----------------------------------------------------------------------
!
!     calculate the (spherical) kinetic energy density 
!     tau = sum_nwf 0.5*|\nabla phi|^2     in Cartesian coordinate
!     = 0.5*sum_nwf occ(nwf)/4pi*{[d(chi/r)/dr]^2 + (chi/r)^2*ll*(ll+1)/r^2}
!     = 0.5*sum_nwf occ(nwf)/4pi*[(1/r*dchi-chi/r^2)^2+(chi/r)^2*ll(ll+1)/r^2)]
!       in spherical coordinate
!     temporary - must be rewritten in a cleaner way
!
      implicit none
! input
      integer ndm, mesh, nwf,ll(nwf)
      real*8 oc(nwf), chi(ndm,2,nwf),r2(ndm),r(ndm),dx
! output
      real*8 tau(ndm,2)
! local
      logical tspin
      integer i,n
      real*8  z,full,half,rho_tot, fourpi,                              &
     &     ocup,ocdn,temp,nll, corr
      parameter ( fourpi = 4.d0*3.141592653589793d+00 )
      real*8 dchi(ndm)
      tspin = .false.
!
      tau (:,:)=0.0d0
!
      do n=1,nwf
         half=2*ll(n)+1
         full=half+half
         nll=ll(n)*(ll(n)+1)
         call deriv5pt(mesh,dx,r,chi(1,1,n),dchi)
         if( oc(n).le.half) then
            ocup=oc(n)
            ocdn=0.d0
         else
            tspin = .true.
            ocup=half
            ocdn=oc(n)-half
         end if
!
         do i=1,mesh
!     kinetic energy density
            temp=(dchi(i)-chi(i,1,n)/r(i))**2.0d0                         &
     &           + chi(i,1,n)**2.0d0/r2(i)*nll
            tau(i,1)=tau(i,1)+ocup*temp
            tau(i,2)=tau(i,2)+ocdn*temp
         end do
      end do
!
      do i=1,mesh
         tau(i,:)=tau(i,:)/fourpi*0.5d0/r2(i)
      end do
      return
      end
!
      subroutine deriv5pt(mesh,dx,r,v,dv)
!
! numerical derivative using 5-point formula - error: O(dx^5)
! Assumes a logarithmic grid  r_i = r_0 exp((i-1)dx)
!
      implicit none
! input
      integer mesh
      real*8 v(mesh), r(mesh), dx
! output:  dv = dv/dr = (1/r) dv/dx
      real*8 dv(mesh)
!
      integer i
!
!
      dv(1)=(-25.d0*v(1)+48.d0*v(2)-36.d0*v(3)+16.d0*v(4)               &
     &       -3.d0*v(5))/(12.d0*dx*r(1))
      dv(2)=(-3.d0*v(1)-10.d0*v(2)+18.d0*v(3)- 6.d0*v(4)                &
     &       +     v(5))/(12.d0*dx*r(2))
!
      do i=3,mesh-2
         dv(i)=(v(i-2)-8.d0*v(i-1)+8.d0*v(i+1)-v(i+2))/(12.d0*dx*r(i))
      end do
!
      dv(mesh-1)=( 3.d0*v(mesh)+10.d0*v(mesh-1)-18.d0*v(mesh-2)+        &
     &             6.d0*v(mesh-3)-     v(mesh-4))/(12.d0*dx*r(mesh-1))
      dv(mesh)=( 25.d0*v(mesh)  -48.d0*v(mesh-1)+36.d0*v(mesh-2)-       &
     &           16.d0*v(mesh-3)+3.d0*v(mesh-4))/(12.d0*dx*r(mesh))
!
      return
      end
