!
!--------------------------------------------------------------------------
      subroutine compute_phius(lam,ik,ns,xc,iflag)
!--------------------------------------------------------------------------
!
!     This routine computes the phi functions by pseusizing the
!     all_electron chi functions. In input it receives, the point
!     ik where the cut is done, the angualar momentum lam,
!     and the correspondence with the all eletron wavefunction
!
!
!
use ld1inc
  implicit none

      real(kind=dp) :: &
               fae,    & ! the value of the all-electron function
               f1ae,   & ! its first derivative
               xc(8),  & ! the coefficients of the fit
               f2ae      ! the second derivative

      integer :: &
               ik, &     ! the point corresponding to rc
               ns, &     ! the function to pseudize
               iflag,&   ! if 1 print
               iok,  &   ! if 0 there are no problem
               lam      ! the angualar momentum


      real(kind=dp) :: &
               f1aep1,f1aem1,jnor, &  ! auxilairy quantities
               bm(2),  &              ! the derivative of the bessel
               ff,     &              ! contain deltah corrections
               fact(2), &             ! factor of normalization
               j1(ndm,8)             ! the bessel functions
     
      real(kind=dp) :: &
            deriv_7pts, deriv2_7pts


      integer :: &
               n, &        ! counter on mesh points
               nc         ! counter on bessel

!
!    compute first and second derivative
!
      ff=1-dx**2/48.d0
!      fae=(psipsus(ik+1,ns)+psipsus(ik,ns))*0.5d0
!      f1aep1=psipsus(ik+1,ns)*ff/sqr(ik+1)
!      f1ae=psipsus(ik,ns)*(-12.d0+10.d0*ff)/sqr(ik)
!      f1aem1=psipsus(ik-1,ns)*ff/sqr(ik-1)
!      f2ae=(f1aep1+f1ae+f1aem1)/dx**2
!      f1ae=(psipsus(ik+1,ns)-psipsus(ik,ns))/(r(ik+1)-r(ik))

      fae=psipsus(ik,ns)
      f1ae=deriv_7pts(psipsus(1,ns),ik,r(ik),dx)
      f2ae=deriv2_7pts(psipsus(1,ns),ik,r(ik),dx)

!
!    find the q_i of the bessel functions
!      
      call find_qi(f1ae/fae,xc(4),ik,lam,2,1,iok)
      if (iok.ne.0) &
             call errore('compute_phius','problems with find_qi',1)
!
!    compute the functions
!
      do nc=1,2
         call sph_besr(ik+5,r,xc(3+nc),lam,j1(1,nc))
         fact(nc)=psipsus(ik,ns)/j1(ik,nc)
         do n=1,ik+5
            j1(n,nc)=j1(n,nc)*fact(nc)
         enddo
      enddo
!
!    compute the second derivative and impose continuity of zero, 
!    first and second derivative
!
       
      do nc=1,2
!          f1aep1=j1(ik+1,nc)*ff/sqr(ik+1)
!          f1ae=j1(ik,nc)*(-12.d0+10.d0*ff)/sqr(ik)
!          f1aem1=j1(ik-1,nc)*ff/sqr(ik-1)
!          bm(nc)=(f1aep1+f1ae+f1aem1)/dx**2
         bm(nc)=deriv2_7pts(j1(1,nc),ik,r(ik),dx)
      enddo

      xc(2)=(f2ae-bm(1))/(bm(2)-bm(1))
      xc(1)=1.d0-xc(2)
      if (iflag.eq.1) then
         write(6,110) els(ns),rcutus(ns),2.d0*xc(5)**2
110      format (5x, ' Wfc-us ',a3,' rcutus=',f6.3, &
                '  Estimated cut-off energy= ', f8.2,' Ry')
      endif
!
!    define the phis function
!
      do n=1,ik
         phis(n,ns)=xc(1)*j1(n,1)+xc(2)*j1(n,2)
      enddo

      do n=ik+1,mesh
         phis(n,ns)=psipsus(n,ns)
      enddo

      do nc=1,2
         xc(nc)=xc(nc)*fact(nc)
      enddo
      xc(3)=0.d0
      xc(6)=0.d0
      xc(7)=0.d0
      xc(8)=0.d0

      return
      end
