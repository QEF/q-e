!
!---------------------------------------------------------------
      subroutine starting_potential &
         (ndm,mesh,zval,zed,nwf,oc,nn,ll,r,enl, &
              vxt,vpot,enne,nspin)
!---------------------------------------------------------------
!
! starting potential: generalized thomas-fermi atomic potential
! it is assumed that the effective charge seen by the reference
! electron cannot be smaller than 1 (far from the core)
!
      implicit none
      integer, parameter :: dp = kind(1.d0)
      integer :: nwf, nn(nwf), ll(nwf), ndm, mesh, n, i, nspin
      real(kind=dp) :: r(ndm), vpot(ndm,2), vxt(ndm), enl(nwf), oc(nwf), &
             zed, zval, zz, zen, enne, t,x, vext
      real(kind=dp), parameter :: e2 = 2.d0
      external vext
!
      enne = 0.d0
      zz = max(zed,zval)
      do  n=1,nwf
         enne = enne + oc(n)
         zen= 0.0d0
         do  i=1,nwf
            if(nn(i).lt.nn(n)) zen=zen+oc(i)
            if(nn(i).eq.nn(n).and.ll(i).le.ll(n)) zen=zen+oc(i)
         end do
         zen = max(zz-zen+1.0d0,1.0d0)
         enl(n) =-(zen/nn(n))**2
      end do 
!
      do  i=1,mesh
         vxt(i)=vext(r(i))
         x =r(i)*enne**(1.d0/3.d0)/0.885d0
         t= zz/(1.+sqrt(x)*(0.02747d0-x*(0.1486d0-0.007298d0*x)) &
                       + x*(1.243d0+x*(0.2302d0+0.006944d0*x)))
         t = max(1.0d0,t)
         vpot(i,1) = -e2*t/r(i) + vxt(i)
      enddo
!
      if (nspin.eq.2) then
         do i=1,mesh
            vpot(i,2)=vpot(i,1)
         enddo
      endif
!
      return
      end
