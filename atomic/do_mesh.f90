!
!---------------------------------------------------------------
      subroutine do_mesh(rmax,zmesh,xmin,dx,ibound,ndm, &
                         mesh,r,r2,rab,sqr)
!---------------------------------------------------------------
!
      use kinds, only : DP
      implicit none
      integer::  ibound, ndm, mesh, i
      real(kind=dp) :: zmesh, xmin, xmax, x, dx, rmax, &
             r(ndm), r2(ndm), rab(mesh), sqr(ndm)
!      
      xmax=log(rmax*zmesh)
      mesh=(xmax-xmin)/dx+1
!
!      mesh must be odd for simpson integration. 
!
      mesh=2*(mesh/2)+1
      if(mesh+1 > ndm) &
         call errore('do_mesh','ndm is too small',1)
      if(ibound == 1) xmin=xmax-dx*(mesh-1)
!
      do i=1,mesh
         x=xmin+dble(i-1)*dx
         r(i)=exp(x)/zmesh
         r2(i)=r(i)*r(i)
         rab(i)=r(i)*dx
         sqr(i)=sqrt(r(i))
      end do
!
      return
      end
