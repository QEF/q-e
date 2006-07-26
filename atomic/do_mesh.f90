!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
      subroutine do_mesh(rmax,zmesh,xmin,dx,ibound,ndm, &
                         mesh,r,r2,rab,sqr)
!---------------------------------------------------------------
!
      use kinds, only : DP
      implicit none
      integer::  ibound, ndm, mesh, i
      real(DP) :: zmesh, xmin, xmax, x, dx, rmax, &
             r(ndm), r2(ndm), rab(ndm), sqr(ndm)
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
      !!! r(0)=0.0_dp
      do i=1,mesh
         x=xmin+DBLE(i-1)*dx
         r(i)=exp(x)/zmesh
         rab(i)=r(i)*dx
         !!! r(i)=exp(xmin)*(exp(i*dx)-1.0_dp)/zmesh
         !!! rab(i)=(r(i)+exp(xmin)/zmesh)*dx
         r2(i)=r(i)*r(i)
         sqr(i)=sqrt(r(i))
      end do
!
      return
      end subroutine do_mesh
