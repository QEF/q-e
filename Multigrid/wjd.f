c* ///////////////////////////////////////////////////////////////////////////
c* @file    wjd.f
c* @author  Michael Holst
c* @brief   Weighted Jacobi iteration.
c* @version $Id: wjd.f,v 1.1 2008-06-11 10:47:38 degironc Exp $
c* @attention
c* @verbatim
c*
c* PMG -- Parallel algebraic MultiGrid
c* Copyright (c) 1994-2006.  Michael Holst.
c*
c* Michael Holst <mholst@math.ucsd.edu>
c* University of California, San Diego
c* Department of Mathematics, 5739 AP&M
c* 9500 Gilman Drive, Dept. 0112
c* La Jolla, CA 92093-0112 USA                                                  
c* http://math.ucsd.edu/~mholst
c*
c* This file is part of PMG.
c*
c* PMG is free software; you can redistribute it and/or modify
c* it under the terms of the GNU General Public License as published by
c* the Free Software Foundation; either version 2 of the License, or
c* (at your option) any later version.
c*
c* PMG is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c* GNU General Public License for more details.
c*
c* You should have received a copy of the GNU General Public License
c* along with PMG; if not, write to the Free Software
c* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
c*
c* Linking PMG statically or dynamically with other modules is making a
c* combined work based on PMG. Thus, the terms and conditions of the GNU
c* General Public License cover the whole combination.
c* 
c* SPECIAL GPL EXCEPTION
c* In addition, as a special exception, the copyright holders of PMG
c* give you permission to combine the PMG program with free software
c* programs and libraries that are released under the GNU LGPL or with
c* code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision.
c* Such combined software may be linked with PMG and redistributed together 
c* in original or modified form as mere aggregation without requirement that 
c* the entire work be under the scope of the GNU General Public License.
c* This special exception permission is also extended to any software listed
c* in the SPECIAL GPL EXCEPTION clauses by the FEtk and APBS libraries.
c* 
c* Note that people who make modified versions of PMG are not obligated
c* to grant this special exception for their modified versions; it is
c* their choice whether to do so. The GNU General Public License gives
c* permission to release a modified version without this exception; this
c* exception also makes it possible to release a modified version which
c* carries forward this exception.
c*
c* @endverbatim
c* ///////////////////////////////////////////////////////////////////////////

      subroutine wjac(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2   itmax,iters,errtol,omega,iresid,iadjoint,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    call the fast diagonal iterative method.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          numdia
      integer          whichbc(3)
      double precision omega,errtol
      double precision rpc(*),ac(nx*ny*nz,*),cc(nx,ny,nz),fc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** do in one step ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call wjac7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),
     3      x,w1,w2,r,
     4      itmax,iters,errtol,omega,iresid,iadjoint,whichbc)
      elseif (numdia .eq. 27) then
         call wjac27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3      ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4      x,w1,w2,r,
     5      itmax,iters,errtol,omega,iresid,iadjoint)
      else
         print*,'% WJAC: invalid stencil type given...'
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine wjac7(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint,whichbc)
c* *********************************************************************
c* purpose: 
c*
c*   fast 7 diagonal weighted jacobi routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k
      integer          whichbc(3)
      integer          im,ip,jm,jp,km,kp
      double precision omega,errtol,fac
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
      double precision xdot
c*
cmdir 0 0
c*
c*    *** do the jacobi iteration itmax times ***
      do 30 iters = 1, itmax
c*
c*       *** do it ***
cmdir 3 1
         w1 = 0.d0

         do 10 k=2-whichbc(3),nz-1+whichbc(3)
         km = mod( k - 1 + nz - 1, nz ) + 1
         kp = mod( k - 1 + nz + 1, nz ) + 1
         do 11 j=2-whichbc(2),ny-1+whichbc(2)
            jm = mod( j - 1 + ny - 1, ny ) + 1
            jp = mod( j - 1 + ny + 1, ny ) + 1
            do 12 i=2-whichbc(1),nx-1+whichbc(1)
               im = mod( i - 1 + nx - 1, nx ) + 1
               ip = mod( i - 1 + nx + 1, nx ) + 1
                  w1(i,j,k) = (fc(i,j,k)+(
     2               +  oN(i,j,k)        * x(i,jp,k)
     3               +  oN(i,jm,k)      * x(i,jm,k)
     4               +  oE(i,j,k)        * x(ip,j,k)
     5               +  oE(im,j,k)      * x(im,j,k)
     6               +  uC(i,j,km)      * x(i,j,km)
     7               +  uC(i,j,k)        * x(i,j,kp)
     8               )) /  (oC(i,j,k) + cc(i,j,k))
 12            continue
 11         continue
 10      continue
c*
c*       *** copy temp back to solution, with over-relaxation ***
         fac = 1.0d0 - omega
cmdir 3 1
         do 20 k=2-whichbc(3),nz-1+whichbc(3)
cmdir 3 2
            do 21 j=2-whichbc(2),ny-1+whichbc(2)
cmdir 3 3
               do 22 i=2-whichbc(1),nx-1+whichbc(1)
                  x(i,j,k) = omega*w1(i,j,k) + fac*x(i,j,k)
 22            continue
 21         continue
 20      continue
c*
c*       *** main loop ***
 30   continue
c*
c*    *** if specified, return the new residual as well ***
      if (iresid .eq. 1) then
         call mresid7_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r,
     2                   whichbc)
         !write(*,*) "wjac residue", xdot(nx,ny,nz,r,r), itmax
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine wjac27(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*    fast 27 diagonal weighted jacobi routine.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k
      double precision omega,errtol,fac
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
      double precision tmpO,tmpU,tmpD
c*
cmdir 0 0
c*
c*    *** do the gauss-seidel iteration itmax times ***
      do 30 iters = 1, itmax
c*
c*       *** do all of the points ***
cmdir 3 1
         do 10 k=2,nz-1
cmdir 3 2
            do 11 j=2,ny-1
cmdir 3 3
               do 12 i=2,nx-1
                  tmpO =
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  oNE(i,j,k)       * x(i+1,j+1,k)
     7               +  oNW(i,j,k)       * x(i-1,j+1,k)
     8               +  oNW(i+1,j-1,k)   * x(i+1,j-1,k)
     9               +  oNE(i-1,j-1,k)   * x(i-1,j-1,k)
                  tmpU =
     2               +  uC(i,j,k)        * x(i,j,k+1)
     3               +  uN(i,j,k)        * x(i,j+1,k+1)
     4               +  uS(i,j,k)        * x(i,j-1,k+1)
     5               +  uE(i,j,k)        * x(i+1,j,k+1)
     6               +  uW(i,j,k)        * x(i-1,j,k+1)
     7               +  uNE(i,j,k)       * x(i+1,j+1,k+1)
     8               +  uNW(i,j,k)       * x(i-1,j+1,k+1)
     9               +  uSE(i,j,k)       * x(i+1,j-1,k+1)
     9               +  uSW(i,j,k)       * x(i-1,j-1,k+1)
                  tmpD =
     2               +  uC(i,j,k-1)      * x(i,j,k-1)
     3               +  uS(i,j+1,k-1)    * x(i,j+1,k-1)
     4               +  uN(i,j-1,k-1)    * x(i,j-1,k-1)
     5               +  uW(i+1,j,k-1)    * x(i+1,j,k-1)
     6               +  uE(i-1,j,k-1)    * x(i-1,j,k-1)
     7               +  uSW(i+1,j+1,k-1) * x(i+1,j+1,k-1)
     8               +  uSE(i-1,j+1,k-1) * x(i-1,j+1,k-1)
     9               +  uNW(i+1,j-1,k-1) * x(i+1,j-1,k-1)
     9               +  uNE(i-1,j-1,k-1) * x(i-1,j-1,k-1)
                  w1(i,j,k) = (fc(i,j,k)+(tmpO + tmpU + tmpD
     9               )) /  (oC(i,j,k) + cc(i,j,k))
 12            continue
 11         continue
 10      continue
c*
c*       *** copy temp back to solution, with over-relaxation ***
         fac = 1.0d0 - omega
cmdir 3 1
         do 20 k=2,nz-1
cmdir 3 2
            do 21 j=2,ny-1
cmdir 3 3
               do 22 i=2,nx-1
                  x(i,j,k) = omega*w1(i,j,k) + fac*x(i,j,k)
 22            continue
 21         continue
 20      continue
c*
c*       *** main loop ***
 30   continue
c*
c*    *** if specified, return the new residual as well ***
      if (iresid .eq. 1) then
         call mresid27_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,
     2      oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3      x,r)
      endif
c*
c*    *** return and end ***
      return
      end

