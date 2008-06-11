c* ///////////////////////////////////////////////////////////////////////////
c* @file    pded.f
c* @author  Michael Holst
c* @brief   PDE definition file.
c* @version $Id: pded.f,v 1.1 2008-06-11 10:47:38 degironc Exp $
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

c* ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
c*
c* DIFFERENTIAL EQUATION: Poisson's
c*
c* BOUNDARY CONDITIONS:
c*
c*     East   Face (xmin):  Dirichlet, homogeneous
c*     West   Face (xmax):  Dirichlet, homogeneous
c*     North  Face (ymin):  Dirichlet, homogeneous
c*     South  Face (ymax):  Dirichlet, homogeneous
c*     Top    Face (zmin):  Dirichlet, homogeneous
c*     Bottom Face (zmax):  Dirichlet, homogeneous
c*
c* MESH:                  hx  = (xmax-xmin) / (nx-1)
c*                        hy  = (ymax-ymin) / (ny-1)
c*                        hz  = (zmax-zmin) / (nz-1)
c*                        xi = xmin + (i-1) * hx,  i=1,...,nx
c*                        yi = ymin + (j-1) * hy,  j=1,...,ny
c*                        zi = zmin + (k-1) * hk,  k=1,...,nz
c*
c* CHOSEN TRUE SOLUTION:  u(x,y,z) = Sin(n pi x)*Sin(m pi y)*Sin(l pi z)
c*
c* author:  michael holst
c* ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
c*
      function c_scal(coef,u,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the nonlinearity (scalar version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision coef,u,c_scal
      integer          ipkey
c*
c*    *** create the nonlinear term ***
      c_scal = coef * u
c*
c*    *** end it ***
      return
      end
      function dc_scal(coef,u,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the derivative of the nonlinearity (scalar version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      double precision coef,u,dc_scal
      integer          ipkey
c*
c*    *** create the nonlinear term ***
      dc_scal = coef
c*
c*    *** end it ***
      return
      end
      subroutine c_vec(coef,uin,uout,nx,ny,nz,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the nonlinearity (vector version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,ipkey
      double precision uin(*),uout(*),coef(*)
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            uout(i) = coef(i) * uin(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         uout(i) = coef(i) * uin(i)
 20   continue
c*
c*    *** end it ***
      return
      end
      subroutine dc_vec(coef,uin,uout,nx,ny,nz,ipkey)
c* *********************************************************************
c* purpose:
c*
c*    define the derivative of the nonlinearity (vector version)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nx,ny,nz,ipkey
      double precision uin(*),uout(*),coef(*)
      integer          n,i,ii,nproc,ipara,ivect
      parameter        (nproc=1)
c*
cmdir 0 0
c*
c*    *** find parallel loops (ipara), remainder (ivect) ***
      n     = nx * ny * nz
      ipara = n / nproc
      ivect = mod(n,nproc)
c*
c*    *** do parallel loops ***
cmdir 2 1
      do 10 ii = 1, nproc
cmdir 2 2
         do 11 i = 1+(ipara*(ii-1)), ipara*ii
            uout(i) = coef(i)
 11      continue
 10   continue
c*
c*    *** do vector loops ***
cmdir 1 1
      do 20 i = ipara*nproc+1, n
         uout(i) = coef(i)
 20   continue
c*
c*    *** end it ***
      return
      end
      subroutine fillco(iparm,rparm,nx,ny,nz,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,
     3   source,v,eps,kap,mr1,mr2,mr3,d1,d2,d3)
c* *********************************************************************
c* purpose:
c*
c*   this subroutine defines poisson's equation on the unit cube.
c*   boundary conditions are zero dirichlet.
c*
c* details:
c*
c*   this routine sets up the coefficients of the following
c*   three-dimensional, 2nd order linear elliptic partial 
c*   differential equation:
c*
c*      lu = f, u in omega
c*       u = g, u on boundary of omega
c*
c*   where
c*      omega = [xmin,xmax]x[ymin,ymax]x[zmin,zmax]
c*
c*   and l is taken to be in the following form:
c*
c*      lu = - \nabla \cdot (a \nabla u) + c u
c*
c*   the coefficients must be supplied as follows:
c*
c*      a1cf(,,,):  at x-half grid points, hence array is:  (nx-1)*ny*nz
c*      a2cf(,,,):  at y-half grid points, hence array is:  nx*(ny-1)*nz
c*      a3cf(,,,):  at z-half grid points, hence array is:  nx*ny*(nz-1)
c*                  (we index all three as:  nx*ny*nz however)
c*
c*      ccf:        at grid points, hence array is:  nx*ny*nz
c*      fcf:        at grid points, hence array is:  nx*ny*nz
c*      u:          must be set to have appropriate boundary values
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iparm(*),nx,ny,nz
      double precision rparm(*)
      double precision a1cf(nx,ny,nz),a2cf(nx,ny,nz),a3cf(nx,ny,nz)
      double precision ccf(nx,ny,nz),fcf(nx,ny,nz),tcf(nx,ny,nz)
      double precision xf(nx),yf(ny),zf(nz)
      double precision gxcf(ny,nz,*),gycf(nx,nz,*),gzcf(nx,ny,*)
      integer          mr1, mr2, mr3
      double precision d1, d2, d3
      double precision source(mr1,mr2,mr3), v(mr1,mr2,mr3),
     2                 eps(mr1,mr2,mr3,3), kap(mr1,mr2,mr3)
c*
c*    *** variable declarations ***
      integer          i,j,k,iinfo
      double precision hx,hy,hz,xmin,xmax,ymin,ymax,zmin,zmax
      double precision hxo2,hyo2,hzo2
c*
c*    *** some parameters ***
      double precision zn,zm,zl,pi
      parameter        (zn=1.0d0,zm=1.0d0,zl=1.0d0)
      pi               = 4.0d0 * datan(1.0d0)
c*
cmdir 0 0
c*
c*    *** info messages ***
      iinfo = iparm(12)
c*
c* *********************************************************************
c* *** definition of the domain
c* *********************************************************************
c*
c*    *** set the correct boundary ranges ***
      xmin = 0.0d0
      xmax = d1 * dble( mr1 - 1 )
      ymin = 0.0d0
      ymax = d2 * dble( mr2 - 1 )
      zmin = 0.0d0
      zmax = d3 * dble( mr3 - 1 )
c*
c*    *** store in rparm array ***
      rparm(3)  = xmin
      rparm(4)  = xmax
      rparm(5)  = ymin
      rparm(6)  = ymax
      rparm(7)  = zmin
      rparm(8)  = zmax
c*
c* *********************************************************************
c* *** boundary value definitions
c* *********************************************************************
c*
c*    *** determine mesh widths ***
      hx  = (xmax-xmin) / dble(nx-1)
      hy  = (ymax-ymin) / dble(ny-1)
      hz  = (zmax-zmin) / dble(nz-1)
      hxo2= hx / 2.0d0
      hyo2= hy / 2.0d0
      hzo2= hz / 2.0d0
c*
c* *********************************************************************
c* *** pde coefficient value definitions
c* *********************************************************************
c*
c*    *** fill coefficient arrays ***
cmdir 3 1
      do 10 k = 1, nz
         zf(k) = zmin + dble(k-1) * hz
cmdir 3 2
         do 11 j = 1, ny
            yf(j) = ymin + dble(j-1) * hy
cmdir 3 3
            do 12 i = 1, nx
               xf(i) = xmin + dble(i-1) * hx
c*
c*             *** the diagonal tensor (2nd derivative) entries ***
               a1cf(i,j,k) = eps(i,j,k,1)
               a2cf(i,j,k) = eps(i,j,k,2)
               a3cf(i,j,k) = eps(i,j,k,3)
c*
c*             *** the scalar (0th derivative) entry ***
               ccf(i,j,k) = kap(i,j,k)
c*
c*             *** the rhs entry ***
CZZZ           fcf(i,j,k) = ((zn*zn+zm*zm+zl*zl)*pi*pi
CZZZ 2         *dsin(zn*pi*xf(i))*dsin(zm*pi*yf(j))*dsin(zl*pi*zf(k)))
               fcf(i,j,k) = source(i,j,k)
c*
c*             *** the chosen true solution ***
               tcf(i,j,k) = v(i,j,k)
c     2         (dsin(zn*pi*xf(i))*dsin(zm*pi*yf(j))*dsin(zl*pi*zf(k)))
 12         continue
 11      continue
 10   continue
c*
c*    *** the (i=1) and (i=nx) boundaries (dirichlet) ***
cmdir 2 1
      do 50 k = 1, nz
cmdir 2 2
         do 51 j = 1, ny
            gxcf(j,k,1) = tcf(1 ,j,k)
            gxcf(j,k,2) = tcf(nx,j,k)
            gxcf(j,k,3) = 0.0d0
            gxcf(j,k,4) = 0.0d0
 51      continue
 50   continue
c*
c*    *** the (j=1) and (j=ny) boundaries (dirichlet) ***
cmdir 2 1
      do 60 k = 1, nz
cmdir 2 2
         do 61 i = 1, nx
            gycf(i,k,1) = tcf(i,1 ,k)
            gycf(i,k,2) = tcf(i,ny,k)
            gycf(i,k,3) = 0.0d0
            gycf(i,k,4) = 0.0d0
 61      continue
 60   continue
c*
c*    *** the (k=1) and (k=nz) boundaries (dirichlet) ***
cmdir 2 1
      do 70 j = 1, ny
cmdir 2 2
         do 71 i = 1, nx
            gzcf(i,j,1) = tcf(i,j,1 )
            gzcf(i,j,2) = tcf(i,j,nz)
            gzcf(i,j,3) = 0.0d0
            gzcf(i,j,4) = 0.0d0
 71      continue
 70   continue
c*
c*    *** end it ***
      return
      end

