c* ///////////////////////////////////////////////////////////////////////////
c* @file    buildAd.f
c* @author  Michael Holst
c* @brief   PDE Discretization routines.
c* @version $Id: buildAd.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

      subroutine buildA (nx,ny,nz,ipkey,mgdisc,numdia,
     2   ipc,rpc,ac,cc,fc,
     3   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c* *********************************************************************
c* purpose: 
c*
c*    break the matrix data-structure into diagonals
c*    and then call the matrix build routine.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,ipkey,mgdisc,numdia
      double precision rpc(*),ac(nx*ny*nz,*),cc(*),fc(*)
      double precision xf(*),yf(*),zf(*),gxcf(*),gycf(*),gzcf(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      if (mgdisc .eq. 0) then
!         print*,'% BUILDA:   (BOX) performing discretization...'
         call buildA_fv (nx,ny,nz,ipkey,numdia,
     2      ipc,rpc,ac(1,1),cc,fc,
     3      ac(1,2),ac(1,3),ac(1,4),
     4      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
      elseif (mgdisc .eq. 1) then
!         print*,'% BUILDA:   (FEM) performing discretization...'
         call buildA_fe (nx,ny,nz,ipkey,numdia,
     2      ipc,rpc,ac(1,1),cc,fc,
     3      ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     4      ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     5      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
      else
!         print*,'% BUILDA: invalid discretization requested...'
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine buildA_fv (nx,ny,nz,ipkey,numdia,
     2   ipc,rpc,oC,cc,fc,
     3   oE,oN,uC,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c* *********************************************************************
c* purpose:  
c*
c*    box method (finite volume) discretization of a 3d pde on a
c*    tensor product (axi-parallel) three-dimensional mesh.
c*
c*    this subroutine discretizes the elliptic boundary value problem:
c*
c*          lu = f, u in omega
c*           u = g, u on boundary of omega
c*
c*    the multigrid code requires the operator in the form:
c*
c*         - \nabla \cdot (a \nabla u) + b \cdot u + c u = f
c*
c*    or: 
c*
c*        lu = (a11 u_x)_x + (a22 u_y)_y + (a33 u_z)_z
c*             + b1 u_x + b2 u_y + b3 u_z + c u
c*   
c*    here, we consider only the case: b=b1=b2=b3=0.
c*    then we assume:
c*
c*    the tensor a=diag(a11,a22,a33) has components which are
c*    then scalar functions, a11(x,y,z),a22(x,y,z),a33(x,y,z)
c*    and the functions c(x,y,z) and f(x,y,z) are also scalar.
c*    functions.  All are allowed to be possibly discontinuous on
c*    omega (the discontinuities must be along grid lines on fine grid).
c*    the boundary function g(x,y,z) is smooth on boundary of omega.
c*
c*    we will take the following conventions:
c*    (north,south,east,west refers to x-y plane.  
c*    up/down refers to the z-axis)
c*      (u(x+h_x,y,z) = u^+   (east neighbor of u(x,y,z))
c*      (u(x-h_x,y,z) = u^-   (west neighbor of u(x,y,z))
c*      (u(x,y+h_y,z) = u_+   (north neighbor of u(x,y,z))
c*      (u(x,y-h_y,z) = u_-   (south neighbor of u(x,y,z))
c*      (u(x,y,z+h_z) = u.+   (up neighbor of u(x,y,z))
c*      (u(x,y,z-h_z) = u.-   (down neighbor u(x,y,z))
c*
c*    below, we will denote:  
c*          u(x+h_x,y,z)       as u^+
c*          u(x+(1/2)h_x,y,z)  as u^++
c*          u(x-h_x,y,z)       as u^+
c*          u(x-(1/2)h_x,y,z)  as u^--
c*    and similarly for u_-,u_--,u_+,u_++,u.-,u.--,u.+,u.++.
c*
c*    we use the 3d analogue of the box method (see varga, pg. 191)
c*    which results in the following difference scheme:
c*
c*    u            : [ + (a11^++ + a11^--) * (h_y*h_z/h_x) 
c*                     + (a22_++ + a22_--) * (h_x*h_z/h_y) 
c*                     + (a33.++ + a33.--) * (h_x*h_y/h_z) 
c*                     + c * (h_x*h_y*h_z) ] u
c*    u^+ (e nbr)  : [ - (a11^++) * (h_y*h_z/h_x) ]
c*    u^_ (w nbr)  : [ - (a11^--) * (h_y*h_z/h_x) ] 
c*    u_+ (n nbr)  : [ - (a22_++) * (h_x*h_z/h_y) ] 
c*    u_- (s nbr)  : [ - (a22_--) * (h_x*h_z/h_y) ] 
c*    u.+ (u nbr)  : [ - (a33.++) * (h_x*h_y/h_z) ] 
c*    u.- (d nbr)  : [ - (a33.--) * (h_x*h_y/h_z) ] 
c*    f            : [ h_x*h_y*h_z ] f
c*
c*    note: fast way to do a conditional: we wish to set (a=coef), 
c*          unless we are on the (z=0)-boundary, which occurs with k=1:
c*
c*             coef = ...etc...
c*             ike  = min0(1,iabs(k-1))
c*             a (index1) = (ike)*coef
c*             b(index2) = b(index2) - (1-ike)*coef*bnd_data
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k,ike,jke,kke
      integer          nxm1,nym1,nzm1,ipkey,numdia
      double precision hx,hy,hz
      double precision hxm1,hym1,hzm1,coef_fc
      double precision bc_cond_e,bc_cond_w,bc_cond_n,bc_cond_s
      double precision bc_cond_u,bc_cond_d
      double precision coef_oE,coef_oN,coef_uC
      double precision coef_oEm1,coef_oNm1,coef_uCm1
      double precision diag
      double precision a1cf(nx,ny,nz),a2cf(nx,ny,nz),a3cf(nx,ny,nz)
      double precision ccf(nx,ny,nz),fcf(nx,ny,nz)
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz),oC(nx,ny,nz)
      double precision gxcf(ny,nz,*),gycf(nx,nz,*),gzcf(nx,ny,*)
      double precision xf(*),yf(*),zf(*)
      integer          im, ip, jm, jp, km, kp
c*
cmdir 0 0
c*
c*    *** save the problem key with this operator ***
      ipc(10) = ipkey
c*
c*    *** note how many nonzeros in this discretization stencil
      ipc(11) = 7
      ipc(12) = 1
      numdia  = 4
c*
c*    *** define n and determine number of mesh points ***
      nxm1    = nx - 1
      nym1    = ny - 1
      nzm1    = nz - 1
c*
c*    *** determine diag scale factor ***
c*    *** (would like something close to ones on the main diagonal) ***
      diag = 1.0d0
c*
c* *********************************************************************
c* *** interior points ***
c* *********************************************************************
c*
c*    *** build the operator ***
cmdir 3 1
      do 10 k = 1, nz
         km = mod( k - 1 + nz - 1, nz ) + 1
         kp = mod( k - 1 + nz + 1, nz ) + 1
c         hzm1 = zf(k)   - zf(km)
c         hz   = zf(kp) - zf(k)
          hzm1 = zf(2) - zf(1)
          hz = hzm1
cmdir 3 2
         do 11 j = 1, ny
            jm = mod( j - 1 + ny - 1, ny ) + 1
            jp = mod( j - 1 + ny + 1, ny ) + 1
c            hym1 = yf(j)   - yf(jm)
c            hy   = yf(jp) - yf(j)
             hym1 = yf(2) - yf(1)
             hy = hym1
cmdir 3 3
            do 12 i = 1, nx
               im = mod( i - 1 + nx - 1, nx ) + 1
               ip = mod( i - 1 + nx + 1, nx ) + 1
c               hxm1 = xf(i)   - xf(im)
c               hx   = xf(ip) - xf(i)
                hxm1 = xf(2) - xf(1)
                hx = hxm1
c*
c*             *** some coefficients ***
               coef_oE   = diag * (hym1+hy)*(hzm1+hz)/(4.*hx)
               coef_oEm1 = diag * (hym1+hy)*(hzm1+hz)/(4.*hxm1)
               coef_oN   = diag * (hxm1+hx)*(hzm1+hz)/(4.*hy)
               coef_oNm1 = diag * (hxm1+hx)*(hzm1+hz)/(4.*hym1)
               coef_uC   = diag * (hxm1+hx)*(hym1+hy)/(4.*hz)
               coef_uCm1 = diag * (hxm1+hx)*(hym1+hy)/(4.*hzm1)
               coef_fc   = diag * (hxm1+hx)*(hym1+hy)*(hzm1+hz)/8.
c*
c*             *** the coefficient and source function ***
               fc(i,j,k) = coef_fc * fcf(i,j,k)
               cc(i,j,k) = coef_fc * ccf(i,j,k)
c*
c*             *** the diagonal for matvecs and smoothings ***
               oC(i,j,k)  = 
     2            coef_oE   * a1cf(i,j,k)
     3          + coef_oEm1 * a1cf(im,j,k)
     4          + coef_oN   * a2cf(i,j,k)
     5          + coef_oNm1 * a2cf(i,jm,k)
     6          + coef_uC   * a3cf(i,j,k)
     7          + coef_uCm1 * a3cf(i,j,km)
c*
c*             *** east neighbor ***
c               ike  = min0(1,iabs(i-nxm1))
               ike = 1
               oE(i,j,k) = ike    *coef_oE*a1cf(i,j,k)
c               bc_cond_e  = (1-ike)*coef_oE*a1cf(i,j,k)*gxcf(j,k,2)
c               fc(i,j,k)  = fc(i,j,k) + bc_cond_e
c*
c*             *** north neighbor ***
c               jke  = min0(1,iabs(j-nym1))
               jke = 1
               oN(i,j,k) = jke    *coef_oN*a2cf(i,j,k)
c               bc_cond_n  = (1-jke)*coef_oN*a2cf(i,j,k)*gycf(i,k,2)
c               fc(i,j,k)  = fc(i,j,k) + bc_cond_n
c*
c*             *** up neighbor ***
c               kke  = min0(1,iabs(k-nzm1))
               kke = 1
               uC(i,j,k) = kke    *coef_uC*a3cf(i,j,k)
c               bc_cond_u  = (1-kke)*coef_uC*a3cf(i,j,k)*gzcf(i,j,2)
c               fc(i,j,k)  = fc(i,j,k) + bc_cond_u
c*
c*             *** west neighbor (just handle b.c.) ***
c               ike  = min0(1,iabs(i-2))
c               ike = 1
c               bc_cond_w  = (1-ike)*coef_oEm1*a1cf(im,j,k)*gxcf(j,k,1)
c               fc(i,j,k)  = fc(i,j,k) + bc_cond_w
c*
c*             *** south neighbor (just handle b.c.) ***
c               jke  = min0(1,iabs(j-2))
c               jke = 1
c               bc_cond_s  = (1-jke)*coef_oNm1*a2cf(i,jm,k)*gycf(i,k,1)
c               fc(i,j,k)  = fc(i,j,k) + bc_cond_s
c*
c*             *** down neighbor (just handle b.c.) ***
c               kke  = min0(1,iabs(k-2))
c               kke = 1
c               bc_cond_d  = (1-kke)*coef_uCm1*a3cf(i,j,km)*gzcf(i,j,1)
c               fc(i,j,k)  = fc(i,j,k) + bc_cond_d
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildA_fe (nx,ny,nz,ipkey,numdia,
     2   ipc,rpc,oC,cc,fc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf)
c* *********************************************************************
c* purpose:  
c*
c*    finite element method discretization of a 3d pde on a
c*    tensor product (axi-parallel) three-dimensional mesh.
c*
c*    KEY RESTRICTION: the coefficients in the pde below must
c*                     be piecewise constant in the elements
c*                     for this discretization to be formally
c*                     correct.
c*
c*    this subroutine discretizes the elliptic boundary value problem:
c*
c*          lu = f, u in omega
c*           u = g, u on boundary of omega
c*
c*    the multigrid code requires the operator in the form:
c*
c*         - \nabla \cdot (a \nabla u) + b \cdot u + c u = f
c*
c*    or: 
c*
c*        lu = (a11 u_x)_x + (a22 u_y)_y + (a33 u_z)_z
c*             + b1 u_x + b2 u_y + b3 u_z + c u
c*   
c*    here, we consider only the case: b=b1=b2=b3=0.
c*    then we assume:
c*
c*    the tensor a=diag(a11,a22,a33) has components which are
c*    then scalar functions, a11(x,y,z),a22(x,y,z),a33(x,y,z)
c*    and the functions c(x,y,z) and f(x,y,z) are also scalar.
c*    functions.  All are allowed to be possibly discontinuous on
c*    omega (the discontinuities must be along grid lines on fine grid).
c*    the boundary function g(x,y,z) is smooth on boundary of omega.
c*
c*    we will take the following conventions:
c*    (north,south,east,west refers to x-y plane.  
c*    up/down refers to the z-axis)
c*      (u(x+h_x,y,z) = u^+   (east neighbor of u(x,y,z))
c*      (u(x-h_x,y,z) = u^-   (west neighbor of u(x,y,z))
c*      (u(x,y+h_y,z) = u_+   (north neighbor of u(x,y,z))
c*      (u(x,y-h_y,z) = u_-   (south neighbor of u(x,y,z))
c*      (u(x,y,z+h_z) = u.+   (up neighbor of u(x,y,z))
c*      (u(x,y,z-h_z) = u.-   (down neighbor u(x,y,z))
c*
c*    below, we will denote:  
c*          u(x+h_x,y,z)       as u^+
c*          u(x+(1/2)h_x,y,z)  as u^++
c*          u(x-h_x,y,z)       as u^+
c*          u(x-(1/2)h_x,y,z)  as u^--
c*    and similarly for u_-,u_--,u_+,u_++,u.-,u.--,u.+,u.++.
c*
c*    we use trilinear basis functions and hexahedral elements
c*    to perform this standard finite element discretization, 
c*    which results in the following difference scheme:
c*
c*    u            : [ + (a11^++ + a11^--) * (h_y*h_z/h_x) 
c*                     + (a22_++ + a22_--) * (h_x*h_z/h_y) 
c*                     + (a33.++ + a33.--) * (h_x*h_y/h_z) 
c*                     + c * (h_x*h_y*h_z) ] u
c*    u^+ (e nbr)  : [ - (a11^++) * (h_y*h_z/h_x) ]
c*    u^_ (w nbr)  : [ - (a11^--) * (h_y*h_z/h_x) ] 
c*    u_+ (n nbr)  : [ - (a22_++) * (h_x*h_z/h_y) ] 
c*    u_- (s nbr)  : [ - (a22_--) * (h_x*h_z/h_y) ] 
c*    u.+ (u nbr)  : [ - (a33.++) * (h_x*h_y/h_z) ] 
c*    u.- (d nbr)  : [ - (a33.--) * (h_x*h_y/h_z) ] 
c*    f            : [ h_x*h_y*h_z ] f
c*
c*    note: fast way to do a conditional: we wish to set (a=coef), 
c*          unless we are on the (z=0)-boundary, which occurs with k=1:
c*
c*             coef = ...etc...
c*             ike  = min0(1,iabs(k-1))
c*             a (index1) = (ike)*coef
c*             b(index2) = b(index2) - (1-ike)*coef*bnd_data
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nx,ny,nz,i,j,k,ike,jke,kke
      integer          nxm1,nym1,nzm1,ipkey,numdia
      double precision hx,hy,hz
      double precision hxm1,hym1,hzm1,coef_fc,coef_uC
      double precision diag
      double precision a1cf(nx,ny,nz),a2cf(nx,ny,nz),a3cf(nx,ny,nz)
      double precision ccf(nx,ny,nz),fcf(nx,ny,nz)
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision oNE(nx,ny,nz),oNW(nx,ny,nz),uE(nx,ny,nz)
      double precision uW(nx,ny,nz),uN(nx,ny,nz),uS(nx,ny,nz)
      double precision uNE(nx,ny,nz),uNW(nx,ny,nz),uSE(nx,ny,nz)
      double precision uSW(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz),oC(nx,ny,nz)
      double precision gxcf(ny,nz,*),gycf(nx,nz,*),gzcf(nx,ny,*)
      double precision xf(*),yf(*),zf(*)
c*
cmdir 0 0
c*
c*    *** save the problem key with this operator ***
      ipc(10) = ipkey
c*
c*    *** note how many nonzeros in this discretization stencil
      ipc(11) = 27
      ipc(12) = 1
      numdia  = 14
c*
c*    *** define n and determine number of mesh points ***
      nxm1    = nx - 1
      nym1    = ny - 1
      nzm1    = nz - 1
c*
c*    *** determine diag scale factor ***
c*    *** (would like something close to ones on the main diagonal) ***
      diag = 1.0d0
c*
c* *********************************************************************
c* *** interior points ***
c* *********************************************************************
c*
c*    *** build the operator ***
cmdir 3 1
      do 10 k = 2, nz-1
         hzm1 = zf(k)   - zf(k-1)
         hz   = zf(k+1) - zf(k)
cmdir 3 2
         do 11 j = 2, ny-1
            hym1 = yf(j)   - yf(j-1)
            hy   = yf(j+1) - yf(j)
cmdir 3 3
            do 12 i = 2, nx-1
               hxm1 = xf(i)   - xf(i-1)
               hx   = xf(i+1) - xf(i)
c*
c*             *** some coefficients ***
               coef_uC = hx * hy / (2.0d0 * hz)
               coef_fc = hx * hy * hz
c*
c*             *** the coefficient and source function ***
               fc(i,j,k) = coef_fc * fcf(i,j,k)
               cc(i,j,k) = coef_fc * ccf(i,j,k)
c*
c*             *** the diagonal for matvecs and smoothings ***
               oC(i,j,k)  = coef_uC * 32.0d0 / 6.0d0
c*
c*             *** east neighbor ***
               ike  = min0(1,iabs(i-nxm1))
               oE(i,j,k) = 0
c*
c*             *** north neighbor ***
               jke  = min0(1,iabs(j-nym1))
               oN(i,j,k)  = 0
c*
c*             *** up neighbor ***
               kke  = min0(1,iabs(k-nzm1))
               uC(i,j,k)  = 0
c*
c*             *** north-east neighbor ***
               ike = min0(1,iabs(i-nxm1)) * min0(1,iabs(j-nym1))
               oNE(i,j,k) = ike    *coef_uC * 2.0d0 / 6.0d0
c*
c*             *** north-west neighbor ***
               ike = min0(1,iabs(i-2)) * min0(1,iabs(j-nym1))
               oNW(i,j,k) = ike    *coef_uC * 2.0d0 / 6.0d0
c*
c*             *** up-east neighbor ***
               ike = min0(1,iabs(i-nxm1)) * min0(1,iabs(k-nzm1))
               uE(i,j,k)  = ike    *coef_uC * 2.0d0 / 6.0d0
c*
c*             *** up-west neighbor ***
               ike = min0(1,iabs(i-2)) * min0(1,iabs(k-nzm1))
               uW(i,j,k)  = ike    *coef_uC * 2.0d0 / 6.0d0
c*
c*             *** up-north neighbor ***
               ike = min0(1,iabs(j-nym1)) * min0(1,iabs(k-nzm1))
               uN(i,j,k)  = ike    *coef_uC * 2.0d0 / 6.0d0
c*
c*             *** up-south neighbor ***
               ike = min0(1,iabs(j-2)) * min0(1,iabs(k-nzm1))
               uS(i,j,k)  = ike    *coef_uC * 2.0d0 / 6.0d0
c*
c*             *** up-north-east neighbor ***
               ike = min0(1,iabs(i-nxm1)) * min0(1,iabs(j-nym1)) 
     2             * min0(1,iabs(k-nzm1))
               uNE(i,j,k)  = ike    *coef_uC * 1.0d0 / 6.0d0
c*
c*             *** up-north-west neighbor ***
               ike = min0(1,iabs(i-2)) * min0(1,iabs(j-nym1)) 
     2             * min0(1,iabs(k-nzm1))
               uNW(i,j,k)  = ike    *coef_uC * 1.0d0 / 6.0d0
c*
c*             *** up-south-east neighbor ***
               ike = min0(1,iabs(i-nxm1)) * min0(1,iabs(j-2)) 
     2             * min0(1,iabs(k-nzm1))
               uSE(i,j,k)  = ike    *coef_uC * 1.0d0 / 6.0d0
c*
c*             *** up-south-west neighbor ***
               ike = min0(1,iabs(i-2)) * min0(1,iabs(j-2)) 
     2             * min0(1,iabs(k-nzm1))
               uSW(i,j,k)  = ike    *coef_uC * 1.0d0 / 6.0d0
c*
 12   continue
 11   continue
 10   continue
c*
c*    *** return and end ***
      return
      end
