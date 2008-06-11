c* ///////////////////////////////////////////////////////////////////////////
c* @file    ninterpd.f
c* @author  Michael Holst
c* @brief   Nonlinear version of operator-based interpolation.
c* @version $Id: ninterpd.f,v 1.1 2008-06-11 10:47:38 degironc Exp $
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

      subroutine ninterp(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,pc,
     2   ipc,rpc,ac,cc,fc)
c* *********************************************************************
c* purpose:
c*
c*    driver for nonlinear operator-based interpolation.
c*
c* author:  michael holst
c* *********************************************************************
      integer          nxf,nyf,nzf,nxc,nyc,nzc,numdia
      double precision xin(*),xout(*),pc(nxc*nyc*nzc,*)
      integer          ipc(*)
      double precision rpc(*),ac(nxf*nyf*nzf,*),cc(*),fc(*)
c*
c*    *** doit ***
      numdia = ipc(11)
      if (numdia .eq. 7) then
         call ninterp7(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,
     2      pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),
     3      pc(1,7),pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),
     4      pc(1,13),pc(1,14),pc(1,15),pc(1,16),pc(1,17),pc(1,18),
     5      pc(1,19),pc(1,20),pc(1,21),pc(1,22),pc(1,23),pc(1,24),
     6      pc(1,25),pc(1,26),pc(1,27),
     7      ipc,rpc,ac(1,1),ac(1,2),ac(1,3),ac(1,4),cc,fc)
      elseif (numdia .eq. 27) then
         call ninterp27(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,
     2      pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),
     3      pc(1,7),pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),
     4      pc(1,13),pc(1,14),pc(1,15),pc(1,16),pc(1,17),pc(1,18),
     5      pc(1,19),pc(1,20),pc(1,21),pc(1,22),pc(1,23),pc(1,24),
     6      pc(1,25),pc(1,26),pc(1,27),
     7      ipc,rpc,ac(1,1),ac(1,2),ac(1,3),ac(1,4),
     8      ac(1,5),ac(1,6),ac(1,7),ac(1,8),ac(1,9),ac(1,10),ac(1,11),
     9      ac(1,12),ac(1,13),ac(1,14),cc,fc)
      else
         print*,'% NINTERP: invalid stencil type given...'
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine ninterp7(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,
     2   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     3   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     4   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW,
     5   ipc,rpc,oC,oE,oN,uC,cc,fc)
c* *********************************************************************
c* purpose:
c*
c*    7-diagonal nonlinear operator-based interpolation.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxc,nyc,nzc,nxf,nyf,nzf,i,j,k,jj,ii,kk
      double precision xin(nxc,nyc,nzc),xout(nxf,nyf,nzf)
      double precision oPC(nxc,nyc,nzc),oPN(nxc,nyc,nzc)
      double precision oPS(nxc,nyc,nzc),oPE(nxc,nyc,nzc)
      double precision oPW(nxc,nyc,nzc),oPNE(nxc,nyc,nzc)
      double precision oPNW(nxc,nyc,nzc),oPSE(nxc,nyc,nzc)
      double precision oPSW(nxc,nyc,nzc),uPC(nxc,nyc,nzc)
      double precision uPN(nxc,nyc,nzc),uPS(nxc,nyc,nzc)
      double precision uPE(nxc,nyc,nzc),uPW(nxc,nyc,nzc)
      double precision uPNE(nxc,nyc,nzc),uPNW(nxc,nyc,nzc)
      double precision uPSE(nxc,nyc,nzc),uPSW(nxc,nyc,nzc)
      double precision dPC(nxc,nyc,nzc),dPN(nxc,nyc,nzc)
      double precision dPS(nxc,nyc,nzc),dPE(nxc,nyc,nzc)
      double precision dPW(nxc,nyc,nzc),dPNE(nxc,nyc,nzc)
      double precision dPNW(nxc,nyc,nzc),dPSE(nxc,nyc,nzc)
      double precision dPSW(nxc,nyc,nzc)
      integer          ipc(*)
      double precision rpc(*),cc(nxf,nyf,nzf),fc(nxf,nyf,nzf)
      double precision oC(nxf,nyf,nzf)
      double precision oE(nxf,nyf,nzf),oN(nxf,nyf,nzf),uC(nxf,nyf,nzf)
      integer          i1,j1,k1,key
      double precision c_scal,dc_scal
      double precision UU,AA,BB,DD,FF,zNN,DzNN,zntol,zndenom,znerror
      double precision change
      integer          nitmax,niters,ifail_tol,ipkey
c*
c* *********************************************************************
c* setup 
c* *********************************************************************
c*
cmdir 0 0
c*
c*    *** nonlinear iteration tolerance and itmax ***
      nitmax    = 10
      zntol     = 1.0e-5
      ifail_tol = 0
      ipkey     = ipc(10)
      key       = 0
c*
c*    *** verify correctness of the input boundary points ***
      call fbound00(nxc,nyc,nzc,xin)
c*
c*    *** doit ***
cmdir 3 1
      do 10 k = 1, nzf-2, 2
         kk = (k-1) / 2 + 1
cmdir 3 2
         do 11 j = 1, nyf-2, 2
            jj = (j-1) / 2 + 1
cmdir 3 3
            do 12 i = 1, nxf-2, 2
               ii = (i-1) / 2 + 1
c*
c* *********************************************************************
c* type 1 -- fine grid points common to a coarse grid point
c* *********************************************************************
c*
c*             *** copy coinciding points from coarse grid to fine grid ***
               xout(i,j,k) = xin(ii,jj,kk)
c*
c* *********************************************************************
c* type 2 -- fine grid points common to a coarse grid plane
c* *********************************************************************
c*
c* *********************************************************************
c*             *** fine grid pts common only to y-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on x-row) ***
               xout(i+1,j,k) = oPE(ii,jj,kk)   * xin(ii,jj,kk) 
     2                       + oPW(ii+1,jj,kk) * xin(ii+1,jj,kk)
               i1 = i+1
               j1 = j
               k1 = k
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - uC(i1,j1,k1) 
     3               - uC(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1))
     4               - oN(i1,j1,k1) 
     5               - oN(i1,max0(1,j1-1),k1) * dble(min0(1,j1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 23
               if (niters .gt. nitmax) goto 22
 21            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 23
                  if (niters .gt. nitmax) goto 22
c*
c*             *** loop ***
               goto 21
c*
c*             *** tolerance not reached ***
 22            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 23            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to x-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a y-row) ***
               xout(i,j+1,k) = oPN(ii,jj,kk)   * xin(ii,jj,kk)
     2                       + oPS(ii,jj+1,kk) * xin(ii,jj+1,kk) 
               i1 = i
               j1 = j+1
               k1 = k
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - oE(i1,j1,k1) 
     3               - oE(max0(1,i1-1),j1,k1) * dble(min0(1,i1-1))
     4               - uC(i1,j1,k1) 
     5               - uC(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 33
               if (niters .gt. nitmax) goto 32
 31            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 33
                  if (niters .gt. nitmax) goto 32
c*
c*             *** loop ***
               goto 31
c*
c*             *** tolerance not reached ***
 32            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 33            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to x-y planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a z-row) ***
               xout(i,j,k+1) = uPC(ii,jj,kk)   * xin(ii,jj,kk)
     2                       + dPC(ii,jj,kk+1) * xin(ii,jj,kk+1) 
               i1 = i
               j1 = j
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1)
     2               - oE(i1,j1,k1)
     3               - oE(max0(1,i1-1),j1,k1) * dble(min0(1,i1-1))
     4               - oN(i1,j1,k1)
     5               - oN(i1,max0(1,j1-1),k1) * dble(min0(1,j1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 43
               if (niters .gt. nitmax) goto 42
 41            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 43
                  if (niters .gt. nitmax) goto 42
c*
c*             *** loop ***
               goto 41
c*
c*             *** tolerance not reached ***
 42            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 43            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c* type 3 -- fine grid points common to a coarse grid line
c* *********************************************************************
c*
c* *********************************************************************
c*             *** fine grid pts common only to z planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xy-plane***
               xout(i+1,j+1,k) = oPNE(ii,jj,kk)     * xin(ii,jj,kk)   
     2                         + oPNW(ii+1,jj,kk)   * xin(ii+1,jj,kk)
     3                         + oPSE(ii,jj+1,kk)   * xin(ii,jj+1,kk)
     4                         + oPSW(ii+1,jj+1,kk) * xin(ii+1,jj+1,kk)
               i1 = i+1
               j1 = j+1
               k1 = k
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - uC(i1,j1,k1) 
     3               - uC(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 53
               if (niters .gt. nitmax) goto 52
 51            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 53
                  if (niters .gt. nitmax) goto 52
c*
c*             *** loop ***
               goto 51
c*
c*             *** tolerance not reached ***
 52            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 53            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to y planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xz-plane***
               xout(i+1,j,k+1) = uPE(ii,jj,kk)     * xin(ii,jj,kk)
     2                         + uPW(ii+1,jj,kk)   * xin(ii+1,jj,kk)
     3                         + dPE(ii,jj,kk+1)   * xin(ii,jj,kk+1)
     4                         + dPW(ii+1,jj,kk+1) * xin(ii+1,jj,kk+1)
               i1 = i+1
               j1 = j
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - oN(i1,j1,k1) 
     3               - oN(i1,max0(1,j1-1),k1) * dble(min0(1,j1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 63
               if (niters .gt. nitmax) goto 62
 61            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 63
                  if (niters .gt. nitmax) goto 62
c*
c*             *** loop ***
               goto 61
c*
c*             *** tolerance not reached ***
 62            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 63            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to x planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the yz-plane***
               xout(i,j+1,k+1) = uPN(ii,jj,kk)     * xin(ii,jj,kk)
     2                         + uPS(ii,jj+1,kk)   * xin(ii,jj+1,kk)
     3                         + dPN(ii,jj,kk+1)   * xin(ii,jj,kk+1)
     4                         + dPS(ii,jj+1,kk+1) * xin(ii,jj+1,kk+1)
               i1 = i
               j1 = j+1
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - oE(i1,j1,k1) 
     3               - oE(max0(1,i1-1),j1,k1) * dble(min0(1,i1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 73
               if (niters .gt. nitmax) goto 72
 71            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 73
                  if (niters .gt. nitmax) goto 72
c*
c*             *** loop ***
               goto 71
c*
c*             *** tolerance not reached ***
 72            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 73            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c* type 4 -- fine grid points not common to coarse grid pts/lines/planes
c* *********************************************************************
c*
c*             *** completely interior points ***
               xout(i+1,j+1,k+1) = 
     2            + uPNE(ii,jj,kk)       * xin(ii,jj,kk)
     3            + uPNW(ii+1,jj,kk)     * xin(ii+1,jj,kk)
     4            + uPSE(ii,jj+1,kk)     * xin(ii,jj+1,kk)
     5            + uPSW(ii+1,jj+1,kk)   * xin(ii+1,jj+1,kk)
     6            + dPNE(ii,jj,kk+1)     * xin(ii,jj,kk+1)
     7            + dPNW(ii+1,jj,kk+1)   * xin(ii+1,jj,kk+1)
     8            + dPSE(ii,jj+1,kk+1)   * xin(ii,jj+1,kk+1)
     9            + dPSW(ii+1,jj+1,kk+1) * xin(ii+1,jj+1,kk+1)
               i1 = i+1
               j1 = j+1
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = oC(i1,j1,k1)
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 83
               if (niters .gt. nitmax) goto 82
 81            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 83
                  if (niters .gt. nitmax) goto 82
c*
c*             *** loop ***
               goto 81
c*
c*             *** tolerance not reached ***
 82            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 83            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** verify correctness of the output boundary points ***
      call fbound00(nxf,nyf,nzf,xout)
c*
c*    *** messages ***
      if (ifail_tol .gt. 0) then
         print*,'% NINTERP7: 1d newton tolerance failures: ',ifail_tol
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine ninterp27(nxc,nyc,nzc,nxf,nyf,nzf,xin,xout,
     2   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     3   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     4   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW,
     5   ipc,rpc,
     6   oC,oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     7   cc,fc)
c* *********************************************************************
c* purpose:
c*
c*    27-diagonal nonlinear operator-based interpolation.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxc,nyc,nzc,nxf,nyf,nzf,i,j,k,jj,ii,kk
      double precision xin(nxc,nyc,nzc),xout(nxf,nyf,nzf)
      double precision oPC(nxc,nyc,nzc),oPN(nxc,nyc,nzc)
      double precision oPS(nxc,nyc,nzc),oPE(nxc,nyc,nzc)
      double precision oPW(nxc,nyc,nzc),oPNE(nxc,nyc,nzc)
      double precision oPNW(nxc,nyc,nzc),oPSE(nxc,nyc,nzc)
      double precision oPSW(nxc,nyc,nzc),uPC(nxc,nyc,nzc)
      double precision uPN(nxc,nyc,nzc),uPS(nxc,nyc,nzc)
      double precision uPE(nxc,nyc,nzc),uPW(nxc,nyc,nzc)
      double precision uPNE(nxc,nyc,nzc),uPNW(nxc,nyc,nzc)
      double precision uPSE(nxc,nyc,nzc),uPSW(nxc,nyc,nzc)
      double precision dPC(nxc,nyc,nzc),dPN(nxc,nyc,nzc)
      double precision dPS(nxc,nyc,nzc),dPE(nxc,nyc,nzc)
      double precision dPW(nxc,nyc,nzc),dPNE(nxc,nyc,nzc)
      double precision dPNW(nxc,nyc,nzc),dPSE(nxc,nyc,nzc)
      double precision dPSW(nxc,nyc,nzc)
      integer          ipc(*)
      double precision rpc(*),cc(nxf,nyf,nzf),fc(nxf,nyf,nzf)
      double precision oC(nxf,nyf,nzf),oE(nxf,nyf,nzf)
      double precision oN(nxf,nyf,nzf),uC(nxf,nyf,nzf)
      double precision oNE(nxf,nyf,nzf),oNW(nxf,nyf,nzf)
      double precision uE(nxf,nyf,nzf),uW(nxf,nyf,nzf)
      double precision uN(nxf,nyf,nzf),uS(nxf,nyf,nzf)
      double precision uNE(nxf,nyf,nzf),uNW(nxf,nyf,nzf)
      double precision uSE(nxf,nyf,nzf),uSW(nxf,nyf,nzf)
      integer          i1,j1,k1,key
      double precision c_scal,dc_scal
      double precision UU,AA,BB,DD,FF,zNN,DzNN,zntol,zndenom,znerror
      double precision change
      integer          nitmax,niters,ifail_tol,ipkey
c*
c* *********************************************************************
c* setup 
c* *********************************************************************
c*
cmdir 0 0
c*
c*    *** nonlinear iteration tolerance and itmax ***
      nitmax    = 10
      zntol     = 1.0e-5
      ifail_tol = 0
      ipkey     = ipc(10)
      key       = 0
c*
c*    *** verify correctness of the input boundary points ***
      call fbound00(nxc,nyc,nzc,xin)
c*
c*    *** doit ***
cmdir 3 1
      do 10 k = 1, nzf-2, 2
         kk = (k-1) / 2 + 1
cmdir 3 2
         do 11 j = 1, nyf-2, 2
            jj = (j-1) / 2 + 1
cmdir 3 3
            do 12 i = 1, nxf-2, 2
               ii = (i-1) / 2 + 1
c*
c* *********************************************************************
c* type 1 -- fine grid points common to a coarse grid point
c* *********************************************************************
c*
c*             *** copy coinciding points from coarse grid to fine grid ***
               xout(i,j,k) = xin(ii,jj,kk)
c*
c* *********************************************************************
c* type 2 -- fine grid points common to a coarse grid plane
c* *********************************************************************
c*
c* *********************************************************************
c*             *** fine grid pts common only to y-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on x-row) ***
               xout(i+1,j,k) = oPE(ii,jj,kk)   * xin(ii,jj,kk) 
     2                       + oPW(ii+1,jj,kk) * xin(ii+1,jj,kk)
               i1 = i+1
               j1 = j
               k1 = k
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - uC(i1,j1,k1) 
     3               - uC(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1))
     4               - oN(i1,j1,k1) 
     5               - oN(i1,max0(1,j1-1),k1) * dble(min0(1,j1-1))
     6               - uN(i1,j1,k1) 
     7               - uS(i1,j1,k1) 
     8               - uN(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1))
     9               - uS(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 23
               if (niters .gt. nitmax) goto 22
 21            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 23
                  if (niters .gt. nitmax) goto 22
c*
c*             *** loop ***
               goto 21
c*
c*             *** tolerance not reached ***
 22            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 23            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to x-z planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a y-row) ***
               xout(i,j+1,k) = oPN(ii,jj,kk)   * xin(ii,jj,kk)
     2                       + oPS(ii,jj+1,kk) * xin(ii,jj+1,kk) 
               i1 = i
               j1 = j+1
               k1 = k
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - oE(i1,j1,k1) 
     3               - oE(max0(1,i1-1),j1,k1) * dble(min0(1,i1-1))
     4               - uC(i1,j1,k1) 
     5               - uC(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1))
     6               - uE(i1,j1,k1) 
     7               - uW(i1,j1,k1) 
     8               - uE(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1))
     9               - uW(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 33
               if (niters .gt. nitmax) goto 32
 31            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 33
                  if (niters .gt. nitmax) goto 32
c*
c*             *** loop ***
               goto 31
c*
c*             *** tolerance not reached ***
 32            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 33            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to x-y planes on coarse grid ***
c*             *** (intermediate pts between 2 grid points on a z-row) ***
               xout(i,j,k+1) = uPC(ii,jj,kk)   * xin(ii,jj,kk)
     2                       + dPC(ii,jj,kk+1) * xin(ii,jj,kk+1) 
               i1 = i
               j1 = j
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1)
     2               - oE(i1,j1,k1)
     3               - oE(max0(1,i1-1),j1,k1) * dble(min0(1,i1-1))
     4               - oN(i1,j1,k1)
     5               - oN(i1,max0(1,j1-1),k1) * dble(min0(1,j1-1))
     6               - oNE(i1,j1,k1)
     7               - oNW(i1,j1,k1)
     8               - oNE(i1-1,j1-1,k1) * dble(min0(1,i1-1)) 
     8                                   * dble(min0(1,j1-1))
     9               - oNW(i1+1,j1-1,k1) * dble(min0(1,j1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 43
               if (niters .gt. nitmax) goto 42
 41            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 43
                  if (niters .gt. nitmax) goto 42
c*
c*             *** loop ***
               goto 41
c*
c*             *** tolerance not reached ***
 42            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 43            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c* type 3 -- fine grid points common to a coarse grid line
c* *********************************************************************
c*
c* *********************************************************************
c*             *** fine grid pts common only to z planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xy-plane***
               xout(i+1,j+1,k) = oPNE(ii,jj,kk)     * xin(ii,jj,kk)   
     2                         + oPNW(ii+1,jj,kk)   * xin(ii+1,jj,kk)
     3                         + oPSE(ii,jj+1,kk)   * xin(ii,jj+1,kk)
     4                         + oPSW(ii+1,jj+1,kk) * xin(ii+1,jj+1,kk)
               i1 = i+1
               j1 = j+1
               k1 = k
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - uC(i1,j1,k1) 
     3               - uC(i1,j1,max0(1,k1-1)) * dble(min0(1,k1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 53
               if (niters .gt. nitmax) goto 52
 51            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 53
                  if (niters .gt. nitmax) goto 52
c*
c*             *** loop ***
               goto 51
c*
c*             *** tolerance not reached ***
 52            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 53            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to y planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the xz-plane***
               xout(i+1,j,k+1) = uPE(ii,jj,kk)     * xin(ii,jj,kk)
     2                         + uPW(ii+1,jj,kk)   * xin(ii+1,jj,kk)
     3                         + dPE(ii,jj,kk+1)   * xin(ii,jj,kk+1)
     4                         + dPW(ii+1,jj,kk+1) * xin(ii+1,jj,kk+1)
               i1 = i+1
               j1 = j
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - oN(i1,j1,k1) 
     3               - oN(i1,max0(1,j1-1),k1) * dble(min0(1,j1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 63
               if (niters .gt. nitmax) goto 62
 61            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 63
                  if (niters .gt. nitmax) goto 62
c*
c*             *** loop ***
               goto 61
c*
c*             *** tolerance not reached ***
 62            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 63            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c*             *** fine grid pts common only to x planes on coarse grid ***
c*             *** (intermediate pts between 4 grid pts on the yz-plane***
               xout(i,j+1,k+1) = uPN(ii,jj,kk)     * xin(ii,jj,kk)
     2                         + uPS(ii,jj+1,kk)   * xin(ii,jj+1,kk)
     3                         + dPN(ii,jj,kk+1)   * xin(ii,jj,kk+1)
     4                         + dPS(ii,jj+1,kk+1) * xin(ii,jj+1,kk+1)
               i1 = i
               j1 = j+1
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = (  oC(i1,j1,k1) 
     2               - oE(i1,j1,k1) 
     3               - oE(max0(1,i1-1),j1,k1) * dble(min0(1,i1-1)) )
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 73
               if (niters .gt. nitmax) goto 72
 71            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 73
                  if (niters .gt. nitmax) goto 72
c*
c*             *** loop ***
               goto 71
c*
c*             *** tolerance not reached ***
 72            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 73            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c* *********************************************************************
c* type 4 -- fine grid points not common to coarse grid pts/lines/planes
c* *********************************************************************
c*
c*             *** completely interior points ***
               xout(i+1,j+1,k+1) = 
     2            + uPNE(ii,jj,kk)       * xin(ii,jj,kk)
     3            + uPNW(ii+1,jj,kk)     * xin(ii+1,jj,kk)
     4            + uPSE(ii,jj+1,kk)     * xin(ii,jj+1,kk)
     5            + uPSW(ii+1,jj+1,kk)   * xin(ii+1,jj+1,kk)
     6            + dPNE(ii,jj,kk+1)     * xin(ii,jj,kk+1)
     7            + dPNW(ii+1,jj,kk+1)   * xin(ii+1,jj,kk+1)
     8            + dPSE(ii,jj+1,kk+1)   * xin(ii,jj+1,kk+1)
     9            + dPSW(ii+1,jj+1,kk+1) * xin(ii+1,jj+1,kk+1)
               i1 = i+1
               j1 = j+1
               k1 = k+1
c*
c*             *** setup the 1d equation ***
               DD = oC(i1,j1,k1)
               if (DD .eq. 0.0d0) then
                  DD = 1.0d0
               else
                  DD = 1.0d0 / DD
               endif
               UU = xout(i1,j1,k1)
               AA = 1.0d0
               BB = cc(i1,j1,k1)
               FF = key * fc(i1,j1,k1) + xout(i1,j1,k1)
c*
c*             *** evaluate residual of 1d system
               zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*             *** setup for newton's method
               niters = 0
               zndenom = zNN
               if (zndenom .eq. 0.0d0) zndenom = 1.0d0
               znerror = zNN / zndenom
c*
c*             *** the 1d newton's method ***
               if (znerror .lt. zntol) goto 83
               if (niters .gt. nitmax) goto 82
 81            continue
                  niters = niters + 1
c*
c*                *** construct jacobian matrix of NN ***
                  DzNN = AA + DD * dc_scal(BB,UU,ipkey)
c*
c*                *** solve the linear system ***
                  change = - zNN / DzNN
c*
c*                *** update the solution ***
                  UU = UU + change
c*
c*                *** evaluate residual of 1d system ***
                  zNN = AA * UU + DD * c_scal(BB,UU,ipkey) - FF
c*
c*                *** compute error ***
                  znerror = zNN / zndenom
c*
c*                *** tolerance and itmax check ***
                  if (znerror .lt. zntol) goto 83
                  if (niters .gt. nitmax) goto 82
c*
c*             *** loop ***
               goto 81
c*
c*             *** tolerance not reached ***
 82            continue
               ifail_tol = ifail_tol + 1
c*
c*             *** tolerance reached ***
 83            continue
c*
c*             *** newton's method complete -- update solution value ***
               xout(i1,j1,k1) = UU
c*
c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** verify correctness of the output boundary points ***
      call fbound00(nxf,nyf,nzf,xout)
c*
c*    *** messages ***
      if (ifail_tol .gt. 0) then
         print*,'% NINTERP27: 1d newton tolerance failures: ',ifail_tol
      endif
c*
c*    *** return and end ***
      return
      end

