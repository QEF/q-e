c* ///////////////////////////////////////////////////////////////////////////
c* @file    buildPd.f
c* @author  Michael Holst
c* @brief   Routines to build the normal and opertor-based prolongations.
c* @version $Id: buildPd.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

      subroutine buildP (nxf,nyf,nzf,nxc,nyc,nzc,mgprol,
     2   ipc,rpc,pc,ac,xf,yf,zf)
c* *********************************************************************
c* purpose: 
c*
c*    form the prolongation operator, which is one of:
c*
c*       (1) trilinear interpolation
c*       (2) standard operator-based prolongation
c*       (3) a modified operator-based prolongation
c*
c*    we differentiate between 7 and 27 point fine grid matrix in
c*    each case for efficiency in building the prolongation operator.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc,mgprol,numdia
      double precision rpc(*),pc(nxc*nyc*nzc,*),ac(nxf*nyf*nzf,*)
      double precision xf(*),yf(*),zf(*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      if (mgprol .eq. 0) then
!         print*,'% BUILDP:   (TRILINEAR) building prolongation...'
         call buildP_trilin (nxf,nyf,nzf,nxc,nyc,nzc,pc(1,1),xf,yf,zf)
      elseif (mgprol .eq. 1) then
         numdia = ipc(11)
         if (numdia .eq. 7) then
!            print*,'% BUILDP:   (OPERATOR-7) building prolongation...'
            call buildP_op7 (nxf,nyf,nzf,nxc,nyc,nzc,
     2         ipc,rpc,ac(1,1),pc(1,1))
         elseif (numdia .eq. 27) then
!           print*,'% BUILDP:   (OPERATOR-27) building prolongation...'
            call buildP_op27 (nxf,nyf,nzf,nxc,nyc,nzc,
     2         ipc,rpc,ac(1,1),pc(1,1))
         else
            print*,'% BUILDP: invalid stencil type give...'
         endif
      elseif (mgprol .eq. 2) then
         numdia = ipc(11)
         if (numdia .eq. 7) then
!            print*,'% BUILDP: (MODIFY-OP-7) building prolongation...'
            call buildP_modop7 (nxf,nyf,nzf,nxc,nyc,nzc,
     2         ipc,rpc,ac(1,1),pc(1,1))
         elseif (numdia .eq. 27) then
!            print*,'% BUILDP: (MODIFY-OP-27) building prolongation...'
            call buildP_modop27 (nxf,nyf,nzf,nxc,nyc,nzc,
     2         ipc,rpc,ac(1,1),pc(1,1))
         else
            print*,'% BUILDP: invalid stencil type give...'
         endif
      else
         print*,'% BUILDP: invalid prolongation requested...'
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine buildP_trilin (nxf,nyf,nzf,nxc,nyc,nzc,pc,xf,yf,zf)
c* *********************************************************************
c* purpose: 
c*
c*    call the routine to form the trilinear prolongation operator.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nxc,nyc,nzc
      double precision pc(nxc*nyc*nzc,*),xf(*),yf(*),zf(*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      call buildPb_trilin (nxf,nyf,nzf,nxc,nyc,nzc,

     2   pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),pc(1,7),
     3   pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),pc(1,13),pc(1,14),
     4   pc(1,15),pc(1,16),pc(1,17),pc(1,18),pc(1,19),pc(1,20),pc(1,21),
     5   pc(1,22),pc(1,23),pc(1,24),pc(1,25),pc(1,26),pc(1,27),

     6   xf,yf,zf)
c*
c*    *** return and end ***
      return
      end
      subroutine buildPb_trilin (nxf,nyf,nzf,nxc,nyc,nzc,
     4   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     5   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     6   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW,
     7   xf,yf,zf)
c* *********************************************************************
c* purpose:  
c*
c*    form the trilinear prolongation operator.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nxc,nyc,nzc,i,j,k

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

      double precision xf(*),yf(*),zf(*)
      double precision won,half,quarter,eighth
c*
c*    *** interpolation stencil ***
      won     = 1.0d0
      half    = 1.0d0 / 2.0d0
      quarter = 1.0d0 / 4.0d0
      eighth  = 1.0d0 / 8.0d0
cmdir 3 1
      do 10 k = 1, nzc
cmdir 3 2
         do 11 j = 1, nyc
cmdir 3 3
            do 12 i = 1, nxc

               oPC(i,j,k)  = won

               oPN(i,j,k)  = half
               oPS(i,j,k)  = half
               oPE(i,j,k)  = half
               oPW(i,j,k)  = half
               uPC(i,j,k)  = half
               dPC(i,j,k)  = half

               oPNE(i,j,k) = quarter
               oPNW(i,j,k) = quarter
               oPSE(i,j,k) = quarter
               oPSW(i,j,k) = quarter
               dPN(i,j,k)  = quarter
               dPS(i,j,k)  = quarter
               dPE(i,j,k)  = quarter
               dPW(i,j,k)  = quarter
               uPN(i,j,k)  = quarter
               uPS(i,j,k)  = quarter
               uPE(i,j,k)  = quarter
               uPW(i,j,k)  = quarter

               dPNE(i,j,k) = eighth
               dPNW(i,j,k) = eighth
               dPSE(i,j,k) = eighth
               dPSW(i,j,k) = eighth
               uPNE(i,j,k) = eighth
               uPNW(i,j,k) = eighth
               uPSE(i,j,k) = eighth
               uPSW(i,j,k) = eighth

 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildP_op7 (nxf,nyf,nzf,nxc,nyc,nzc,ipc,rpc,ac,pc)
c* *********************************************************************
c* purpose: 
c*
c*    standard 7-pt operator-based prologation.
c*
c*    call the routine to form the prolongation operator from a
c*    7 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc
      double precision rpc(*),ac(nxf*nyf*nzf,*),pc(nxc*nyc*nzc,*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      call buildPb_op7 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,ac(1,1),

     3   ac(1,2),ac(1,3),ac(1,4),

     5   pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),pc(1,7),
     6   pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),pc(1,13),pc(1,14),
     7   pc(1,15),pc(1,16),pc(1,17),pc(1,18),pc(1,19),pc(1,20),pc(1,21),
     8   pc(1,22),pc(1,23),pc(1,24),pc(1,25),pc(1,26),pc(1,27))
c*
c*    *** return and end ***
      return
      end
      subroutine buildPb_op7 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,oC,
     3   oE,oN,uC,
     4   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     5   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     6   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW)
c* *********************************************************************
c* purpose:  
c*
c*    standard 7-pt operator-based prologation.
c*
c*    form the prolongation operator from a 7 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc,i,j,k,ii,jj,kk
      integer          im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
      integer          iim1,iip1,jjm1,jjp1,kkm1,kkp1
      double precision rpc(*),oC(nxf,nyf,nzf)

      double precision oE(nxf,nyf,nzf),oN(nxf,nyf,nzf),uC(nxf,nyf,nzf)

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

      double precision won,half,quarter,eighth
c*
c*    *** interpolation stencil ***
      won     = 1.0d0
      half    = 1.0d0 / 2.0d0
      quarter = 1.0d0 / 4.0d0
      eighth  = 1.0d0 / 8.0d0
      do 10 kk = 2, nzc-1
         k = 2 * kk - 1
         do 11 jj = 2, nyc-1
            j = 2 * jj - 1
            do 12 ii = 2, nxc-1
               i = 2 * ii - 1
c*
c*             *** index computations ***
               im1 = i-1
               ip1 = i+1
               im2 = i-2
               ip2 = i+2
               jm1 = j-1
               jp1 = j+1
               jm2 = j-2
               jp2 = j+2
               km1 = k-1
               kp1 = k+1
               km2 = k-2
               kp2 = k+2
               iim1 = ii-1
               iip1 = ii+1
               jjm1 = jj-1
               jjp1 = jj+1
               kkm1 = kk-1
               kkp1 = kk+1
c*
c* ********************************************************************
c* *** > oPC;
c* ********************************************************************

      oPC(ii,jj,kk) = won

c* ********************************************************************
c* *** > oPN;
c* ********************************************************************

      oPN(ii,jj,kk) = 
     2  oN(i,j,k)
     2  /  ( oC(i,jp1,k) 
     2     - oE(im1,jp1,k) 
     2     - oE(i,jp1,k) 
     2     - uC(i,jp1,km1) 
     2     - uC(i,jp1,k))

c* ********************************************************************
c* *** > oPS;
c* ********************************************************************

      oPS(ii,jj,kk) = 
     2  oN(i,jm1,k)
     2  /  ( oC(i,jm1,k) 
     2     - oE(im1,jm1,k) 
     2     - oE(i,jm1,k) 
     2     - uC(i,jm1,km1)
     2     - uC(i,jm1,k))

c* ********************************************************************
c* *** > oPE;
c* ********************************************************************

      oPE(ii,jj,kk) = 
     2  oE(i,j,k)
     2  /  ( oC(ip1,j,k) 
     2     - uC(ip1,j,km1) 
     2     - uC(ip1,j,k) 
     2     - oN(ip1,j,k)
     2     - oN(ip1,jm1,k))

c* ********************************************************************
c* *** > oPW;
c* ********************************************************************

      oPW(ii,jj,kk) = 
     2  oE(im1,j,k)
     2  /  ( oC(im1,j,k) 
     2     - uC(im1,j,km1) 
     2     - uC(im1,j,k) 
     2     - oN(im1,j,k)
     2     - oN(im1,jm1,k))

c* ********************************************************************
c* *** > oPNE;
c* ********************************************************************

      oPNE(ii,jj,kk) = 
     2   (oN(ip1,j,k) * oPE(ii,jj,kk)
     2  + oE(i,jp1,k) * oPN(ii,jj,kk))
     2  /  ( oC(ip1,jp1,k) 
     2     - uC(ip1,jp1,km1) 
     2     - uC(ip1,jp1,k))

c* ********************************************************************
c* *** > oPNW;
c* ********************************************************************

      oPNW(ii,jj,kk) = 
     2   (oN(im1,j,k) * oPW(ii,jj,kk) 
     2  + oE(im1,jp1,k) * oPN(ii,jj,kk))
     2  /  ( oC(im1,jp1,k) 
     2     - uC(im1,jp1,km1) 
     2     - uC(im1,jp1,k))

c* ********************************************************************
c* *** > oPSE;
c* ********************************************************************

      oPSE(ii,jj,kk) = 
     2    (oN(ip1,jm1,k) * oPE(ii,jj,kk) 
     2  +  oE(i,jm1,k) * oPS(ii,jj,kk))
     2  /  ( oC(ip1,jm1,k) 
     2     - uC(ip1,jm1,km1) 
     2     - uC(ip1,jm1,k))

c* ********************************************************************
c* *** > oPSW;
c* ********************************************************************

      oPSW(ii,jj,kk) = 
     2    (oN(im1,jm1,k) * oPW(ii,jj,kk)
     2  +  oE(im1,jm1,k) * oPS(ii,jj,kk))
     2  /  ( oC(im1,jm1,k) 
     2     - uC(im1,jm1,km1) 
     2     - uC(im1,jm1,k))

c* ********************************************************************
c* *** > dPC;
c* ********************************************************************

      dPC(ii,jj,kk) = 
     2  uC(i,j,km1)
     2  /  ( oC(i,j,km1) 
     2     - oN(i,j,km1) 
     2     - oN(i,jm1,km1) 
     2     - oE(im1,j,km1)
     2     - oE(i,j,km1))

c* ********************************************************************
c* *** > dPN;
c* ********************************************************************

      dPN(ii,jj,kk) = 
     2    (oN(i,j,km1) * dPC(ii,jj,kk) 
     2  +  uC(i,jp1,km1) * oPN(ii,jj,kk))
     2  /  ( oC(i,jp1,km1) 
     2     - oE(im1,jp1,km1) 
     2     - oE(i,jp1,km1))

c* ********************************************************************
c* *** > dPS;
c* ********************************************************************

      dPS(ii,jj,kk) = 
     2    (oN(i,jm1,km1) * dPC(ii,jj,kk)
     2  +  uC(i,jm1,km1) * oPS(ii,jj,kk))
     2  /  ( oC(i,jm1,km1) 
     2     - oE(im1,jm1,km1) 
     2     - oE(i,jm1,km1))

c* ********************************************************************
c* *** > dPE;
c* ********************************************************************

      dPE(ii,jj,kk) = 
     2   (uC(ip1,j,km1) * oPE(ii,jj,kk)
     2 +  oE(i,j,km1) * dPC(ii,jj,kk))
     2 /  ( oC(ip1,j,km1) 
     2    - oN(ip1,j,km1) 
     2    - oN(ip1,jm1,km1))

c* ********************************************************************
c* *** > dPW;
c* ********************************************************************

      dPW(ii,jj,kk) = 
     2    (uC(im1,j,km1) * oPW(ii,jj,kk)
     2  +  oE(im1,j,km1) * dPC(ii,jj,kk))
     2  /  ( oC(im1,j,km1) 
     2     - oN(im1,j,km1) 
     2     - oN(im1,jm1,km1))

c* ********************************************************************
c* *** > dPNE;
c* ********************************************************************

      dPNE(ii,jj,kk) = 
     2    (uC(ip1,jp1,km1) * oPNE(ii,jj,kk)
     2   + oE(i,jp1,km1) * dPN(ii,jj,kk)
     2   + oN(ip1,j,km1) * dPE(ii,jj,kk))
     2   /  oC(ip1,jp1,km1)

c* ********************************************************************
c* *** > dPNW;
c* ********************************************************************

      dPNW(ii,jj,kk) = 
     2   (uC(im1,jp1,km1) * oPNW(ii,jj,kk)
     2  + oE(im1,jp1,km1) * dPN(ii,jj,kk)
     2  + oN(im1,j,km1) * dPW(ii,jj,kk))
     2  /  oC(im1,jp1,km1)

c* ********************************************************************
c* *** > dPSE;
c* ********************************************************************

      dPSE(ii,jj,kk) = 
     2   (uC(ip1,jm1,km1) * oPSE(ii,jj,kk)
     2  + oE(i,jm1,km1) * dPS(ii,jj,kk)
     2  + oN(ip1,jm1,km1) * dPE(ii,jj,kk))
     2  /  oC(ip1,jm1,km1)

c* ********************************************************************
c* *** > dPSW;
c* ********************************************************************

      dPSW(ii,jj,kk) = 
     2    (uC(im1,jm1,km1) * oPSW(ii,jj,kk)
     2   + oE(im1,jm1,km1) * dPS(ii,jj,kk)
     2   + oN(im1,jm1,km1) * dPW(ii,jj,kk))
     2   /  oC(im1,jm1,km1)

c* ********************************************************************
c* *** > uPC;
c* ********************************************************************

      uPC(ii,jj,kk) = 
     2  uC(i,j,k)
     2  /  ( oC(i,j,kp1) 
     2     - oN(i,j,kp1) 
     2     - oN(i,jm1,kp1) 
     2     - oE(im1,j,kp1)
     2     - oE(i,j,kp1))

c* ********************************************************************
c* *** > uPN;
c* ********************************************************************

      uPN(ii,jj,kk) = 
     2   (oN(i,j,kp1) * uPC(ii,jj,kk) 
     2 +  uC(i,jp1,k) * oPN(ii,jj,kk))
     2 /  ( oC(i,jp1,kp1) 
     2    - oE(im1,jp1,kp1) 
     2    - oE(i,jp1,kp1))

c* ********************************************************************
c* *** > uPS;
c* ********************************************************************

      uPS(ii,jj,kk) = 
     2   (oN(i,jm1,kp1) * uPC(ii,jj,kk)
     2 +  uC(i,jm1,k) * oPS(ii,jj,kk))
     2 /  ( oC(i,jm1,kp1) 
     2    - oE(im1,jm1,kp1) 
     2    - oE(i,jm1,kp1))

c* ********************************************************************
c* *** > uPE;
c* ********************************************************************

      uPE(ii,jj,kk) = 
     2   (uC(ip1,j,k) * oPE(ii,jj,kk)
     2 +  oE(i,j,kp1) * uPC(ii,jj,kk))
     2 /  ( oC(ip1,j,kp1) 
     2    - oN(ip1,j,kp1) 
     2    - oN(ip1,jm1,kp1))

c* ********************************************************************
c* *** > uPW;
c* ********************************************************************

      uPW(ii,jj,kk) = 
     2    (uC(im1,j,k) * oPW(ii,jj,kk)
     2 +  oE(im1,j,kp1) * uPC(ii,jj,kk))
     2 /  ( oC(im1,j,kp1) 
     2    - oN(im1,j,kp1) 
     2    - oN(im1,jm1,kp1))

c* ********************************************************************
c* *** > uPNE;
c* ********************************************************************

      uPNE(ii,jj,kk) = 
     2    (uC(ip1,jp1,k) * oPNE(ii,jj,kk)
     2  +  oE(i,jp1,kp1) * uPN(ii,jj,kk)
     2  +  oN(ip1,j,kp1) * uPE(ii,jj,kk))
     2  /  oC(ip1,jp1,kp1)

c* ********************************************************************
c* *** > uPNW;
c* ********************************************************************

      uPNW(ii,jj,kk) = 
     2   (uC(im1,jp1,k) * oPNW(ii,jj,kk)
     2  + oE(im1,jp1,kp1) * uPN(ii,jj,kk)
     2  + oN(im1,j,kp1) * uPW(ii,jj,kk))
     2  /  oC(im1,jp1,kp1)

c* ********************************************************************
c* *** > uPSE;
c* ********************************************************************

      uPSE(ii,jj,kk) = 
     2   (uC(ip1,jm1,k) * oPSE(ii,jj,kk)
     2  + oE(i,jm1,kp1) * uPS(ii,jj,kk)
     2  + oN(ip1,jm1,kp1) * uPE(ii,jj,kk))
     2  /  oC(ip1,jm1,kp1)

c* ********************************************************************
c* *** > uPSW;
c* ********************************************************************

      uPSW(ii,jj,kk) = 
     2   (uC(im1,jm1,k) * oPSW(ii,jj,kk)
     2  + oE(im1,jm1,kp1) * uPS(ii,jj,kk)
     2  + oN(im1,jm1,kp1) * uPW(ii,jj,kk))
     2  /  oC(im1,jm1,kp1)

c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildP_op27 (nxf,nyf,nzf,nxc,nyc,nzc,ipc,rpc,ac,pc)
c* *********************************************************************
c* purpose: 
c*
c*    standard 27-pt operator-based prologation.
c*
c*    call the routine to form the prolongation operator from a 
c*    27 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc
      double precision rpc(*),ac(nxf*nyf*nzf,*),pc(nxc*nyc*nzc,*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      call buildPb_op27 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,ac(1,1),

     3   ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     4   ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),

     5   pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),pc(1,7),
     6   pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),pc(1,13),pc(1,14),
     7   pc(1,15),pc(1,16),pc(1,17),pc(1,18),pc(1,19),pc(1,20),pc(1,21),
     8   pc(1,22),pc(1,23),pc(1,24),pc(1,25),pc(1,26),pc(1,27))
c*
c*    *** return and end ***
      return
      end
      subroutine buildPb_op27 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,oC,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     5   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     6   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW)
c* *********************************************************************
c* purpose:  
c*
c*    standard 27-pt operator-based prologation.
c*
c*    form the prolongation operator from a 27 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc,i,j,k,ii,jj,kk
      integer          im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
      integer          iim1,iip1,jjm1,jjp1,kkm1,kkp1
      double precision rpc(*),oC(nxf,nyf,nzf)

      double precision oE(nxf,nyf,nzf),oN(nxf,nyf,nzf)
      double precision uC(nxf,nyf,nzf),oNE(nxf,nyf,nzf)
      double precision oNW(nxf,nyf,nzf),uE(nxf,nyf,nzf)
      double precision uW(nxf,nyf,nzf),uN(nxf,nyf,nzf)
      double precision uS(nxf,nyf,nzf),uNE(nxf,nyf,nzf)
      double precision uNW(nxf,nyf,nzf),uSE(nxf,nyf,nzf)
      double precision uSW(nxf,nyf,nzf)

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

      double precision won,half,quarter,eighth
c*
c*    *** interpolation stencil ***
      won     = 1.0d0
      half    = 1.0d0 / 2.0d0
      quarter = 1.0d0 / 4.0d0
      eighth  = 1.0d0 / 8.0d0
      do 10 kk = 2, nzc-1
         k = 2 * kk - 1
         do 11 jj = 2, nyc-1
            j = 2 * jj - 1
            do 12 ii = 2, nxc-1
               i = 2 * ii - 1
c*
c*             *** index computations ***
               im1 = i-1
               ip1 = i+1
               im2 = i-2
               ip2 = i+2
               jm1 = j-1
               jp1 = j+1
               jm2 = j-2
               jp2 = j+2
               km1 = k-1
               kp1 = k+1
               km2 = k-2
               kp2 = k+2
               iim1 = ii-1
               iip1 = ii+1
               jjm1 = jj-1
               jjp1 = jj+1
               kkm1 = kk-1
               kkp1 = kk+1
c*
c* ********************************************************************
c* *** > oPC;
c* ********************************************************************

      oPC(ii,jj,kk) = won

c* ********************************************************************
c* *** > oPN;
c* ********************************************************************

      oPN(ii,jj,kk) = 
     2   (uNE(im1,j,km1) + uN(i,j,km1) + uNW(ip1,j,km1)
     2  + oNE(im1,j,k) + oN(i,j,k) + oNW(ip1,j,k) 
     2  + uSW(i,jp1,k) + uS(i,jp1,k) + uSE(i,jp1,k))  
     2  /  ( oC(i,jp1,k) 
     2     - oE(im1,jp1,k) - oE(i,jp1,k) 
     2     - uC(i,jp1,km1) - uE(im1,jp1,km1) - uW(ip1,jp1,km1) 
     2     - uC(i,jp1,k) - uW(i,jp1,k) - uE(i,jp1,k))

c* ********************************************************************
c* *** > oPS;
c* ********************************************************************

      oPS(ii,jj,kk) = 
     2   (uSE(im1,j,km1) + uS(i,j,km1) + uSW(ip1,j,km1)
     2  + oNW(i,jm1,k) + oN(i,jm1,k) + oNE(i,jm1,k) 
     2  + uNW(i,jm1,k) + uN(i,jm1,k) + uNE(i,jm1,k))  
     2  /  ( oC(i,jm1,k) 
     2     - oE(im1,jm1,k) - oE(i,jm1,k) 
     2     - uC(i,jm1,km1) - uE(im1,jm1,km1) - uW(ip1,jm1,km1) 
     2     - uC(i,jm1,k) - uW(i,jm1,k) - uE(i,jm1,k))

c* ********************************************************************
c* *** > oPE;
c* ********************************************************************

      oPE(ii,jj,kk) = 
     2   (uSE(i,jp1,km1) + oNW(ip1,j,k) + uNW(ip1,j,k)
     2  + uE(i,j,km1) + oE(i,j,k) + uW(ip1,j,k) 
     2  + uNE(i,jm1,km1) + oNE(i,jm1,k) + uSW(ip1,j,k))  
     2  /  ( oC(ip1,j,k) 
     2     - uC(ip1,j,km1) - uC(ip1,j,k) 
     2     - oN(ip1,j,k) - uS(ip1,jp1,km1) - uN(ip1,j,k) 
     2     - oN(ip1,jm1,k) - uN(ip1,jm1,km1) - uS(ip1,j,k))

c* ********************************************************************
c* *** > oPW;
c* ********************************************************************

      oPW(ii,jj,kk) = 
     2   (uSW(i,jp1,km1) + oNE(im1,j,k) + uNE(im1,j,k)
     2  + uW(i,j,km1) + oE(im1,j,k) + uE(im1,j,k)
     2  + uNW(i,jm1,km1) + oNW(i,jm1,k) + uSE(im1,j,k))  
     2  /  ( oC(im1,j,k) 
     2     - uC(im1,j,km1) - uC(im1,j,k) 
     2     - oN(im1,j,k) - uS(im1,jp1,km1) - uN(im1,j,k) 
     2     - oN(im1,jm1,k) - uN(im1,jm1,km1) - uS(im1,j,k))

c* ********************************************************************
c* *** > oPNE;
c* ********************************************************************

      oPNE(ii,jj,kk) = 
     2   (uNE(i,j,km1) + oNE(i,j,k) + uSW(ip1,jp1,k) 
     2 + (uN(ip1,j,km1) + oN(ip1,j,k) + uS(ip1,jp1,k))*oPE(ii,jj,kk)
     2 + (uE(i,jp1,km1) + oE(i,jp1,k) + uW(ip1,jp1,k))*oPN(ii,jj,kk))
     2 /  ( oC(ip1,jp1,k) 
     2    - uC(ip1,jp1,km1) - uC(ip1,jp1,k))

c* ********************************************************************
c* *** > oPNW;
c* ********************************************************************

      oPNW(ii,jj,kk) = 
     2   (uNW(i,j,km1) + oNW(i,j,k) + uSE(im1,jp1,k) 
     2 + (uN(im1,j,km1) + oN(im1,j,k) + uS(im1,jp1,k))*oPW(ii,jj,kk) 
     2 + (uW(i,jp1,km1) + oE(im1,jp1,k) + uE(im1,jp1,k))*oPN(ii,jj,kk))
     2 /  ( oC(im1,jp1,k) 
     2    - uC(im1,jp1,km1) - uC(im1,jp1,k))

c* ********************************************************************
c* *** > oPSE;
c* ********************************************************************

      oPSE(ii,jj,kk) = 
     2   (uSE(i,j,km1) + oNW(ip1,jm1,k) + uNW(ip1,jm1,k)
     2 + (uS(ip1,j,km1) + oN(ip1,jm1,k) + uN(ip1,jm1,k))*oPE(ii,jj,kk) 
     2 + (uE(i,jm1,km1) + oE(i,jm1,k) + uW(ip1,jm1,k))*oPS(ii,jj,kk))
     2 /  ( oC(ip1,jm1,k) 
     2    - uC(ip1,jm1,km1) - uC(ip1,jm1,k))

c* ********************************************************************
c* *** > oPSW;
c* ********************************************************************

      oPSW(ii,jj,kk) = 
     2   (uSW(i,j,km1) + oNE(im1,jm1,k) + uNE(im1,jm1,k)
     2 + (uS(im1,j,km1) + oN(im1,jm1,k) + uN(im1,jm1,k))*oPW(ii,jj,kk)
     2 + (uW(i,jm1,km1) + oE(im1,jm1,k) + uE(im1,jm1,k))*oPS(ii,jj,kk))
     2 /  ( oC(im1,jm1,k) 
     2    - uC(im1,jm1,km1) - uC(im1,jm1,k))

c* ********************************************************************
c* *** > dPC;
c* ********************************************************************

      dPC(ii,jj,kk) = 
     2   (uNW(i,j,km1) + uW(i,j,km1) + uSW(i,j,km1)
     2  + uN(i,j,km1) + uC(i,j,km1) + uS(i,j,km1) + uNE(i,j,km1)
     2  + uE(i,j,km1) + uSE(i,j,km1))  
     2  /  ( oC(i,j,km1) 
     2     - oN(i,j,km1) - oN(i,jm1,km1) 
     2     - oNW(i,j,km1) - oE(im1,j,km1) - oNE(im1,jm1,km1) 
     2     - oNE(i,j,km1) - oE(i,j,km1) - oNW(ip1,jm1,km1))

c* ********************************************************************
c* *** > dPN;
c* ********************************************************************

      dPN(ii,jj,kk) = 
     2   (uSW(i,jp1,km1) + uS(i,jp1,km1) + uSE(i,jp1,km1) 
     2 + (oNE(im1,j,km1) + oN(i,j,km1) + oNW(ip1,j,km1))*dPC(ii,jj,kk) 
     2 + (uW(i,jp1,km1) + uC(i,jp1,km1) + uE(i,jp1,km1))*oPN(ii,jj,kk))
     2 /  ( oC(i,jp1,km1) 
     2    - oE(im1,jp1,km1) - oE(i,jp1,km1))

c* ********************************************************************
c* *** > dPS;
c* ********************************************************************

      dPS(ii,jj,kk) = 
     2   (uNW(i,jm1,km1) + uN(i,jm1,km1) + uNE(i,jm1,km1) 
     2 + (oNW(i,jm1,km1) + oN(i,jm1,km1) + oNE(i,jm1,km1))*dPC(ii,jj,kk)
     2 + (uW(i,jm1,km1) + uC(i,jm1,km1) + uE(i,jm1,km1))*oPS(ii,jj,kk))
     2 /  ( oC(i,jm1,km1) 
     2    - oE(im1,jm1,km1) - oE(i,jm1,km1))

c* ********************************************************************
c* *** > dPE;
c* ********************************************************************

      dPE(ii,jj,kk) = 
     2   (uNW(ip1,j,km1) + uW(ip1,j,km1) + uSW(ip1,j,km1) 
     2 + (uN(ip1,j,km1) + uC(ip1,j,km1) + uS(ip1,j,km1))*oPE(ii,jj,kk)
     2 + (oNW(ip1,j,km1) + oE(i,j,km1) + oNE(i,jm1,km1))*dPC(ii,jj,kk))
     2 /  ( oC(ip1,j,km1) 
     2    - oN(ip1,j,km1) - oN(ip1,jm1,km1))

c* ********************************************************************
c* *** > dPW;
c* ********************************************************************

      dPW(ii,jj,kk) = 
     2  (uNE(im1,j,km1) + uE(im1,j,km1) + uSE(im1,j,km1) 
     2 +(uN(im1,j,km1) + uC(im1,j,km1) + uS(im1,j,km1))*oPW(ii,jj,kk)
     2 +(oNE(im1,j,km1) + oE(im1,j,km1) + oNW(i,jm1,km1))*dPC(ii,jj,kk))
     2 /  ( oC(im1,j,km1) 
     2    - oN(im1,j,km1) - oN(im1,jm1,km1))

c* ********************************************************************
c* *** > dPNE;
c* ********************************************************************

      dPNE(ii,jj,kk) = 
     2   (uSW(ip1,jp1,km1) 
     2  + uW(ip1,jp1,km1) * oPN(ii,jj,kk)
     2  + uS(ip1,jp1,km1) * oPE(ii,jj,kk) 
     2  + uC(ip1,jp1,km1) * oPNE(ii,jj,kk)
     2  + oNE(i,j,km1) * dPC(ii,jj,kk)
     2  + oE(i,jp1,km1) * dPN(ii,jj,kk)
     2  + oN(ip1,j,km1) * dPE(ii,jj,kk))
     2  /  oC(ip1,jp1,km1)

c* ********************************************************************
c* *** > dPNW;
c* ********************************************************************

      dPNW(ii,jj,kk) = 
     2   (uSE(im1,jp1,km1) 
     2  + uE(im1,jp1,km1) * oPN(ii,jj,kk)
     2  + uS(im1,jp1,km1) * oPW(ii,jj,kk)
     2  + uC(im1,jp1,km1) * oPNW(ii,jj,kk)
     2  + oNW(i,j,km1) * dPC(ii,jj,kk)
     2  + oE(im1,jp1,km1) * dPN(ii,jj,kk)
     2  + oN(im1,j,km1) * dPW(ii,jj,kk))
     2  /  oC(im1,jp1,km1)

c* ********************************************************************
c* *** > dPSE;
c* ********************************************************************

      dPSE(ii,jj,kk) = 
     2    (uNW(ip1,jm1,km1)
     2  + uW(ip1,jm1,km1) * oPS(ii,jj,kk)
     2  + uN(ip1,jm1,km1) * oPE(ii,jj,kk)
     2  + uC(ip1,jm1,km1) * oPSE(ii,jj,kk)
     2  + oNW(ip1,jm1,km1) * dPC(ii,jj,kk)
     2  + oE(i,jm1,km1) * dPS(ii,jj,kk)
     2  + oN(ip1,jm1,km1) * dPE(ii,jj,kk))
     2  /  oC(ip1,jm1,km1)

c* ********************************************************************
c* *** > dPSW;
c* ********************************************************************

      dPSW(ii,jj,kk) = 
     2    (uNE(im1,jm1,km1)
     2   + uE(im1,jm1,km1) * oPS(ii,jj,kk)
     2   + uN(im1,jm1,km1) * oPW(ii,jj,kk)
     2   + uC(im1,jm1,km1) * oPSW(ii,jj,kk)
     2   + oNE(im1,jm1,km1) * dPC(ii,jj,kk)
     2   + oE(im1,jm1,km1) * dPS(ii,jj,kk)
     2   + oN(im1,jm1,km1) * dPW(ii,jj,kk))
     2   /  oC(im1,jm1,km1)

c* ********************************************************************
c* *** > uPC;
c* ********************************************************************

      uPC(ii,jj,kk) = 
     2   (uSE(im1,jp1,k) + uE(im1,j,k) + uNE(im1,jm1,k)
     2  + uS(i,jp1,k) + uC(i,j,k) + uN(i,jm1,k) 
     2  + uSW(ip1,jp1,k) + uW(ip1,j,k) + uNW(ip1,jm1,k))  
     2  /  ( oC(i,j,kp1) 
     2     - oN(i,j,kp1) - oN(i,jm1,kp1) 
     2     - oNW(i,j,kp1) - oE(im1,j,kp1) - oNE(im1,jm1,kp1) 
     2     - oNE(i,j,kp1) - oE(i,j,kp1) - oNW(ip1,jm1,kp1))

c* ********************************************************************
c* *** > uPN;
c* ********************************************************************

      uPN(ii,jj,kk) = 
     2   (uNE(im1,j,k) + uN(i,j,k) + uNW(ip1,j,k) 
     2 + (oNE(im1,j,kp1) + oN(i,j,kp1) + oNW(ip1,j,kp1))*uPC(ii,jj,kk) 
     2 + (uE(im1,jp1,k) + uC(i,jp1,k) + uW(ip1,jp1,k))*oPN(ii,jj,kk))
     2 /  ( oC(i,jp1,kp1) 
     2    - oE(im1,jp1,kp1) - oE(i,jp1,kp1))

c* ********************************************************************
c* *** > uPS;
c* ********************************************************************

      uPS(ii,jj,kk) = 
     2   (uSE(im1,j,k) + uS(i,j,k) + uSW(ip1,j,k) 
     2 + (oNW(i,jm1,kp1) + oN(i,jm1,kp1) + oNE(i,jm1,kp1))*uPC(ii,jj,kk)
     2 + (uE(im1,jm1,k) + uC(i,jm1,k) + uW(ip1,jm1,k))*oPS(ii,jj,kk))
     2 /  ( oC(i,jm1,kp1) 
     2    - oE(im1,jm1,kp1) - oE(i,jm1,kp1))

c* ********************************************************************
c* *** > uPE;
c* ********************************************************************

      uPE(ii,jj,kk) = 
     2   (uSE(i,jp1,k) + uS(ip1,jp1,k) + uNE(i,jm1,k) 
     2 + (uS(ip1,jp1,k) + uC(ip1,j,k) + uN(ip1,jm1,k))*oPE(ii,jj,kk)
     2 + (oNW(ip1,j,kp1) + oE(i,j,kp1) + oNE(i,jm1,kp1))*uPC(ii,jj,kk))
     2 /  ( oC(ip1,j,kp1) 
     2    - oN(ip1,j,kp1) - oN(ip1,jm1,kp1))

c* ********************************************************************
c* *** > uPW;
c* ********************************************************************

      uPW(ii,jj,kk) = 
     2  (uSW(i,jp1,k) + uW(i,j,k) + uNW(i,jm1,k) 
     2 +(uS(im1,jp1,k) + uC(im1,j,k) + uN(im1,jm1,k))*oPW(ii,jj,kk)
     2 +(oNE(im1,j,kp1) + oE(im1,j,kp1) + oNW(i,jm1,kp1))*uPC(ii,jj,kk))
     2 /  ( oC(im1,j,kp1) 
     2    - oN(im1,j,kp1) - oN(im1,jm1,kp1))

c* ********************************************************************
c* *** > uPNE;
c* ********************************************************************

      uPNE(ii,jj,kk) = 
     2   (uNE(i,j,k) 
     2  + uE(i,jp1,k) * oPN(ii,jj,kk)
     2  + uN(ip1,j,k) * oPE(ii,jj,kk)
     2  + uC(ip1,jp1,k) * oPNE(ii,jj,kk)
     2  + oNE(i,j,kp1) * uPC(ii,jj,kk)
     2  + oE(i,jp1,kp1) * uPN(ii,jj,kk)
     2  + oN(ip1,j,kp1) * uPE(ii,jj,kk))
     2  /  oC(ip1,jp1,kp1)

c* ********************************************************************
c* *** > uPNW;
c* ********************************************************************

      uPNW(ii,jj,kk) = 
     2   (uNW(i,j,k) 
     2  + uW(i,jp1,k) * oPN(ii,jj,kk)
     2  + uN(im1,j,k) * oPW(ii,jj,kk)
     2  + uC(im1,jp1,k) * oPNW(ii,jj,kk)
     2  + oNW(i,j,kp1) * uPC(ii,jj,kk)
     2  + oE(im1,jp1,kp1) * uPN(ii,jj,kk)
     2  + oN(im1,j,kp1) * uPW(ii,jj,kk))
     2  /  oC(im1,jp1,kp1)

c* ********************************************************************
c* *** > uPSE;
c* ********************************************************************

      uPSE(ii,jj,kk) = 
     2   (uSE(i,j,k) 
     2  + uE(i,jm1,k) * oPS(ii,jj,kk)
     2  + uS(ip1,j,k) * oPE(ii,jj,kk)
     2  + uC(ip1,jm1,k) * oPSE(ii,jj,kk)
     2  + oNW(ip1,jm1,kp1) * uPC(ii,jj,kk)
     2  + oE(i,jm1,kp1) * uPS(ii,jj,kk)
     2  + oN(ip1,jm1,kp1) * uPE(ii,jj,kk))
     2  /  oC(ip1,jm1,kp1)

c* ********************************************************************
c* *** > uPSW;
c* ********************************************************************

      uPSW(ii,jj,kk) = 
     2   (uSW(i,j,k) 
     2  + uW(i,jm1,k) * oPS(ii,jj,kk)
     2  + uS(im1,j,k) * oPW(ii,jj,kk)
     2  + uC(im1,jm1,k) * oPSW(ii,jj,kk)
     2  + oNE(im1,jm1,kp1) * uPC(ii,jj,kk)
     2  + oE(im1,jm1,kp1) * uPS(ii,jj,kk)
     2  + oN(im1,jm1,kp1) * uPW(ii,jj,kk))
     2  /  oC(im1,jm1,kp1)

c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildP_modop7 (nxf,nyf,nzf,nxc,nyc,nzc,ipc,rpc,ac,pc)
c* *********************************************************************
c* purpose: 
c*
c*    a modified 7-pt operator-based prologation.
c*
c*    call the routine to form the prolongation operator from a
c*    7 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc
      double precision rpc(*),ac(nxf*nyf*nzf,*),pc(nxc*nyc*nzc,*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      call buildPb_modop7 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,ac(1,1),

     3   ac(1,2),ac(1,3),ac(1,4),

     5   pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),pc(1,7),
     6   pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),pc(1,13),pc(1,14),
     7   pc(1,15),pc(1,16),pc(1,17),pc(1,18),pc(1,19),pc(1,20),pc(1,21),
     8   pc(1,22),pc(1,23),pc(1,24),pc(1,25),pc(1,26),pc(1,27))
c*
c*    *** return and end ***
      return
      end
      subroutine buildPb_modop7 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,oC,
     3   oE,oN,uC,
     4   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     5   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     6   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW)
c* *********************************************************************
c* purpose:  
c*
c*    a modified 7-pt operator-based prologation.
c*
c*    form the prolongation operator from a 7 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc,i,j,k,ii,jj,kk
      integer          im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
      integer          iim1,iip1,jjm1,jjp1,kkm1,kkp1
      double precision rpc(*),oC(nxf,nyf,nzf)

      double precision oE(nxf,nyf,nzf),oN(nxf,nyf,nzf),uC(nxf,nyf,nzf)

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

      double precision won,half,quarter,eighth

      double precision oCijk
      double precision oCim1jk,oCip1jk
      double precision oCijm1k,oCijp1k
      double precision oCim1jm1k,oCip1jm1k
      double precision oCim1jp1k,oCip1jp1k
      double precision oCijkm1
      double precision oCim1jkm1,oCip1jkm1
      double precision oCijm1km1,oCijp1km1
      double precision oCim1jm1km1,oCip1jm1km1
      double precision oCim1jp1km1,oCip1jp1km1
      double precision oCijkp1
      double precision oCim1jkp1,oCip1jkp1
      double precision oCijm1kp1,oCijp1kp1
      double precision oCim1jm1kp1,oCip1jm1kp1
      double precision oCim1jp1kp1,oCip1jp1kp1
c*
c*    *** interpolation stencil ***
      won     = 1.0d0
      half    = 1.0d0 / 2.0d0
      quarter = 1.0d0 / 4.0d0
      eighth  = 1.0d0 / 8.0d0
      do 10 kk = 2, nzc-1
         k = 2 * kk - 1
         do 11 jj = 2, nyc-1
            j = 2 * jj - 1
            do 12 ii = 2, nxc-1
               i = 2 * ii - 1
c*
c*             *** index computations ***
               im1 = i-1
               ip1 = i+1
               im2 = i-2
               ip2 = i+2
               jm1 = j-1
               jp1 = j+1
               jm2 = j-2
               jp2 = j+2
               km1 = k-1
               kp1 = k+1
               km2 = k-2
               kp2 = k+2
               iim1 = ii-1
               iip1 = ii+1
               jjm1 = jj-1
               jjp1 = jj+1
               kkm1 = kk-1
               kkp1 = kk+1
c*
c*             *** fix the operator ***
               oCijk =
     2     + dabs(oE(im1,j,k))
     2     + dabs(oE(i,j,k))
     2     + dabs(oN(i,jm1,k))
     2     + dabs(oN(i,j,k))
     2     + dabs(uC(i,j,km1))
     2     + dabs(uC(i,j,k))
               oCijp1k =
     2     + dabs(oE(im1,jp1,k))
     2     + dabs(oE(i,jp1,k))
     2     + dabs(oN(i,j,k))
     2     + dabs(oN(i,jp1,k))
     2     + dabs(uC(i,jp1,km1))
     2     + dabs(uC(i,jp1,k))
               oCijm1k =
     2     + dabs(oE(im1,jm1,k))
     2     + dabs(oE(i,jm1,k))
     2     + dabs(oN(i,jm2,k))
     2     + dabs(oN(i,jm1,k))
     2     + dabs(uC(i,jm1,km1))
     2     + dabs(uC(i,jm1,k))
               oCip1jk =
     2     + dabs(oE(i,j,k))
     2     + dabs(oE(ip1,j,k))
     2     + dabs(oN(ip1,jm1,k))
     2     + dabs(oN(ip1,j,k))
     2     + dabs(uC(ip1,j,km1))
     2     + dabs(uC(ip1,j,k))
               oCim1jk =
     2     + dabs(oE(im2,j,k))
     2     + dabs(oE(im1,j,k))
     2     + dabs(oN(im1,jm1,k))
     2     + dabs(oN(im1,j,k))
     2     + dabs(uC(im1,j,km1))
     2     + dabs(uC(im1,j,k))
               oCip1jp1k =
     2     + dabs(oE(i,jp1,k))
     2     + dabs(oE(ip1,jp1,k))
     2     + dabs(oN(ip1,j,k))
     2     + dabs(oN(ip1,jp1,k))
     2     + dabs(uC(ip1,jp1,km1))
     2     + dabs(uC(ip1,jp1,k))
               oCim1jp1k =
     2     + dabs(oE(im1,jp1,k))
     2     + dabs(oE(i,jp1,k))
     2     + dabs(oN(i,j,k))
     2     + dabs(oN(i,jp1,k))
     2     + dabs(uC(i,jp1,km1))
     2     + dabs(uC(i,jp1,k))
               oCip1jm1k =
     2     + dabs(oE(i,jm1,k))
     2     + dabs(oE(ip1,jm1,k))
     2     + dabs(oN(ip1,jm2,k))
     2     + dabs(oN(ip1,jm1,k))
     2     + dabs(uC(ip1,jm1,km1))
     2     + dabs(uC(ip1,jm1,k))
               oCim1jm1k =
     2     + dabs(oE(im1,jm1,k))
     2     + dabs(oE(i,jm1,k))
     2     + dabs(oN(i,jm2,k))
     2     + dabs(oN(i,jm1,k))
     2     + dabs(uC(i,jm1,km1))
     2     + dabs(uC(i,jm1,k))

c* ********************************************************************

               oCijkm1 =
     2     + dabs(oE(im1,j,km1))
     2     + dabs(oE(i,j,km1))
     2     + dabs(oN(i,jm1,km1))
     2     + dabs(oN(i,j,km1))
     2     + dabs(uC(i,j,km2))
     2     + dabs(uC(i,j,km1))
               oCijp1km1 =
     2     + dabs(oE(im1,jp1,km1))
     2     + dabs(oE(i,jp1,km1))
     2     + dabs(oN(i,j,km1))
     2     + dabs(oN(i,jp1,km1))
     2     + dabs(uC(i,jp1,km2))
     2     + dabs(uC(i,jp1,km1))
               oCijm1km1 =
     2     + dabs(oE(im1,jm1,km1))
     2     + dabs(oE(i,jm1,km1))
     2     + dabs(oN(i,jm2,km1))
     2     + dabs(oN(i,jm1,km1))
     2     + dabs(uC(i,jm1,km2))
     2     + dabs(uC(i,jm1,km1))
               oCip1jkm1 =
     2     + dabs(oE(i,j,km1))
     2     + dabs(oE(ip1,j,km1))
     2     + dabs(oN(ip1,jm1,km1))
     2     + dabs(oN(ip1,j,km1))
     2     + dabs(uC(ip1,j,km2))
     2     + dabs(uC(ip1,j,km1))
               oCim1jkm1 =
     2     + dabs(oE(im2,j,km1))
     2     + dabs(oE(im1,j,km1))
     2     + dabs(oN(im1,jm1,km1))
     2     + dabs(oN(im1,j,km1))
     2     + dabs(uC(im1,j,km2))
     2     + dabs(uC(im1,j,km1))
               oCip1jp1km1 =
     2     + dabs(oE(i,jp1,km1))
     2     + dabs(oE(ip1,jp1,km1))
     2     + dabs(oN(ip1,j,km1))
     2     + dabs(oN(ip1,jp1,km1))
     2     + dabs(uC(ip1,jp1,km2))
     2     + dabs(uC(ip1,jp1,km1))
               oCim1jp1km1 =
     2     + dabs(oE(im1,jp1,km1))
     2     + dabs(oE(i,jp1,km1))
     2     + dabs(oN(i,j,km1))
     2     + dabs(oN(i,jp1,km1))
     2     + dabs(uC(i,jp1,km2))
     2     + dabs(uC(i,jp1,km1))
               oCip1jm1km1 =
     2     + dabs(oE(i,jm1,km1))
     2     + dabs(oE(ip1,jm1,km1))
     2     + dabs(oN(ip1,jm2,km1))
     2     + dabs(oN(ip1,jm1,km1))
     2     + dabs(uC(ip1,jm1,km2))
     2     + dabs(uC(ip1,jm1,km1))
               oCim1jm1km1 =
     2     + dabs(oE(im1,jm1,km1))
     2     + dabs(oE(i,jm1,km1))
     2     + dabs(oN(i,jm2,km1))
     2     + dabs(oN(i,jm1,km1))
     2     + dabs(uC(i,jm1,km2))
     2     + dabs(uC(i,jm1,km1))

c* ********************************************************************

               oCijkp1 =
     2     + dabs(oE(im1,j,kp1))
     2     + dabs(oE(i,j,kp1))
     2     + dabs(oN(i,jm1,kp1))
     2     + dabs(oN(i,j,kp1))
     2     + dabs(uC(i,j,k))
     2     + dabs(uC(i,j,kp1))
               oCijp1kp1 =
     2     + dabs(oE(im1,jp1,kp1))
     2     + dabs(oE(i,jp1,kp1))
     2     + dabs(oN(i,j,kp1))
     2     + dabs(oN(i,jp1,kp1))
     2     + dabs(uC(i,jp1,k))
     2     + dabs(uC(i,jp1,kp1))
               oCijm1kp1 =
     2     + dabs(oE(im1,jm1,kp1))
     2     + dabs(oE(i,jm1,kp1))
     2     + dabs(oN(i,jm2,kp1))
     2     + dabs(oN(i,jm1,kp1))
     2     + dabs(uC(i,jm1,k))
     2     + dabs(uC(i,jm1,kp1))
               oCip1jkp1 =
     2     + dabs(oE(i,j,kp1))
     2     + dabs(oE(ip1,j,kp1))
     2     + dabs(oN(ip1,jm1,kp1))
     2     + dabs(oN(ip1,j,kp1))
     2     + dabs(uC(ip1,j,k))
     2     + dabs(uC(ip1,j,kp1))
               oCim1jkp1 =
     2     + dabs(oE(im2,j,kp1))
     2     + dabs(oE(im1,j,kp1))
     2     + dabs(oN(im1,jm1,kp1))
     2     + dabs(oN(im1,j,kp1))
     2     + dabs(uC(im1,j,k))
     2     + dabs(uC(im1,j,kp1))
               oCip1jp1kp1 =
     2     + dabs(oE(i,jp1,kp1))
     2     + dabs(oE(ip1,jp1,kp1))
     2     + dabs(oN(ip1,j,kp1))
     2     + dabs(oN(ip1,jp1,kp1))
     2     + dabs(uC(ip1,jp1,k))
     2     + dabs(uC(ip1,jp1,kp1))
               oCim1jp1kp1 =
     2     + dabs(oE(im1,jp1,kp1))
     2     + dabs(oE(i,jp1,kp1))
     2     + dabs(oN(i,j,kp1))
     2     + dabs(oN(i,jp1,kp1))
     2     + dabs(uC(i,jp1,k))
     2     + dabs(uC(i,jp1,kp1))
               oCip1jm1kp1 =
     2     + dabs(oE(i,jm1,kp1))
     2     + dabs(oE(ip1,jm1,kp1))
     2     + dabs(oN(ip1,jm2,kp1))
     2     + dabs(oN(ip1,jm1,kp1))
     2     + dabs(uC(ip1,jm1,k))
     2     + dabs(uC(ip1,jm1,kp1))
               oCim1jm1kp1 =
     2     + dabs(oE(im1,jm1,kp1))
     2     + dabs(oE(i,jm1,kp1))
     2     + dabs(oN(i,jm2,kp1))
     2     + dabs(oN(i,jm1,kp1))
     2     + dabs(uC(i,jm1,k))
     2     + dabs(uC(i,jm1,kp1))

c* ********************************************************************
c* *** > oPC;
c* ********************************************************************

      oPC(ii,jj,kk) = won

c* ********************************************************************
c* *** > oPN;
c* ********************************************************************

      oPN(ii,jj,kk) = 
     2  dabs(oN(i,j,k))
     2  /  ( oCijp1k 
     2     - dabs(oE(im1,jp1,k))
     2     - dabs(oE(i,jp1,k))
     2     - dabs(uC(i,jp1,km1))
     2     - dabs(uC(i,jp1,k)))

c* ********************************************************************
c* *** > oPS;
c* ********************************************************************

      oPS(ii,jj,kk) = 
     2  dabs(oN(i,jm1,k))
     2  /  ( oCijm1k 
     2     - dabs(oE(im1,jm1,k))
     2     - dabs(oE(i,jm1,k))
     2     - dabs(uC(i,jm1,km1))
     2     - dabs(uC(i,jm1,k)))

c* ********************************************************************
c* *** > oPE;
c* ********************************************************************

      oPE(ii,jj,kk) = 
     2  dabs(oE(i,j,k))
     2  /  ( oCip1jk 
     2     - dabs(uC(ip1,j,km1))
     2     - dabs(uC(ip1,j,k))
     2     - dabs(oN(ip1,j,k))
     2     - dabs(oN(ip1,jm1,k)))

c* ********************************************************************
c* *** > oPW;
c* ********************************************************************

      oPW(ii,jj,kk) = 
     2  dabs(oE(im1,j,k))
     2  /  ( oCim1jk 
     2     - dabs(uC(im1,j,km1))
     2     - dabs(uC(im1,j,k))
     2     - dabs(oN(im1,j,k))
     2     - dabs(oN(im1,jm1,k)))

c* ********************************************************************
c* *** > oPNE;
c* ********************************************************************

      oPNE(ii,jj,kk) = 
     2   (dabs(oN(ip1,j,k)) * oPE(ii,jj,kk)
     2  + dabs(oE(i,jp1,k)) * oPN(ii,jj,kk))
     2  /  ( oCip1jp1k 
     2     - dabs(uC(ip1,jp1,km1))
     2     - dabs(uC(ip1,jp1,k)))

c* ********************************************************************
c* *** > oPNW;
c* ********************************************************************

      oPNW(ii,jj,kk) = 
     2   (dabs(oN(im1,j,k)) * oPW(ii,jj,kk) 
     2  + dabs(oE(im1,jp1,k)) * oPN(ii,jj,kk))
     2  /  ( oCim1jp1k 
     2     - dabs(uC(im1,jp1,km1)) 
     2     - dabs(uC(im1,jp1,k)))

c* ********************************************************************
c* *** > oPSE;
c* ********************************************************************

      oPSE(ii,jj,kk) = 
     2    (dabs(oN(ip1,jm1,k)) * oPE(ii,jj,kk) 
     2  +  dabs(oE(i,jm1,k)) * oPS(ii,jj,kk))
     2  /  ( oCip1jm1k 
     2     - dabs(uC(ip1,jm1,km1)) 
     2     - dabs(uC(ip1,jm1,k)))

c* ********************************************************************
c* *** > oPSW;
c* ********************************************************************

      oPSW(ii,jj,kk) = 
     2    (dabs(oN(im1,jm1,k)) * oPW(ii,jj,kk)
     2  +  dabs(oE(im1,jm1,k)) * oPS(ii,jj,kk))
     2  /  ( oCim1jm1k 
     2     - dabs(uC(im1,jm1,km1)) 
     2     - dabs(uC(im1,jm1,k)))

c* ********************************************************************
c* *** > dPC;
c* ********************************************************************

      dPC(ii,jj,kk) = 
     2  dabs(uC(i,j,km1))
     2  /  ( oCijkm1 
     2     - dabs(oN(i,j,km1))
     2     - dabs(oN(i,jm1,km1))
     2     - dabs(oE(im1,j,km1))
     2     - dabs(oE(i,j,km1)))

c* ********************************************************************
c* *** > dPN;
c* ********************************************************************

      dPN(ii,jj,kk) = 
     2    (dabs(oN(i,j,km1)) * dPC(ii,jj,kk) 
     2  +  dabs(uC(i,jp1,km1)) * oPN(ii,jj,kk))
     2  /  ( oCijp1km1 
     2     - dabs(oE(im1,jp1,km1))
     2     - dabs(oE(i,jp1,km1)))

c* ********************************************************************
c* *** > dPS;
c* ********************************************************************

      dPS(ii,jj,kk) = 
     2    (dabs(oN(i,jm1,km1)) * dPC(ii,jj,kk)
     2  +  dabs(uC(i,jm1,km1)) * oPS(ii,jj,kk))
     2  /  ( oCijm1km1 
     2     - dabs(oE(im1,jm1,km1))
     2     - dabs(oE(i,jm1,km1)))

c* ********************************************************************
c* *** > dPE;
c* ********************************************************************

      dPE(ii,jj,kk) = 
     2   (dabs(uC(ip1,j,km1)) * oPE(ii,jj,kk)
     2 +  dabs(oE(i,j,km1)) * dPC(ii,jj,kk))
     2 /  ( oCip1jkm1 
     2    - dabs(oN(ip1,j,km1))
     2    - dabs(oN(ip1,jm1,km1)))

c* ********************************************************************
c* *** > dPW;
c* ********************************************************************

      dPW(ii,jj,kk) = 
     2    (dabs(uC(im1,j,km1)) * oPW(ii,jj,kk)
     2  +  dabs(oE(im1,j,km1)) * dPC(ii,jj,kk))
     2  /  ( oCim1jkm1 
     2     - dabs(oN(im1,j,km1)) 
     2     - dabs(oN(im1,jm1,km1)))

c* ********************************************************************
c* *** > dPNE;
c* ********************************************************************

      dPNE(ii,jj,kk) = 
     2    (dabs(uC(ip1,jp1,km1)) * oPNE(ii,jj,kk)
     2   + dabs(oE(i,jp1,km1)) * dPN(ii,jj,kk)
     2   + dabs(oN(ip1,j,km1)) * dPE(ii,jj,kk))
     2   /  oCip1jp1km1

c* ********************************************************************
c* *** > dPNW;
c* ********************************************************************

      dPNW(ii,jj,kk) = 
     2   (dabs(uC(im1,jp1,km1)) * oPNW(ii,jj,kk)
     2  + dabs(oE(im1,jp1,km1)) * dPN(ii,jj,kk)
     2  + dabs(oN(im1,j,km1)) * dPW(ii,jj,kk))
     2  /  oCim1jp1km1

c* ********************************************************************
c* *** > dPSE;
c* ********************************************************************

      dPSE(ii,jj,kk) = 
     2   (dabs(uC(ip1,jm1,km1)) * oPSE(ii,jj,kk)
     2  + dabs(oE(i,jm1,km1)) * dPS(ii,jj,kk)
     2  + dabs(oN(ip1,jm1,km1)) * dPE(ii,jj,kk))
     2  /  oCip1jm1km1

c* ********************************************************************
c* *** > dPSW;
c* ********************************************************************

      dPSW(ii,jj,kk) = 
     2    (dabs(uC(im1,jm1,km1)) * oPSW(ii,jj,kk)
     2   + dabs(oE(im1,jm1,km1)) * dPS(ii,jj,kk)
     2   + dabs(oN(im1,jm1,km1)) * dPW(ii,jj,kk))
     2   /  oCim1jm1km1

c* ********************************************************************
c* *** > uPC;
c* ********************************************************************

      uPC(ii,jj,kk) = 
     2  dabs(uC(i,j,k))
     2  /  ( oCijkp1 
     2     - dabs(oN(i,j,kp1))
     2     - dabs(oN(i,jm1,kp1))
     2     - dabs(oE(im1,j,kp1))
     2     - dabs(oE(i,j,kp1)))

c* ********************************************************************
c* *** > uPN;
c* ********************************************************************

      uPN(ii,jj,kk) = 
     2   (dabs(oN(i,j,kp1)) * uPC(ii,jj,kk) 
     2 +  dabs(uC(i,jp1,k)) * oPN(ii,jj,kk))
     2 /  ( oCijp1kp1 
     2    - dabs(oE(im1,jp1,kp1)) 
     2    - dabs(oE(i,jp1,kp1)))

c* ********************************************************************
c* *** > uPS;
c* ********************************************************************

      uPS(ii,jj,kk) = 
     2   (dabs(oN(i,jm1,kp1)) * uPC(ii,jj,kk)
     2 +  dabs(uC(i,jm1,k)) * oPS(ii,jj,kk))
     2 /  ( oCijm1kp1 
     2    - dabs(oE(im1,jm1,kp1)) 
     2    - dabs(oE(i,jm1,kp1)))

c* ********************************************************************
c* *** > uPE;
c* ********************************************************************

      uPE(ii,jj,kk) = 
     2   (dabs(uC(ip1,j,k)) * oPE(ii,jj,kk)
     2 +  dabs(oE(i,j,kp1)) * uPC(ii,jj,kk))
     2 /  ( oCip1jkp1 
     2    - dabs(oN(ip1,j,kp1))
     2    - dabs(oN(ip1,jm1,kp1)))

c* ********************************************************************
c* *** > uPW;
c* ********************************************************************

      uPW(ii,jj,kk) = 
     2    (dabs(uC(im1,j,k)) * oPW(ii,jj,kk)
     2 +   dabs(oE(im1,j,kp1)) * uPC(ii,jj,kk))
     2 /  ( oCim1jkp1
     2    - dabs(oN(im1,j,kp1)) 
     2    - dabs(oN(im1,jm1,kp1)))

c* ********************************************************************
c* *** > uPNE;
c* ********************************************************************

      uPNE(ii,jj,kk) = 
     2    (dabs(uC(ip1,jp1,k)) * oPNE(ii,jj,kk)
     2  +  dabs(oE(i,jp1,kp1)) * uPN(ii,jj,kk)
     2  +  dabs(oN(ip1,j,kp1)) * uPE(ii,jj,kk))
     2  /  oCip1jp1kp1

c* ********************************************************************
c* *** > uPNW;
c* ********************************************************************

      uPNW(ii,jj,kk) = 
     2   (dabs(uC(im1,jp1,k)) * oPNW(ii,jj,kk)
     2  + dabs(oE(im1,jp1,kp1)) * uPN(ii,jj,kk)
     2  + dabs(oN(im1,j,kp1)) * uPW(ii,jj,kk))
     2  /  oCim1jp1kp1

c* ********************************************************************
c* *** > uPSE;
c* ********************************************************************

      uPSE(ii,jj,kk) = 
     2   (dabs(uC(ip1,jm1,k)) * oPSE(ii,jj,kk)
     2  + dabs(oE(i,jm1,kp1)) * uPS(ii,jj,kk)
     2  + dabs(oN(ip1,jm1,kp1)) * uPE(ii,jj,kk))
     2  /  oCip1jm1kp1

c* ********************************************************************
c* *** > uPSW;
c* ********************************************************************

      uPSW(ii,jj,kk) = 
     2   (dabs(uC(im1,jm1,k)) * oPSW(ii,jj,kk)
     2  + dabs(oE(im1,jm1,kp1)) * uPS(ii,jj,kk)
     2  + dabs(oN(im1,jm1,kp1)) * uPW(ii,jj,kk))
     2  /  oCim1jm1kp1

c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildP_modop27 (nxf,nyf,nzf,nxc,nyc,nzc,ipc,rpc,ac,pc)
c* *********************************************************************
c* purpose: 
c*
c*    a modified 27-pt operator-based prologation.
c*
c*    call the routine to form the prolongation operator from a 
c*    27 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc
      double precision rpc(*),ac(nxf*nyf*nzf,*),pc(nxc*nyc*nzc,*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      call buildPb_modop27 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,ac(1,1),

     3   ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     4   ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),

     5   pc(1,1),pc(1,2),pc(1,3),pc(1,4),pc(1,5),pc(1,6),pc(1,7),
     6   pc(1,8),pc(1,9),pc(1,10),pc(1,11),pc(1,12),pc(1,13),pc(1,14),
     7   pc(1,15),pc(1,16),pc(1,17),pc(1,18),pc(1,19),pc(1,20),pc(1,21),
     8   pc(1,22),pc(1,23),pc(1,24),pc(1,25),pc(1,26),pc(1,27))
c*
c*    *** return and end ***
      return
      end
      subroutine buildPb_modop27 (nxf,nyf,nzf,nxc,nyc,nzc,
     2   ipc,rpc,oC,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     5   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     6   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW)
c* *********************************************************************
c* purpose:  
c*
c*    a modified 27-pt operator-based prologation.
c*
c*    form the prolongation operator from a 27 diagonal fine grid matrix.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),nxf,nyf,nzf,nxc,nyc,nzc,i,j,k,ii,jj,kk
      integer          im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
      integer          iim1,iip1,jjm1,jjp1,kkm1,kkp1
      double precision rpc(*),oC(nxf,nyf,nzf)

      double precision oE(nxf,nyf,nzf),oN(nxf,nyf,nzf)
      double precision uC(nxf,nyf,nzf),oNE(nxf,nyf,nzf)
      double precision oNW(nxf,nyf,nzf),uE(nxf,nyf,nzf)
      double precision uW(nxf,nyf,nzf),uN(nxf,nyf,nzf)
      double precision uS(nxf,nyf,nzf),uNE(nxf,nyf,nzf)
      double precision uNW(nxf,nyf,nzf),uSE(nxf,nyf,nzf)
      double precision uSW(nxf,nyf,nzf)

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

      double precision won,half,quarter,eighth
c*
c*    *** interpolation stencil ***
      won     = 1.0d0
      half    = 1.0d0 / 2.0d0
      quarter = 1.0d0 / 4.0d0
      eighth  = 1.0d0 / 8.0d0
      do 10 kk = 2, nzc-1
         k = 2 * kk - 1
         do 11 jj = 2, nyc-1
            j = 2 * jj - 1
            do 12 ii = 2, nxc-1
               i = 2 * ii - 1
c*
c*             *** index computations ***
               im1 = i-1
               ip1 = i+1
               im2 = i-2
               ip2 = i+2
               jm1 = j-1
               jp1 = j+1
               jm2 = j-2
               jp2 = j+2
               km1 = k-1
               kp1 = k+1
               km2 = k-2
               kp2 = k+2
               iim1 = ii-1
               iip1 = ii+1
               jjm1 = jj-1
               jjp1 = jj+1
               kkm1 = kk-1
               kkp1 = kk+1
c*
c* ********************************************************************
c* *** > oPC;
c* ********************************************************************

      oPC(ii,jj,kk) = won

c* ********************************************************************
c* *** > oPN;
c* ********************************************************************

      oPN(ii,jj,kk) = 
     2   (uNE(im1,j,km1) + uN(i,j,km1) + uNW(ip1,j,km1)
     2  + oNE(im1,j,k) + oN(i,j,k) + oNW(ip1,j,k) 
     2  + uSW(i,jp1,k) + uS(i,jp1,k) + uSE(i,jp1,k))  
     2  /  ( oC(i,jp1,k) 
     2     - oE(im1,jp1,k) - oE(i,jp1,k) 
     2     - uC(i,jp1,km1) - uE(im1,jp1,km1) - uW(ip1,jp1,km1) 
     2     - uC(i,jp1,k) - uW(i,jp1,k) - uE(i,jp1,k))

c* ********************************************************************
c* *** > oPS;
c* ********************************************************************

      oPS(ii,jj,kk) = 
     2   (uSE(im1,j,km1) + uS(i,j,km1) + uSW(ip1,j,km1)
     2  + oNW(i,jm1,k) + oN(i,jm1,k) + oNE(i,jm1,k) 
     2  + uNW(i,jm1,k) + uN(i,jm1,k) + uNE(i,jm1,k))  
     2  /  ( oC(i,jm1,k) 
     2     - oE(im1,jm1,k) - oE(i,jm1,k) 
     2     - uC(i,jm1,km1) - uE(im1,jm1,km1) - uW(ip1,jm1,km1) 
     2     - uC(i,jm1,k) - uW(i,jm1,k) - uE(i,jm1,k))

c* ********************************************************************
c* *** > oPE;
c* ********************************************************************

      oPE(ii,jj,kk) = 
     2   (uSE(i,jp1,km1) + oNW(ip1,j,k) + uNW(ip1,j,k)
     2  + uE(i,j,km1) + oE(i,j,k) + uW(ip1,j,k) 
     2  + uNE(i,jm1,km1) + oNE(i,jm1,k) + uSW(ip1,j,k))  
     2  /  ( oC(ip1,j,k) 
     2     - uC(ip1,j,km1) - uC(ip1,j,k) 
     2     - oN(ip1,j,k) - uS(ip1,jp1,km1) - uN(ip1,j,k) 
     2     - oN(ip1,jm1,k) - uN(ip1,jm1,km1) - uS(ip1,j,k))

c* ********************************************************************
c* *** > oPW;
c* ********************************************************************

      oPW(ii,jj,kk) = 
     2   (uSW(i,jp1,km1) + oNE(im1,j,k) + uNE(im1,j,k)
     2  + uW(i,j,km1) + oE(im1,j,k) + uE(im1,j,k)
     2  + uNW(i,jm1,km1) + oNW(i,jm1,k) + uSE(im1,j,k))  
     2  /  ( oC(im1,j,k) 
     2     - uC(im1,j,km1) - uC(im1,j,k) 
     2     - oN(im1,j,k) - uS(im1,jp1,km1) - uN(im1,j,k) 
     2     - oN(im1,jm1,k) - uN(im1,jm1,km1) - uS(im1,j,k))

c* ********************************************************************
c* *** > oPNE;
c* ********************************************************************

      oPNE(ii,jj,kk) = 
     2   (uNE(i,j,km1) + oNE(i,j,k) + uSW(ip1,jp1,k) 
     2 + (uN(ip1,j,km1) + oN(ip1,j,k) + uS(ip1,jp1,k))*oPE(ii,jj,kk)
     2 + (uE(i,jp1,km1) + oE(i,jp1,k) + uW(ip1,jp1,k))*oPN(ii,jj,kk))
     2 /  ( oC(ip1,jp1,k) 
     2    - uC(ip1,jp1,km1) - uC(ip1,jp1,k))

c* ********************************************************************
c* *** > oPNW;
c* ********************************************************************

      oPNW(ii,jj,kk) = 
     2   (uNW(i,j,km1) + oNW(i,j,k) + uSE(im1,jp1,k) 
     2 + (uN(im1,j,km1) + oN(im1,j,k) + uS(im1,jp1,k))*oPW(ii,jj,kk) 
     2 + (uW(i,jp1,km1) + oE(im1,jp1,k) + uE(im1,jp1,k))*oPN(ii,jj,kk))
     2 /  ( oC(im1,jp1,k) 
     2    - uC(im1,jp1,km1) - uC(im1,jp1,k))

c* ********************************************************************
c* *** > oPSE;
c* ********************************************************************

      oPSE(ii,jj,kk) = 
     2   (uSE(i,j,km1) + oNW(ip1,jm1,k) + uNW(ip1,jm1,k)
     2 + (uS(ip1,j,km1) + oN(ip1,jm1,k) + uN(ip1,jm1,k))*oPE(ii,jj,kk) 
     2 + (uE(i,jm1,km1) + oE(i,jm1,k) + uW(ip1,jm1,k))*oPS(ii,jj,kk))
     2 /  ( oC(ip1,jm1,k) 
     2    - uC(ip1,jm1,km1) - uC(ip1,jm1,k))

c* ********************************************************************
c* *** > oPSW;
c* ********************************************************************

      oPSW(ii,jj,kk) = 
     2   (uSW(i,j,km1) + oNE(im1,jm1,k) + uNE(im1,jm1,k)
     2 + (uS(im1,j,km1) + oN(im1,jm1,k) + uN(im1,jm1,k))*oPW(ii,jj,kk)
     2 + (uW(i,jm1,km1) + oE(im1,jm1,k) + uE(im1,jm1,k))*oPS(ii,jj,kk))
     2 /  ( oC(im1,jm1,k) 
     2    - uC(im1,jm1,km1) - uC(im1,jm1,k))

c* ********************************************************************
c* *** > dPC;
c* ********************************************************************

      dPC(ii,jj,kk) = 
     2   (uNW(i,j,km1) + uW(i,j,km1) + uSW(i,j,km1)
     2  + uN(i,j,km1) + uC(i,j,km1) + uS(i,j,km1) + uNE(i,j,km1)
     2  + uE(i,j,km1) + uSE(i,j,km1))  
     2  /  ( oC(i,j,km1) 
     2     - oN(i,j,km1) - oN(i,jm1,km1) 
     2     - oNW(i,j,km1) - oE(im1,j,km1) - oNE(im1,jm1,km1) 
     2     - oNE(i,j,km1) - oE(i,j,km1) - oNW(ip1,jm1,km1))

c* ********************************************************************
c* *** > dPN;
c* ********************************************************************

      dPN(ii,jj,kk) = 
     2   (uSW(i,jp1,km1) + uS(i,jp1,km1) + uSE(i,jp1,km1) 
     2 + (oNE(im1,j,km1) + oN(i,j,km1) + oNW(ip1,j,km1))*dPC(ii,jj,kk) 
     2 + (uW(i,jp1,km1) + uC(i,jp1,km1) + uE(i,jp1,km1))*oPN(ii,jj,kk))
     2 /  ( oC(i,jp1,km1) 
     2    - oE(im1,jp1,km1) - oE(i,jp1,km1))

c* ********************************************************************
c* *** > dPS;
c* ********************************************************************

      dPS(ii,jj,kk) = 
     2   (uNW(i,jm1,km1) + uN(i,jm1,km1) + uNE(i,jm1,km1) 
     2 + (oNW(i,jm1,km1) + oN(i,jm1,km1) + oNE(i,jm1,km1))*dPC(ii,jj,kk)
     2 + (uW(i,jm1,km1) + uC(i,jm1,km1) + uE(i,jm1,km1))*oPS(ii,jj,kk))
     2 /  ( oC(i,jm1,km1) 
     2    - oE(im1,jm1,km1) - oE(i,jm1,km1))

c* ********************************************************************
c* *** > dPE;
c* ********************************************************************

      dPE(ii,jj,kk) = 
     2   (uNW(ip1,j,km1) + uW(ip1,j,km1) + uSW(ip1,j,km1) 
     2 + (uN(ip1,j,km1) + uC(ip1,j,km1) + uS(ip1,j,km1))*oPE(ii,jj,kk)
     2 + (oNW(ip1,j,km1) + oE(i,j,km1) + oNE(i,jm1,km1))*dPC(ii,jj,kk))
     2 /  ( oC(ip1,j,km1) 
     2    - oN(ip1,j,km1) - oN(ip1,jm1,km1))

c* ********************************************************************
c* *** > dPW;
c* ********************************************************************

      dPW(ii,jj,kk) = 
     2  (uNE(im1,j,km1) + uE(im1,j,km1) + uSE(im1,j,km1) 
     2 +(uN(im1,j,km1) + uC(im1,j,km1) + uS(im1,j,km1))*oPW(ii,jj,kk)
     2 +(oNE(im1,j,km1) + oE(im1,j,km1) + oNW(i,jm1,km1))*dPC(ii,jj,kk))
     2 /  ( oC(im1,j,km1) 
     2    - oN(im1,j,km1) - oN(im1,jm1,km1))

c* ********************************************************************
c* *** > dPNE;
c* ********************************************************************

      dPNE(ii,jj,kk) = 
     2   (uSW(ip1,jp1,km1) 
     2  + uW(ip1,jp1,km1) * oPN(ii,jj,kk)
     2  + uS(ip1,jp1,km1) * oPE(ii,jj,kk) 
     2  + uC(ip1,jp1,km1) * oPNE(ii,jj,kk)
     2  + oNE(i,j,km1) * dPC(ii,jj,kk)
     2  + oE(i,jp1,km1) * dPN(ii,jj,kk)
     2  + oN(ip1,j,km1) * dPE(ii,jj,kk))
     2  /  oC(ip1,jp1,km1)

c* ********************************************************************
c* *** > dPNW;
c* ********************************************************************

      dPNW(ii,jj,kk) = 
     2   (uSE(im1,jp1,km1) 
     2  + uE(im1,jp1,km1) * oPN(ii,jj,kk)
     2  + uS(im1,jp1,km1) * oPW(ii,jj,kk)
     2  + uC(im1,jp1,km1) * oPNW(ii,jj,kk)
     2  + oNW(i,j,km1) * dPC(ii,jj,kk)
     2  + oE(im1,jp1,km1) * dPN(ii,jj,kk)
     2  + oN(im1,j,km1) * dPW(ii,jj,kk))
     2  /  oC(im1,jp1,km1)

c* ********************************************************************
c* *** > dPSE;
c* ********************************************************************

      dPSE(ii,jj,kk) = 
     2    (uNW(ip1,jm1,km1)
     2  + uW(ip1,jm1,km1) * oPS(ii,jj,kk)
     2  + uN(ip1,jm1,km1) * oPE(ii,jj,kk)
     2  + uC(ip1,jm1,km1) * oPSE(ii,jj,kk)
     2  + oNW(ip1,jm1,km1) * dPC(ii,jj,kk)
     2  + oE(i,jm1,km1) * dPS(ii,jj,kk)
     2  + oN(ip1,jm1,km1) * dPE(ii,jj,kk))
     2  /  oC(ip1,jm1,km1)

c* ********************************************************************
c* *** > dPSW;
c* ********************************************************************

      dPSW(ii,jj,kk) = 
     2    (uNE(im1,jm1,km1)
     2   + uE(im1,jm1,km1) * oPS(ii,jj,kk)
     2   + uN(im1,jm1,km1) * oPW(ii,jj,kk)
     2   + uC(im1,jm1,km1) * oPSW(ii,jj,kk)
     2   + oNE(im1,jm1,km1) * dPC(ii,jj,kk)
     2   + oE(im1,jm1,km1) * dPS(ii,jj,kk)
     2   + oN(im1,jm1,km1) * dPW(ii,jj,kk))
     2   /  oC(im1,jm1,km1)

c* ********************************************************************
c* *** > uPC;
c* ********************************************************************

      uPC(ii,jj,kk) = 
     2   (uSE(im1,jp1,k) + uE(im1,j,k) + uNE(im1,jm1,k)
     2  + uS(i,jp1,k) + uC(i,j,k) + uN(i,jm1,k) 
     2  + uSW(ip1,jp1,k) + uW(ip1,j,k) + uNW(ip1,jm1,k))  
     2  /  ( oC(i,j,kp1) 
     2     - oN(i,j,kp1) - oN(i,jm1,kp1) 
     2     - oNW(i,j,kp1) - oE(im1,j,kp1) - oNE(im1,jm1,kp1) 
     2     - oNE(i,j,kp1) - oE(i,j,kp1) - oNW(ip1,jm1,kp1))

c* ********************************************************************
c* *** > uPN;
c* ********************************************************************

      uPN(ii,jj,kk) = 
     2   (uNE(im1,j,k) + uN(i,j,k) + uNW(ip1,j,k) 
     2 + (oNE(im1,j,kp1) + oN(i,j,kp1) + oNW(ip1,j,kp1))*uPC(ii,jj,kk) 
     2 + (uE(im1,jp1,k) + uC(i,jp1,k) + uW(ip1,jp1,k))*oPN(ii,jj,kk))
     2 /  ( oC(i,jp1,kp1) 
     2    - oE(im1,jp1,kp1) - oE(i,jp1,kp1))

c* ********************************************************************
c* *** > uPS;
c* ********************************************************************

      uPS(ii,jj,kk) = 
     2   (uSE(im1,j,k) + uS(i,j,k) + uSW(ip1,j,k) 
     2 + (oNW(i,jm1,kp1) + oN(i,jm1,kp1) + oNE(i,jm1,kp1))*uPC(ii,jj,kk)
     2 + (uE(im1,jm1,k) + uC(i,jm1,k) + uW(ip1,jm1,k))*oPS(ii,jj,kk))
     2 /  ( oC(i,jm1,kp1) 
     2    - oE(im1,jm1,kp1) - oE(i,jm1,kp1))

c* ********************************************************************
c* *** > uPE;
c* ********************************************************************

      uPE(ii,jj,kk) = 
     2   (uSE(i,jp1,k) + uS(ip1,jp1,k) + uNE(i,jm1,k) 
     2 + (uS(ip1,jp1,k) + uC(ip1,j,k) + uN(ip1,jm1,k))*oPE(ii,jj,kk)
     2 + (oNW(ip1,j,kp1) + oE(i,j,kp1) + oNE(i,jm1,kp1))*uPC(ii,jj,kk))
     2 /  ( oC(ip1,j,kp1) 
     2    - oN(ip1,j,kp1) - oN(ip1,jm1,kp1))

c* ********************************************************************
c* *** > uPW;
c* ********************************************************************

      uPW(ii,jj,kk) = 
     2  (uSW(i,jp1,k) + uW(i,j,k) + uNW(i,jm1,k) 
     2 +(uS(im1,jp1,k) + uC(im1,j,k) + uN(im1,jm1,k))*oPW(ii,jj,kk)
     2 +(oNE(im1,j,kp1) + oE(im1,j,kp1) + oNW(i,jm1,kp1))*uPC(ii,jj,kk))
     2 /  ( oC(im1,j,kp1) 
     2    - oN(im1,j,kp1) - oN(im1,jm1,kp1))

c* ********************************************************************
c* *** > uPNE;
c* ********************************************************************

      uPNE(ii,jj,kk) = 
     2   (uNE(i,j,k) 
     2  + uE(i,jp1,k) * oPN(ii,jj,kk)
     2  + uN(ip1,j,k) * oPE(ii,jj,kk)
     2  + uC(ip1,jp1,k) * oPNE(ii,jj,kk)
     2  + oNE(i,j,kp1) * uPC(ii,jj,kk)
     2  + oE(i,jp1,kp1) * uPN(ii,jj,kk)
     2  + oN(ip1,j,kp1) * uPE(ii,jj,kk))
     2  /  oC(ip1,jp1,kp1)

c* ********************************************************************
c* *** > uPNW;
c* ********************************************************************

      uPNW(ii,jj,kk) = 
     2   (uNW(i,j,k) 
     2  + uW(i,jp1,k) * oPN(ii,jj,kk)
     2  + uN(im1,j,k) * oPW(ii,jj,kk)
     2  + uC(im1,jp1,k) * oPNW(ii,jj,kk)
     2  + oNW(i,j,kp1) * uPC(ii,jj,kk)
     2  + oE(im1,jp1,kp1) * uPN(ii,jj,kk)
     2  + oN(im1,j,kp1) * uPW(ii,jj,kk))
     2  /  oC(im1,jp1,kp1)

c* ********************************************************************
c* *** > uPSE;
c* ********************************************************************

      uPSE(ii,jj,kk) = 
     2   (uSE(i,j,k) 
     2  + uE(i,jm1,k) * oPS(ii,jj,kk)
     2  + uS(ip1,j,k) * oPE(ii,jj,kk)
     2  + uC(ip1,jm1,k) * oPSE(ii,jj,kk)
     2  + oNW(ip1,jm1,kp1) * uPC(ii,jj,kk)
     2  + oE(i,jm1,kp1) * uPS(ii,jj,kk)
     2  + oN(ip1,jm1,kp1) * uPE(ii,jj,kk))
     2  /  oC(ip1,jm1,kp1)

c* ********************************************************************
c* *** > uPSW;
c* ********************************************************************

      uPSW(ii,jj,kk) = 
     2   (uSW(i,j,k) 
     2  + uW(i,jm1,k) * oPS(ii,jj,kk)
     2  + uS(im1,j,k) * oPW(ii,jj,kk)
     2  + uC(im1,jm1,k) * oPSW(ii,jj,kk)
     2  + oNE(im1,jm1,kp1) * uPC(ii,jj,kk)
     2  + oE(im1,jm1,kp1) * uPS(ii,jj,kk)
     2  + oN(im1,jm1,kp1) * uPW(ii,jj,kk))
     2  /  oC(im1,jm1,kp1)

c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
