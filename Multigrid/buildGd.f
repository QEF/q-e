c* ///////////////////////////////////////////////////////////////////////////
c* @file    buildGd.f
c* @author  Michael Holst
c* @brief   The heart of PMG: MAPLE-generated Galerkin product expressions.
c* @version $Id: buildGd.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

      subroutine buildG (nxf,nyf,nzf,nxc,nyc,nzc,numdia,pcFF,acFF,ac)
c* *********************************************************************
c* purpose: 
c*
c*    form the galerkin coarse grid matrix from the fine grid
c*    matrix and the prolongation operator.
c*
c*    differentiate between 1, 7, and 27 diagonal fine grid matrices
c*    for efficiency.
c*    
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nxc,nyc,nzc,numdia
      double precision pcFF(nxc*nyc*nzc,*)
      double precision acFF(nxf*nyf*nzf,*),ac(nxc*nyc*nzc,*)
c*
cmdir 0 0
c*
c*    *** call the build routine ***
      if (numdia .eq. 1) then

         call buildG_1 (nxf,nyf,nzf,nxc,nyc,nzc,
    
     1      pcFF(1,1),pcFF(1,2),pcFF(1,3),pcFF(1,4),pcFF(1,5),
     2      pcFF(1,6),pcFF(1,7),pcFF(1,8),pcFF(1,9),pcFF(1,10),
     3      pcFF(1,11),pcFF(1,12),pcFF(1,13),pcFF(1,14),pcFF(1,15),
     4      pcFF(1,16),pcFF(1,17),pcFF(1,18),pcFF(1,19),pcFF(1,20),
     5      pcFF(1,21),pcFF(1,22),pcFF(1,23),pcFF(1,24),pcFF(1,25),
     6      pcFF(1,26),pcFF(1,27),
    
     1      acFF(1,1),
    
     1      ac(1,1),ac(1,2),ac(1,3),ac(1,4),
     2      ac(1,5),ac(1,6),ac(1,7),ac(1,8),ac(1,9),ac(1,10),ac(1,11),
     3      ac(1,12),ac(1,13),ac(1,14))

      elseif (numdia .eq. 7) then

         call buildG_7 (nxf,nyf,nzf,nxc,nyc,nzc,
    
     1      pcFF(1,1),pcFF(1,2),pcFF(1,3),pcFF(1,4),pcFF(1,5),
     2      pcFF(1,6),pcFF(1,7),pcFF(1,8),pcFF(1,9),pcFF(1,10),
     3      pcFF(1,11),pcFF(1,12),pcFF(1,13),pcFF(1,14),pcFF(1,15),
     4      pcFF(1,16),pcFF(1,17),pcFF(1,18),pcFF(1,19),pcFF(1,20),
     5      pcFF(1,21),pcFF(1,22),pcFF(1,23),pcFF(1,24),pcFF(1,25),
     6      pcFF(1,26),pcFF(1,27),
    
     1      acFF(1,1),acFF(1,2),acFF(1,3),acFF(1,4),
    
     1      ac(1,1),ac(1,2),ac(1,3),ac(1,4),
     2      ac(1,5),ac(1,6),ac(1,7),ac(1,8),ac(1,9),ac(1,10),ac(1,11),
     3      ac(1,12),ac(1,13),ac(1,14))

      elseif (numdia .eq. 27) then

         call buildG_27 (nxf,nyf,nzf,nxc,nyc,nzc,
    
     1      pcFF(1,1),pcFF(1,2),pcFF(1,3),pcFF(1,4),pcFF(1,5),
     2      pcFF(1,6),pcFF(1,7),pcFF(1,8),pcFF(1,9),pcFF(1,10),
     3      pcFF(1,11),pcFF(1,12),pcFF(1,13),pcFF(1,14),pcFF(1,15),
     4      pcFF(1,16),pcFF(1,17),pcFF(1,18),pcFF(1,19),pcFF(1,20),
     5      pcFF(1,21),pcFF(1,22),pcFF(1,23),pcFF(1,24),pcFF(1,25),
     6      pcFF(1,26),pcFF(1,27),
    
     1      acFF(1,1),acFF(1,2),acFF(1,3),acFF(1,4),
     2      acFF(1,5),acFF(1,6),acFF(1,7),acFF(1,8),acFF(1,9),
     3      acFF(1,10),acFF(1,11),acFF(1,12),acFF(1,13),acFF(1,14),
    
     1      ac(1,1),ac(1,2),ac(1,3),ac(1,4),
     2      ac(1,5),ac(1,6),ac(1,7),ac(1,8),ac(1,9),ac(1,10),ac(1,11),
     3      ac(1,12),ac(1,13),ac(1,14))

      else
         print*,'% BUILDG: invalid stencil type give...'
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine buildG_1 (nxf,nyf,nzf,nx,ny,nz,
     2   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     3   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     4   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW,
     5   oC,
     6   XoC,XoE,XoN,XuC,XoNE,XoNW,XuE,XuW,XuN,XuS,XuNE,XuNW,XuSE,XuSW)
c* ********************************************************************
c* purpose: compute a 27-point galerkin coarse grid matrix from 
c*          a 1-point (i.e., diagonal) fine grid matrix.
c*
c* expressions for the galerkin coarse grid stencil XA in terms of
c* the fine grid matrix stencil A and the interpolation operator
c* stencil P.  these stencils have the form:
c* 
c*    XA := array([
c*
c*      matrix([
c*         [ -XdNW(i,j,k), -XdN(i,j,k), -XdNE(i,j,k) ],
c*         [ -XdW(i,j,k),  -XdC(i,j,k), -XdE(i,j,k)  ],
c*         [ -XdSW(i,j,k), -XdS(i,j,k), -XdSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -XoNW(i,j,k), -XoN(i,j,k), -XoNE(i,j,k) ],
c*         [ -XoW(i,j,k),   XoC(i,j,k), -XoE(i,j,k)  ],
c*         [ -XoSW(i,j,k), -XoS(i,j,k), -XoSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -XuNW(i,j,k), -XuN(i,j,k), -XuNE(i,j,k) ],
c*         [ -XuW(i,j,k),  -XuC(i,j,k), -XuE(i,j,k)  ],
c*         [ -XuSW(i,j,k), -XuS(i,j,k), -XuSE(i,j,k) ] 
c*      ])
c*    ]):
c*
c*    A := array([
c*
c*      matrix([
c*         [  0,           0,          0          ],
c*         [  0,           0,          0          ],
c*         [  0,           0,          0          ]
c*      ]),
c*
c*      matrix([
c*         [  0,           0,          0          ],
c*         [  0,           oC(i,j,k),  0          ],
c*         [  0,           0,          0          ]
c*      ]),
c*
c*      matrix([
c*         [  0,           0,          0          ],
c*         [  0,           0,          0          ],
c*         [  0,           0,          0          ]
c*      ])
c*
c*   P := array([
c*
c*      matrix([
c*         [ dPNW(i,j,k), dPN(i,j,k), dPNE(i,j,k) ],
c*         [ dPW(i,j,k),  dPC(i,j,k), dPE(i,j,k)  ],
c*         [ dPSW(i,j,k), dPS(i,j,k), dPSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ oPNW(i,j,k), oPN(i,j,k), oPNE(i,j,k) ],
c*         [ oPW(i,j,k),  oPC(i,j,k), oPE(i,j,k)  ],
c*         [ oPSW(i,j,k), oPS(i,j,k), oPSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ uPNW(i,j,k), uPN(i,j,k), uPNE(i,j,k) ],
c*         [ uPW(i,j,k),  uPC(i,j,k), uPE(i,j,k)  ],
c*         [ uPSW(i,j,k), uPS(i,j,k), uPSE(i,j,k) ] 
c*      ])
c*    ]):
c*
c* author:  michael holst
c* ********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k,ii,jj,kk
      integer          im1,ip1,jm1,jp1,km1,kp1
      integer          iim1,iip1,jjm1,jjp1,kkm1,kkp1
      integer          nxm1,nym1,nzm1,nxf,nyf,nzf

      double precision TMP1_XoC,TMP2_XoC,TMP3_XoC

      double precision oC(nxf,nyf,nzf)

      double precision XoC(nx,ny,nz),XoE(nx,ny,nz)
      double precision XoN(nx,ny,nz),XuC(nx,ny,nz)
      double precision XoNE(nx,ny,nz),XoNW(nx,ny,nz)
      double precision XuE(nx,ny,nz),XuW(nx,ny,nz)
      double precision XuN(nx,ny,nz),XuS(nx,ny,nz)
      double precision XuNE(nx,ny,nz),XuNW(nx,ny,nz)
      double precision XuSE(nx,ny,nz),XuSW(nx,ny,nz)

      double precision oPC(nx,ny,nz),oPN(nx,ny,nz),oPS(nx,ny,nz)
      double precision oPE(nx,ny,nz),oPW(nx,ny,nz),oPNE(nx,ny,nz)
      double precision oPNW(nx,ny,nz),oPSE(nx,ny,nz),oPSW(nx,ny,nz)
      double precision uPC(nx,ny,nz),uPN(nx,ny,nz),uPS(nx,ny,nz)
      double precision uPE(nx,ny,nz),uPW(nx,ny,nz),uPNE(nx,ny,nz)
      double precision uPNW(nx,ny,nz),uPSE(nx,ny,nz),uPSW(nx,ny,nz)
      double precision dPC(nx,ny,nz),dPN(nx,ny,nz),dPS(nx,ny,nz)
      double precision dPE(nx,ny,nz),dPW(nx,ny,nz),dPNE(nx,ny,nz)
      double precision dPNW(nx,ny,nz),dPSE(nx,ny,nz),dPSW(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** define n and determine number of mesh points ***
      nxm1    = nx - 1
      nym1    = ny - 1
      nzm1    = nz - 1
c*
c*    *** build the operator ***
cmdir 3 1
      do 10 kk = 2, nz-1
         k = 2 * kk - 1
cmdir 3 2
         do 11 jj = 2, ny-1
            j = 2 * jj - 1
cmdir 3 3
            do 12 ii = 2, nx-1
               i = 2 * ii - 1
c*
c*             *** index computations ***
               im1 = i-1
               ip1 = i+1
               jm1 = j-1
               jp1 = j+1
               km1 = k-1
               kp1 = k+1
               iim1 = ii-1
               iip1 = ii+1
               jjm1 = jj-1
               jjp1 = jj+1
               kkm1 = kk-1
               kkp1 = kk+1
c*
c* ********************************************************************
c* *** > oC;
c* ********************************************************************

c* ***XoC(ii,jj,kk) =
      TMP1_XoC =
     2     uPS(ii,jj,kk)  * uPS(ii,jj,kk)  * oC(i,jm1,kp1) 
     2   + dPSW(ii,jj,kk) * dPSW(ii,jj,kk) * oC(im1,jm1,km1)
     2   + oPSW(ii,jj,kk) * oPSW(ii,jj,kk) * oC(im1,jm1,k)
     2   + uPSW(ii,jj,kk) * uPSW(ii,jj,kk) * oC(im1,jm1,kp1)
     2   + dPW(ii,jj,kk)  * dPW(ii,jj,kk)  * oC(im1,j,km1) 
     2   + oPW(ii,jj,kk)  * oPW(ii,jj,kk)  * oC(im1,j,k)
     2   + uPNW(ii,jj,kk) * uPNW(ii,jj,kk) * oC(im1,jp1,kp1)
     2   + dPS(ii,jj,kk)  * dPS(ii,jj,kk)  * oC(i,jm1,km1) 
     2   + oPS(ii,jj,kk)  * oPS(ii,jj,kk)  * oC(i,jm1,k)

      TMP2_XoC =
     2   + dPC(ii,jj,kk)  * dPC(ii,jj,kk)  * oC(i,j,km1) 
     2   + oPC(ii,jj,kk)  * oPC(ii,jj,kk)  * oC(i,j,k)
     2   + uPC(ii,jj,kk)  * uPC(ii,jj,kk)  * oC(i,j,kp1) 
     2   + dPN(ii,jj,kk)  * dPN(ii,jj,kk)  * oC(i,jp1,km1)
     2   + oPN(ii,jj,kk)  * oPN(ii,jj,kk)  * oC(i,jp1,k) 
     2   + uPW(ii,jj,kk)  * uPW(ii,jj,kk)  * oC(im1,j,kp1)
     2   + dPNW(ii,jj,kk) * dPNW(ii,jj,kk) * oC(im1,jp1,km1)
     2   + oPNW(ii,jj,kk) * oPNW(ii,jj,kk) * oC(im1,jp1,k) 
     2   + oPE(ii,jj,kk)  * oPE(ii,jj,kk)  * oC(ip1,j,k)

      TMP3_XoC =
     2   + uPE(ii,jj,kk)  * uPE(ii,jj,kk)  * oC(ip1,j,kp1)
     2   + dPNE(ii,jj,kk) * dPNE(ii,jj,kk) * oC(ip1,jp1,km1)
     2   + oPNE(ii,jj,kk) * oPNE(ii,jj,kk) * oC(ip1,jp1,k)
     2   + uPNE(ii,jj,kk) * uPNE(ii,jj,kk) * oC(ip1,jp1,kp1)
     2   + uPN(ii,jj,kk)  * uPN(ii,jj,kk)  * oC(i,jp1,kp1)
     2   + dPSE(ii,jj,kk) * dPSE(ii,jj,kk) * oC(ip1,jm1,km1)
     2   + oPSE(ii,jj,kk) * oPSE(ii,jj,kk) * oC(ip1,jm1,k)
     2   + uPSE(ii,jj,kk) * uPSE(ii,jj,kk) * oC(ip1,jm1,kp1)
     2   + dPE(ii,jj,kk)  * dPE(ii,jj,kk)  * oC(ip1,j,km1)

      XoC(ii,jj,kk) = TMP1_XoC + TMP2_XoC + TMP3_XoC

c* ********************************************************************
c* *** > oE;
c* ********************************************************************

      XoE(ii,jj,kk) =
     2   - dPSE(ii,jj,kk) * oC(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2   - oPSE(ii,jj,kk) * oC(ip1,jm1,k)   * oPSW(iip1,jj,kk)
     2   - uPSE(ii,jj,kk) * oC(ip1,jm1,kp1) * uPSW(iip1,jj,kk)
     2   - dPE(ii,jj,kk)  * oC(ip1,j,km1)   * dPW(iip1,jj,kk)
     2   - oPE(ii,jj,kk)  * oC(ip1,j,k)     * oPW(iip1,jj,kk)
     2   - uPE(ii,jj,kk)  * oC(ip1,j,kp1)   * uPW(iip1,jj,kk)
     2   - dPNE(ii,jj,kk) * oC(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2   - oPNE(ii,jj,kk) * oC(ip1,jp1,k)   * oPNW(iip1,jj,kk)
     2   - uPNE(ii,jj,kk) * oC(ip1,jp1,kp1) * uPNW(iip1,jj,kk)

c* ********************************************************************
c* *** > oN;
c* ********************************************************************

      XoN(ii,jj,kk) =
     2   - dPNW(ii,jj,kk) * oC(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2   - oPNW(ii,jj,kk) * oC(im1,jp1,k)   * oPSW(ii,jjp1,kk)
     2   - uPNW(ii,jj,kk) * oC(im1,jp1,kp1) * uPSW(ii,jjp1,kk)
     2   - dPN(ii,jj,kk)  * oC(i,jp1,km1)   * dPS(ii,jjp1,kk)
     2   - oPN(ii,jj,kk)  * oC(i,jp1,k)     * oPS(ii,jjp1,kk)
     2   - uPN(ii,jj,kk)  * oC(i,jp1,kp1)   * uPS(ii,jjp1,kk)
     2   - dPNE(ii,jj,kk) * oC(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2   - oPNE(ii,jj,kk) * oC(ip1,jp1, k)  * oPSE(ii,jjp1,kk)
     2   - uPNE(ii,jj,kk) * oC(ip1,jp1,kp1) * uPSE(ii,jjp1,kk)

c* ********************************************************************
c* *** > uC;
c* ********************************************************************

      XuC(ii,jj,kk) =
     2   - dPSW(ii,jj,kkp1) * oC(im1,jm1,kp1) * uPSW(ii,jj,kk)
     2   - dPW(ii,jj,kkp1)  * oC(im1,j,kp1)   * uPW(ii,jj,kk)
     2   - dPNW(ii,jj,kkp1) * oC(im1,jp1,kp1) * uPNW(ii,jj,kk)
     2   - dPS(ii,jj,kkp1)  * oC(i,jm1,kp1)   * uPS(ii,jj,kk)
     2   - dPC(ii,jj,kkp1)  * oC(i,j,kp1)     * uPC(ii,jj,kk)
     2   - dPN(ii,jj,kkp1)  * oC(i,jp1,kp1)   * uPN(ii,jj,kk)
     2   - dPSE(ii,jj,kkp1) * oC(ip1,jm1,kp1) * uPSE(ii,jj,kk)
     2   - dPE(ii,jj,kkp1)  * oC(ip1,j,kp1)   * uPE(ii,jj,kk)
     2   - dPNE(ii,jj,kkp1) * oC(ip1,jp1,kp1) * uPNE(ii,jj,kk)

c* ********************************************************************
c* *** > oNE;
c* ********************************************************************

      XoNE(ii,jj,kk) =
     2    - dPNE(ii,jj,kk) * oC(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2    - oPNE(ii,jj,kk) * oC(ip1,jp1,k)   * oPSW(iip1,jjp1,kk)
     2    - uPNE(ii,jj,kk) * oC(ip1,jp1,kp1) * uPSW(iip1,jjp1,kk)


c* ********************************************************************
c* *** > oNW;
c* ********************************************************************

      XoNW(ii,jj,kk) =
     2    - dPNW(ii,jj,kk) * oC(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2    - oPNW(ii,jj,kk) * oC(im1,jp1,k)   * oPSE(iim1,jjp1,kk)
     2    - uPNW(ii,jj,kk) * oC(im1,jp1,kp1) * uPSE(iim1,jjp1,kk)

c* ********************************************************************
c* *** > uE;
c* ********************************************************************

      XuE(ii,jj,kk) =
     2   - uPSE(ii,jj,kk) * oC(ip1,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2   - uPE(ii,jj,kk)  * oC(ip1,j,kp1)   * dPW(iip1,jj,kkp1)
     2   - uPNE(ii,jj,kk) * oC(ip1,jp1,kp1) * dPNW(iip1,jj,kkp1)

c* ********************************************************************
c* *** > uW;
c* ********************************************************************

      XuW(ii,jj,kk) =
     2   - uPSW(ii,jj,kk) * oC(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2   - uPW(ii,jj,kk)  * oC(im1,j,kp1)   * dPE(iim1,jj,kkp1)
     2   - uPNW(ii,jj,kk) * oC(im1,jp1,kp1) * dPNE(iim1,jj,kkp1)

c* ********************************************************************
c* *** > uN;
c* ********************************************************************

      XuN(ii,jj,kk) =
     2   - uPNW(ii,jj,kk) * oC(im1,jp1,kp1) * dPSW(ii,jjp1,kkp1)
     2   - uPN(ii,jj,kk)  * oC(i,jp1,kp1)   * dPS(ii,jjp1,kkp1)
     2   - uPNE(ii,jj,kk) * oC(ip1,jp1,kp1) * dPSE(ii,jjp1,kkp1)

c* ********************************************************************
c* *** > uS;
c* ********************************************************************

      XuS(ii,jj,kk) =
     2   - uPSW(ii,jj,kk) * oC(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2   - uPS(ii,jj,kk)  * oC(i,jm1,kp1)   * dPN(ii,jjm1,kkp1)
     2   - uPSE(ii,jj,kk) * oC(ip1,jm1,kp1) * dPNE(ii,jjm1,kkp1)

c* ********************************************************************
c* *** > uNE;
c* ********************************************************************

      XuNE(ii,jj,kk) =
     2    - uPNE(ii,jj,kk) * oC(ip1,jp1,kp1) * dPSW(iip1,jjp1,kkp1)

c* ********************************************************************
c* *** > uNW;
c* ********************************************************************

      XuNW(ii,jj,kk) =
     2    - uPNW(ii,jj,kk) * oC(im1,jp1,kp1) * dPSE(iim1,jjp1,kkp1)

c* ********************************************************************
c* *** > uSE;
c* ********************************************************************

      XuSE(ii,jj,kk) =
     2    - uPSE(ii,jj,kk) * oC(ip1,jm1,kp1) * dPNW(iip1,jjm1,kkp1)

c* ********************************************************************
c* *** > uSW;
c* ********************************************************************

      XuSW(ii,jj,kk) =
     2    - uPSW(ii,jj,kk) * oC(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)

c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildG_7 (nxf,nyf,nzf,nx,ny,nz,
     2   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     3   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     4   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW,
     5   oC,oE,oN,uC,
     6   XoC,XoE,XoN,XuC,XoNE,XoNW,XuE,XuW,XuN,XuS,XuNE,XuNW,XuSE,XuSW)
c* ********************************************************************
c* purpose: compute a 27-point galerkin coarse grid matrix from 
c*          a 7-point fine grid matrix.
c*
c* expressions for the galerkin coarse grid stencil XA in terms of
c* the fine grid matrix stencil A and the interpolation operator
c* stencil P.  these stencils have the form:
c* 
c*    XA := array([
c*
c*      matrix([
c*         [ -XdNW(i,j,k), -XdN(i,j,k), -XdNE(i,j,k) ],
c*         [ -XdW(i,j,k),  -XdC(i,j,k), -XdE(i,j,k)  ],
c*         [ -XdSW(i,j,k), -XdS(i,j,k), -XdSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -XoNW(i,j,k), -XoN(i,j,k), -XoNE(i,j,k) ],
c*         [ -XoW(i,j,k),   XoC(i,j,k), -XoE(i,j,k)  ],
c*         [ -XoSW(i,j,k), -XoS(i,j,k), -XoSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -XuNW(i,j,k), -XuN(i,j,k), -XuNE(i,j,k) ],
c*         [ -XuW(i,j,k),  -XuC(i,j,k), -XuE(i,j,k)  ],
c*         [ -XuSW(i,j,k), -XuS(i,j,k), -XuSE(i,j,k) ] 
c*      ])
c*    ]):
c*
c*    A := array([
c*
c*      matrix([
c*         [  0,           0,          0          ],
c*         [  0,          -dC(i,j,k),  0          ],
c*         [  0,           0,          0          ]
c*      ]),
c*
c*      matrix([
c*         [  0,          -oN(i,j,k),  0          ],
c*         [ -oW(i,j,k),   oC(i,j,k), -oE(i,j,k)  ],
c*         [  0,          -oS(i,j,k),  0          ]
c*      ]),
c*
c*      matrix([
c*         [  0,           0,          0          ],
c*         [  0,          -uC(i,j,k),  0          ],
c*         [  0,           0,          0          ]
c*      ])
c*
c*   P := array([
c*
c*      matrix([
c*         [ dPNW(i,j,k), dPN(i,j,k), dPNE(i,j,k) ],
c*         [ dPW(i,j,k),  dPC(i,j,k), dPE(i,j,k)  ],
c*         [ dPSW(i,j,k), dPS(i,j,k), dPSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ oPNW(i,j,k), oPN(i,j,k), oPNE(i,j,k) ],
c*         [ oPW(i,j,k),  oPC(i,j,k), oPE(i,j,k)  ],
c*         [ oPSW(i,j,k), oPS(i,j,k), oPSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ uPNW(i,j,k), uPN(i,j,k), uPNE(i,j,k) ],
c*         [ uPW(i,j,k),  uPC(i,j,k), uPE(i,j,k)  ],
c*         [ uPSW(i,j,k), uPS(i,j,k), uPSE(i,j,k) ] 
c*      ])
c*    ]):
c*
c* in addition, A is assumed to be symmetric so that:
c*
c*    oS  := proc(x,y,z) RETURN( oN(x,y-1,z) ): end:
c*    oW  := proc(x,y,z) RETURN( oE(x-1,y,z) ): end:
c*    dC  := proc(x,y,z) RETURN( uC(x,y,z-1) ): end:
c*
c* author:  michael holst
c* ********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k,ii,jj,kk
      integer          im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
      integer          iim1,iip1,jjm1,jjp1,kkm1,kkp1
      integer          nxm1,nym1,nzm1,nxf,nyf,nzf

      double precision TMP1_XoC,TMP2_XoC,TMP3_XoC,TMP4_XoC
      double precision TMP5_XoC,TMP6_XoC,TMP7_XoC,TMP8_XoC
      double precision TMP9_XoC
      double precision TMP1_XoE,TMP2_XoE,TMP3_XoE,TMP4_XoE
      double precision TMP1_XoN,TMP2_XoN,TMP3_XoN,TMP4_XoN
      double precision TMP1_XuC,TMP2_XuC,TMP3_XuC,TMP4_XuC

      double precision oC(nxf,nyf,nzf),oE(nxf,nyf,nzf)
      double precision oN(nxf,nyf,nzf),uC(nxf,nyf,nzf)

      double precision XoC(nx,ny,nz),XoE(nx,ny,nz)
      double precision XoN(nx,ny,nz),XuC(nx,ny,nz)
      double precision XoNE(nx,ny,nz),XoNW(nx,ny,nz)
      double precision XuE(nx,ny,nz),XuW(nx,ny,nz)
      double precision XuN(nx,ny,nz),XuS(nx,ny,nz)
      double precision XuNE(nx,ny,nz),XuNW(nx,ny,nz)
      double precision XuSE(nx,ny,nz),XuSW(nx,ny,nz)

      double precision oPC(nx,ny,nz),oPN(nx,ny,nz),oPS(nx,ny,nz)
      double precision oPE(nx,ny,nz),oPW(nx,ny,nz),oPNE(nx,ny,nz)
      double precision oPNW(nx,ny,nz),oPSE(nx,ny,nz),oPSW(nx,ny,nz)
      double precision uPC(nx,ny,nz),uPN(nx,ny,nz),uPS(nx,ny,nz)
      double precision uPE(nx,ny,nz),uPW(nx,ny,nz),uPNE(nx,ny,nz)
      double precision uPNW(nx,ny,nz),uPSE(nx,ny,nz),uPSW(nx,ny,nz)
      double precision dPC(nx,ny,nz),dPN(nx,ny,nz),dPS(nx,ny,nz)
      double precision dPE(nx,ny,nz),dPW(nx,ny,nz),dPNE(nx,ny,nz)
      double precision dPNW(nx,ny,nz),dPSE(nx,ny,nz),dPSW(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** define n and determine number of mesh points ***
      nxm1    = nx - 1
      nym1    = ny - 1
      nzm1    = nz - 1
c*
c*    *** build the operator ***
cmdir 3 1
      do 10 kk = 2, nz-1
         k = 2 * kk - 1
cmdir 3 2
         do 11 jj = 2, ny-1
            j = 2 * jj - 1
cmdir 3 3
            do 12 ii = 2, nx-1
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
c* *** > oC;
c* ********************************************************************

c* ***XoC(ii,jj,kk) =
      TMP1_XoC =
     2     dPSW(ii,jj,kk) * ( oC(im1,jm1,km1) * dPSW(ii,jj,kk)
     2                      - uC(im1,jm1,km1) * oPSW(ii,jj,kk)
     2                      - oN(im1,jm1,km1) * dPW(ii,jj,kk)
     2                      - oE(im1,jm1,km1) * dPS(ii,jj,kk)) 

     2   + oPNE(ii,jj,kk) * (- oE(i,jp1,k) * oPN(ii,jj,kk) 
     2                       - oN(ip1,j,k) * oPE(ii,jj,kk)
     2                       - uC(ip1,jp1,km1) * dPNE(ii,jj,kk)
     2                       + oC(ip1,jp1,k) * oPNE(ii,jj,kk) 
     2                       - uC(ip1,jp1,k) * uPNE(ii,jj,kk))

     2   + dPNE(ii,jj,kk) * (- oE(i,jp1,km1) * dPN(ii,jj,kk)
     2                       - oN(ip1,j,km1) * dPE(ii,jj,kk) 
     2                       + oC(ip1,jp1,km1) * dPNE(ii,jj,kk)
     2                       - uC(ip1,jp1,km1) * oPNE(ii,jj,kk)) 

      TMP2_XoC =
     2   + dPSE(ii,jj,kk) * (- oE(i,jm1,km1) * dPS(ii,jj,kk) 
     2                       + oC(ip1,jm1,km1) * dPSE(ii,jj,kk)
     2                       - uC(ip1,jm1,km1) * oPSE(ii,jj,kk)
     2                       - oN(ip1,jm1,km1) * dPE(ii,jj,kk)) 

     2   + uPE(ii,jj,kk) * (- oE(i,j,kp1) * uPC(ii,jj,kk) 
     2                      - oN(ip1,jm1,kp1) * uPSE(ii,jj,kk)
     2                      - uC(ip1,j,k) * oPE(ii,jj,kk) 
     2                      + oC(ip1,j,kp1) * uPE(ii,jj,kk)
     2                      - oN(ip1,j,kp1) * uPNE(ii,jj,kk)) 

     2   + oPE(ii,jj,kk) * (- oE(i,j,k) * oPC(ii,jj,kk) 
     2                      - oN(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                      - uC(ip1,j,km1) * dPE(ii,jj,kk) 
     2                      + oC(ip1,j,k) * oPE(ii,jj,kk)
     2                      - uC(ip1,j,k) * uPE(ii,jj,kk) 
     2                      - oN(ip1,j,k) * oPNE(ii,jj,kk))

      TMP3_XoC =
     2   + dPE(ii,jj,kk) * (- oE(i,j,km1) * dPC(ii,jj,kk)
     2                      - oN(ip1,jm1,km1) * dPSE(ii,jj,kk) 
     2                      + oC(ip1,j,km1) * dPE(ii,jj,kk)
     2                      - uC(ip1,j,km1) * oPE(ii,jj,kk) 
     2                      - oN(ip1,j,km1) * dPNE(ii,jj,kk))

     2   + uPSE(ii,jj,kk) * (- oE(i,jm1,kp1) * uPS(ii,jj,kk)
     2                       - uC(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                       + oC(ip1,jm1,kp1) * uPSE(ii,jj,kk)
     2                       - oN(ip1,jm1,kp1) * uPE(ii,jj,kk))

     2   + uPNE(ii,jj,kk) * (- oE(i,jp1,kp1) * uPN(ii,jj,kk) 
     2                       - oN(ip1,j,kp1) * uPE(ii,jj,kk)
     2                       - uC(ip1,jp1,k) * oPNE(ii,jj,kk)
     2                       + oC(ip1,jp1,kp1) * uPNE(ii,jj,kk))

      TMP4_XoC =
     2   + oPS(ii,jj,kk) * (- oE(im1,jm1,k) * oPSW(ii,jj,kk) 
     2                      - uC(i,jm1,km1) * dPS(ii,jj,kk)
     2                      + oC(i,jm1,k) * oPS(ii,jj,kk) 
     2                      - uC(i,jm1,k) * uPS(ii,jj,kk)
     2                      - oN(i,jm1,k) * oPC(ii,jj,kk) 
     2                      - oE(i,jm1,k) * oPSE(ii,jj,kk))

     2   + dPS(ii,jj,kk) * (- oE(im1,jm1,km1) * dPSW(ii,jj,kk)
     2                      + oC(i,jm1,km1) * dPS(ii,jj,kk) 
     2                      - uC(i,jm1,km1) * oPS(ii,jj,kk)
     2                      - oN(i,jm1,km1) * dPC(ii,jj,kk) 
     2                      - oE(i,jm1,km1) * dPSE(ii,jj,kk))

     2   + oPSE(ii,jj,kk) * (- oE(i,jm1,k) * oPS(ii,jj,kk)
     2                       - uC(ip1,jm1,km1) * dPSE(ii,jj,kk)
     2                       + oC(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                       - uC(ip1,jm1,k) * uPSE(ii,jj,kk)
     2                       - oN(ip1,jm1,k) * oPE(ii,jj,kk))

      TMP5_XoC =
     2   + dPN(ii,jj,kk) * (- oE(im1,jp1,km1) * dPNW(ii,jj,kk) 
     2                      - oN(i,j,km1) * dPC(ii,jj,kk)
     2                      + oC(i,jp1,km1) * dPN(ii,jj,kk) 
     2                      - uC(i,jp1,km1) * oPN(ii,jj,kk)
     2                      - oE(i,jp1,km1) * dPNE(ii,jj,kk))

     2   + uPC(ii,jj,kk) * (- oE(im1,j,kp1) * uPW(ii,jj,kk) 
     2                      - oN(i,jm1,kp1) * uPS(ii,jj,kk)
     2                      - uC(i,j,k) * oPC(ii,jj,kk) 
     2                      + oC(i,j,kp1) * uPC(ii,jj,kk)
     2                      - oN(i,j,kp1) * uPN(ii,jj,kk) 
     2                      - oE(i,j,kp1) * uPE(ii,jj,kk))

     2   + oPC(ii,jj,kk) * (- oE(im1,j,k) * oPW(ii,jj,kk) 
     2                      - oN(i,jm1,k) * oPS(ii,jj,kk)
     2                      - uC(i,j,km1) * dPC(ii,jj,kk) 
     2                      + oC(i,j,k) * oPC(ii,jj,kk)
     2                      - uC(i,j,k) * uPC(ii,jj,kk) 
     2                      - oN(i,j,k) * oPN(ii,jj,kk)
     2                      - oE(i,j,k) * oPE(ii,jj,kk))

      TMP6_XoC =
     2   + dPC(ii,jj,kk) * (- oE(im1,j,km1) * dPW(ii,jj,kk) 
     2                      - oN(i,jm1,km1) * dPS(ii,jj,kk)
     2                      + oC(i,j,km1) * dPC(ii,jj,kk) 
     2                      - uC(i,j,km1) * oPC(ii,jj,kk)
     2                      - oN(i,j,km1) * dPN(ii,jj,kk)
     2                      - oE(i,j,km1) * dPE(ii,jj,kk))

     2   + uPS(ii,jj,kk) * (- oE(im1,jm1,kp1) * uPSW(ii,jj,kk)
     2                      - uC(i,jm1,k) * oPS(ii,jj,kk) 
     2                      + oC(i,jm1,kp1) * uPS(ii,jj,kk)
     2                      - oN(i,jm1,kp1) * uPC(ii,jj,kk) 
     2                      - oE(i,jm1,kp1) * uPSE(ii,jj,kk))

     2   + uPNW(ii,jj,kk) * (- oN(im1,j,kp1) * uPW(ii,jj,kk)
     2                       - uC(im1,jp1,k) * oPNW(ii,jj,kk)
     2                       + oC(im1,jp1,kp1) * uPNW(ii,jj,kk)
     2                       - oE(im1,jp1,kp1) * uPN(ii,jj,kk))

      TMP7_XoC =
     2   + dPW(ii,jj,kk) * (- oN(im1,jm1,km1) * dPSW(ii,jj,kk) 
     2                      + oC(im1,j,km1) * dPW(ii,jj,kk)
     2                      - uC(im1,j,km1) * oPW(ii,jj,kk) 
     2                      - oN(im1,j,km1) * dPNW(ii,jj,kk)
     2                      - oE(im1,j,km1) * dPC(ii,jj,kk))

     2   + uPSW(ii,jj,kk) * (- uC(im1,jm1,k) * oPSW(ii,jj,kk) 
     2                       + oC(im1,jm1,kp1) * uPSW(ii,jj,kk)
     2                       - oN(im1,jm1,kp1) * uPW(ii,jj,kk)
     2                       - oE(im1,jm1,kp1) * uPS(ii,jj,kk))

     2   + oPSW(ii,jj,kk) * (- uC(im1,jm1,km1) * dPSW(ii,jj,kk) 
     2                       + oC(im1,jm1,k) * oPSW(ii,jj,kk)
     2                       - uC(im1,jm1,k) * uPSW(ii,jj,kk)
     2                       - oN(im1,jm1,k) * oPW(ii,jj,kk)
     2                       - oE(im1,jm1,k) * oPS(ii,jj,kk))

      TMP8_XoC =
     2   + oPNW(ii,jj,kk) * (- oN(im1,j,k) * oPW(ii,jj,kk) 
     2                       - uC(im1,jp1,km1) * dPNW(ii,jj,kk)
     2                       + oC(im1,jp1,k) * oPNW(ii,jj,kk)
     2                       - uC(im1,jp1,k) * uPNW(ii,jj,kk)
     2                       - oE(im1,jp1,k) * oPN(ii,jj,kk))

     2   + dPNW(ii,jj,kk) * (- oN(im1,j,km1) * dPW(ii,jj,kk) 
     2                       + oC(im1,jp1,km1) * dPNW(ii,jj,kk)
     2                       - uC(im1,jp1,km1) * oPNW(ii,jj,kk)
     2                       - oE(im1,jp1,km1) * dPN(ii,jj,kk))

     2   + oPW(ii,jj,kk) * (- oN(im1,jm1,k) * oPSW(ii,jj,kk) 
     2                      - uC(im1,j,km1) * dPW(ii,jj,kk)
     2                      + oC(im1,j,k) * oPW(ii,jj,kk)
     2                      - uC(im1,j,k) * uPW(ii,jj,kk)
     2                      - oN(im1,j,k) * oPNW(ii,jj,kk)
     2                      - oE(im1,j,k) * oPC(ii,jj,kk))

      TMP9_XoC =
     2   + uPW(ii,jj,kk) * (- oN(im1,jm1,kp1) * uPSW(ii,jj,kk)
     2                      - uC(im1,j,k) * oPW(ii,jj,kk)
     2                      + oC(im1,j,kp1) * uPW(ii,jj,kk)
     2                      - oN(im1,j,kp1) * uPNW(ii,jj,kk)
     2                      - oE(im1,j,kp1) * uPC(ii,jj,kk))

     2   + uPN(ii,jj,kk) * (- oE(im1,jp1,kp1) * uPNW(ii,jj,kk)
     2                      - oN(i,j,kp1) * uPC(ii,jj,kk) 
     2                      - uC(i,jp1,k) * oPN(ii,jj,kk)
     2                      + oC(i,jp1,kp1) * uPN(ii,jj,kk)
     2                      - oE(i,jp1,kp1) * uPNE(ii,jj,kk))

     2   + oPN(ii,jj,kk) * (- oE(im1,jp1,k) * oPNW(ii,jj,kk)
     2                      - oN(i,j,k) * oPC(ii,jj,kk) 
     2                      - uC(i,jp1,km1) * dPN(ii,jj,kk)
     2                      + oC(i,jp1,k) * oPN(ii,jj,kk)
     2                      - uC(i,jp1,k) * uPN(ii,jj,kk)
     2                      - oE(i,jp1,k) * oPNE(ii,jj,kk))

      XoC(ii,jj,kk) = TMP1_XoC + TMP2_XoC + TMP3_XoC + TMP4_XoC
     2   + TMP5_XoC + TMP6_XoC + TMP7_XoC + TMP8_XoC + TMP9_XoC

c* ********************************************************************
c* *** > oE;
c* ********************************************************************

c* ***XoE(ii,jj,kk) =
      TMP1_XoE =
     2     dPS(ii,jj,kk) * oE(i,jm1,km1) * dPSW(iip1,jj,kk)
     2   + oPS(ii,jj,kk) * oE(i,jm1,k) * oPSW(iip1,jj,kk)
     2   + uPS(ii,jj,kk) * oE(i,jm1,kp1) * uPSW(iip1,jj,kk)
     2   + dPC(ii,jj,kk) * oE(i,j,km1) * dPW(iip1,jj,kk)
     2   + oPC(ii,jj,kk) * oE(i,j,k) * oPW(iip1,jj,kk)
     2   + uPC(ii,jj,kk) * oE(i,j,kp1) * uPW(iip1,jj,kk)
     2   + dPN(ii,jj,kk) * oE(i,jp1,km1) * dPNW(iip1,jj,kk)
     2   + oPN(ii,jj,kk) * oE(i,jp1,k) * oPNW(iip1,jj,kk)
     2   + uPN(ii,jj,kk) * oE(i,jp1,kp1) * uPNW(iip1,jj,kk) 

     2   - dPSE(ii,jj,kk) * (  oC(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                       - uC(ip1,jm1,km1) * oPSW(iip1,jj,kk)
     2                       - oN(ip1,jm1,km1) * dPW(iip1,jj,kk)
     2                       - oE(ip1,jm1,km1) * dPS(iip1,jj,kk))

      TMP2_XoE =
     2   - oPSE(ii,jj,kk) * (- uC(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                       + oC(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                       - uC(ip1,jm1,k) * uPSW(iip1,jj,kk)
     2                       - oN(ip1,jm1,k) * oPW(iip1,jj,kk)
     2                       - oE(ip1,jm1,k) * oPS(iip1,jj,kk))

     2   - uPSE(ii,jj,kk) * (- uC(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                       + oC(ip1,jm1,kp1) * uPSW(iip1,jj,kk)
     2                       - oN(ip1,jm1,kp1) * uPW(iip1,jj,kk)
     2                       - oE(ip1,jm1,kp1) * uPS(iip1,jj,kk))

     2   - dPE(ii,jj,kk) * (- oN(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                      + oC(ip1,j,km1) * dPW(iip1,jj,kk)
     2                      - uC(ip1,j,km1) * oPW(iip1,jj,kk)
     2                      - oN(ip1,j,km1) * dPNW(iip1,jj,kk)
     2                      - oE(ip1,j,km1) * dPC(iip1,jj,kk))

      TMP3_XoE =
     2   - oPE(ii,jj,kk) * (- oN(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                      - uC(ip1,j,km1) * dPW(iip1,jj,kk)
     2                      + oC(ip1,j,k) * oPW(iip1,jj,kk)
     2                      - uC(ip1,j,k) * uPW(iip1,jj,kk)
     2                      - oN(ip1,j,k) * oPNW(iip1,jj,kk)
     2                      - oE(ip1,j,k) * oPC(iip1,jj,kk))

     2   - uPE(ii,jj,kk) * (- oN(ip1,jm1,kp1) * uPSW(iip1,jj,kk)
     2                      - uC(ip1,j,k) * oPW(iip1,jj,kk) 
     2                      + oC(ip1,j,kp1) * uPW(iip1,jj,kk)
     2                      - oN(ip1,j,kp1) * uPNW(iip1,jj,kk)
     2                      - oE(ip1,j,kp1) * uPC(iip1,jj,kk))

     2   - dPNE(ii,jj,kk) * (- oN(ip1,j,km1) * dPW(iip1,jj,kk)
     2                       + oC(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2                       - uC(ip1,jp1,km1) * oPNW(iip1,jj,kk)
     2                       - oE(ip1,jp1,km1) * dPN(iip1,jj,kk))

      TMP4_XoE =
     2   - oPNE(ii,jj,kk) * (- oN(ip1,j,k) * oPW(iip1,jj,kk)
     2                       - uC(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2                       + oC(ip1,jp1,k) * oPNW(iip1,jj,kk)
     2                       - uC(ip1,jp1,k) * uPNW(iip1,jj,kk)
     2                       - oE(ip1,jp1,k) * oPN(iip1,jj,kk))

     2   - uPNE(ii,jj,kk) * (- oN(ip1,j,kp1) * uPW(iip1,jj,kk)
     2                       - uC(ip1,jp1,k) * oPNW(iip1,jj,kk)
     2                       + oC(ip1,jp1,kp1) * uPNW(iip1,jj,kk)
     2                       - oE(ip1,jp1,kp1) * uPN(iip1,jj,kk))

      XoE(ii,jj,kk) = TMP1_XoE + TMP2_XoE + TMP3_XoE + TMP4_XoE

c* ********************************************************************
c* *** > oN;
c* ********************************************************************

c* ***XoN(ii,jj,kk) =
      TMP1_XoN =
     2     dPW(ii,jj,kk) * oN(im1,j,km1) * dPSW(ii,jjp1,kk)
     2   + oPW(ii,jj,kk) * oN(im1,j,k) * oPSW(ii,jjp1,kk)
     2   + uPW(ii,jj,kk) * oN(im1,j,kp1) * uPSW(ii,jjp1,kk) 

     2   - dPNW(ii,jj,kk) * (  oC(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                       - uC(im1,jp1,km1) * oPSW(ii,jjp1,kk)
     2                       - oN(im1,jp1,km1) * dPW(ii,jjp1,kk)
     2                       - oE(im1,jp1,km1) * dPS(ii,jjp1,kk))

     2   - oPNW(ii,jj,kk) * (- uC(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                       + oC(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                       - uC(im1,jp1,k) * uPSW(ii,jjp1,kk)
     2                       - oN(im1,jp1,k) * oPW(ii,jjp1,kk)
     2                       - oE(im1,jp1,k) * oPS(ii,jjp1,kk))

      TMP2_XoN =
     2   - uPNW(ii,jj,kk) * (- uC(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                       + oC(im1,jp1,kp1) * uPSW(ii,jjp1,kk)
     2                       - oN(im1,jp1,kp1) * uPW(ii,jjp1,kk)
     2                       - oE(im1,jp1,kp1) * uPS(ii,jjp1,kk))

     2   + dPC(ii,jj,kk) * oN(i,j,km1) * dPS(ii,jjp1,kk)
     2   + oPC(ii,jj,kk) * oN(i,j,k) * oPS(ii,jjp1,kk)
     2   + uPC(ii,jj,kk) * oN(i,j,kp1) * uPS(ii,jjp1,kk) 

     2   - dPN(ii,jj,kk) * (- oE(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                      + oC(i,jp1,km1) * dPS(ii,jjp1,kk)
     2                      - uC(i,jp1,km1) * oPS(ii,jjp1,kk)
     2                      - oN(i,jp1,km1) * dPC(ii,jjp1,kk)
     2                      - oE(i,jp1,km1) * dPSE(ii,jjp1,kk))

      TMP3_XoN =
     2   - oPN(ii,jj,kk) * (- oE(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                      - uC(i,jp1,km1) * dPS(ii,jjp1,kk)
     2                      + oC(i,jp1,k) * oPS(ii,jjp1,kk)
     2                      - uC(i,jp1,k) * uPS(ii,jjp1,kk) 
     2                      - oN(i,jp1,k) * oPC(ii,jjp1,kk)
     2                      - oE(i,jp1,k) * oPSE(ii,jjp1,kk)) 

     2   - uPN(ii,jj,kk) * (- oE(im1,jp1,kp1) * uPSW(ii,jjp1,kk)
     2                      - uC(i,jp1,k) * oPS(ii,jjp1,kk) 
     2                      + oC(i,jp1,kp1) * uPS(ii,jjp1,kk)
     2                      - oN(i,jp1,kp1) * uPC(ii,jjp1,kk)
     2                      - oE(i,jp1,kp1) * uPSE(ii,jjp1,kk))

     2   + dPE(ii,jj,kk) * oN(ip1,j,km1) * dPSE(ii,jjp1,kk)
     2   + oPE(ii,jj,kk) * oN(ip1,j,k) * oPSE(ii,jjp1,kk)
     2   + uPE(ii,jj,kk) * oN(ip1,j,kp1) * uPSE(ii,jjp1,kk) 

      TMP4_XoN =
     2   - dPNE(ii,jj,kk) * (- oE(i,jp1,km1) * dPS(ii,jjp1,kk)
     2                       + oC(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2                       - uC(ip1,jp1,km1) * oPSE(ii,jjp1,kk)
     2                       - oN(ip1,jp1,km1) * dPE(ii,jjp1,kk))

     2   - oPNE(ii,jj,kk) * (- oE(i,jp1,k) * oPS(ii,jjp1,kk)
     2                       - uC(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2                       + oC(ip1,jp1,k) * oPSE(ii,jjp1,kk)
     2                       - uC(ip1,jp1,k) * uPSE(ii,jjp1,kk)
     2                       - oN(ip1,jp1,k) * oPE(ii,jjp1,kk))

     2   - uPNE(ii,jj,kk) * (- oE(i,jp1,kp1) * uPS(ii,jjp1,kk)
     2                       - uC(ip1,jp1,k) * oPSE(ii,jjp1,kk)
     2                       + oC(ip1,jp1,kp1) * uPSE(ii,jjp1,kk)
     2                       - oN(ip1,jp1,kp1) * uPE(ii,jjp1,kk))

      XoN(ii,jj,kk) = TMP1_XoN + TMP2_XoN + TMP3_XoN + TMP4_XoN

c* ********************************************************************
c* *** > uC;
c* ********************************************************************

c* ***XuC(ii,jj,kk) =
      TMP1_XuC =
     2     oPSW(ii,jj,kk) * uC(im1,jm1,k) * dPSW(ii,jj,kkp1) 

     2   - uPSW(ii,jj,kk) * (  oC(im1,jm1,kp1) * dPSW(ii,jj,kkp1)
     2                       - uC(im1,jm1,kp1) * oPSW(ii,jj,kkp1)
     2                       - oN(im1,jm1,kp1) * dPW(ii,jj,kkp1)
     2                       - oE(im1,jm1,kp1) * dPS(ii,jj,kkp1))

     2   + oPW(ii,jj,kk) * uC(im1,j,k) * dPW(ii,jj,kkp1) 

     2   - uPW(ii,jj,kk) * (- oN(im1,jm1,kp1) * dPSW(ii,jj,kkp1)
     2                      + oC(im1,j,kp1) * dPW(ii,jj,kkp1)
     2                      - uC(im1,j,kp1) * oPW(ii,jj,kkp1)
     2                      - oN(im1,j,kp1) * dPNW(ii,jj,kkp1)
     2                      - oE(im1,j,kp1) * dPC(ii,jj,kkp1))

     2   + oPNW(ii,jj,kk) * uC(im1,jp1,k) * dPNW(ii,jj,kkp1) 

      TMP2_XuC =
     2   - uPNW(ii,jj,kk) * (- oN(im1,j,kp1) * dPW(ii,jj,kkp1)
     2                       + oC(im1,jp1,kp1) * dPNW(ii,jj,kkp1)
     2                       - uC(im1,jp1,kp1) * oPNW(ii,jj,kkp1)
     2                       - oE(im1,jp1,kp1) * dPN(ii,jj,kkp1))

     2   + oPS(ii,jj,kk) * uC(i,jm1,k) * dPS(ii,jj,kkp1) 

     2   - uPS(ii,jj,kk) * (- oE(im1,jm1,kp1) * dPSW(ii,jj,kkp1)
     2                      + oC(i,jm1,kp1) * dPS(ii,jj,kkp1)
     2                      - uC(i,jm1,kp1) * oPS(ii,jj,kkp1)
     2                      - oN(i,jm1,kp1) * dPC(ii,jj,kkp1)
     2                      - oE(i,jm1,kp1) * dPSE(ii,jj,kkp1))

     2   + oPC(ii,jj,kk) * uC(i,j,k) * dPC(ii,jj,kkp1)

     2   - uPC(ii,jj,kk) * (- oE(im1,j,kp1) * dPW(ii,jj,kkp1)
     2                      - oN(i,jm1,kp1) * dPS(ii,jj,kkp1)
     2                      + oC(i,j,kp1) * dPC(ii,jj,kkp1)
     2                      - uC(i,j,kp1) * oPC(ii,jj,kkp1) 
     2                      - oN(i,j,kp1) * dPN(ii,jj,kkp1)
     2                      - oE(i,j,kp1) * dPE(ii,jj,kkp1))

      TMP3_XuC =
     2   + oPN(ii,jj,kk) * uC(i,jp1,k) * dPN(ii,jj,kkp1) 

     2   - uPN(ii,jj,kk) * (- oE(im1,jp1,kp1) * dPNW(ii,jj,kkp1)
     2                      - oN(i,j,kp1) * dPC(ii,jj,kkp1) 
     2                      + oC(i,jp1,kp1) * dPN(ii,jj,kkp1)
     2                      - uC(i,jp1,kp1) * oPN(ii,jj,kkp1)
     2                      - oE(i,jp1,kp1) * dPNE(ii,jj,kkp1))

     2   + oPSE(ii,jj,kk) * uC(ip1,jm1,k) * dPSE(ii,jj,kkp1) 

     2   - uPSE(ii,jj,kk) * (- oE(i,jm1,kp1) * dPS(ii,jj,kkp1)
     2                       + oC(ip1,jm1,kp1) * dPSE(ii,jj,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPSE(ii,jj,kkp1)
     2                       - oN(ip1,jm1,kp1) * dPE(ii,jj,kkp1))

      TMP4_XuC =
     2   + oPE(ii,jj,kk) * uC(ip1,j,k) * dPE(ii,jj,kkp1) 

     2   - uPE(ii,jj,kk) * (- oE(i,j,kp1) * dPC(ii,jj,kkp1)
     2                      - oN(ip1,jm1,kp1) * dPSE(ii,jj,kkp1)
     2                      + oC(ip1,j,kp1) * dPE(ii,jj,kkp1)
     2                      - uC(ip1,j,kp1) * oPE(ii,jj,kkp1)
     2                      - oN(ip1,j,kp1) * dPNE(ii,jj,kkp1))

     2   + oPNE(ii,jj,kk) * uC(ip1,jp1,k) * dPNE(ii,jj,kkp1)

     2   - uPNE(ii,jj,kk) * (- oE(i,jp1,kp1) * dPN(ii,jj,kkp1)
     2                       - oN(ip1,j,kp1) * dPE(ii,jj,kkp1)
     2                       + oC(ip1,jp1,kp1) * dPNE(ii,jj,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPNE(ii,jj,kkp1))

      XuC(ii,jj,kk) = TMP1_XuC + TMP2_XuC + TMP3_XuC + TMP4_XuC

c* ********************************************************************
c* *** > oNE;
c* ********************************************************************

      XoNE(ii,jj,kk) =
     2     dPN(ii,jj,kk) * oE(i,jp1,km1) * dPSW(iip1,jjp1,kk)
     2   + oPN(ii,jj,kk) * oE(i,jp1,k) * oPSW(iip1,jjp1,kk)
     2   + uPN(ii,jj,kk) * oE(i,jp1,kp1) * uPSW(iip1,jjp1,kk)
     2   + dPE(ii,jj,kk) * oN(ip1,j,km1) * dPSW(iip1,jjp1,kk)
     2   + oPE(ii,jj,kk) * oN(ip1,j,k) * oPSW(iip1,jjp1,kk)
     2   + uPE(ii,jj,kk) * oN(ip1,j,kp1) * uPSW(iip1,jjp1,kk)

     2   - dPNE(ii,jj,kk) * (  oC(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                       - uC(ip1,jp1,km1) * oPSW(iip1,jjp1,kk)
     2                       - oN(ip1,jp1,km1) * dPW(iip1,jjp1,kk)
     2                       - oE(ip1,jp1,km1) * dPS(iip1,jjp1,kk))

     2   - oPNE(ii,jj,kk) * (- uC(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                       + oC(ip1,jp1,k) * oPSW(iip1,jjp1,kk)
     2                       - uC(ip1,jp1,k) * uPSW(iip1,jjp1,kk)
     2                       - oN(ip1,jp1,k) * oPW(iip1,jjp1,kk)
     2                       - oE(ip1,jp1,k) * oPS(iip1,jjp1,kk))

     2   - uPNE(ii,jj,kk) * (- uC(ip1,jp1,k) * oPSW(iip1,jjp1,kk)
     2                       + oC(ip1,jp1,kp1) * uPSW(iip1,jjp1,kk)
     2                       - oN(ip1,jp1,kp1) * uPW(iip1,jjp1,kk)
     2                       - oE(ip1,jp1,kp1) * uPS(iip1,jjp1,kk))

c* ********************************************************************
c* *** > oNW;
c* ********************************************************************

      XoNW(ii,jj,kk) =
     2     dPW(ii,jj,kk) * oN(im1,j,km1) * dPSE(iim1,jjp1,kk)
     2   + oPW(ii,jj,kk) * oN(im1,j,k) * oPSE(iim1,jjp1,kk)
     2   + uPW(ii,jj,kk) * oN(im1,j,kp1) * uPSE(iim1,jjp1,kk)

     2   - dPNW(ii,jj,kk) * (- oE(im2,jp1,km1) * dPS(iim1,jjp1,kk)
     2                       + oC(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                       - uC(im1,jp1,km1) * oPSE(iim1,jjp1,kk)
     2                       - oN(im1,jp1,km1) * dPE(iim1,jjp1,kk)) 

     2   - oPNW(ii,jj,kk) * (- oE(im2,jp1,k) * oPS(iim1,jjp1,kk)
     2                       - uC(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                       + oC(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                       - uC(im1,jp1,k) * uPSE(iim1,jjp1,kk)
     2                       - oN(im1,jp1,k) * oPE(iim1,jjp1,kk))

     2   - uPNW(ii,jj,kk) * (- oE(im2,jp1,kp1) * uPS(iim1,jjp1,kk)
     2                       - uC(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                       + oC(im1,jp1,kp1) * uPSE(iim1,jjp1,kk)
     2                       - oN(im1,jp1,kp1) * uPE(iim1,jjp1,kk))

     2   + dPN(ii,jj,kk) * oE(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2   + oPN(ii,jj,kk) * oE(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2   + uPN(ii,jj,kk) * oE(im1,jp1,kp1) * uPSE(iim1,jjp1,kk)

c* ********************************************************************
c* *** > uE;
c* ********************************************************************

      XuE(ii,jj,kk) =
     2     uPS(ii,jj,kk) * oE(i,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2   + uPC(ii,jj,kk) * oE(i,j,kp1) * dPW(iip1,jj,kkp1)
     2   + uPN(ii,jj,kk) * oE(i,jp1,kp1) * dPNW(iip1,jj,kkp1)
     2   + oPSE(ii,jj,kk) * uC(ip1,jm1,k) * dPSW(iip1,jj,kkp1)

     2   - uPSE(ii,jj,kk) * (  oC(ip1,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPSW(iip1,jj,kkp1)
     2                       - oN(ip1,jm1,kp1) * dPW(iip1,jj,kkp1)
     2                       - oE(ip1,jm1,kp1) * dPS(iip1,jj,kkp1))

     2   + oPE(ii,jj,kk) * uC(ip1,j,k) * dPW(iip1,jj,kkp1) 

     2   - uPE(ii,jj,kk) * (- oN(ip1,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2                      + oC(ip1,j,kp1) * dPW(iip1,jj,kkp1)
     2                      - uC(ip1,j,kp1) * oPW(iip1,jj,kkp1)
     2                      - oN(ip1,j,kp1) * dPNW(iip1,jj,kkp1)
     2                      - oE(ip1,j,kp1) * dPC(iip1,jj,kkp1))

     2   + oPNE(ii,jj,kk) * uC(ip1,jp1,k) * dPNW(iip1,jj,kkp1)

     2   - uPNE(ii,jj,kk) * (- oN(ip1,j,kp1) * dPW(iip1,jj,kkp1)
     2                       + oC(ip1,jp1,kp1) * dPNW(iip1,jj,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPNW(iip1,jj,kkp1)
     2                       - oE(ip1,jp1,kp1) * dPN(iip1,jj,kkp1))

c* ********************************************************************
c* *** > uW;
c* ********************************************************************

      XuW(ii,jj,kk) =
     2     oPSW(ii,jj,kk) * uC(im1,jm1,k) * dPSE(iim1,jj,kkp1) 

     2   - uPSW(ii,jj,kk) * (- oE(im2,jm1,kp1) * dPS(iim1,jj,kkp1)
     2                       + oC(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2                       - uC(im1,jm1,kp1) * oPSE(iim1,jj,kkp1)
     2                       - oN(im1,jm1,kp1) * dPE(iim1,jj,kkp1))

     2   + oPW(ii,jj,kk) * uC(im1,j,k) * dPE(iim1,jj,kkp1) 

     2   - uPW(ii,jj,kk) * (- oE(im2,j,kp1) * dPC(iim1,jj,kkp1)
     2                      - oN(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2                      + oC(im1,j,kp1) * dPE(iim1,jj,kkp1)
     2                      - uC(im1,j,kp1) * oPE(iim1,jj,kkp1)
     2                      - oN(im1,j,kp1) * dPNE(iim1,jj,kkp1))

     2   + oPNW(ii,jj,kk) * uC(im1,jp1,k) * dPNE(iim1,jj,kkp1)

     2   - uPNW(ii,jj,kk) * (- oE(im2,jp1,kp1) * dPN(iim1,jj,kkp1)
     2                       - oN(im1,j,kp1) * dPE(iim1,jj,kkp1)
     2                       + oC(im1,jp1,kp1) * dPNE(iim1,jj,kkp1)
     2                       - uC(im1,jp1,kp1) * oPNE(iim1,jj,kkp1))

     2   + uPS(ii,jj,kk) * oE(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2   + uPC(ii,jj,kk) * oE(im1,j,kp1) * dPE(iim1,jj,kkp1)
     2   + uPN(ii,jj,kk) * oE(im1,jp1,kp1) * dPNE(iim1,jj,kkp1)

c* ********************************************************************
c* *** > uN;
c* ********************************************************************

      XuN(ii,jj,kk) =
     2     uPW(ii,jj,kk) * oN(im1,j,kp1) * dPSW(ii,jjp1,kkp1)
     2   + oPNW(ii,jj,kk) * uC(im1,jp1,k) * dPSW(ii,jjp1,kkp1) 

     2   - uPNW(ii,jj,kk) * (  oC(im1,jp1,kp1) * dPSW(ii,jjp1,kkp1)
     2                       - uC(im1,jp1,kp1) * oPSW(ii,jjp1,kkp1)
     2                       - oN(im1,jp1,kp1) * dPW(ii,jjp1,kkp1)
     2                       - oE(im1,jp1,kp1) * dPS(ii,jjp1,kkp1))

     2   + uPC(ii,jj,kk) * oN(i,j,kp1) * dPS(ii,jjp1,kkp1)
     2   + oPN(ii,jj,kk) * uC(i,jp1,k) * dPS(ii,jjp1,kkp1) 

     2   - uPN(ii,jj,kk) * (- oE(im1,jp1,kp1) * dPSW(ii,jjp1,kkp1)
     2                      + oC(i,jp1,kp1) * dPS(ii,jjp1,kkp1)
     2                      - uC(i,jp1,kp1) * oPS(ii,jjp1,kkp1)
     2                      - oN(i,jp1,kp1) * dPC(ii,jjp1,kkp1)
     2                      - oE(i,jp1,kp1) * dPSE(ii,jjp1,kkp1))

     2   + uPE(ii,jj,kk) * oN(ip1,j,kp1) * dPSE(ii,jjp1,kkp1)
     2   + oPNE(ii,jj,kk) * uC(ip1,jp1,k) * dPSE(ii,jjp1,kkp1)

     2   - uPNE(ii,jj,kk) * (- oE(i,jp1,kp1) * dPS(ii,jjp1,kkp1)
     2                       + oC(ip1,jp1,kp1) * dPSE(ii,jjp1,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPSE(ii,jjp1,kkp1)
     2                       - oN(ip1,jp1,kp1) * dPE(ii,jjp1,kkp1))

c* ********************************************************************
c* *** > uS;
c* ********************************************************************

      XuS(ii,jj,kk) =
     2     oPSW(ii,jj,kk) * uC(im1,jm1,k) * dPNW(ii,jjm1,kkp1) 

     2   - uPSW(ii,jj,kk) * (- oN(im1,jm2,kp1) * dPW(ii,jjm1,kkp1)
     2                       + oC(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2                       - uC(im1,jm1,kp1) * oPNW(ii,jjm1,kkp1)
     2                       - oE(im1,jm1,kp1) * dPN(ii,jjm1,kkp1))

     2   + uPW(ii,jj,kk) * oN(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2   + oPS(ii,jj,kk) * uC(i,jm1,k) * dPN(ii,jjm1,kkp1) 

     2   - uPS(ii,jj,kk) * (- oE(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2                      - oN(i,jm2,kp1) * dPC(ii,jjm1,kkp1)
     2                      + oC(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                      - uC(i,jm1,kp1) * oPN(ii,jjm1,kkp1)
     2                      - oE(i,jm1,kp1) * dPNE(ii,jjm1,kkp1))

     2   + uPC(ii,jj,kk) * oN(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2   + oPSE(ii,jj,kk) * uC(ip1,jm1,k) * dPNE(ii,jjm1,kkp1)

     2   - uPSE(ii,jj,kk) * (- oE(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                       - oN(ip1,jm2,kp1) * dPE(ii,jjm1,kkp1)
     2                       + oC(ip1,jm1,kp1) * dPNE(ii,jjm1,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPNE(ii,jjm1,kkp1))

     2   + uPE(ii,jj,kk) * oN(ip1,jm1,kp1) * dPNE(ii,jjm1,kkp1)

c* ********************************************************************
c* *** > uNE;
c* ********************************************************************

      XuNE(ii,jj,kk) =
     2     uPN(ii,jj,kk) * oE(i,jp1,kp1) * dPSW(iip1,jjp1,kkp1)
     2   + uPE(ii,jj,kk) * oN(ip1,j,kp1) * dPSW(iip1,jjp1,kkp1)
     2   + oPNE(ii,jj,kk) * uC(ip1,jp1,k) * dPSW(iip1,jjp1,kkp1)

     2   - uPNE(ii,jj,kk) * (  oC(ip1,jp1,kp1) * dPSW(iip1,jjp1,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPSW(iip1,jjp1,kkp1)
     2                       - oN(ip1,jp1,kp1) * dPW(iip1,jjp1,kkp1)
     2                       - oE(ip1,jp1,kp1) * dPS(iip1,jjp1,kkp1))

c* ********************************************************************
c* *** > uNW;
c* ********************************************************************

      XuNW(ii,jj,kk) =
     2     uPW(ii,jj,kk) * oN(im1,j,kp1) * dPSE(iim1,jjp1,kkp1)
     2   + oPNW(ii,jj,kk) * uC(im1,jp1,k) * dPSE(iim1,jjp1,kkp1)

     2   - uPNW(ii,jj,kk) * (- oE(im2,jp1,kp1) * dPS(iim1,jjp1,kkp1)
     2                       + oC(im1,jp1,kp1) * dPSE(iim1,jjp1,kkp1)
     2                       - uC(im1,jp1,kp1) * oPSE(iim1,jjp1,kkp1)
     2                       - oN(im1,jp1,kp1) * dPE(iim1,jjp1,kkp1))

     2   + uPN(ii,jj,kk) * oE(im1,jp1,kp1) * dPSE(iim1,jjp1,kkp1)

c* ********************************************************************
c* *** > uSE;
c* ********************************************************************

      XuSE(ii,jj,kk) =
     2     uPS(ii,jj,kk) * oE(i,jm1,kp1) * dPNW(iip1,jjm1,kkp1)
     2   + oPSE(ii,jj,kk) * uC(ip1,jm1,k) * dPNW(iip1,jjm1,kkp1)

     2   - uPSE(ii,jj,kk) * (- oN(ip1,jm2,kp1) * dPW(iip1,jjm1,kkp1)
     2                       + oC(ip1,jm1,kp1) * dPNW(iip1,jjm1,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPNW(iip1,jjm1,kkp1)
     2                       - oE(ip1,jm1,kp1) * dPN(iip1,jjm1,kkp1))

     2   + uPE(ii,jj,kk) * oN(ip1,jm1,kp1) * dPNW(iip1,jjm1,kkp1)

c* ********************************************************************
c* *** > uSW;
c* ********************************************************************

      XuSW(ii,jj,kk) =
     2     oPSW(ii,jj,kk) * uC(im1,jm1,k) * dPNE(iim1,jjm1,kkp1)

     2   - uPSW(ii,jj,kk) * (- oE(im2,jm1,kp1) * dPN(iim1,jjm1,kkp1)
     2                       - oN(im1,jm2,kp1) * dPE(iim1,jjm1,kkp1)
     2                       + oC(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)
     2                       - uC(im1,jm1,kp1) * oPNE(iim1,jjm1,kkp1))

     2   + uPW(ii,jj,kk) * oN(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)
     2   + uPS(ii,jj,kk) * oE(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)

c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildG_27 (nxf,nyf,nzf,nx,ny,nz,
     2   oPC,oPN,oPS,oPE,oPW,oPNE,oPNW,oPSE,oPSW,
     3   uPC,uPN,uPS,uPE,uPW,uPNE,uPNW,uPSE,uPSW,
     4   dPC,dPN,dPS,dPE,dPW,dPNE,dPNW,dPSE,dPSW,
     5   oC,oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     6   XoC,XoE,XoN,XuC,XoNE,XoNW,XuE,XuW,XuN,XuS,XuNE,XuNW,XuSE,XuSW)
c* ********************************************************************
c* purpose: compute a 27-point galerkin coarse grid matrix from 
c*          a 27-point fine grid matrix.
c*
c* expressions for the galerkin coarse grid stencil XA in terms of
c* the fine grid matrix stencil A and the interpolation operator
c* stencil P.  these stencils have the form:
c* 
c*    XA := array([
c*
c*      matrix([
c*         [ -XdNW(i,j,k), -XdN(i,j,k), -XdNE(i,j,k) ],
c*         [ -XdW(i,j,k),  -XdC(i,j,k), -XdE(i,j,k)  ],
c*         [ -XdSW(i,j,k), -XdS(i,j,k), -XdSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -XoNW(i,j,k), -XoN(i,j,k), -XoNE(i,j,k) ],
c*         [ -XoW(i,j,k),   XoC(i,j,k), -XoE(i,j,k)  ],
c*         [ -XoSW(i,j,k), -XoS(i,j,k), -XoSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -XuNW(i,j,k), -XuN(i,j,k), -XuNE(i,j,k) ],
c*         [ -XuW(i,j,k),  -XuC(i,j,k), -XuE(i,j,k)  ],
c*         [ -XuSW(i,j,k), -XuS(i,j,k), -XuSE(i,j,k) ] 
c*      ])
c*    ]):
c*
c*    A := array([
c*
c*      matrix([
c*         [ -dNW(i,j,k), -dN(i,j,k), -dNE(i,j,k) ],
c*         [ -dW(i,j,k),  -dC(i,j,k), -dE(i,j,k)  ],
c*         [ -dSW(i,j,k), -dS(i,j,k), -dSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -oNW(i,j,k), -oN(i,j,k), -oNE(i,j,k) ],
c*         [ -oW(i,j,k),   oC(i,j,k), -oE(i,j,k)  ],
c*         [ -oSW(i,j,k), -oS(i,j,k), -oSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ -uNW(i,j,k), -uN(i,j,k), -uNE(i,j,k) ],
c*         [ -uW(i,j,k),  -uC(i,j,k), -uE(i,j,k)  ],
c*         [ -uSW(i,j,k), -uS(i,j,k), -uSE(i,j,k) ] 
c*      ])
c*    ]):
c*
c*   P := array([
c*
c*      matrix([
c*         [ dPNW(i,j,k), dPN(i,j,k), dPNE(i,j,k) ],
c*         [ dPW(i,j,k),  dPC(i,j,k), dPE(i,j,k)  ],
c*         [ dPSW(i,j,k), dPS(i,j,k), dPSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ oPNW(i,j,k), oPN(i,j,k), oPNE(i,j,k) ],
c*         [ oPW(i,j,k),  oPC(i,j,k), oPE(i,j,k)  ],
c*         [ oPSW(i,j,k), oPS(i,j,k), oPSE(i,j,k) ] 
c*      ]),
c*
c*      matrix([
c*         [ uPNW(i,j,k), uPN(i,j,k), uPNE(i,j,k) ],
c*         [ uPW(i,j,k),  uPC(i,j,k), uPE(i,j,k)  ],
c*         [ uPSW(i,j,k), uPS(i,j,k), uPSE(i,j,k) ] 
c*      ])
c*    ]):
c*
c* in addition, A is assumed to be symmetric so that:
c*
c*    oS  := proc(x,y,z) RETURN( oN(x,y-1,z) ): end:
c*    oW  := proc(x,y,z) RETURN( oE(x-1,y,z) ): end:
c*    oSE := proc(x,y,z) RETURN( oNW(x+1,y-1,z) ): end:
c*    oSW := proc(x,y,z) RETURN( oNE(x-1,y-1,z) ): end:
c* 
c*    dC  := proc(x,y,z) RETURN( uC(x,y,z-1) ): end:
c*    dW  := proc(x,y,z) RETURN( uE(x-1,y,z-1) ): end:
c*    dE  := proc(x,y,z) RETURN( uW(x+1,y,z-1) ): end:
c* 
c*    dN  := proc(x,y,z) RETURN( uS(x,y+1,z-1) ): end:
c*    dNW := proc(x,y,z) RETURN( uSE(x-1,y+1,z-1) ): end:
c*    dNE := proc(x,y,z) RETURN( uSW(x+1,y+1,z-1) ): end:
c* 
c*    dS  := proc(x,y,z) RETURN( uN(x,y-1,z-1) ): end:
c*    dSW := proc(x,y,z) RETURN( uNE(x-1,y-1,z-1) ): end:
c*    dSE := proc(x,y,z) RETURN( uNW(x+1,y-1,z-1) ): end:
c* 
c* author:  michael holst
c* ********************************************************************
      implicit         none
      integer          nx,ny,nz,i,j,k,ii,jj,kk
      integer          nxm1,nym1,nzm1,nxf,nyf,nzf
      integer          im1,ip1,im2,ip2,jm1,jp1,jm2,jp2,km1,kp1,km2,kp2
      integer          iim1,iip1,jjm1,jjp1,kkm1,kkp1

      double precision TMP1_XoC,TMP2_XoC,TMP3_XoC,TMP4_XoC
      double precision TMP5_XoC,TMP6_XoC,TMP7_XoC,TMP8_XoC
      double precision TMP9_XoC,TMP10_XoC,TMP11_XoC,TMP12_XoC
      double precision TMP13_XoC,TMP14_XoC,TMP15_XoC,TMP16_XoC
      double precision TMP17_XoC,TMP18_XoC,TMP19_XoC,TMP20_XoC
      double precision TMP21_XoC,TMP22_XoC,TMP23_XoC,TMP24_XoC
      double precision TMP25_XoC,TMP26_XoC,TMP27_XoC
      double precision TMP1_XoE,TMP2_XoE,TMP3_XoE,TMP4_XoE
      double precision TMP5_XoE,TMP6_XoE,TMP7_XoE,TMP8_XoE
      double precision TMP9_XoE,TMP10_XoE,TMP11_XoE,TMP12_XoE
      double precision TMP1_XoN,TMP2_XoN,TMP3_XoN,TMP4_XoN
      double precision TMP5_XoN,TMP6_XoN,TMP7_XoN,TMP8_XoN
      double precision TMP9_XoN,TMP10_XoN,TMP11_XoN,TMP12_XoN
      double precision TMP1_XuC,TMP2_XuC,TMP3_XuC,TMP4_XuC
      double precision TMP5_XuC,TMP6_XuC,TMP7_XuC,TMP8_XuC
      double precision TMP9_XuC,TMP10_XuC,TMP11_XuC,TMP12_XuC
      double precision TMP1_XoNE,TMP2_XoNE,TMP3_XoNE,TMP4_XoNE
      double precision TMP5_XoNE,TMP6_XoNE
      double precision TMP1_XoNW,TMP2_XoNW,TMP3_XoNW,TMP4_XoNW
      double precision TMP5_XoNW,TMP6_XoNW
      double precision TMP1_XuE,TMP2_XuE,TMP3_XuE,TMP4_XuE
      double precision TMP5_XuE,TMP6_XuE
      double precision TMP1_XuW,TMP2_XuW,TMP3_XuW,TMP4_XuW
      double precision TMP5_XuW,TMP6_XuW
      double precision TMP1_XuN,TMP2_XuN,TMP3_XuN,TMP4_XuN
      double precision TMP5_XuN,TMP6_XuN
      double precision TMP1_XuS,TMP2_XuS,TMP3_XuS,TMP4_XuS
      double precision TMP5_XuS,TMP6_XuS
      double precision TMP1_XuNE,TMP2_XuNE,TMP1_XuNW,TMP2_XuNW
      double precision TMP1_XuSE,TMP2_XuSE,TMP1_XuSW,TMP2_XuSW

      double precision oC(nxf,nyf,nzf),oE(nxf,nyf,nzf)
      double precision oN(nxf,nyf,nzf),uC(nxf,nyf,nzf)
      double precision oNE(nxf,nyf,nzf),oNW(nxf,nyf,nzf)
      double precision uE(nxf,nyf,nzf),uW(nxf,nyf,nzf)
      double precision uN(nxf,nyf,nzf),uS(nxf,nyf,nzf)
      double precision uNE(nxf,nyf,nzf),uNW(nxf,nyf,nzf)
      double precision uSE(nxf,nyf,nzf),uSW(nxf,nyf,nzf)

      double precision XoC(nx,ny,nz),XoE(nx,ny,nz)
      double precision XoN(nx,ny,nz),XuC(nx,ny,nz)
      double precision XoNE(nx,ny,nz),XoNW(nx,ny,nz)
      double precision XuE(nx,ny,nz),XuW(nx,ny,nz)
      double precision XuN(nx,ny,nz),XuS(nx,ny,nz)
      double precision XuNE(nx,ny,nz),XuNW(nx,ny,nz)
      double precision XuSE(nx,ny,nz),XuSW(nx,ny,nz)

      double precision oPC(nx,ny,nz),oPN(nx,ny,nz),oPS(nx,ny,nz)
      double precision oPE(nx,ny,nz),oPW(nx,ny,nz),oPNE(nx,ny,nz)
      double precision oPNW(nx,ny,nz),oPSE(nx,ny,nz),oPSW(nx,ny,nz)
      double precision uPC(nx,ny,nz),uPN(nx,ny,nz),uPS(nx,ny,nz)
      double precision uPE(nx,ny,nz),uPW(nx,ny,nz),uPNE(nx,ny,nz)
      double precision uPNW(nx,ny,nz),uPSE(nx,ny,nz),uPSW(nx,ny,nz)
      double precision dPC(nx,ny,nz),dPN(nx,ny,nz),dPS(nx,ny,nz)
      double precision dPE(nx,ny,nz),dPW(nx,ny,nz),dPNE(nx,ny,nz)
      double precision dPNW(nx,ny,nz),dPSE(nx,ny,nz),dPSW(nx,ny,nz)
c*
cmdir 0 0
c*
c*    *** define n and determine number of mesh points ***
      nxm1    = nx - 1
      nym1    = ny - 1
      nzm1    = nz - 1
c*
c*    *** build the operator ***
cmdir 3 1
      do 10 kk = 2, nz-1
         k = 2 * kk - 1
cmdir 3 2
         do 11 jj = 2, ny-1
            j = 2 * jj - 1
cmdir 3 3
            do 12 ii = 2, nx-1
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
c* *** > oC;
c* ********************************************************************

c* ***XoC(ii,jj,kk) =
      TMP1_XoC =
     2     oPN(ii,jj,kk) * (- uNE(im1,j,km1) * dPW(ii,jj,kk)
     2                      - oNE(im1,j,k) * oPW(ii,jj,kk)
     2                      - uSW(i,jp1,k) * uPW(ii,jj,kk)
     2                      - uE(im1,jp1,km1) * dPNW(ii,jj,kk)
     2                      - oE(im1,jp1,k) * oPNW(ii,jj,kk)
     2                      - uW(i,jp1,k) * uPNW(ii,jj,kk)
     2                      - uN(i,j,km1) * dPC(ii,jj,kk)
     2                      - oN(i,j,k) * oPC(ii,jj,kk)
     2                      - uS(i,jp1,k) * uPC(ii,jj,kk)
     2                      - uC(i,jp1,km1) * dPN(ii,jj,kk)
     2                      + oC(i,jp1,k) * oPN(ii,jj,kk)
     2                      - uC(i,jp1,k) * uPN(ii,jj,kk)
     2                      - uNW(ip1,j,km1) * dPE(ii,jj,kk)
     2                      - oNW(ip1,j,k) * oPE(ii,jj,kk)
     2                      - uSE(i,jp1,k) * uPE(ii,jj,kk)
     2                      - uW(ip1,jp1,km1) * dPNE(ii,jj,kk)
     2                      - oE(i,jp1,k) * oPNE(ii,jj,kk)
     2                      - uE(i,jp1,k) * uPNE(ii,jj,kk)) 

      TMP2_XoC =
     2   + dPN(ii,jj,kk) * (- oNE(im1,j,km1) * dPW(ii,jj,kk)
     2                      - uSW(i,jp1,km1) * oPW(ii,jj,kk)
     2                      - oE(im1,jp1,km1) * dPNW(ii,jj,kk)
     2                      - uW(i,jp1,km1) * oPNW(ii,jj,kk) 
     2                      - oN(i,j,km1) * dPC(ii,jj,kk)
     2                      - uS(i,jp1,km1) * oPC(ii,jj,kk) 
     2                      + oC(i,jp1,km1) * dPN(ii,jj,kk)
     2                      - uC(i,jp1,km1) * oPN(ii,jj,kk) 
     2                      - oNW(ip1,j,km1) * dPE(ii,jj,kk)
     2                      - uSE(i,jp1,km1) * oPE(ii,jj,kk) 
     2                      - oE(i,jp1,km1) * dPNE(ii,jj,kk)
     2                      - uE(i,jp1,km1) * oPNE(ii,jj,kk)) 

      TMP3_XoC =
     2   + dPC(ii,jj,kk) * (- oNE(im1,jm1,km1) * dPSW(ii,jj,kk) 
     2                      - uSW(i,j,km1) * oPSW(ii,jj,kk)
     2                      - oE(im1,j,km1) * dPW(ii,jj,kk) 
     2                      - uW(i,j,km1) * oPW(ii,jj,kk)
     2                      - oNW(i,j,km1) * dPNW(ii,jj,kk) 
     2                      - uNW(i,j,km1) * oPNW(ii,jj,kk)
     2                      - oN(i,jm1,km1) * dPS(ii,jj,kk) 
     2                      - uS(i,j,km1) * oPS(ii,jj,kk)
     2                      + oC(i,j,km1) * dPC(ii,jj,kk) 
     2                      - uC(i,j,km1) * oPC(ii,jj,kk)
     2                      - oN(i,j,km1) * dPN(ii,jj,kk) 
     2                      - uN(i,j,km1) * oPN(ii,jj,kk)
     2                      - oNW(ip1,jm1,km1) * dPSE(ii,jj,kk) 
     2                      - uSE(i,j,km1) * oPSE(ii,jj,kk)
     2                      - oE(i,j,km1) * dPE(ii,jj,kk) 
     2                      - uE(i,j,km1) * oPE(ii,jj,kk)
     2                      - oNE(i,j,km1) * dPNE(ii,jj,kk) 
     2                      - uNE(i,j,km1) * oPNE(ii,jj,kk)) 

      TMP4_XoC =
     2   + uPC(ii,jj,kk) * (- uNE(im1,jm1,k) * oPSW(ii,jj,kk)
     2                      - oNE(im1,jm1,kp1) * uPSW(ii,jj,kk) 
     2                      - uE(im1,j,k) * oPW(ii,jj,kk)
     2                      - oE(im1,j,kp1) * uPW(ii,jj,kk) 
     2                      - uSE(im1,jp1,k) * oPNW(ii,jj,kk)
     2                      - oNW(i,j,kp1) * uPNW(ii,jj,kk) 
     2                      - uN(i,jm1,k) * oPS(ii,jj,kk)
     2                      - oN(i,jm1,kp1) * uPS(ii,jj,kk) 
     2                      - uC(i,j,k) * oPC(ii,jj,kk)
     2                      + oC(i,j,kp1) * uPC(ii,jj,kk) 
     2                      - uS(i,jp1,k) * oPN(ii,jj,kk)
     2                      - oN(i,j,kp1) * uPN(ii,jj,kk) 
     2                      - uNW(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                      - oNW(ip1,jm1,kp1) * uPSE(ii,jj,kk) 
     2                      - uW(ip1,j,k) * oPE(ii,jj,kk)
     2                      - oE(i,j,kp1) * uPE(ii,jj,kk) 
     2                      - uSW(ip1,jp1,k) * oPNE(ii,jj,kk)
     2                      - oNE(i,j,kp1) * uPNE(ii,jj,kk)) 

      TMP5_XoC =
     2   + oPC(ii,jj,kk) * (- uW(ip1,j,km1) * dPE(ii,jj,kk) 
     2                      - oE(im1,j,k) * oPW(ii,jj,kk)
     2                      - uSE(im1,jp1,km1) * dPNW(ii,jj,kk)
     2                      - uNE(im1,jm1,km1) * dPSW(ii,jj,kk)
     2                      - uN(i,jm1,km1) * dPS(ii,jj,kk) 
     2                      - oNE(im1,jm1,k) * oPSW(ii,jj,kk)
     2                      - uE(im1,j,km1) * dPW(ii,jj,kk) 
     2                      - oNW(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                      - uC(i,j,km1) * dPC(ii,jj,kk) 
     2                      - uNW(ip1,jm1,km1) * dPSE(ii,jj,kk) 
     2                      - uSW(ip1,jp1,km1) * dPNE(ii,jj,kk)
     2 - uS(i,jp1,km1) * dPN(ii,jj,kk) - oN(i,jm1,k) * oPS(ii,jj,kk)
     2 - uNE(i,j,k) * uPNE(ii,jj,kk) - oNE(i,j,k) * oPNE(ii,jj,kk)
     2 - uE(i,j,k) * uPE(ii,jj,kk) - uSE(i,j,k) * uPSE(ii,jj,kk)
     2 - oN(i,j,k) * oPN(ii,jj,kk) - oE(i,j,k) * oPE(ii,jj,kk)
     2 - uS(i,j,k) * uPS(ii,jj,kk) + oC(i,j,k) * oPC(ii,jj,kk)
     2 - uSW(i,j,k) * uPSW(ii,jj,kk) - uN(i,j,k) * uPN(ii,jj,kk)
     2 - uC(i,j,k) * uPC(ii,jj,kk) - uW(i,j,k) * uPW(ii,jj,kk)
     2 - oNW(i,j,k) * oPNW(ii,jj,kk) - uNW(i,j,k) * uPNW(ii,jj,kk)) 

      TMP6_XoC =
     2   + uPS(ii,jj,kk) * (- uE(im1,jm1,k) * oPSW(ii,jj,kk)
     2                      - oE(im1,jm1,kp1) * uPSW(ii,jj,kk) 
     2                      - uSE(im1,j,k) * oPW(ii,jj,kk)
     2                      - oNW(i,jm1,kp1) * uPW(ii,jj,kk) 
     2                      - uC(i,jm1,k) * oPS(ii,jj,kk)
     2                      + oC(i,jm1,kp1) * uPS(ii,jj,kk) 
     2                      - uS(i,j,k) * oPC(ii,jj,kk)
     2                      - oN(i,jm1,kp1) * uPC(ii,jj,kk) 
     2                      - uW(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                      - oE(i,jm1,kp1) * uPSE(ii,jj,kk) 
     2                      - uSW(ip1,j,k) * oPE(ii,jj,kk)
     2                      - oNE(i,jm1,kp1) * uPE(ii,jj,kk)) 

      TMP7_XoC =
     2   + oPS(ii,jj,kk) * (- uE(im1,jm1,km1) * dPSW(ii,jj,kk) 
     2                      - oE(im1,jm1,k) * oPSW(ii,jj,kk)
     2                      - uW(i,jm1,k) * uPSW(ii,jj,kk) 
     2                      - uSE(im1,j,km1) * dPW(ii,jj,kk)
     2                      - oNW(i,jm1,k) * oPW(ii,jj,kk) 
     2                      - uNW(i,jm1,k) * uPW(ii,jj,kk)
     2                      - uC(i,jm1,km1) * dPS(ii,jj,kk) 
     2                      + oC(i,jm1,k) * oPS(ii,jj,kk)
     2                      - uC(i,jm1,k) * uPS(ii,jj,kk) 
     2                      - uS(i,j,km1) * dPC(ii,jj,kk)
     2                      - oN(i,jm1,k) * oPC(ii,jj,kk) 
     2                      - uN(i,jm1,k) * uPC(ii,jj,kk)
     2                      - uW(ip1,jm1,km1) * dPSE(ii,jj,kk) 
     2                      - oE(i,jm1,k) * oPSE(ii,jj,kk)
     2                      - uE(i,jm1,k) * uPSE(ii,jj,kk) 
     2                      - uSW(ip1,j,km1) * dPE(ii,jj,kk)
     2                      - oNE(i,jm1,k) * oPE(ii,jj,kk) 
     2                      - uNE(i,jm1,k) * uPE(ii,jj,kk)) 

      TMP8_XoC =
     2   + dPS(ii,jj,kk) * (- oE(im1,jm1,km1) * dPSW(ii,jj,kk)
     2                      - uW(i,jm1,km1) * oPSW(ii,jj,kk) 
     2                      - oNW(i,jm1,km1) * dPW(ii,jj,kk)
     2                      - uNW(i,jm1,km1) * oPW(ii,jj,kk) 
     2                      + oC(i,jm1,km1) * dPS(ii,jj,kk)
     2                      - uC(i,jm1,km1) * oPS(ii,jj,kk) 
     2                      - oN(i,jm1,km1) * dPC(ii,jj,kk)
     2                      - uN(i,jm1,km1) * oPC(ii,jj,kk) 
     2                      - oE(i,jm1,km1) * dPSE(ii,jj,kk)
     2                      - uE(i,jm1,km1) * oPSE(ii,jj,kk) 
     2                      - oNE(i,jm1,km1) * dPE(ii,jj,kk)
     2                      - uNE(i,jm1,km1) * oPE(ii,jj,kk)) 

      TMP9_XoC =
     2   + uPNW(ii,jj,kk) * (- uN(im1,j,k) * oPW(ii,jj,kk) 
     2                       - oN(im1,j,kp1) * uPW(ii,jj,kk)
     2                       - uC(im1,jp1,k) * oPNW(ii,jj,kk)
     2                       + oC(im1,jp1,kp1) * uPNW(ii,jj,kk) 
     2                       - uNW(i,j,k) * oPC(ii,jj,kk)
     2                       - oNW(i,j,kp1) * uPC(ii,jj,kk) 
     2                       - uW(i,jp1,k) * oPN(ii,jj,kk)
     2                       - oE(im1,jp1,kp1) * uPN(ii,jj,kk)) 

      TMP10_XoC =
     2   + oPNW(ii,jj,kk) * (- uN(im1,j,km1) * dPW(ii,jj,kk) 
     2                       - oN(im1,j,k) * oPW(ii,jj,kk)
     2                       - uS(im1,jp1,k) * uPW(ii,jj,kk) 
     2                       - uC(im1,jp1,km1) * dPNW(ii,jj,kk)
     2                       + oC(im1,jp1,k) * oPNW(ii,jj,kk) 
     2                       - uC(im1,jp1,k) * uPNW(ii,jj,kk)
     2                       - uNW(i,j,km1) * dPC(ii,jj,kk) 
     2                       - oNW(i,j,k) * oPC(ii,jj,kk)
     2                       - uSE(im1,jp1,k) * uPC(ii,jj,kk) 
     2                       - uW(i,jp1,km1) * dPN(ii,jj,kk)
     2                       - oE(im1,jp1,k) * oPN(ii,jj,kk) 
     2                       - uE(im1,jp1,k) * uPN(ii,jj,kk)) 

      TMP11_XoC =
     2   + uPW(ii,jj,kk) * (- uN(im1,jm1,k) * oPSW(ii,jj,kk)
     2                      - oN(im1,jm1,kp1) * uPSW(ii,jj,kk) 
     2                      - uC(im1,j,k) * oPW(ii,jj,kk)
     2                      + oC(im1,j,kp1) * uPW(ii,jj,kk) 
     2                      - uS(im1,jp1,k) * oPNW(ii,jj,kk)
     2                      - oN(im1,j,kp1) * uPNW(ii,jj,kk) 
     2                      - uNW(i,jm1,k) * oPS(ii,jj,kk)
     2                      - oNW(i,jm1,kp1) * uPS(ii,jj,kk) 
     2                      - uW(i,j,k) * oPC(ii,jj,kk)
     2                      - oE(im1,j,kp1) * uPC(ii,jj,kk) 
     2                      - uSW(i,jp1,k) * oPN(ii,jj,kk)
     2                      - oNE(im1,j,kp1) * uPN(ii,jj,kk)) 

      TMP12_XoC =
     2   + dPNW(ii,jj,kk) * (- oN(im1,j,km1) * dPW(ii,jj,kk) 
     2                       - uS(im1,jp1,km1) * oPW(ii,jj,kk)
     2                       + oC(im1,jp1,km1) * dPNW(ii,jj,kk)
     2                       - uC(im1,jp1,km1) * oPNW(ii,jj,kk) 
     2                       - oNW(i,j,km1) * dPC(ii,jj,kk)
     2                       - uSE(im1,jp1,km1) * oPC(ii,jj,kk)
     2                       - oE(im1,jp1,km1) * dPN(ii,jj,kk)
     2                       - uE(im1,jp1,km1) * oPN(ii,jj,kk)) 

      TMP13_XoC =
     2   + oPW(ii,jj,kk) * (- uN(im1,jm1,km1) * dPSW(ii,jj,kk) 
     2                      - oN(im1,jm1,k) * oPSW(ii,jj,kk)
     2                      - uS(im1,j,k) * uPSW(ii,jj,kk) 
     2                      - uC(im1,j,km1) * dPW(ii,jj,kk)
     2                      + oC(im1,j,k) * oPW(ii,jj,kk) 
     2                      - uC(im1,j,k) * uPW(ii,jj,kk)
     2                      - uS(im1,jp1,km1) * dPNW(ii,jj,kk) 
     2                      - oN(im1,j,k) * oPNW(ii,jj,kk)
     2                      - uN(im1,j,k) * uPNW(ii,jj,kk) 
     2                      - uNW(i,jm1,km1) * dPS(ii,jj,kk)
     2                      - oNW(i,jm1,k) * oPS(ii,jj,kk) 
     2                      - uSE(im1,j,k) * uPS(ii,jj,kk)
     2                      - uW(i,j,km1) * dPC(ii,jj,kk) 
     2                      - oE(im1,j,k) * oPC(ii,jj,kk)
     2                      - uE(im1,j,k) * uPC(ii,jj,kk) 
     2                      - uSW(i,jp1,km1) * dPN(ii,jj,kk)
     2                      - oNE(im1,j,k) * oPN(ii,jj,kk) 
     2                      - uNE(im1,j,k) * uPN(ii,jj,kk)) 

      TMP14_XoC =
     2   + uPSW(ii,jj,kk) * (- uC(im1,jm1,k) * oPSW(ii,jj,kk)
     2                       + oC(im1,jm1,kp1) * uPSW(ii,jj,kk) 
     2                       - uS(im1,j,k) * oPW(ii,jj,kk)
     2                       - oN(im1,jm1,kp1) * uPW(ii,jj,kk) 
     2                       - uW(i,jm1,k) * oPS(ii,jj,kk)
     2                       - oE(im1,jm1,kp1) * uPS(ii,jj,kk) 
     2                       - uSW(i,j,k) * oPC(ii,jj,kk)
     2                       - oNE(im1,jm1,kp1) * uPC(ii,jj,kk)) 

      TMP15_XoC =
     2   + oPSW(ii,jj,kk) * (- uC(im1,jm1,km1) * dPSW(ii,jj,kk) 
     2                       + oC(im1,jm1,k) * oPSW(ii,jj,kk)
     2                       - uC(im1,jm1,k) * uPSW(ii,jj,kk) 
     2                       - uS(im1,j,km1) * dPW(ii,jj,kk)
     2                       - oN(im1,jm1,k) * oPW(ii,jj,kk) 
     2                       - uN(im1,jm1,k) * uPW(ii,jj,kk)
     2                       - uW(i,jm1,km1) * dPS(ii,jj,kk) 
     2                       - oE(im1,jm1,k) * oPS(ii,jj,kk)
     2                       - uE(im1,jm1,k) * uPS(ii,jj,kk) 
     2                       - uSW(i,j,km1) * dPC(ii,jj,kk)
     2                       - oNE(im1,jm1,k) * oPC(ii,jj,kk) 
     2                       - uNE(im1,jm1,k) * uPC(ii,jj,kk))

      TMP16_XoC =
     2   + dPW(ii,jj,kk) * (- oN(im1,jm1,km1) * dPSW(ii,jj,kk)
     2                      - uS(im1,j,km1) * oPSW(ii,jj,kk) 
     2                      + oC(im1,j,km1) * dPW(ii,jj,kk)
     2                      - uC(im1,j,km1) * oPW(ii,jj,kk) 
     2                      - oN(im1,j,km1) * dPNW(ii,jj,kk)
     2                      - uN(im1,j,km1) * oPNW(ii,jj,kk) 
     2                      - oNW(i,jm1,km1) * dPS(ii,jj,kk)
     2                      - uSE(im1,j,km1) * oPS(ii,jj,kk) 
     2                      - oE(im1,j,km1) * dPC(ii,jj,kk)
     2                      - uE(im1,j,km1) * oPC(ii,jj,kk) 
     2                      - oNE(im1,j,km1) * dPN(ii,jj,kk)
     2                      - uNE(im1,j,km1) * oPN(ii,jj,kk)) 

      TMP17_XoC =
     2   + uPNE(ii,jj,kk) * (- uNE(i,j,k) * oPC(ii,jj,kk) 
     2                       - oNE(i,j,kp1) * uPC(ii,jj,kk)
     2                       - uE(i,jp1,k) * oPN(ii,jj,kk) 
     2                       - oE(i,jp1,kp1) * uPN(ii,jj,kk)
     2                       - uN(ip1,j,k) * oPE(ii,jj,kk) 
     2                       - oN(ip1,j,kp1) * uPE(ii,jj,kk)
     2                       - uC(ip1,jp1,k) * oPNE(ii,jj,kk)
     2                       + oC(ip1,jp1,kp1) * uPNE(ii,jj,kk)) 

      TMP18_XoC =
     2   + uPE(ii,jj,kk) * (- uNE(i,jm1,k) * oPS(ii,jj,kk) 
     2                      - oNE(i,jm1,kp1) * uPS(ii,jj,kk)
     2                      - uE(i,j,k) * oPC(ii,jj,kk) 
     2                      - oE(i,j,kp1) * uPC(ii,jj,kk)
     2                      - uSE(i,jp1,k) * oPN(ii,jj,kk) 
     2                      - oNW(ip1,j,kp1) * uPN(ii,jj,kk)
     2                      - uN(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                      - oN(ip1,jm1,kp1) * uPSE(ii,jj,kk) 
     2                      - uC(ip1,j,k) * oPE(ii,jj,kk)
     2                      + oC(ip1,j,kp1) * uPE(ii,jj,kk) 
     2                      - uS(ip1,jp1,k) * oPNE(ii,jj,kk)
     2                      - oN(ip1,j,kp1) * uPNE(ii,jj,kk)) 

      TMP19_XoC =
     2   + dPNE(ii,jj,kk) * (- oNE(i,j,km1) * dPC(ii,jj,kk) 
     2                       - uSW(ip1,jp1,km1) * oPC(ii,jj,kk)
     2                       - oE(i,jp1,km1) * dPN(ii,jj,kk) 
     2                       - uW(ip1,jp1,km1) * oPN(ii,jj,kk)
     2                       - oN(ip1,j,km1) * dPE(ii,jj,kk) 
     2                       - uS(ip1,jp1,km1) * oPE(ii,jj,kk)
     2                       + oC(ip1,jp1,km1) * dPNE(ii,jj,kk)
     2                       - uC(ip1,jp1,km1) * oPNE(ii,jj,kk)) 

      TMP20_XoC =
     2   + oPNE(ii,jj,kk) * (- uNE(i,j,km1) * dPC(ii,jj,kk) 
     2                       - oNE(i,j,k) * oPC(ii,jj,kk)
     2                       - uSW(ip1,jp1,k) * uPC(ii,jj,kk) 
     2                       - uE(i,jp1,km1) * dPN(ii,jj,kk)
     2                       - oE(i,jp1,k) * oPN(ii,jj,kk) 
     2                       - uW(ip1,jp1,k) * uPN(ii,jj,kk)
     2                       - uN(ip1,j,km1) * dPE(ii,jj,kk) 
     2                       - oN(ip1,j,k) * oPE(ii,jj,kk)
     2                       - uS(ip1,jp1,k) * uPE(ii,jj,kk) 
     2                       - uC(ip1,jp1,km1) * dPNE(ii,jj,kk)
     2                       + oC(ip1,jp1,k) * oPNE(ii,jj,kk) 
     2                       - uC(ip1,jp1,k) * uPNE(ii,jj,kk))

      TMP21_XoC =
     2   + oPSE(ii,jj,kk) * (- uE(i,jm1,km1) * dPS(ii,jj,kk)
     2                       - oE(i,jm1,k) * oPS(ii,jj,kk) 
     2                       - uW(ip1,jm1,k) * uPS(ii,jj,kk)
     2                       - uSE(i,j,km1) * dPC(ii,jj,kk) 
     2                       - oNW(ip1,jm1,k) * oPC(ii,jj,kk)
     2                       - uNW(ip1,jm1,k) * uPC(ii,jj,kk)
     2                       - uC(ip1,jm1,km1) * dPSE(ii,jj,kk)
     2                       + oC(ip1,jm1,k) * oPSE(ii,jj,kk) 
     2                       - uC(ip1,jm1,k) * uPSE(ii,jj,kk)
     2                       - uS(ip1,j,km1) * dPE(ii,jj,kk) 
     2                       - oN(ip1,jm1,k) * oPE(ii,jj,kk)
     2                       - uN(ip1,jm1,k) * uPE(ii,jj,kk)) 

      TMP22_XoC =
     2   + dPSE(ii,jj,kk) * (- oE(i,jm1,km1) * dPS(ii,jj,kk) 
     2                       - uW(ip1,jm1,km1) * oPS(ii,jj,kk)
     2                       - oNW(ip1,jm1,km1) * dPC(ii,jj,kk)
     2                       - uNW(ip1,jm1,km1) * oPC(ii,jj,kk)
     2                       + oC(ip1,jm1,km1) * dPSE(ii,jj,kk)
     2                       - uC(ip1,jm1,km1) * oPSE(ii,jj,kk)
     2                       - oN(ip1,jm1,km1) * dPE(ii,jj,kk)
     2                       - uN(ip1,jm1,km1) * oPE(ii,jj,kk)) 

      TMP23_XoC =
     2   + uPSE(ii,jj,kk) * (- uE(i,jm1,k) * oPS(ii,jj,kk) 
     2                       - oE(i,jm1,kp1) * uPS(ii,jj,kk)
     2                       - uSE(i,j,k) * oPC(ii,jj,kk) 
     2                       - oNW(ip1,jm1,kp1) * uPC(ii,jj,kk)
     2                       - uC(ip1,jm1,k) * oPSE(ii,jj,kk)
     2                       + oC(ip1,jm1,kp1) * uPSE(ii,jj,kk) 
     2                       - uS(ip1,j,k) * oPE(ii,jj,kk)
     2                       - oN(ip1,jm1,kp1) * uPE(ii,jj,kk)) 

      TMP24_XoC =
     2   + oPE(ii,jj,kk) * (- uNE(i,jm1,km1) * dPS(ii,jj,kk) 
     2                      - oNE(i,jm1,k) * oPS(ii,jj,kk)
     2                      - uSW(ip1,j,k) * uPS(ii,jj,kk) 
     2                      - uE(i,j,km1) * dPC(ii,jj,kk)
     2                      - oE(i,j,k) * oPC(ii,jj,kk) 
     2                      - uW(ip1,j,k) * uPC(ii,jj,kk)
     2                      - uSE(i,jp1,km1) * dPN(ii,jj,kk) 
     2                      - oNW(ip1,j,k) * oPN(ii,jj,kk)
     2                      - uNW(ip1,j,k) * uPN(ii,jj,kk) 
     2                      - uN(ip1,jm1,km1) * dPSE(ii,jj,kk)
     2                      - oN(ip1,jm1,k) * oPSE(ii,jj,kk) 
     2                      - uS(ip1,j,k) * uPSE(ii,jj,kk)
     2                      - uC(ip1,j,km1) * dPE(ii,jj,kk) 
     2                      + oC(ip1,j,k) * oPE(ii,jj,kk)
     2                      - uC(ip1,j,k) * uPE(ii,jj,kk) 
     2                      - uS(ip1,jp1,km1) * dPNE(ii,jj,kk)
     2                      - oN(ip1,j,k) * oPNE(ii,jj,kk) 
     2                      - uN(ip1,j,k) * uPNE(ii,jj,kk)) 

      TMP25_XoC =
     2   + dPE(ii,jj,kk) * (- oNE(i,jm1,km1) * dPS(ii,jj,kk)
     2                      - uSW(ip1,j,km1) * oPS(ii,jj,kk) 
     2                      - oE(i,j,km1) * dPC(ii,jj,kk)
     2                      - uW(ip1,j,km1) * oPC(ii,jj,kk) 
     2                      - oNW(ip1,j,km1) * dPN(ii,jj,kk)
     2                      - uNW(ip1,j,km1) * oPN(ii,jj,kk)
     2                      - oN(ip1,jm1,km1) * dPSE(ii,jj,kk)
     2                      - uS(ip1,j,km1) * oPSE(ii,jj,kk) 
     2                      + oC(ip1,j,km1) * dPE(ii,jj,kk)
     2                      - uC(ip1,j,km1) * oPE(ii,jj,kk) 
     2                      - oN(ip1,j,km1) * dPNE(ii,jj,kk)
     2                      - uN(ip1,j,km1) * oPNE(ii,jj,kk)) 

      TMP26_XoC =
     2   + uPN(ii,jj,kk) * (- uNE(im1,j,k) * oPW(ii,jj,kk) 
     2                      - oNE(im1,j,kp1) * uPW(ii,jj,kk)
     2                      - uE(im1,jp1,k) * oPNW(ii,jj,kk)
     2                      - oE(im1,jp1,kp1) * uPNW(ii,jj,kk) 
     2                      - uN(i,j,k) * oPC(ii,jj,kk)
     2                      - oN(i,j,kp1) * uPC(ii,jj,kk) 
     2                      - uC(i,jp1,k) * oPN(ii,jj,kk)
     2                      + oC(i,jp1,kp1) * uPN(ii,jj,kk) 
     2                      - uNW(ip1,j,k) * oPE(ii,jj,kk)
     2                      - oNW(ip1,j,kp1) * uPE(ii,jj,kk) 
     2                      - uW(ip1,jp1,k) * oPNE(ii,jj,kk)
     2                      - oE(i,jp1,kp1) * uPNE(ii,jj,kk)) 

      TMP27_XoC =
     2   + dPSW(ii,jj,kk) * ( oC(im1,jm1,km1) * dPSW(ii,jj,kk)
     2                      - uC(im1,jm1,km1) * oPSW(ii,jj,kk)
     2                      - oN(im1,jm1,km1) * dPW(ii,jj,kk)
     2                      - uN(im1,jm1,km1) * oPW(ii,jj,kk)
     2                      - oE(im1,jm1,km1) * dPS(ii,jj,kk)
     2                      - uE(im1,jm1,km1) * oPS(ii,jj,kk)
     2                      - oNE(im1,jm1,km1) * dPC(ii,jj,kk)
     2                      - uNE(im1,jm1,km1) * oPC(ii,jj,kk))

      XoC(ii,jj,kk) = TMP1_XoC + TMP2_XoC + TMP3_XoC + TMP4_XoC
     2   + TMP5_XoC + TMP6_XoC + TMP7_XoC + TMP8_XoC + TMP9_XoC 
     3   + TMP10_XoC + TMP11_XoC + TMP12_XoC + TMP13_XoC + TMP14_XoC 
     4   + TMP15_XoC + TMP16_XoC + TMP17_XoC + TMP18_XoC + TMP19_XoC 
     5   + TMP20_XoC + TMP21_XoC + TMP22_XoC + TMP23_XoC + TMP24_XoC 
     6   + TMP25_XoC + TMP26_XoC + TMP27_XoC

c* ********************************************************************
c* *** > oE;
c* ********************************************************************

c* ***XoE(ii,jj,kk) =
      TMP1_XoE =
     2   - dPS(ii,jj,kk) * (- oE(i,jm1,km1) * dPSW(iip1,jj,kk)
     2                      - uE(i,jm1,km1) * oPSW(iip1,jj,kk)
     2                      - oNE(i,jm1,km1) * dPW(iip1,jj,kk)
     2                      - uNE(i,jm1,km1) * oPW(iip1,jj,kk)) 

     2   - oPS(ii,jj,kk) * (- uW(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                      - oE(i,jm1,k) * oPSW(iip1,jj,kk) 
     2                      - uE(i,jm1,k) * uPSW(iip1,jj,kk)
     2                      - uSW(ip1,j,km1) * dPW(iip1,jj,kk)
     2                      - oNE(i,jm1,k) * oPW(iip1,jj,kk) 
     2                      - uNE(i,jm1,k) * uPW(iip1,jj,kk)) 

     2   - uPS(ii,jj,kk) * (- uW(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                      - oE(i,jm1,kp1) * uPSW(iip1,jj,kk)
     2                      - uSW(ip1,j,k) * oPW(iip1,jj,kk)
     2                      - oNE(i,jm1,kp1) * uPW(iip1,jj,kk)) 

      TMP2_XoE =
     2   - dPC(ii,jj,kk) * (- oNW(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                      - uSE(i,j,km1) * oPSW(iip1,jj,kk) 
     2                      - oE(i,j,km1) * dPW(iip1,jj,kk)
     2                      - uE(i,j,km1) * oPW(iip1,jj,kk) 
     2                      - oNE(i,j,km1) * dPNW(iip1,jj,kk)
     2                      - uNE(i,j,km1) * oPNW(iip1,jj,kk)) 

     2   - oPC(ii,jj,kk) * (- uNW(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                      - oNW(ip1,jm1,k) * oPSW(iip1,jj,kk) 
     2                      - uSE(i,j,k) * uPSW(iip1,jj,kk)
     2                      - uW(ip1,j,km1) * dPW(iip1,jj,kk) 
     2                      - oE(i,j,k) * oPW(iip1,jj,kk)
     2                      - uE(i,j,k) * uPW(iip1,jj,kk)
     2                      - uSW(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2                      - oNE(i,j,k) * oPNW(iip1,jj,kk) 
     2                      - uNE(i,j,k) * uPNW(iip1,jj,kk)) 

      TMP3_XoE =
     2   - uPC(ii,jj,kk) * (- uNW(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                      - oNW(ip1,jm1,kp1) * uPSW(iip1,jj,kk)
     2                      - uW(ip1,j,k) * oPW(iip1,jj,kk) 
     2                      - oE(i,j,kp1) * uPW(iip1,jj,kk)
     2                      - uSW(ip1,jp1,k) * oPNW(iip1,jj,kk)
     2                      - oNE(i,j,kp1) * uPNW(iip1,jj,kk)) 

     2   - dPN(ii,jj,kk) * (- oNW(ip1,j,km1) * dPW(iip1,jj,kk)
     2                      - uSE(i,jp1,km1) * oPW(iip1,jj,kk)
     2                      - oE(i,jp1,km1) * dPNW(iip1,jj,kk)
     2                      - uE(i,jp1,km1) * oPNW(iip1,jj,kk)) 

     2   - oPN(ii,jj,kk) * (- uNW(ip1,j,km1) * dPW(iip1,jj,kk) 
     2                      - oNW(ip1,j,k) * oPW(iip1,jj,kk)
     2                      - uSE(i,jp1,k) * uPW(iip1,jj,kk)
     2                      - uW(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2                      - oE(i,jp1,k) * oPNW(iip1,jj,kk) 
     2                      - uE(i,jp1,k) * uPNW(iip1,jj,kk)) 

      TMP4_XoE =
     2   - uPN(ii,jj,kk) * (- uNW(ip1,j,k) * oPW(iip1,jj,kk)
     2                      - oNW(ip1,j,kp1) * uPW(iip1,jj,kk)
     2                      - uW(ip1,jp1,k) * oPNW(iip1,jj,kk)
     2                      - oE(i,jp1,kp1) * uPNW(iip1,jj,kk)) 

     2   - dPSE(ii,jj,kk) * (  oC(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                       - uC(ip1,jm1,km1) * oPSW(iip1,jj,kk)
     2                       - oN(ip1,jm1,km1) * dPW(iip1,jj,kk)
     2                       - uN(ip1,jm1,km1) * oPW(iip1,jj,kk)
     2                       - oE(ip1,jm1,km1) * dPS(iip1,jj,kk)
     2                       - uE(ip1,jm1,km1) * oPS(iip1,jj,kk)
     2                       - oNE(ip1,jm1,km1) * dPC(iip1,jj,kk)
     2                       - uNE(ip1,jm1,km1) * oPC(iip1,jj,kk)) 

      TMP5_XoE =
     2   - oPSE(ii,jj,kk) * (- uC(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                       + oC(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                       - uC(ip1,jm1,k) * uPSW(iip1,jj,kk)
     2                       - uS(ip1,j,km1) * dPW(iip1,jj,kk)
     2                       - oN(ip1,jm1,k) * oPW(iip1,jj,kk)
     2                       - uN(ip1,jm1,k) * uPW(iip1,jj,kk)
     2                       - uW(ip2,jm1,km1) * dPS(iip1,jj,kk)
     2                       - oE(ip1,jm1,k) * oPS(iip1,jj,kk)
     2                       - uE(ip1,jm1,k) * uPS(iip1,jj,kk)
     2                       - uSW(ip2,j,km1) * dPC(iip1,jj,kk)
     2                       - oNE(ip1,jm1,k) * oPC(iip1,jj,kk)
     2                       - uNE(ip1,jm1,k) * uPC(iip1,jj,kk)) 

      TMP6_XoE =
     2   - uPSE(ii,jj,kk) * (- uC(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                       + oC(ip1,jm1,kp1) * uPSW(iip1,jj,kk)
     2                       - uS(ip1,j,k) * oPW(iip1,jj,kk)
     2                       - oN(ip1,jm1,kp1) * uPW(iip1,jj,kk)
     2                       - uW(ip2,jm1,k) * oPS(iip1,jj,kk)
     2                       - oE(ip1,jm1,kp1) * uPS(iip1,jj,kk)
     2                       - uSW(ip2,j,k) * oPC(iip1,jj,kk)
     2                       - oNE(ip1,jm1,kp1) * uPC(iip1,jj,kk)) 

      TMP7_XoE =
     2   - dPE(ii,jj,kk) * (- oN(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                      - uS(ip1,j,km1) * oPSW(iip1,jj,kk)
     2                      + oC(ip1,j,km1) * dPW(iip1,jj,kk)
     2                      - uC(ip1,j,km1) * oPW(iip1,jj,kk)
     2                      - oN(ip1,j,km1) * dPNW(iip1,jj,kk)
     2                      - uN(ip1,j,km1) * oPNW(iip1,jj,kk)
     2                      - oNW(ip2,jm1,km1) * dPS(iip1,jj,kk)
     2                      - uSE(ip1,j,km1) * oPS(iip1,jj,kk)
     2                      - oE(ip1,j,km1) * dPC(iip1,jj,kk)
     2                      - uE(ip1,j,km1) * oPC(iip1,jj,kk)
     2                      - oNE(ip1,j,km1) * dPN(iip1,jj,kk)
     2                      - uNE(ip1,j,km1) * oPN(iip1,jj,kk)) 

      TMP8_XoE =
     2   - oPE(ii,jj,kk) * (- uN(ip1,jm1,km1) * dPSW(iip1,jj,kk)
     2                      - oN(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                      - uS(ip1,j,k) * uPSW(iip1,jj,kk) 
     2                      - uC(ip1,j,km1) * dPW(iip1,jj,kk)
     2                      + oC(ip1,j,k) * oPW(iip1,jj,kk) 
     2                      - uC(ip1,j,k) * uPW(iip1,jj,kk)
     2                      - uS(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2                      - oN(ip1,j,k) * oPNW(iip1,jj,kk) 
     2                      - uN(ip1,j,k) * uPNW(iip1,jj,kk)
     2                      - uNW(ip2,jm1,km1) * dPS(iip1,jj,kk)
     2                      - oNW(ip2,jm1,k) * oPS(iip1,jj,kk)
     2                      - uSE(ip1,j,k) * uPS(iip1,jj,kk) 
     2                      - uW(ip2,j,km1) * dPC(iip1,jj,kk)
     2                      - oE(ip1,j,k) * oPC(iip1,jj,kk) 
     2                      - uE(ip1,j,k) * uPC(iip1,jj,kk)
     2                      - uSW(ip2,jp1,km1) * dPN(iip1,jj,kk)
     2                      - oNE(ip1,j,k) * oPN(iip1,jj,kk) 
     2                      - uNE(ip1,j,k) * uPN(iip1,jj,kk)) 

      TMP9_XoE =
     2   - uPE(ii,jj,kk) * (- uN(ip1,jm1,k) * oPSW(iip1,jj,kk)
     2                      - oN(ip1,jm1,kp1) * uPSW(iip1,jj,kk)
     2                      - uC(ip1,j,k) * oPW(iip1,jj,kk) 
     2                      + oC(ip1,j,kp1) * uPW(iip1,jj,kk)
     2                      - uS(ip1,jp1,k) * oPNW(iip1,jj,kk)
     2                      - oN(ip1,j,kp1) * uPNW(iip1,jj,kk)
     2                      - uNW(ip2,jm1,k) * oPS(iip1,jj,kk)
     2                      - oNW(ip2,jm1,kp1) * uPS(iip1,jj,kk)
     2                      - uW(ip2,j,k) * oPC(iip1,jj,kk) 
     2                      - oE(ip1,j,kp1) * uPC(iip1,jj,kk)
     2                      - uSW(ip2,jp1,k) * oPN(iip1,jj,kk)
     2                      - oNE(ip1,j,kp1) * uPN(iip1,jj,kk)) 

      TMP10_XoE =
     2   - dPNE(ii,jj,kk) * (- oN(ip1,j,km1) * dPW(iip1,jj,kk)
     2                       - uS(ip1,jp1,km1) * oPW(iip1,jj,kk)
     2                       + oC(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2                       - uC(ip1,jp1,km1) * oPNW(iip1,jj,kk)
     2                       - oNW(ip2,j,km1) * dPC(iip1,jj,kk)
     2                       - uSE(ip1,jp1,km1) * oPC(iip1,jj,kk)
     2                       - oE(ip1,jp1,km1) * dPN(iip1,jj,kk)
     2                       - uE(ip1,jp1,km1) * oPN(iip1,jj,kk))

      TMP11_XoE =
     2   - oPNE(ii,jj,kk) * (- uN(ip1,j,km1) * dPW(iip1,jj,kk) 
     2                       - oN(ip1,j,k) * oPW(iip1,jj,kk)
     2                       - uS(ip1,jp1,k) * uPW(iip1,jj,kk)
     2                       - uC(ip1,jp1,km1) * dPNW(iip1,jj,kk)
     2                       + oC(ip1,jp1,k) * oPNW(iip1,jj,kk)
     2                       - uC(ip1,jp1,k) * uPNW(iip1,jj,kk)
     2                       - uNW(ip2,j,km1) * dPC(iip1,jj,kk)
     2                       - oNW(ip2,j,k) * oPC(iip1,jj,kk)
     2                       - uSE(ip1,jp1,k) * uPC(iip1,jj,kk)
     2                       - uW(ip2,jp1,km1) * dPN(iip1,jj,kk)
     2                       - oE(ip1,jp1,k) * oPN(iip1,jj,kk)
     2                       - uE(ip1,jp1,k) * uPN(iip1,jj,kk)) 

      TMP12_XoE =
     2   - uPNE(ii,jj,kk) * (- uN(ip1,j,k) * oPW(iip1,jj,kk) 
     2                       - oN(ip1,j,kp1) * uPW(iip1,jj,kk)
     2                       - uC(ip1,jp1,k) * oPNW(iip1,jj,kk)
     2                       + oC(ip1,jp1,kp1) * uPNW(iip1,jj,kk)
     2                       - uNW(ip2,j,k) * oPC(iip1,jj,kk)
     2                       - oNW(ip2,j,kp1) * uPC(iip1,jj,kk)
     2                       - uW(ip2,jp1,k) * oPN(iip1,jj,kk)
     2                       - oE(ip1,jp1,kp1) * uPN(iip1,jj,kk))

      XoE(ii,jj,kk) = TMP1_XoE + TMP2_XoE + TMP3_XoE + TMP4_XoE
     2   + TMP5_XoE + TMP6_XoE + TMP7_XoE + TMP8_XoE + TMP9_XoE
     3   + TMP10_XoE + TMP11_XoE + TMP12_XoE

c* ********************************************************************
c* *** > oN;
c* ********************************************************************

c* ***XoN(ii,jj,kk) =
      TMP1_XoN =
     2   - dPW(ii,jj,kk) * (- oN(im1,j,km1) * dPSW(ii,jjp1,kk)
     2                      - uN(im1,j,km1) * oPSW(ii,jjp1,kk)
     2                      - oNE(im1,j,km1) * dPS(ii,jjp1,kk)
     2                      - uNE(im1,j,km1) * oPS(ii,jjp1,kk)) 

     2   - oPW(ii,jj,kk) * (- uS(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                      - oN(im1,j,k) * oPSW(ii,jjp1,kk) 
     2                      - uN(im1,j,k) * uPSW(ii,jjp1,kk)
     2                      - uSW(i,jp1,km1) * dPS(ii,jjp1,kk)
     2                      - oNE(im1,j,k) * oPS(ii,jjp1,kk) 
     2                      - uNE(im1,j,k) * uPS(ii,jjp1,kk)) 

     2   - uPW(ii,jj,kk) * (- uS(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                      - oN(im1,j,kp1) * uPSW(ii,jjp1,kk)
     2                      - uSW(i,jp1,k) * oPS(ii,jjp1,kk)
     2                      - oNE(im1,j,kp1) * uPS(ii,jjp1,kk)) 

      TMP2_XoN =
     2   - dPNW(ii,jj,kk) * (  oC(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                       - uC(im1,jp1,km1) * oPSW(ii,jjp1,kk)
     2                       - oN(im1,jp1,km1) * dPW(ii,jjp1,kk)
     2                       - uN(im1,jp1,km1) * oPW(ii,jjp1,kk)
     2                       - oE(im1,jp1,km1) * dPS(ii,jjp1,kk)
     2                       - uE(im1,jp1,km1) * oPS(ii,jjp1,kk)
     2                       - oNE(im1,jp1,km1) * dPC(ii,jjp1,kk)
     2                       - uNE(im1,jp1,km1) * oPC(ii,jjp1,kk)) 

      TMP3_XoN =
     2   - oPNW(ii,jj,kk) * (- uC(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                       + oC(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                       - uC(im1,jp1,k) * uPSW(ii,jjp1,kk)
     2                       - uS(im1,jp2,km1) * dPW(ii,jjp1,kk)
     2                       - oN(im1,jp1,k) * oPW(ii,jjp1,kk)
     2                       - uN(im1,jp1,k) * uPW(ii,jjp1,kk)
     2                       - uW(i,jp1,km1) * dPS(ii,jjp1,kk)
     2                       - oE(im1,jp1,k) * oPS(ii,jjp1,kk)
     2                       - uE(im1,jp1,k) * uPS(ii,jjp1,kk)
     2                       - uSW(i,jp2,km1) * dPC(ii,jjp1,kk)
     2                       - oNE(im1,jp1,k) * oPC(ii,jjp1,kk)
     2                       - uNE(im1,jp1,k) * uPC(ii,jjp1,kk)) 

      TMP4_XoN =
     2   - uPNW(ii,jj,kk) * (- uC(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                       + oC(im1,jp1,kp1) * uPSW(ii,jjp1,kk)
     2                       - uS(im1,jp2,k) * oPW(ii,jjp1,kk)
     2                       - oN(im1,jp1,kp1) * uPW(ii,jjp1,kk)
     2                       - uW(i,jp1,k) * oPS(ii,jjp1,kk)
     2                       - oE(im1,jp1,kp1) * uPS(ii,jjp1,kk)
     2                       - uSW(i,jp2,k) * oPC(ii,jjp1,kk)
     2                       - oNE(im1,jp1,kp1) * uPC(ii,jjp1,kk)) 

     2   - dPC(ii,jj,kk) * (- oNW(i,j,km1) * dPSW(ii,jjp1,kk) 
     2                      - uNW(i,j,km1) * oPSW(ii,jjp1,kk)
     2                      - oN(i,j,km1) * dPS(ii,jjp1,kk) 
     2                      - uN(i,j,km1) * oPS(ii,jjp1,kk)
     2                      - oNE(i,j,km1) * dPSE(ii,jjp1,kk) 
     2                      - uNE(i,j,km1) * oPSE(ii,jjp1,kk))

      TMP5_XoN =
     2   - oPC(ii,jj,kk) * (- uSE(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                      - oNW(i,j,k) * oPSW(ii,jjp1,kk) 
     2                      - uNW(i,j,k) * uPSW(ii,jjp1,kk)
     2                      - uS(i,jp1,km1) * dPS(ii,jjp1,kk) 
     2                      - oN(i,j,k) * oPS(ii,jjp1,kk)
     2                      - uN(i,j,k) * uPS(ii,jjp1,kk)
     2                      - uSW(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2                      - oNE(i,j,k) * oPSE(ii,jjp1,kk) 
     2                      - uNE(i,j,k) * uPSE(ii,jjp1,kk)) 

     2   - uPC(ii,jj,kk) * (- uSE(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                      - oNW(i,j,kp1) * uPSW(ii,jjp1,kk) 
     2                      - uS(i,jp1,k) * oPS(ii,jjp1,kk)
     2                      - oN(i,j,kp1) * uPS(ii,jjp1,kk)
     2                      - uSW(ip1,jp1,k) * oPSE(ii,jjp1,kk)
     2                      - oNE(i,j,kp1) * uPSE(ii,jjp1,kk)) 

      TMP6_XoN =
     2   - dPN(ii,jj,kk) * (- oE(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                      - uW(i,jp1,km1) * oPSW(ii,jjp1,kk)
     2                      - oNW(i,jp1,km1) * dPW(ii,jjp1,kk)
     2                      - uNW(i,jp1,km1) * oPW(ii,jjp1,kk)
     2                      + oC(i,jp1,km1) * dPS(ii,jjp1,kk)
     2                      - uC(i,jp1,km1) * oPS(ii,jjp1,kk)
     2                      - oN(i,jp1,km1) * dPC(ii,jjp1,kk)
     2                      - uN(i,jp1,km1) * oPC(ii,jjp1,kk)
     2                      - oE(i,jp1,km1) * dPSE(ii,jjp1,kk)
     2                      - uE(i,jp1,km1) * oPSE(ii,jjp1,kk)
     2                      - oNE(i,jp1,km1) * dPE(ii,jjp1,kk)
     2                      - uNE(i,jp1,km1) * oPE(ii,jjp1,kk)) 

      TMP7_XoN =
     2   - oPN(ii,jj,kk) * (- uE(im1,jp1,km1) * dPSW(ii,jjp1,kk)
     2                      - oE(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                      - uW(i,jp1,k) * uPSW(ii,jjp1,kk)
     2                      - uSE(im1,jp2,km1) * dPW(ii,jjp1,kk)
     2                      - oNW(i,jp1,k) * oPW(ii,jjp1,kk) 
     2                      - uNW(i,jp1,k) * uPW(ii,jjp1,kk)
     2                      - uC(i,jp1,km1) * dPS(ii,jjp1,kk) 
     2                      + oC(i,jp1,k) * oPS(ii,jjp1,kk)
     2                      - uC(i,jp1,k) * uPS(ii,jjp1,kk) 
     2                      - uS(i,jp2,km1) * dPC(ii,jjp1,kk)
     2                      - oN(i,jp1,k) * oPC(ii,jjp1,kk) 
     2                      - uN(i,jp1,k) * uPC(ii,jjp1,kk)
     2                      - uW(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2                      - oE(i,jp1,k) * oPSE(ii,jjp1,kk) 
     2                      - uE(i,jp1,k) * uPSE(ii,jjp1,kk)
     2                      - uSW(ip1,jp2,km1) * dPE(ii,jjp1,kk)
     2                      - oNE(i,jp1,k) * oPE(ii,jjp1,kk) 
     2                      - uNE(i,jp1,k) * uPE(ii,jjp1,kk)) 

      TMP8_XoN =
     2   - uPN(ii,jj,kk) * (- uE(im1,jp1,k) * oPSW(ii,jjp1,kk)
     2                      - oE(im1,jp1,kp1) * uPSW(ii,jjp1,kk)
     2                      - uSE(im1,jp2,k) * oPW(ii,jjp1,kk)
     2                      - oNW(i,jp1,kp1) * uPW(ii,jjp1,kk) 
     2                      - uC(i,jp1,k) * oPS(ii,jjp1,kk)
     2                      + oC(i,jp1,kp1) * uPS(ii,jjp1,kk) 
     2                      - uS(i,jp2,k) * oPC(ii,jjp1,kk)
     2                      - oN(i,jp1,kp1) * uPC(ii,jjp1,kk)
     2                      - uW(ip1,jp1,k) * oPSE(ii,jjp1,kk)
     2                      - oE(i,jp1,kp1) * uPSE(ii,jjp1,kk)
     2                      - uSW(ip1,jp2,k) * oPE(ii,jjp1,kk)
     2                      - oNE(i,jp1,kp1) * uPE(ii,jjp1,kk)) 

     2   - dPE(ii,jj,kk) * (- oNW(ip1,j,km1) * dPS(ii,jjp1,kk)
     2                      - uNW(ip1,j,km1) * oPS(ii,jjp1,kk)
     2                      - oN(ip1,j,km1) * dPSE(ii,jjp1,kk)
     2                      - uN(ip1,j,km1) * oPSE(ii,jjp1,kk)) 

      TMP9_XoN =
     2   - oPE(ii,jj,kk) * (- uSE(i,jp1,km1) * dPS(ii,jjp1,kk) 
     2                      - oNW(ip1,j,k) * oPS(ii,jjp1,kk)
     2                      - uNW(ip1,j,k) * uPS(ii,jjp1,kk)
     2                      - uS(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2                      - oN(ip1,j,k) * oPSE(ii,jjp1,kk) 
     2                      - uN(ip1,j,k) * uPSE(ii,jjp1,kk)) 

     2   - uPE(ii,jj,kk) * (- uSE(i,jp1,k) * oPS(ii,jjp1,kk)
     2                      - oNW(ip1,j,kp1) * uPS(ii,jjp1,kk)
     2                      - uS(ip1,jp1,k) * oPSE(ii,jjp1,kk)
     2                      - oN(ip1,j,kp1) * uPSE(ii,jjp1,kk)) 

      TMP10_XoN =
     2   - dPNE(ii,jj,kk) * (- oE(i,jp1,km1) * dPS(ii,jjp1,kk)
     2                       - uW(ip1,jp1,km1) * oPS(ii,jjp1,kk)
     2                       - oNW(ip1,jp1,km1) * dPC(ii,jjp1,kk)
     2                       - uNW(ip1,jp1,km1) * oPC(ii,jjp1,kk)
     2                       + oC(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2                       - uC(ip1,jp1,km1) * oPSE(ii,jjp1,kk)
     2                       - oN(ip1,jp1,km1) * dPE(ii,jjp1,kk)
     2                       - uN(ip1,jp1,km1) * oPE(ii,jjp1,kk))

      TMP11_XoN =
     2    - oPNE(ii,jj,kk) * (- uE(i,jp1,km1) * dPS(ii,jjp1,kk) 
     2                        - oE(i,jp1,k) * oPS(ii,jjp1,kk)
     2                        - uW(ip1,jp1,k) * uPS(ii,jjp1,kk)
     2                        - uSE(i,jp2,km1) * dPC(ii,jjp1,kk)
     2                        - oNW(ip1,jp1,k) * oPC(ii,jjp1,kk)
     2                        - uNW(ip1,jp1,k) * uPC(ii,jjp1,kk)
     2                        - uC(ip1,jp1,km1) * dPSE(ii,jjp1,kk)
     2                        + oC(ip1,jp1,k) * oPSE(ii,jjp1,kk)
     2                        - uC(ip1,jp1,k) * uPSE(ii,jjp1,kk)
     2                        - uS(ip1,jp2,km1) * dPE(ii,jjp1,kk)
     2                        - oN(ip1,jp1,k) * oPE(ii,jjp1,kk)
     2                        - uN(ip1,jp1,k) * uPE(ii,jjp1,kk)) 

      TMP12_XoN =
     2   - uPNE(ii,jj,kk) * (- uE(i,jp1,k) * oPS(ii,jjp1,kk) 
     2                       - oE(i,jp1,kp1) * uPS(ii,jjp1,kk)
     2                       - uSE(i,jp2,k) * oPC(ii,jjp1,kk)
     2                       - oNW(ip1,jp1,kp1) * uPC(ii,jjp1,kk)
     2                       - uC(ip1,jp1,k) * oPSE(ii,jjp1,kk)
     2                       + oC(ip1,jp1,kp1) * uPSE(ii,jjp1,kk)
     2                       - uS(ip1,jp2,k) * oPE(ii,jjp1,kk)
     2                       - oN(ip1,jp1,kp1) * uPE(ii,jjp1,kk))

      XoN(ii,jj,kk) = TMP1_XoN + TMP2_XoN + TMP3_XoN + TMP4_XoN
     2   + TMP5_XoN + TMP6_XoN + TMP7_XoN + TMP8_XoN + TMP9_XoN
     3   + TMP10_XoN + TMP11_XoN + TMP12_XoN

c* ********************************************************************
c* *** > uC;
c* ********************************************************************

c* ***XuC(ii,jj,kk) =
      TMP1_XuC =
     2   - oPSW(ii,jj,kk) * (- uC(im1,jm1,k) * dPSW(ii,jj,kkp1)
     2                       - uN(im1,jm1,k) * dPW(ii,jj,kkp1)
     2                       - uE(im1,jm1,k) * dPS(ii,jj,kkp1)
     2                       - uNE(im1,jm1,k) * dPC(ii,jj,kkp1)) 

     2   - uPSW(ii,jj,kk) * (  oC(im1,jm1,kp1) * dPSW(ii,jj,kkp1)
     2                       - uC(im1,jm1,kp1) * oPSW(ii,jj,kkp1)
     2                       - oN(im1,jm1,kp1) * dPW(ii,jj,kkp1)
     2                       - uN(im1,jm1,kp1) * oPW(ii,jj,kkp1)
     2                       - oE(im1,jm1,kp1) * dPS(ii,jj,kkp1)
     2                       - uE(im1,jm1,kp1) * oPS(ii,jj,kkp1)
     2                       - oNE(im1,jm1,kp1) * dPC(ii,jj,kkp1)
     2                       - uNE(im1,jm1,kp1) * oPC(ii,jj,kkp1)) 

      TMP2_XuC =
     2   - oPW(ii,jj,kk) * (- uS(im1,j,k) * dPSW(ii,jj,kkp1) 
     2                      - uC(im1,j,k) * dPW(ii,jj,kkp1)
     2                      - uN(im1,j,k) * dPNW(ii,jj,kkp1) 
     2                      - uSE(im1,j,k) * dPS(ii,jj,kkp1)
     2                      - uE(im1,j,k) * dPC(ii,jj,kkp1) 
     2                      - uNE(im1,j,k) * dPN(ii,jj,kkp1)) 

      TMP3_XuC =
     2   - uPW(ii,jj,kk) * (- oN(im1,jm1,kp1) * dPSW(ii,jj,kkp1)
     2                      - uS(im1,j,kp1) * oPSW(ii,jj,kkp1)
     2                      + oC(im1,j,kp1) * dPW(ii,jj,kkp1)
     2                      - uC(im1,j,kp1) * oPW(ii,jj,kkp1)
     2                      - oN(im1,j,kp1) * dPNW(ii,jj,kkp1)
     2                      - uN(im1,j,kp1) * oPNW(ii,jj,kkp1)
     2                      - oNW(i,jm1,kp1) * dPS(ii,jj,kkp1)
     2                      - uSE(im1,j,kp1) * oPS(ii,jj,kkp1)
     2                      - oE(im1,j,kp1) * dPC(ii,jj,kkp1)
     2                      - uE(im1,j,kp1) * oPC(ii,jj,kkp1)
     2                      - oNE(im1,j,kp1) * dPN(ii,jj,kkp1)
     2                      - uNE(im1,j,kp1) * oPN(ii,jj,kkp1)) 

     2   - oPNW(ii,jj,kk) * (- uS(im1,jp1,k) * dPW(ii,jj,kkp1)
     2                       - uC(im1,jp1,k) * dPNW(ii,jj,kkp1)
     2                       - uSE(im1,jp1,k) * dPC(ii,jj,kkp1)
     2                       - uE(im1,jp1,k) * dPN(ii,jj,kkp1)) 

      TMP4_XuC =
     2   - uPNW(ii,jj,kk) * (- oN(im1,j,kp1) * dPW(ii,jj,kkp1)
     2                       - uS(im1,jp1,kp1) * oPW(ii,jj,kkp1)
     2                       + oC(im1,jp1,kp1) * dPNW(ii,jj,kkp1)
     2                       - uC(im1,jp1,kp1) * oPNW(ii,jj,kkp1)
     2                       - oNW(i,j,kp1) * dPC(ii,jj,kkp1)
     2                       - uSE(im1,jp1,kp1) * oPC(ii,jj,kkp1)
     2                       - oE(im1,jp1,kp1) * dPN(ii,jj,kkp1)
     2                       - uE(im1,jp1,kp1) * oPN(ii,jj,kkp1))

     2   - oPS(ii,jj,kk) * (- uW(i,jm1,k) * dPSW(ii,jj,kkp1) 
     2                      - uNW(i,jm1,k) * dPW(ii,jj,kkp1)
     2                      - uC(i,jm1,k) * dPS(ii,jj,kkp1) 
     2                      - uN(i,jm1,k) * dPC(ii,jj,kkp1)
     2                      - uE(i,jm1,k) * dPSE(ii,jj,kkp1) 
     2                      - uNE(i,jm1,k) * dPE(ii,jj,kkp1)) 

      TMP5_XuC =
     2   - uPS(ii,jj,kk) * (- oE(im1,jm1,kp1) * dPSW(ii,jj,kkp1)
     2                      - uW(i,jm1,kp1) * oPSW(ii,jj,kkp1)
     2                      - oNW(i,jm1,kp1) * dPW(ii,jj,kkp1)
     2                      - uNW(i,jm1,kp1) * oPW(ii,jj,kkp1)
     2                      + oC(i,jm1,kp1) * dPS(ii,jj,kkp1)
     2                      - uC(i,jm1,kp1) * oPS(ii,jj,kkp1)
     2                      - oN(i,jm1,kp1) * dPC(ii,jj,kkp1)
     2                      - uN(i,jm1,kp1) * oPC(ii,jj,kkp1)
     2                      - oE(i,jm1,kp1) * dPSE(ii,jj,kkp1)
     2                      - uE(i,jm1,kp1) * oPSE(ii,jj,kkp1)
     2                      - oNE(i,jm1,kp1) * dPE(ii,jj,kkp1)
     2                      - uNE(i,jm1,kp1) * oPE(ii,jj,kkp1)) 

      TMP6_XuC =
     2   - oPC(ii,jj,kk) * (- uSW(i,j,k) * dPSW(ii,jj,kkp1) 
     2                      - uW(i,j,k) * dPW(ii,jj,kkp1)
     2                      - uNW(i,j,k) * dPNW(ii,jj,kkp1) 
     2                      - uS(i,j,k) * dPS(ii,jj,kkp1)
     2                      - uC(i,j,k) * dPC(ii,jj,kkp1) 
     2                      - uN(i,j,k) * dPN(ii,jj,kkp1)
     2                      - uSE(i,j,k) * dPSE(ii,jj,kkp1) 
     2                      - uE(i,j,k) * dPE(ii,jj,kkp1)
     2                      - uNE(i,j,k) * dPNE(ii,jj,kkp1)) 

      TMP7_XuC =
     2   - uPC(ii,jj,kk) * (- oNE(im1,jm1,kp1) * dPSW(ii,jj,kkp1)
     2                      - uSW(i,j,kp1) * oPSW(ii,jj,kkp1)
     2                      - oE(im1,j,kp1) * dPW(ii,jj,kkp1) 
     2                      - uW(i,j,kp1) * oPW(ii,jj,kkp1)
     2                      - oNW(i,j,kp1) * dPNW(ii,jj,kkp1) 
     2                      - uNW(i,j,kp1) * oPNW(ii,jj,kkp1)
     2                      - oN(i,jm1,kp1) * dPS(ii,jj,kkp1) 
     2                      - uS(i,j,kp1) * oPS(ii,jj,kkp1)
     2                      + oC(i,j,kp1) * dPC(ii,jj,kkp1) 
     2                      - uC(i,j,kp1) * oPC(ii,jj,kkp1)
     2                      - oN(i,j,kp1) * dPN(ii,jj,kkp1) 
     2                      - uN(i,j,kp1) * oPN(ii,jj,kkp1)
     2                      - oNW(ip1,jm1,kp1) * dPSE(ii,jj,kkp1)
     2                      - uSE(i,j,kp1) * oPSE(ii,jj,kkp1) 
     2                      - oE(i,j,kp1) * dPE(ii,jj,kkp1)
     2                      - uE(i,j,kp1) * oPE(ii,jj,kkp1) 
     2                      - oNE(i,j,kp1) * dPNE(ii,jj,kkp1)
     2                      - uNE(i,j,kp1) * oPNE(ii,jj,kkp1)) 

      TMP8_XuC =
     2   - oPN(ii,jj,kk) * (- uSW(i,jp1,k) * dPW(ii,jj,kkp1) 
     2                      - uW(i,jp1,k) * dPNW(ii,jj,kkp1)
     2                      - uS(i,jp1,k) * dPC(ii,jj,kkp1) 
     2                      - uC(i,jp1,k) * dPN(ii,jj,kkp1)
     2                      - uSE(i,jp1,k) * dPE(ii,jj,kkp1) 
     2                      - uE(i,jp1,k) * dPNE(ii,jj,kkp1)) 

      TMP9_XuC =
     2   - uPN(ii,jj,kk) * (- oNE(im1,j,kp1) * dPW(ii,jj,kkp1)
     2                      - uSW(i,jp1,kp1) * oPW(ii,jj,kkp1)
     2                      - oE(im1,jp1,kp1) * dPNW(ii,jj,kkp1)
     2                      - uW(i,jp1,kp1) * oPNW(ii,jj,kkp1) 
     2                      - oN(i,j,kp1) * dPC(ii,jj,kkp1)
     2                      - uS(i,jp1,kp1) * oPC(ii,jj,kkp1)
     2                      + oC(i,jp1,kp1) * dPN(ii,jj,kkp1)
     2                      - uC(i,jp1,kp1) * oPN(ii,jj,kkp1)
     2                      - oNW(ip1,j,kp1) * dPE(ii,jj,kkp1)
     2                      - uSE(i,jp1,kp1) * oPE(ii,jj,kkp1)
     2                      - oE(i,jp1,kp1) * dPNE(ii,jj,kkp1)
     2                      - uE(i,jp1,kp1) * oPNE(ii,jj,kkp1)) 

     2   - oPSE(ii,jj,kk) * (- uW(ip1,jm1,k) * dPS(ii,jj,kkp1)
     2                       - uNW(ip1,jm1,k) * dPC(ii,jj,kkp1)
     2                       - uC(ip1,jm1,k) * dPSE(ii,jj,kkp1)
     2                       - uN(ip1,jm1,k) * dPE(ii,jj,kkp1)) 

      TMP10_XuC =
     2   - uPSE(ii,jj,kk) * (- oE(i,jm1,kp1) * dPS(ii,jj,kkp1)
     2                       - uW(ip1,jm1,kp1) * oPS(ii,jj,kkp1)
     2                       - oNW(ip1,jm1,kp1) * dPC(ii,jj,kkp1)
     2                       - uNW(ip1,jm1,kp1) * oPC(ii,jj,kkp1)
     2                       + oC(ip1,jm1,kp1) * dPSE(ii,jj,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPSE(ii,jj,kkp1)
     2                       - oN(ip1,jm1,kp1) * dPE(ii,jj,kkp1)
     2                       - uN(ip1,jm1,kp1) * oPE(ii,jj,kkp1))

     2   - oPE(ii,jj,kk) * (- uSW(ip1,j,k) * dPS(ii,jj,kkp1) 
     2                      - uW(ip1,j,k) * dPC(ii,jj,kkp1)
     2                      - uNW(ip1,j,k) * dPN(ii,jj,kkp1) 
     2                      - uS(ip1,j,k) * dPSE(ii,jj,kkp1)
     2                      - uC(ip1,j,k) * dPE(ii,jj,kkp1) 
     2                      - uN(ip1,j,k) * dPNE(ii,jj,kkp1)) 

      TMP11_XuC =
     2   - uPE(ii,jj,kk) * (- oNE(i,jm1,kp1) * dPS(ii,jj,kkp1)
     2                      - uSW(ip1,j,kp1) * oPS(ii,jj,kkp1) 
     2                      - oE(i,j,kp1) * dPC(ii,jj,kkp1)
     2                      - uW(ip1,j,kp1) * oPC(ii,jj,kkp1)
     2                      - oNW(ip1,j,kp1) * dPN(ii,jj,kkp1)
     2                      - uNW(ip1,j,kp1) * oPN(ii,jj,kkp1)
     2                      - oN(ip1,jm1,kp1) * dPSE(ii,jj,kkp1)
     2                      - uS(ip1,j,kp1) * oPSE(ii,jj,kkp1)
     2                      + oC(ip1,j,kp1) * dPE(ii,jj,kkp1)
     2                      - uC(ip1,j,kp1) * oPE(ii,jj,kkp1)
     2                      - oN(ip1,j,kp1) * dPNE(ii,jj,kkp1)
     2                      - uN(ip1,j,kp1) * oPNE(ii,jj,kkp1)) 

      TMP12_XuC =
     2   - oPNE(ii,jj,kk) * (- uSW(ip1,jp1,k) * dPC(ii,jj,kkp1)
     2                       - uW(ip1,jp1,k) * dPN(ii,jj,kkp1)
     2                       - uS(ip1,jp1,k) * dPE(ii,jj,kkp1)
     2                       - uC(ip1,jp1,k) * dPNE(ii,jj,kkp1))

     2   - uPNE(ii,jj,kk) * (- oNE(i,j,kp1) * dPC(ii,jj,kkp1)
     2                       - uSW(ip1,jp1,kp1) * oPC(ii,jj,kkp1)
     2                       - oE(i,jp1,kp1) * dPN(ii,jj,kkp1)
     2                       - uW(ip1,jp1,kp1) * oPN(ii,jj,kkp1)
     2                       - oN(ip1,j,kp1) * dPE(ii,jj,kkp1)
     2                       - uS(ip1,jp1,kp1) * oPE(ii,jj,kkp1)
     2                       + oC(ip1,jp1,kp1) * dPNE(ii,jj,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPNE(ii,jj,kkp1))

      XuC(ii,jj,kk) = TMP1_XuC + TMP2_XuC + TMP3_XuC + TMP4_XuC
     2   + TMP5_XuC + TMP6_XuC + TMP7_XuC + TMP8_XuC + TMP9_XuC
     3   + TMP10_XuC + TMP11_XuC + TMP12_XuC

c* ********************************************************************
c* *** > oNE;
c* ********************************************************************

c* ***XoNE(ii,jj,kk) =
      TMP1_XoNE =
     2   - dPC(ii,jj,kk) * (- oNE(i,j,km1) * dPSW(iip1,jjp1,kk)
     2                      - uNE(i,j,km1) * oPSW(iip1,jjp1,kk)) 

     2   - oPC(ii,jj,kk) * (- uSW(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                      - oNE(i,j,k) * oPSW(iip1,jjp1,kk) 
     2                      - uNE(i,j,k) * uPSW(iip1,jjp1,kk))

     2   - uPC(ii,jj,kk) * (- uSW(ip1,jp1,k) * oPSW(iip1,jjp1,kk)
     2                      - oNE(i,j,kp1) * uPSW(iip1,jjp1,kk)) 

     2   - dPN(ii,jj,kk) * (- oE(i,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                      - uE(i,jp1,km1) * oPSW(iip1,jjp1,kk)
     2                      - oNE(i,jp1,km1) * dPW(iip1,jjp1,kk)
     2                      - uNE(i,jp1,km1) * oPW(iip1,jjp1,kk))

      TMP2_XoNE =
     2   - oPN(ii,jj,kk) * (- uW(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                      - oE(i,jp1,k) * oPSW(iip1,jjp1,kk)
     2                      - uE(i,jp1,k) * uPSW(iip1,jjp1,kk)
     2                      - uSW(ip1,jp2,km1) * dPW(iip1,jjp1,kk)
     2                      - oNE(i,jp1,k) * oPW(iip1,jjp1,kk)
     2                      - uNE(i,jp1,k) * uPW(iip1,jjp1,kk)) 

     2   - uPN(ii,jj,kk) * (- uW(ip1,jp1,k) * oPSW(iip1,jjp1,kk)
     2                      - oE(i,jp1,kp1) * uPSW(iip1,jjp1,kk)
     2                      - uSW(ip1,jp2,k) * oPW(iip1,jjp1,kk)
     2                      - oNE(i,jp1,kp1) * uPW(iip1,jjp1,kk))

     2    - dPE(ii,jj,kk) * (- oN(ip1,j,km1) * dPSW(iip1,jjp1,kk)
     2                       - uN(ip1,j,km1) * oPSW(iip1,jjp1,kk)
     2                       - oNE(ip1,j,km1) * dPS(iip1,jjp1,kk)
     2                       - uNE(ip1,j,km1) * oPS(iip1,jjp1,kk))

      TMP3_XoNE =
     2    - oPE(ii,jj,kk) * (- uS(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                       - oN(ip1,j,k) * oPSW(iip1,jjp1,kk)
     2                       - uN(ip1,j,k) * uPSW(iip1,jjp1,kk)
     2                       - uSW(ip2,jp1,km1) * dPS(iip1,jjp1,kk)
     2                       - oNE(ip1,j,k) * oPS(iip1,jjp1,kk)
     2                       - uNE(ip1,j,k) * uPS(iip1,jjp1,kk)) 

     2   - uPE(ii,jj,kk) * (- uS(ip1,jp1,k) * oPSW(iip1,jjp1,kk)
     2                      - oN(ip1,j,kp1) * uPSW(iip1,jjp1,kk)
     2                      - uSW(ip2,jp1,k) * oPS(iip1,jjp1,kk)
     2                      - oNE(ip1,j,kp1) * uPS(iip1,jjp1,kk))

      TMP4_XoNE =
     2    - dPNE(ii,jj,kk) * (  oC(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                        - uC(ip1,jp1,km1) * oPSW(iip1,jjp1,kk)
     2                        - oN(ip1,jp1,km1) * dPW(iip1,jjp1,kk)
     2                        - uN(ip1,jp1,km1) * oPW(iip1,jjp1,kk)
     2                        - oE(ip1,jp1,km1) * dPS(iip1,jjp1,kk)
     2                        - uE(ip1,jp1,km1) * oPS(iip1,jjp1,kk)
     2                        - oNE(ip1,jp1,km1) * dPC(iip1,jjp1,kk)
     2                        - uNE(ip1,jp1,km1) * oPC(iip1,jjp1,kk)) 

      TMP5_XoNE =
     2   - oPNE(ii,jj,kk) * (- uC(ip1,jp1,km1) * dPSW(iip1,jjp1,kk)
     2                       + oC(ip1,jp1,k) * oPSW(iip1,jjp1,kk)
     2                       - uC(ip1,jp1,k) * uPSW(iip1,jjp1,kk)
     2                       - uS(ip1,jp2,km1) * dPW(iip1,jjp1,kk)
     2                       - oN(ip1,jp1,k) * oPW(iip1,jjp1,kk)
     2                       - uN(ip1,jp1,k) * uPW(iip1,jjp1,kk)
     2                       - uW(ip2,jp1,km1) * dPS(iip1,jjp1,kk)
     2                       - oE(ip1,jp1,k) * oPS(iip1,jjp1,kk)
     2                       - uE(ip1,jp1,k) * uPS(iip1,jjp1,kk)
     2                       - uSW(ip2,jp2,km1) * dPC(iip1,jjp1,kk)
     2                       - oNE(ip1,jp1,k) * oPC(iip1,jjp1,kk)
     2                       - uNE(ip1,jp1,k) * uPC(iip1,jjp1,kk)) 

      TMP6_XoNE =
     2   - uPNE(ii,jj,kk) * (- uC(ip1,jp1,k) * oPSW(iip1,jjp1,kk)
     2                       + oC(ip1,jp1,kp1) * uPSW(iip1,jjp1,kk)
     2                       - uS(ip1,jp2,k) * oPW(iip1,jjp1,kk)
     2                       - oN(ip1,jp1,kp1) * uPW(iip1,jjp1,kk)
     2                       - uW(ip2,jp1,k) * oPS(iip1,jjp1,kk)
     2                       - oE(ip1,jp1,kp1) * uPS(iip1,jjp1,kk)
     2                       - uSW(ip2,jp2,k) * oPC(iip1,jjp1,kk)
     2                       - oNE(ip1,jp1,kp1) * uPC(iip1,jjp1,kk))

      XoNE(ii,jj,kk) = TMP1_XoNE + TMP2_XoNE + TMP3_XoNE + TMP4_XoNE
     2               + TMP5_XoNE + TMP6_XoNE

c* ********************************************************************
c* *** > oNW;
c* ********************************************************************

c* ***XoNW(ii,jj,kk) =
      TMP1_XoNW =
     2   - dPW(ii,jj,kk) * (- oNW(im1,j,km1) * dPS(iim1,jjp1,kk)
     2                      - uNW(im1,j,km1) * oPS(iim1,jjp1,kk)
     2                      - oN(im1,j,km1) * dPSE(iim1,jjp1,kk)
     2                      - uN(im1,j,km1) * oPSE(iim1,jjp1,kk)) 

     2   - oPW(ii,jj,kk) * (- uSE(im2,jp1,km1) * dPS(iim1,jjp1,kk)
     2                      - oNW(im1,j,k) * oPS(iim1,jjp1,kk)
     2                      - uNW(im1,j,k) * uPS(iim1,jjp1,kk)
     2                      - uS(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                      - oN(im1,j,k) * oPSE(iim1,jjp1,kk)
     2                      - uN(im1,j,k) * uPSE(iim1,jjp1,kk)) 

      TMP2_XoNW =
     2   - uPW(ii,jj,kk) * (- uSE(im2,jp1,k) * oPS(iim1,jjp1,kk)
     2                      - oNW(im1,j,kp1) * uPS(iim1,jjp1,kk)
     2                      - uS(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                      - oN(im1,j,kp1) * uPSE(iim1,jjp1,kk)) 

     2   - dPNW(ii,jj,kk) * (- oE(im2,jp1,km1) * dPS(iim1,jjp1,kk)
     2                       - uW(im1,jp1,km1) * oPS(iim1,jjp1,kk)
     2                       - oNW(im1,jp1,km1) * dPC(iim1,jjp1,kk)
     2                       - uNW(im1,jp1,km1) * oPC(iim1,jjp1,kk)
     2                       + oC(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                       - uC(im1,jp1,km1) * oPSE(iim1,jjp1,kk)
     2                       - oN(im1,jp1,km1) * dPE(iim1,jjp1,kk)
     2                       - uN(im1,jp1,km1) * oPE(iim1,jjp1,kk)) 

      TMP3_XoNW =
     2   - oPNW(ii,jj,kk) * (- uE(im2,jp1,km1) * dPS(iim1,jjp1,kk)
     2                       - oE(im2,jp1,k) * oPS(iim1,jjp1,kk)
     2                       - uW(im1,jp1,k) * uPS(iim1,jjp1,kk)
     2                       - uSE(im2,jp2,km1) * dPC(iim1,jjp1,kk)
     2                       - oNW(im1,jp1,k) * oPC(iim1,jjp1,kk)
     2                       - uNW(im1,jp1,k) * uPC(iim1,jjp1,kk)
     2                       - uC(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                       + oC(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                       - uC(im1,jp1,k) * uPSE(iim1,jjp1,kk)
     2                       - uS(im1,jp2,km1) * dPE(iim1,jjp1,kk)
     2                       - oN(im1,jp1,k) * oPE(iim1,jjp1,kk)
     2                       - uN(im1,jp1,k) * uPE(iim1,jjp1,kk)) 

      TMP4_XoNW =
     2   - uPNW(ii,jj,kk) * (- uE(im2,jp1,k) * oPS(iim1,jjp1,kk)
     2                       - oE(im2,jp1,kp1) * uPS(iim1,jjp1,kk)
     2                       - uSE(im2,jp2,k) * oPC(iim1,jjp1,kk)
     2                       - oNW(im1,jp1,kp1) * uPC(iim1,jjp1,kk)
     2                       - uC(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                       + oC(im1,jp1,kp1) * uPSE(iim1,jjp1,kk)
     2                       - uS(im1,jp2,k) * oPE(iim1,jjp1,kk)
     2                       - oN(im1,jp1,kp1) * uPE(iim1,jjp1,kk)) 

     2   - dPC(ii,jj,kk) * (- oNW(i,j,km1) * dPSE(iim1,jjp1,kk)
     2                      - uNW(i,j,km1) * oPSE(iim1,jjp1,kk)) 

      TMP5_XoNW =
     2   - oPC(ii,jj,kk) * (- uSE(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                      - oNW(i,j,k) * oPSE(iim1,jjp1,kk) 
     2                      - uNW(i,j,k) * uPSE(iim1,jjp1,kk))

     2   - uPC(ii,jj,kk) * (- uSE(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                      - oNW(i,j,kp1) * uPSE(iim1,jjp1,kk)) 

     2   - dPN(ii,jj,kk) * (- oE(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                      - uW(i,jp1,km1) * oPSE(iim1,jjp1,kk)
     2                      - oNW(i,jp1,km1) * dPE(iim1,jjp1,kk)
     2                      - uNW(i,jp1,km1) * oPE(iim1,jjp1,kk)) 

      TMP6_XoNW =
     2   - oPN(ii,jj,kk) * (- uE(im1,jp1,km1) * dPSE(iim1,jjp1,kk)
     2                      - oE(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                      - uW(i,jp1,k) * uPSE(iim1,jjp1,kk)
     2                      - uSE(im1,jp2,km1) * dPE(iim1,jjp1,kk)
     2                      - oNW(i,jp1,k) * oPE(iim1,jjp1,kk)
     2                      - uNW(i,jp1,k) * uPE(iim1,jjp1,kk)) 

     2   - uPN(ii,jj,kk) * (- uE(im1,jp1,k) * oPSE(iim1,jjp1,kk)
     2                      - oE(im1,jp1,kp1) * uPSE(iim1,jjp1,kk)
     2                      - uSE(im1,jp2,k) * oPE(iim1,jjp1,kk)
     2                      - oNW(i,jp1,kp1) * uPE(iim1,jjp1,kk))

      XoNW(ii,jj,kk) = TMP1_XoNW + TMP2_XoNW + TMP3_XoNW + TMP4_XoNW
     2               + TMP5_XoNW + TMP6_XoNW

c* ********************************************************************
c* *** > uE;
c* ********************************************************************

c* ***XuE(ii,jj,kk) =
      TMP1_XuE =
     2   - oPS(ii,jj,kk) * (- uE(i,jm1,k) * dPSW(iip1,jj,kkp1) 
     2                      - uNE(i,jm1,k) * dPW(iip1,jj,kkp1)) 

     2   - uPS(ii,jj,kk) * (- oE(i,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2                      - uE(i,jm1,kp1) * oPSW(iip1,jj,kkp1)
     2                      - oNE(i,jm1,kp1) * dPW(iip1,jj,kkp1)
     2                      - uNE(i,jm1,kp1) * oPW(iip1,jj,kkp1)) 

     2   - oPC(ii,jj,kk) * (- uSE(i,j,k) * dPSW(iip1,jj,kkp1) 
     2                      - uE(i,j,k) * dPW(iip1,jj,kkp1)
     2                      - uNE(i,j,k) * dPNW(iip1,jj,kkp1)) 

      TMP2_XuE =
     2   - uPC(ii,jj,kk) * (- oNW(ip1,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2                      - uSE(i,j,kp1) * oPSW(iip1,jj,kkp1)
     2                      - oE(i,j,kp1) * dPW(iip1,jj,kkp1) 
     2                      - uE(i,j,kp1) * oPW(iip1,jj,kkp1)
     2                      - oNE(i,j,kp1) * dPNW(iip1,jj,kkp1)
     2                      - uNE(i,j,kp1) * oPNW(iip1,jj,kkp1)) 

     2   - oPN(ii,jj,kk) * (- uSE(i,jp1,k) * dPW(iip1,jj,kkp1) 
     2                      - uE(i,jp1,k) * dPNW(iip1,jj,kkp1)) 

     2   - uPN(ii,jj,kk) * (- oNW(ip1,j,kp1) * dPW(iip1,jj,kkp1)
     2                      - uSE(i,jp1,kp1) * oPW(iip1,jj,kkp1)
     2                      - oE(i,jp1,kp1) * dPNW(iip1,jj,kkp1)
     2                      - uE(i,jp1,kp1) * oPNW(iip1,jj,kkp1)) 

      TMP3_XuE =
     2   - oPSE(ii,jj,kk) * (- uC(ip1,jm1,k) * dPSW(iip1,jj,kkp1)
     2                       - uN(ip1,jm1,k) * dPW(iip1,jj,kkp1)
     2                       - uE(ip1,jm1,k) * dPS(iip1,jj,kkp1)
     2                       - uNE(ip1,jm1,k) * dPC(iip1,jj,kkp1)) 

     2   - uPSE(ii,jj,kk) * (  oC(ip1,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPSW(iip1,jj,kkp1)
     2                       - oN(ip1,jm1,kp1) * dPW(iip1,jj,kkp1)
     2                       - uN(ip1,jm1,kp1) * oPW(iip1,jj,kkp1)
     2                       - oE(ip1,jm1,kp1) * dPS(iip1,jj,kkp1)
     2                       - uE(ip1,jm1,kp1) * oPS(iip1,jj,kkp1)
     2                       - oNE(ip1,jm1,kp1) * dPC(iip1,jj,kkp1)
     2                       - uNE(ip1,jm1,kp1) * oPC(iip1,jj,kkp1)) 

      TMP4_XuE =
     2   - oPE(ii,jj,kk) * (- uS(ip1,j,k) * dPSW(iip1,jj,kkp1) 
     2                      - uC(ip1,j,k) * dPW(iip1,jj,kkp1)
     2                      - uN(ip1,j,k) * dPNW(iip1,jj,kkp1)
     2                      - uSE(ip1,j,k) * dPS(iip1,jj,kkp1) 
     2                      - uE(ip1,j,k) * dPC(iip1,jj,kkp1)
     2                      - uNE(ip1,j,k) * dPN(iip1,jj,kkp1)) 

      TMP5_XuE =
     2   - uPE(ii,jj,kk) * (- oN(ip1,jm1,kp1) * dPSW(iip1,jj,kkp1)
     2                      - uS(ip1,j,kp1) * oPSW(iip1,jj,kkp1)
     2                      + oC(ip1,j,kp1) * dPW(iip1,jj,kkp1)
     2                      - uC(ip1,j,kp1) * oPW(iip1,jj,kkp1)
     2                      - oN(ip1,j,kp1) * dPNW(iip1,jj,kkp1)
     2                      - uN(ip1,j,kp1) * oPNW(iip1,jj,kkp1)
     2                      - oNW(ip2,jm1,kp1) * dPS(iip1,jj,kkp1)
     2                      - uSE(ip1,j,kp1) * oPS(iip1,jj,kkp1)
     2                      - oE(ip1,j,kp1) * dPC(iip1,jj,kkp1)
     2                      - uE(ip1,j,kp1) * oPC(iip1,jj,kkp1)
     2                      - oNE(ip1,j,kp1) * dPN(iip1,jj,kkp1)
     2                      - uNE(ip1,j,kp1) * oPN(iip1,jj,kkp1)) 

      TMP6_XuE =
     2   - oPNE(ii,jj,kk) * (- uS(ip1,jp1,k) * dPW(iip1,jj,kkp1)
     2                       - uC(ip1,jp1,k) * dPNW(iip1,jj,kkp1)
     2                       - uSE(ip1,jp1,k) * dPC(iip1,jj,kkp1)
     2                       - uE(ip1,jp1,k) * dPN(iip1,jj,kkp1)) 

     2   - uPNE(ii,jj,kk) * (- oN(ip1,j,kp1) * dPW(iip1,jj,kkp1)
     2                       - uS(ip1,jp1,kp1) * oPW(iip1,jj,kkp1)
     2                       + oC(ip1,jp1,kp1) * dPNW(iip1,jj,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPNW(iip1,jj,kkp1)
     2                       - oNW(ip2,j,kp1) * dPC(iip1,jj,kkp1)
     2                       - uSE(ip1,jp1,kp1) * oPC(iip1,jj,kkp1)
     2                       - oE(ip1,jp1,kp1) * dPN(iip1,jj,kkp1)
     2                       - uE(ip1,jp1,kp1) * oPN(iip1,jj,kkp1))

      XuE(ii,jj,kk) = TMP1_XuE + TMP2_XuE + TMP3_XuE + TMP4_XuE
     2              + TMP5_XuE + TMP6_XuE

c* ********************************************************************
c* *** > uW;
c* ********************************************************************

c* ***XuW(ii,jj,kk) =
      TMP1_XuW =
     2   - oPSW(ii,jj,kk) * (- uW(im1,jm1,k) * dPS(iim1,jj,kkp1)
     2                       - uNW(im1,jm1,k) * dPC(iim1,jj,kkp1)
     2                       - uC(im1,jm1,k) * dPSE(iim1,jj,kkp1)
     2                       - uN(im1,jm1,k) * dPE(iim1,jj,kkp1)) 

     2   - uPSW(ii,jj,kk) * (- oE(im2,jm1,kp1) * dPS(iim1,jj,kkp1)
     2                       - uW(im1,jm1,kp1) * oPS(iim1,jj,kkp1)
     2                       - oNW(im1,jm1,kp1) * dPC(iim1,jj,kkp1)
     2                       - uNW(im1,jm1,kp1) * oPC(iim1,jj,kkp1)
     2                       + oC(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2                       - uC(im1,jm1,kp1) * oPSE(iim1,jj,kkp1)
     2                       - oN(im1,jm1,kp1) * dPE(iim1,jj,kkp1)
     2                       - uN(im1,jm1,kp1) * oPE(iim1,jj,kkp1)) 

      TMP2_XuW =
     2   - oPW(ii,jj,kk) * (- uSW(im1,j,k) * dPS(iim1,jj,kkp1) 
     2                      - uW(im1,j,k) * dPC(iim1,jj,kkp1)
     2                      - uNW(im1,j,k) * dPN(iim1,jj,kkp1)
     2                      - uS(im1,j,k) * dPSE(iim1,jj,kkp1) 
     2                      - uC(im1,j,k) * dPE(iim1,jj,kkp1)
     2                      - uN(im1,j,k) * dPNE(iim1,jj,kkp1)) 

      TMP3_XuW =
     2   - uPW(ii,jj,kk) * (- oNE(im2,jm1,kp1) * dPS(iim1,jj,kkp1)
     2                      - uSW(im1,j,kp1) * oPS(iim1,jj,kkp1)
     2                      - oE(im2,j,kp1) * dPC(iim1,jj,kkp1)
     2                      - uW(im1,j,kp1) * oPC(iim1,jj,kkp1)
     2                      - oNW(im1,j,kp1) * dPN(iim1,jj,kkp1)
     2                      - uNW(im1,j,kp1) * oPN(iim1,jj,kkp1)
     2                      - oN(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2                      - uS(im1,j,kp1) * oPSE(iim1,jj,kkp1)
     2                      + oC(im1,j,kp1) * dPE(iim1,jj,kkp1)
     2                      - uC(im1,j,kp1) * oPE(iim1,jj,kkp1)
     2                      - oN(im1,j,kp1) * dPNE(iim1,jj,kkp1)
     2                      - uN(im1,j,kp1) * oPNE(iim1,jj,kkp1)) 

      TMP4_XuW =
     2    - oPNW(ii,jj,kk) * (- uSW(im1,jp1,k) * dPC(iim1,jj,kkp1)
     2                        - uW(im1,jp1,k) * dPN(iim1,jj,kkp1)
     2                        - uS(im1,jp1,k) * dPE(iim1,jj,kkp1)
     2                        - uC(im1,jp1,k) * dPNE(iim1,jj,kkp1)) 

     2    - uPNW(ii,jj,kk) * (- oNE(im2,j,kp1) * dPC(iim1,jj,kkp1)
     2                        - uSW(im1,jp1,kp1) * oPC(iim1,jj,kkp1)
     2                        - oE(im2,jp1,kp1) * dPN(iim1,jj,kkp1)
     2                        - uW(im1,jp1,kp1) * oPN(iim1,jj,kkp1)
     2                        - oN(im1,j,kp1) * dPE(iim1,jj,kkp1)
     2                        - uS(im1,jp1,kp1) * oPE(iim1,jj,kkp1)
     2                        + oC(im1,jp1,kp1) * dPNE(iim1,jj,kkp1)
     2                        - uC(im1,jp1,kp1) * oPNE(iim1,jj,kkp1)) 

      TMP5_XuW =
     2   - oPS(ii,jj,kk) * (- uW(i,jm1,k) * dPSE(iim1,jj,kkp1) 
     2                      - uNW(i,jm1,k) * dPE(iim1,jj,kkp1)) 

     2   - uPS(ii,jj,kk) * (- oE(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2                      - uW(i,jm1,kp1) * oPSE(iim1,jj,kkp1)
     2                      - oNW(i,jm1,kp1) * dPE(iim1,jj,kkp1)
     2                      - uNW(i,jm1,kp1) * oPE(iim1,jj,kkp1)) 

     2   - oPC(ii,jj,kk) * (- uSW(i,j,k) * dPSE(iim1,jj,kkp1) 
     2                      - uW(i,j,k) * dPE(iim1,jj,kkp1)
     2                      - uNW(i,j,k) * dPNE(iim1,jj,kkp1)) 

      TMP6_XuW =
     2   - uPC(ii,jj,kk) * (- oNE(im1,jm1,kp1) * dPSE(iim1,jj,kkp1)
     2                      - uSW(i,j,kp1) * oPSE(iim1,jj,kkp1)
     2                      - oE(im1,j,kp1) * dPE(iim1,jj,kkp1)
     2                      - uW(i,j,kp1) * oPE(iim1,jj,kkp1)
     2                      - oNW(i,j,kp1) * dPNE(iim1,jj,kkp1)
     2                      - uNW(i,j,kp1) * oPNE(iim1,jj,kkp1)) 

     2   - oPN(ii,jj,kk) * (- uSW(i,jp1,k) * dPE(iim1,jj,kkp1) 
     2                      - uW(i,jp1,k) * dPNE(iim1,jj,kkp1)) 

     2   - uPN(ii,jj,kk) * (- oNE(im1,j,kp1) * dPE(iim1,jj,kkp1)
     2                      - uSW(i,jp1,kp1) * oPE(iim1,jj,kkp1)
     2                      - oE(im1,jp1,kp1) * dPNE(iim1,jj,kkp1)
     2                      - uW(i,jp1,kp1) * oPNE(iim1,jj,kkp1))

      XuW(ii,jj,kk) = TMP1_XuW + TMP2_XuW + TMP3_XuW + TMP4_XuW
     2              + TMP5_XuW + TMP6_XuW

c* ********************************************************************
c* *** > uN;
c* ********************************************************************

c* ***XuN(ii,jj,kk) =
      TMP1_XuN =
     2   - oPW(ii,jj,kk) * (- uN(im1,j,k) * dPSW(ii,jjp1,kkp1) 
     2                      - uNE(im1,j,k) * dPS(ii,jjp1,kkp1)) 

     2   - uPW(ii,jj,kk) * (- oN(im1,j,kp1) * dPSW(ii,jjp1,kkp1)
     2                      - uN(im1,j,kp1) * oPSW(ii,jjp1,kkp1)
     2                      - oNE(im1,j,kp1) * dPS(ii,jjp1,kkp1)
     2                      - uNE(im1,j,kp1) * oPS(ii,jjp1,kkp1)) 

     2   - oPNW(ii,jj,kk) * (- uC(im1,jp1,k) * dPSW(ii,jjp1,kkp1)
     2                       - uN(im1,jp1,k) * dPW(ii,jjp1,kkp1)
     2                       - uE(im1,jp1,k) * dPS(ii,jjp1,kkp1)
     2                       - uNE(im1,jp1,k) * dPC(ii,jjp1,kkp1)) 

      TMP2_XuN =
     2   - uPNW(ii,jj,kk) * (  oC(im1,jp1,kp1) * dPSW(ii,jjp1,kkp1)
     2                       - uC(im1,jp1,kp1) * oPSW(ii,jjp1,kkp1)
     2                       - oN(im1,jp1,kp1) * dPW(ii,jjp1,kkp1)
     2                       - uN(im1,jp1,kp1) * oPW(ii,jjp1,kkp1)
     2                       - oE(im1,jp1,kp1) * dPS(ii,jjp1,kkp1)
     2                       - uE(im1,jp1,kp1) * oPS(ii,jjp1,kkp1)
     2                       - oNE(im1,jp1,kp1) * dPC(ii,jjp1,kkp1)
     2                       - uNE(im1,jp1,kp1) * oPC(ii,jjp1,kkp1))

     2    - oPC(ii,jj,kk) * (- uNW(i,j,k) * dPSW(ii,jjp1,kkp1) 
     2                       - uN(i,j,k) * dPS(ii,jjp1,kkp1)
     2                       - uNE(i,j,k) * dPSE(ii,jjp1,kkp1)) 

      TMP3_XuN =
     2   - uPC(ii,jj,kk) * (- oNW(i,j,kp1) * dPSW(ii,jjp1,kkp1)
     2                      - uNW(i,j,kp1) * oPSW(ii,jjp1,kkp1)
     2                      - oN(i,j,kp1) * dPS(ii,jjp1,kkp1) 
     2                      - uN(i,j,kp1) * oPS(ii,jjp1,kkp1)
     2                      - oNE(i,j,kp1) * dPSE(ii,jjp1,kkp1)
     2                      - uNE(i,j,kp1) * oPSE(ii,jjp1,kkp1)) 

     2   - oPN(ii,jj,kk) * (- uW(i,jp1,k) * dPSW(ii,jjp1,kkp1) 
     2                      - uNW(i,jp1,k) * dPW(ii,jjp1,kkp1)
     2                      - uC(i,jp1,k) * dPS(ii,jjp1,kkp1) 
     2                      - uN(i,jp1,k) * dPC(ii,jjp1,kkp1)
     2                      - uE(i,jp1,k) * dPSE(ii,jjp1,kkp1)
     2                      - uNE(i,jp1,k) * dPE(ii,jjp1,kkp1)) 

      TMP4_XuN =
     2   - uPN(ii,jj,kk) * (- oE(im1,jp1,kp1) * dPSW(ii,jjp1,kkp1)
     2                      - uW(i,jp1,kp1) * oPSW(ii,jjp1,kkp1)
     2                      - oNW(i,jp1,kp1) * dPW(ii,jjp1,kkp1)
     2                      - uNW(i,jp1,kp1) * oPW(ii,jjp1,kkp1)
     2                      + oC(i,jp1,kp1) * dPS(ii,jjp1,kkp1)
     2                      - uC(i,jp1,kp1) * oPS(ii,jjp1,kkp1)
     2                      - oN(i,jp1,kp1) * dPC(ii,jjp1,kkp1)
     2                      - uN(i,jp1,kp1) * oPC(ii,jjp1,kkp1)
     2                      - oE(i,jp1,kp1) * dPSE(ii,jjp1,kkp1)
     2                      - uE(i,jp1,kp1) * oPSE(ii,jjp1,kkp1)
     2                      - oNE(i,jp1,kp1) * dPE(ii,jjp1,kkp1)
     2                      - uNE(i,jp1,kp1) * oPE(ii,jjp1,kkp1)) 

      TMP5_XuN =
     2   - oPE(ii,jj,kk) * (- uNW(ip1,j,k) * dPS(ii,jjp1,kkp1) 
     2                      - uN(ip1,j,k) * dPSE(ii,jjp1,kkp1)) 

     2   - uPE(ii,jj,kk) * (- oNW(ip1,j,kp1) * dPS(ii,jjp1,kkp1)
     2                      - uNW(ip1,j,kp1) * oPS(ii,jjp1,kkp1)
     2                      - oN(ip1,j,kp1) * dPSE(ii,jjp1,kkp1)
     2                      - uN(ip1,j,kp1) * oPSE(ii,jjp1,kkp1)) 

     2   - oPNE(ii,jj,kk) * (- uW(ip1,jp1,k) * dPS(ii,jjp1,kkp1)
     2                       - uNW(ip1,jp1,k) * dPC(ii,jjp1,kkp1)
     2                       - uC(ip1,jp1,k) * dPSE(ii,jjp1,kkp1)
     2                       - uN(ip1,jp1,k) * dPE(ii,jjp1,kkp1)) 

      TMP6_XuN =
     2   - uPNE(ii,jj,kk) * (- oE(i,jp1,kp1) * dPS(ii,jjp1,kkp1)
     2                       - uW(ip1,jp1,kp1) * oPS(ii,jjp1,kkp1)
     2                       - oNW(ip1,jp1,kp1) * dPC(ii,jjp1,kkp1)
     2                       - uNW(ip1,jp1,kp1) * oPC(ii,jjp1,kkp1)
     2                       + oC(ip1,jp1,kp1) * dPSE(ii,jjp1,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPSE(ii,jjp1,kkp1)
     2                       - oN(ip1,jp1,kp1) * dPE(ii,jjp1,kkp1)
     2                       - uN(ip1,jp1,kp1) * oPE(ii,jjp1,kkp1))

      XuN(ii,jj,kk) = TMP1_XuN + TMP2_XuN + TMP3_XuN + TMP4_XuN
     2              + TMP5_XuN + TMP6_XuN

c* ********************************************************************
c* *** > uS;
c* ********************************************************************

c* ***XuS(ii,jj,kk) =
      TMP1_XuS =
     2   - oPSW(ii,jj,kk) * (- uS(im1,jm1,k) * dPW(ii,jjm1,kkp1)
     2                       - uC(im1,jm1,k) * dPNW(ii,jjm1,kkp1)
     2                       - uSE(im1,jm1,k) * dPC(ii,jjm1,kkp1)
     2                       - uE(im1,jm1,k) * dPN(ii,jjm1,kkp1)) 

     2   - uPSW(ii,jj,kk) * (- oN(im1,jm2,kp1) * dPW(ii,jjm1,kkp1)
     2                       - uS(im1,jm1,kp1) * oPW(ii,jjm1,kkp1)
     2                       + oC(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2                       - uC(im1,jm1,kp1) * oPNW(ii,jjm1,kkp1)
     2                       - oNW(i,jm2,kp1) * dPC(ii,jjm1,kkp1)
     2                       - uSE(im1,jm1,kp1) * oPC(ii,jjm1,kkp1)
     2                       - oE(im1,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                       - uE(im1,jm1,kp1) * oPN(ii,jjm1,kkp1))

      TMP2_XuS =
     2   - oPW(ii,jj,kk) * (- uS(im1,j,k) * dPNW(ii,jjm1,kkp1) 
     2                      - uSE(im1,j,k) * dPN(ii,jjm1,kkp1)) 

     2   - uPW(ii,jj,kk) * (- oN(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2                      - uS(im1,j,kp1) * oPNW(ii,jjm1,kkp1)
     2                      - oNW(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                      - uSE(im1,j,kp1) * oPN(ii,jjm1,kkp1)) 

     2   - oPS(ii,jj,kk) * (- uSW(i,jm1,k) * dPW(ii,jjm1,kkp1) 
     2                      - uW(i,jm1,k) * dPNW(ii,jjm1,kkp1)
     2                      - uS(i,jm1,k) * dPC(ii,jjm1,kkp1) 
     2                      - uC(i,jm1,k) * dPN(ii,jjm1,kkp1)
     2                      - uSE(i,jm1,k) * dPE(ii,jjm1,kkp1)
     2                      - uE(i,jm1,k) * dPNE(ii,jjm1,kkp1)) 

      TMP3_XuS =
     2   - uPS(ii,jj,kk) * (- oNE(im1,jm2,kp1) * dPW(ii,jjm1,kkp1)
     2                      - uSW(i,jm1,kp1) * oPW(ii,jjm1,kkp1)
     2                      - oE(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2                      - uW(i,jm1,kp1) * oPNW(ii,jjm1,kkp1)
     2                      - oN(i,jm2,kp1) * dPC(ii,jjm1,kkp1)
     2                      - uS(i,jm1,kp1) * oPC(ii,jjm1,kkp1)
     2                      + oC(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                      - uC(i,jm1,kp1) * oPN(ii,jjm1,kkp1)
     2                      - oNW(ip1,jm2,kp1) * dPE(ii,jjm1,kkp1)
     2                      - uSE(i,jm1,kp1) * oPE(ii,jjm1,kkp1)
     2                      - oE(i,jm1,kp1) * dPNE(ii,jjm1,kkp1)
     2                      - uE(i,jm1,kp1) * oPNE(ii,jjm1,kkp1)) 

      TMP4_XuS =
     2   - oPC(ii,jj,kk) * (- uSW(i,j,k) * dPNW(ii,jjm1,kkp1) 
     2                      - uS(i,j,k) * dPN(ii,jjm1,kkp1)
     2                      - uSE(i,j,k) * dPNE(ii,jjm1,kkp1)) 

     2   - uPC(ii,jj,kk) * (- oNE(im1,jm1,kp1) * dPNW(ii,jjm1,kkp1)
     2                      - uSW(i,j,kp1) * oPNW(ii,jjm1,kkp1)
     2                      - oN(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                      - uS(i,j,kp1) * oPN(ii,jjm1,kkp1)
     2                      - oNW(ip1,jm1,kp1) * dPNE(ii,jjm1,kkp1)
     2                      - uSE(i,j,kp1) * oPNE(ii,jjm1,kkp1)) 

      TMP5_XuS =
     2   - oPSE(ii,jj,kk) * (- uSW(ip1,jm1,k) * dPC(ii,jjm1,kkp1)
     2                       - uW(ip1,jm1,k) * dPN(ii,jjm1,kkp1)
     2                       - uS(ip1,jm1,k) * dPE(ii,jjm1,kkp1)
     2                       - uC(ip1,jm1,k) * dPNE(ii,jjm1,kkp1)) 

     2   - uPSE(ii,jj,kk) * (- oNE(i,jm2,kp1) * dPC(ii,jjm1,kkp1)
     2                       - uSW(ip1,jm1,kp1) * oPC(ii,jjm1,kkp1)
     2                       - oE(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                       - uW(ip1,jm1,kp1) * oPN(ii,jjm1,kkp1)
     2                       - oN(ip1,jm2,kp1) * dPE(ii,jjm1,kkp1)
     2                       - uS(ip1,jm1,kp1) * oPE(ii,jjm1,kkp1)
     2                       + oC(ip1,jm1,kp1) * dPNE(ii,jjm1,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPNE(ii,jjm1,kkp1))

      TMP6_XuS =
     2    - oPE(ii,jj,kk) * (- uSW(ip1,j,k) * dPN(ii,jjm1,kkp1) 
     2                       - uS(ip1,j,k) * dPNE(ii,jjm1,kkp1)) 

     2    - uPE(ii,jj,kk) * (- oNE(i,jm1,kp1) * dPN(ii,jjm1,kkp1)
     2                       - uSW(ip1,j,kp1) * oPN(ii,jjm1,kkp1)
     2                       - oN(ip1,jm1,kp1) * dPNE(ii,jjm1,kkp1)
     2                       - uS(ip1,j,kp1) * oPNE(ii,jjm1,kkp1))

      XuS(ii,jj,kk) = TMP1_XuS + TMP2_XuS + TMP3_XuS + TMP4_XuS
     2              + TMP5_XuS + TMP6_XuS

c* ********************************************************************
c* *** > uNE;
c* ********************************************************************

c* ***XuNE(ii,jj,kk) =
      TMP1_XuNE =
     2     oPC(ii,jj,kk) * uNE(i,j,k) * dPSW(iip1,jjp1,kkp1)

     2   - uPC(ii,jj,kk) * (- oNE(i,j,kp1) * dPSW(iip1,jjp1,kkp1)
     2                      - uNE(i,j,kp1) * oPSW(iip1,jjp1,kkp1))

     2   - oPN(ii,jj,kk) * (- uE(i,jp1,k) * dPSW(iip1,jjp1,kkp1)
     2                      - uNE(i,jp1,k) * dPW(iip1,jjp1,kkp1))

     2   - uPN(ii,jj,kk) * (- oE(i,jp1,kp1) * dPSW(iip1,jjp1,kkp1)
     2                      - uE(i,jp1,kp1) * oPSW(iip1,jjp1,kkp1)
     2                      - oNE(i,jp1,kp1) * dPW(iip1,jjp1,kkp1)
     2                      - uNE(i,jp1,kp1) * oPW(iip1,jjp1,kkp1))

     2   - oPE(ii,jj,kk) * (- uN(ip1,j,k) * dPSW(iip1,jjp1,kkp1)
     2                      - uNE(ip1,j,k) * dPS(iip1,jjp1,kkp1))

      TMP2_XuNE =
     2   - uPE(ii,jj,kk) * (- oN(ip1,j,kp1) * dPSW(iip1,jjp1,kkp1)
     2                      - uN(ip1,j,kp1) * oPSW(iip1,jjp1,kkp1)
     2                      - oNE(ip1,j,kp1) * dPS(iip1,jjp1,kkp1)
     2                      - uNE(ip1,j,kp1) * oPS(iip1,jjp1,kkp1))

     2   - oPNE(ii,jj,kk) * (- uC(ip1,jp1,k) * dPSW(iip1,jjp1,kkp1)
     2                       - uN(ip1,jp1,k) * dPW(iip1,jjp1,kkp1)
     2                       - uE(ip1,jp1,k) * dPS(iip1,jjp1,kkp1)
     2                       - uNE(ip1,jp1,k) * dPC(iip1,jjp1,kkp1))

     2   - uPNE(ii,jj,kk) * (  oC(ip1,jp1,kp1) * dPSW(iip1,jjp1,kkp1)
     2                       - uC(ip1,jp1,kp1) * oPSW(iip1,jjp1,kkp1)
     2                       - oN(ip1,jp1,kp1) * dPW(iip1,jjp1,kkp1)
     2                       - uN(ip1,jp1,kp1) * oPW(iip1,jjp1,kkp1)
     2                       - oE(ip1,jp1,kp1) * dPS(iip1,jjp1,kkp1)
     2                       - uE(ip1,jp1,kp1) * oPS(iip1,jjp1,kkp1)
     2                       - oNE(ip1,jp1,kp1) * dPC(iip1,jjp1,kkp1)
     2                       - uNE(ip1,jp1,kp1) * oPC(iip1,jjp1,kkp1))
      XuNE(ii,jj,kk) = TMP1_XuNE + TMP2_XuNE

c* ********************************************************************
c* *** > uNW;
c* ********************************************************************

c* ***XuNW(ii,jj,kk) =
      TMP1_XuNW =
     2   - oPW(ii,jj,kk) * (- uNW(im1,j,k) * dPS(iim1,jjp1,kkp1)
     2                      - uN(im1,j,k) * dPSE(iim1,jjp1,kkp1))

     2   - uPW(ii,jj,kk) * (- oNW(im1,j,kp1) * dPS(iim1,jjp1,kkp1)
     2                      - uNW(im1,j,kp1) * oPS(iim1,jjp1,kkp1)
     2                      - oN(im1,j,kp1) * dPSE(iim1,jjp1,kkp1)
     2                      - uN(im1,j,kp1) * oPSE(iim1,jjp1,kkp1))

     2   - oPNW(ii,jj,kk) * (- uW(im1,jp1,k) * dPS(iim1,jjp1,kkp1)
     2                       - uNW(im1,jp1,k) * dPC(iim1,jjp1,kkp1)
     2                       - uC(im1,jp1,k) * dPSE(iim1,jjp1,kkp1)
     2                       - uN(im1,jp1,k) * dPE(iim1,jjp1,kkp1))

      TMP2_XuNW =
     2   - uPNW(ii,jj,kk) * (- oE(im2,jp1,kp1) * dPS(iim1,jjp1,kkp1)
     2                       - uW(im1,jp1,kp1) * oPS(iim1,jjp1,kkp1)
     2                       - oNW(im1,jp1,kp1) * dPC(iim1,jjp1,kkp1)
     2                       - uNW(im1,jp1,kp1) * oPC(iim1,jjp1,kkp1)
     2                       + oC(im1,jp1,kp1) * dPSE(iim1,jjp1,kkp1)
     2                       - uC(im1,jp1,kp1) * oPSE(iim1,jjp1,kkp1)
     2                       - oN(im1,jp1,kp1) * dPE(iim1,jjp1,kkp1)
     2                       - uN(im1,jp1,kp1) * oPE(iim1,jjp1,kkp1))

     2   + oPC(ii,jj,kk) * uNW(i,j,k) * dPSE(iim1,jjp1,kkp1)

     2   - uPC(ii,jj,kk) * (- oNW(i,j,kp1) * dPSE(iim1,jjp1,kkp1)
     2                      - uNW(i,j,kp1) * oPSE(iim1,jjp1,kkp1))

     2   - oPN(ii,jj,kk) * (- uW(i,jp1,k) * dPSE(iim1,jjp1,kkp1)
     2                      - uNW(i,jp1,k) * dPE(iim1,jjp1,kkp1))

     2   - uPN(ii,jj,kk) * (- oE(im1,jp1,kp1) * dPSE(iim1,jjp1,kkp1)
     2                      - uW(i,jp1,kp1) * oPSE(iim1,jjp1,kkp1)
     2                      - oNW(i,jp1,kp1) * dPE(iim1,jjp1,kkp1)
     2                      - uNW(i,jp1,kp1) * oPE(iim1,jjp1,kkp1))
      XuNW(ii,jj,kk) = TMP1_XuNW + TMP2_XuNW

c* ********************************************************************
c* *** > uSE;
c* ********************************************************************

c* ***XuSE(ii,jj,kk) =
      TMP1_XuSE =
     2   - oPS(ii,jj,kk) * (- uSE(i,jm1,k) * dPW(iip1,jjm1,kkp1)
     2                      - uE(i,jm1,k) * dPNW(iip1,jjm1,kkp1))

     2   - uPS(ii,jj,kk) * (- oNW(ip1,jm2,kp1) * dPW(iip1,jjm1,kkp1)
     2                      - uSE(i,jm1,kp1) * oPW(iip1,jjm1,kkp1) 
     2                      - oE(i,jm1,kp1) * dPNW(iip1,jjm1,kkp1)
     2                      - uE(i,jm1,kp1) * oPNW(iip1,jjm1,kkp1))

     2   + oPC(ii,jj,kk) * uSE(i,j,k) * dPNW(iip1,jjm1,kkp1)

     2   - uPC(ii,jj,kk) * (- oNW(ip1,jm1,kp1) * dPNW(iip1,jjm1,kkp1)
     2                      - uSE(i,j,kp1) * oPNW(iip1,jjm1,kkp1))

      TMP2_XuSE =
     2   - oPSE(ii,jj,kk) * (- uS(ip1,jm1,k) * dPW(iip1,jjm1,kkp1)
     2                       - uC(ip1,jm1,k) * dPNW(iip1,jjm1,kkp1)
     2                       - uSE(ip1,jm1,k) * dPC(iip1,jjm1,kkp1) 
     2                       - uE(ip1,jm1,k) * dPN(iip1,jjm1,kkp1))

     2   - uPSE(ii,jj,kk) * (- oN(ip1,jm2,kp1) * dPW(iip1,jjm1,kkp1)
     2                       - uS(ip1,jm1,kp1) * oPW(iip1,jjm1,kkp1)
     2                       + oC(ip1,jm1,kp1) * dPNW(iip1,jjm1,kkp1)
     2                       - uC(ip1,jm1,kp1) * oPNW(iip1,jjm1,kkp1)
     2                       - oNW(ip2,jm2,kp1) * dPC(iip1,jjm1,kkp1)
     2                       - uSE(ip1,jm1,kp1) * oPC(iip1,jjm1,kkp1)
     2                       - oE(ip1,jm1,kp1) * dPN(iip1,jjm1,kkp1)
     2                       - uE(ip1,jm1,kp1) * oPN(iip1,jjm1,kkp1))

     2   - oPE(ii,jj,kk) * (- uS(ip1,j,k) * dPNW(iip1,jjm1,kkp1)
     2                      - uSE(ip1,j,k) * dPN(iip1,jjm1,kkp1))

     2   - uPE(ii,jj,kk) * (- oN(ip1,jm1,kp1) * dPNW(iip1,jjm1,kkp1)
     2                      - uS(ip1,j,kp1) * oPNW(iip1,jjm1,kkp1)
     2                      - oNW(ip2,jm1,kp1) * dPN(iip1,jjm1,kkp1)
     2                      - uSE(ip1,j,kp1) * oPN(iip1,jjm1,kkp1))
      XuSE(ii,jj,kk) = TMP1_XuSE + TMP2_XuSE

c* ********************************************************************
c* *** > uSW;
c* ********************************************************************

c* ***XuSW(ii,jj,kk) =
      TMP1_XuSW =
     2   - oPSW(ii,jj,kk) * (- uSW(im1,jm1,k) * dPC(iim1,jjm1,kkp1)
     2                       - uW(im1,jm1,k) * dPN(iim1,jjm1,kkp1)
     2                       - uS(im1,jm1,k) * dPE(iim1,jjm1,kkp1)
     2                       - uC(im1,jm1,k) * dPNE(iim1,jjm1,kkp1))

     2   - uPSW(ii,jj,kk) * (- oNE(im2,jm2,kp1) * dPC(iim1,jjm1,kkp1)
     2                       - uSW(im1,jm1,kp1) * oPC(iim1,jjm1,kkp1)
     2                       - oE(im2,jm1,kp1) * dPN(iim1,jjm1,kkp1)
     2                       - uW(im1,jm1,kp1) * oPN(iim1,jjm1,kkp1)
     2                       - oN(im1,jm2,kp1) * dPE(iim1,jjm1,kkp1)
     2                       - uS(im1,jm1,kp1) * oPE(iim1,jjm1,kkp1)
     2                       + oC(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)
     2                       - uC(im1,jm1,kp1) * oPNE(iim1,jjm1,kkp1))

     2   - oPW(ii,jj,kk) * (- uSW(im1,j,k) * dPN(iim1,jjm1,kkp1)
     2                      - uS(im1,j,k) * dPNE(iim1,jjm1,kkp1))

      TMP2_XuSW =
     2   - uPW(ii,jj,kk) * (- oNE(im2,jm1,kp1) * dPN(iim1,jjm1,kkp1)
     2                      - uSW(im1,j,kp1) * oPN(iim1,jjm1,kkp1)
     2                      - oN(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)
     2                      - uS(im1,j,kp1) * oPNE(iim1,jjm1,kkp1))

     2   - oPS(ii,jj,kk) * (- uSW(i,jm1,k) * dPE(iim1,jjm1,kkp1)
     2                      - uW(i,jm1,k) * dPNE(iim1,jjm1,kkp1))

     2   - uPS(ii,jj,kk) * (- oNE(im1,jm2,kp1) * dPE(iim1,jjm1,kkp1)
     2                      - uSW(i,jm1,kp1) * oPE(iim1,jjm1,kkp1)
     2                      - oE(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)
     2                      - uW(i,jm1,kp1) * oPNE(iim1,jjm1,kkp1))

     2   + oPC(ii,jj,kk) * uSW(i,j,k) * dPNE(iim1,jjm1,kkp1)

     2   - uPC(ii,jj,kk) * (- oNE(im1,jm1,kp1) * dPNE(iim1,jjm1,kkp1)
     2                      - uSW(i,j,kp1) * oPNE(iim1,jjm1,kkp1))
      XuSW(ii,jj,kk) = TMP1_XuSW + TMP2_XuSW

c*
c*             *** main loop ***
 12         continue
 11      continue
 10   continue
c*
c*    *** return and end ***
      return
      end
