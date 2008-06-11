c* ///////////////////////////////////////////////////////////////////////////
c* @file    nwjd.f
c* @author  Michael Holst
c* @brief   Nonlinear weighted Jacobi iteration.
c* @version $Id: nwjd.f,v 1.1 2008-06-11 10:47:38 degironc Exp $
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

      subroutine nwjac(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
     2   itmax,iters,errtol,omega,iresid,iadjoint)
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
         call nwjac7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),
     3      x,w1,w2,r,
     4      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (numdia .eq. 27) then
         call nwjac27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3      ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4      x,w1,w2,r,
     5      itmax,iters,errtol,omega,iresid,iadjoint)
      else
         print*,'% NWJAC: invalid stencil type given...'
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nwjac7(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*   nonlinear weighted jacobi-newton iteration (ortega-rheinboldt, page 220)
c*
c* algorithm parameters:
c* ---------------------
c*
c*    U0      = initial guess
c*    Llin    = matrix representing linear part of equations
c*    F       = vector representing rhs of equations
c*    nx      = number of unknowns
c*    itmax   = number of iterations to perform
c*    c_scal  = function representing nonlinear part of equations
c*    dc_scal = function representing derivative of nonlinear 
c*              part of equations
c*    U       = final result
c*
c* algorithm idea:
c* ---------------
c*
c*    the ith nonlinear equation is represented as:
c*
c*       N(UU) = AA * UU + c_scal(UU) - FF
c*
c*    which is solved for UU using newton's method, where:
c* 
c*       AA = ith diagonal entry of Llin
c*       FF = F(i), minus the nondiag of ith eqn of Llin times current U
c*       UU = the unknown, initialized to U(i)
c*
c* algorithm (matlab language):
c* ----------------------------
c*
c* function [U] = wjac_go(U0,Llin,F,nx,itmax);
c*
c* omega = 0.5;
c* U = U0;
c* iters = 0;
c* while (iters < itmax)
c*    iters = iters + 1;
c* 
c*    %%% relax each unknown
c*    U1 = U;
c*    for i = 1:nx
c* 
c*       %%% define the 1d rhs
c*       FF = F(i) - Llin(i,1:nx) * U1 + Llin(i,i) * U1(i);
c*       AA = Llin(i,i);
c*       UU = U1(i);
c*       BB = b(i);
c* 
c*       %%% perform 1d newton's method
c*       %%% note that in the linear case, this is equivalent to:  FF/AA
c*       nitmax = 50;
c*       ntol   = 1.0e-12;
c*       niters = 0;
c*       N      = AA * UU + c_scal(BB,UU) - FF ;
c*       ndenom = N;
c*       if (ndenom == 0) 
c*          ndenom = 1;
c*       end;
c*       nerror = N / ndenom;
c* 
c*       while ((nerror > ntol) & (niters < nitmax))
c*          niters = niters + 1;
c* 
c*          %%% construct jacobian matrix of N
c*          DN = AA + dc_scal(BB,UU);
c* 
c*          %%% solve linear system
c*          change = - N / DN;
c* 
c*          %%% update solution
c*          UU = UU + change;
c* 
c*          %%% construct nonlinear residual
c*          N = AA * UU + c_scal(BB,UU) - FF ;
c* 
c*          %%% compute residual
c*          nerror = N / ndenom;
c*       end
c* 
c*       %%% update the unknown
c*       U(i) = UU;
c*    end
c* 
c*    %%% update the vector of unknowns
c*    U = U1 + omega * (U - U1);
c* end
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k
      double precision omega,errtol,fac
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz),r(nx,ny,nz)
      double precision c_scal,dc_scal
      double precision UU,AA,BB,FF,zNN,DzNN,zntol,zndenom,znerror
      double precision change
      integer          nitmax,niters,ifail_tol,ipkey
c*
cmdir 0 0
c*
c*    *** nonlinear iteration tolerance and itmax ***
      nitmax    = 10
      zntol     = 1.0e-5
      ifail_tol = 0
      ipkey     = ipc(10)
c*
c*    *** do the jacobi iteration itmax times ***
      do 30 iters = 1, itmax
c*
c*       *** do it ***
cmdir 3 1
         do 10 k=2,nz-1
cmdir 3 2
            do 11 j=2,ny-1
cmdir 3 3
               do 12 i=2,nx-1
c*
c*                *** for this unknown, do a 1d newton's method
c*            
c*                *** determine the 1d equation
                  UU = x(i,j,k)
                  AA = oC(i,j,k)
                  BB = cc(i,j,k)
                  FF = fc(i,j,k)
     2               +  oN(i,j,k)        * x(i,j+1,k)
     3               +  oN(i,j-1,k)      * x(i,j-1,k)
     4               +  oE(i,j,k)        * x(i+1,j,k)
     5               +  oE(i-1,j,k)      * x(i-1,j,k)
     6               +  uC(i,j,k-1)      * x(i,j,k-1)
     7               +  uC(i,j,k)        * x(i,j,k+1)
c*
c*                *** evaluate residual of 1d system
                  zNN = AA * UU + c_scal(BB,UU,ipkey) - FF
c*
c*                *** setup for newton's method
                  niters = 0
                  zndenom = zNN
                  if (zndenom .eq. 0.0d0) zndenom = 1.0d0
                  znerror = zNN / zndenom
c*
c*                *** the 1d newton's method ***
                  if (znerror .lt. zntol) goto 17
                  if (niters .gt. nitmax) goto 16
 15               continue
                     niters = niters + 1
c*
c*                   *** construct jacobian matrix of NN ***
                     DzNN = AA + dc_scal(BB,UU,ipkey)
c*
c*                   *** solve the linear system ***
                     change = - zNN / DzNN
c*
c*                   *** update the solution ***
                     UU = UU + change
c*
c*                   *** evaluate residual of 1d system ***
                     zNN = AA * UU + c_scal(BB,UU,ipkey) - FF
c*
c*                   *** compute error ***
                     znerror = zNN / zndenom
c*
c*                   *** tolerance and itmax check ***
                     if (znerror .lt. zntol) goto 17
                     if (niters .gt. nitmax) goto 16
c*
c*                *** loop ***
                  goto 15
c*
c*                *** tolerance not reached ***
 16               continue
                  ifail_tol = ifail_tol + 1
c*
c*                *** tolerance reached ***
 17               continue
c*
c*                *** newton's method complete -- update solution value ***
                  w1(i,j,k) = UU
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
         call nmresid7_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,oE,oN,uC,x,r,w1)
      endif
c*
c*    *** messages ***
      if (ifail_tol .gt. 0) then
         print*,'% NWJAC: 1d newton tolerance failures:    ',ifail_tol
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nwjac27(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*   nonlinear weighted jacobi-newton iteration (ortega-rheinboldt, page 220)
c*
c* algorithm parameters:
c* ---------------------
c*
c*    U0      = initial guess
c*    Llin    = matrix representing linear part of equations
c*    F       = vector representing rhs of equations
c*    nx      = number of unknowns
c*    itmax   = number of iterations to perform
c*    c_scal  = function representing nonlinear part of equations
c*    dc_scal = function representing derivative of nonlinear 
c*              part of equations
c*    U       = final result
c*
c* algorithm idea:
c* ---------------
c*
c*    the ith nonlinear equation is represented as:
c*
c*       N(UU) = AA * UU + c_scal(UU) - FF
c*
c*    which is solved for UU using newton's method, where:
c* 
c*       AA = ith diagonal entry of Llin
c*       FF = F(i), minus the nondiag of ith eqn of Llin times current U
c*       UU = the unknown, initialized to U(i)
c*
c* algorithm (matlab language):
c* ----------------------------
c*
c* function [U] = wjac_go(U0,Llin,F,nx,itmax);
c*
c* omega = 0.5;
c* U = U0;
c* iters = 0;
c* while (iters < itmax)
c*    iters = iters + 1;
c* 
c*    %%% relax each unknown
c*    U1 = U;
c*    for i = 1:nx
c* 
c*       %%% define the 1d rhs
c*       FF = F(i) - Llin(i,1:nx) * U1 + Llin(i,i) * U1(i);
c*       AA = Llin(i,i);
c*       UU = U1(i);
c*       BB = b(i);
c* 
c*       %%% perform 1d newton's method
c*       %%% note that in the linear case, this is equivalent to:  FF/AA
c*       nitmax = 50;
c*       ntol   = 1.0e-12;
c*       niters = 0;
c*       N      = AA * UU + c_scal(BB,UU) - FF ;
c*       ndenom = N;
c*       if (ndenom == 0) 
c*          ndenom = 1;
c*       end;
c*       nerror = N / ndenom;
c* 
c*       while ((nerror > ntol) & (niters < nitmax))
c*          niters = niters + 1;
c* 
c*          %%% construct jacobian matrix of N
c*          DN = AA + dc_scal(BB,UU);
c* 
c*          %%% solve linear system
c*          change = - N / DN;
c* 
c*          %%% update solution
c*          UU = UU + change;
c* 
c*          %%% construct nonlinear residual
c*          N = AA * UU + c_scal(BB,UU) - FF ;
c* 
c*          %%% compute residual
c*          nerror = N / ndenom;
c*       end
c* 
c*       %%% update the unknown
c*       U(i) = UU;
c*    end
c* 
c*    %%% update the vector of unknowns
c*    U = U1 + omega * (U - U1);
c* end
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
      double precision c_scal,dc_scal
      double precision UU,AA,BB,FF,zNN,DzNN,zntol,zndenom,znerror
      double precision change
      integer          nitmax,niters,ifail_tol,ipkey
      double precision tmpO,tmpU,tmpD
c*
cmdir 0 0
c*
c*    *** nonlinear iteration tolerance and itmax ***
      nitmax    = 10
      zntol     = 1.0e-5
      ifail_tol = 0
      ipkey     = ipc(10)
c*
c*    *** do the jacobi iteration itmax times ***
      do 30 iters = 1, itmax
c*
c*       *** do it ***
cmdir 3 1
         do 10 k=2,nz-1
cmdir 3 2
            do 11 j=2,ny-1
cmdir 3 3
               do 12 i=2,nx-1
c*
c*                *** for this unknown, do a 1d newton's method
c*            
c*                *** determine the 1d equation
                  UU = x(i,j,k)
                  AA = oC(i,j,k)
                  BB = cc(i,j,k)
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
                  FF = fc(i,j,k) + tmpO + tmpU + tmpD
c*
c*                *** evaluate residual of 1d system
                  zNN = AA * UU + c_scal(BB,UU,ipkey) - FF
c*
c*                *** setup for newton's method
                  niters = 0
                  zndenom = zNN
                  if (zndenom .eq. 0.0d0) zndenom = 1.0d0
                  znerror = zNN / zndenom
c*
c*                *** the 1d newton's method ***
                  if (znerror .lt. zntol) goto 17
                  if (niters .gt. nitmax) goto 16
 15               continue
                     niters = niters + 1
c*
c*                   *** construct jacobian matrix of NN ***
                     DzNN = AA + dc_scal(BB,UU,ipkey)
c*
c*                   *** solve the linear system ***
                     change = - zNN / DzNN
c*
c*                   *** update the solution ***
                     UU = UU + change
c*
c*                   *** evaluate residual of 1d system ***
                     zNN = AA * UU + c_scal(BB,UU,ipkey) - FF
c*
c*                   *** compute error ***
                     znerror = zNN / zndenom
c*
c*                   *** tolerance and itmax check ***
                     if (znerror .lt. zntol) goto 17
                     if (niters .gt. nitmax) goto 16
c*
c*                *** loop ***
                  goto 15
c*
c*                *** tolerance not reached ***
 16               continue
                  ifail_tol = ifail_tol + 1
c*
c*                *** tolerance reached ***
 17               continue
c*
c*                *** newton's method complete -- update solution value ***
                  w1(i,j,k) = UU
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
         call nmresid27_1s(nx,ny,nz,ipc,rpc,oC,cc,fc,
     2      oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     3      x,r,w1)
      endif
c*
c*    *** messages ***
      if (ifail_tol .gt. 0) then
         print*,'% NWJAC: 1d newton tolerance failures:    ',ifail_tol
      endif
c*
c*    *** return and end ***
      return
      end
