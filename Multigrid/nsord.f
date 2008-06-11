c* ///////////////////////////////////////////////////////////////////////////
c* @file    nsord.f
c* @author  Michael Holst
c* @brief   Nonlinear SOR iteration.
c* @version $Id: nsord.f,v 1.1 2008-06-11 10:47:38 degironc Exp $
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

      subroutine nsor(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,
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
         call nsor7(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),
     3      x,w1,w2,r,
     4      itmax,iters,errtol,omega,iresid,iadjoint)
      elseif (numdia .eq. 27) then
         call nsor27(nx,ny,nz,ipc,rpc,ac(1,1),cc,fc,
     2      ac(1,2),ac(1,3),ac(1,4),ac(1,5),ac(1,6),ac(1,7),ac(1,8),
     3      ac(1,9),ac(1,10),ac(1,11),ac(1,12),ac(1,13),ac(1,14),
     4      x,w1,w2,r,
     5      itmax,iters,errtol,omega,iresid,iadjoint)
      else
         print*,'% NSOR: invalid stencil type given...'
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nsor7(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*    fast 7 diagonal nonlinear red/black sor routine.
c*
c*    note that the interior grid points from 2,...,nx-1, etc.
c*    we then begin coloring with a red point, and do a point red/black
c*    coloring.
c*
c*    the red points are:  
c*
c*       if ((j even) and (k even)) then 
c*          begin row at first point, i=2 (i even), or ioff = 0
c*       else if ((j odd) and (k even)) then
c*          begin row at second point, i=3 (i odd), or ioff = 1
c*       else if ((j even) and (k odd)) then
c*          begin row at second point, i=3 (i odd), or ioff = 1
c*       else if ((j odd) and (k odd)) then
c*          begin row at first point, i=2 (i even), or ioff = 0
c*       endif
c*       then: begin row at:  i=2+ioff
c*
c*    the appropriate ioff function for the red points is then:
c*         ioff = dabs(mod(j,2)-mod(k,2))
c*
c*    the appropriate ioff function for the black points is then:
c*         ioff = 1 - dabs(mod(j,2)-mod(k,2))
c*
c*
c*    alternatively, the red points are:
c*
c*       those whose indices add up to an even number.
c*       to see this, consider that all surrounding points are
c*       only one index different, hence the sum will differ by one.
c*       thus, if a given point has an even sum, then the surrounding
c*       points will have an odd sum.
c*
c*    thus, the black points are:
c*
c*       therefore those whose indices add up to an odd number.
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
c* function [U,iterss,errors] = Zgs(U0,UTRUE,a,b,Llin,F,nx,tol,itmax,stp,plt);
c* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c* %%% nonlinear sor-newton iteration (ortega-rheinboldt, page 219)
c* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c*
c* %%% initial guess
c* iters = 0;
c* error = 1.0d0;
c* U     = U0;
c* Lnon  = c_vec(b,U,nx);
c* N     = Llin * U + Lnon - F ;
c* denom = 1.0d0; if (stp == 1) denom = norm(N); end;
c*
c* %%% loop
c* while ((error > tol) & (iters < itmax))
c*    iters = iters + 1;
c*
c*    %%% save old guys
c*    Lnon_OLD = Lnon;
c*    U_OLD    = U;
c*
c*    %%% relax each unknown
c*    for i = 1:nx
c*
c*       %%% define the 1d rhs
c*       FF = F(i) - Llin(i,1:nx) * U + Llin(i,i) * U(i);
c*       AA = Llin(i,i);
c*       UU = U(i);
c*       BB = b(i);
c*
c*       %%% perform 1d newton's method
c*       %%% note that in the linear case, this is equivalent to:  FF/AA
c*       nitmax = 10;
c*       ntol   = 1.0e-7;
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
c*    %%% new residual %%%
c*    Lnon  = c_vec(b,U,nx);
c*    N     = Llin * U + Lnon - F ;
c*    delta = zeros(U);
c*    delta_s = zeros(U);
c*
c* end
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          i,j,k,ioff
      double precision omega,errtol,fac
      double precision rpc(*),oE(nx,ny,nz),oN(nx,ny,nz),uC(nx,ny,nz)
      double precision fc(nx,ny,nz),oC(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),w1(nx,ny,nz),w2(nx,ny,nz)
      double precision r(nx,ny,nz)
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
c*    *** do the sor iteration itmax times ***
      fac   = 1.0d0 - omega
      do 30 iters = 1, itmax
c*
c*       *** do the red points ***
cmdir 3 1
         do 10 k=2,nz-1
cmdir 3 2
            do 11 j=2,ny-1
               ioff = mod((j+k+2),2)
cmdir 3 3
               do 12 i=2+ioff,nx-1, 2
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
                  x(i,j,k) = x(i,j,k)*fac + omega*UU
 12            continue
 11         continue
 10      continue
c*
c*       *** do the black points ***
cmdir 3 1
         do 20 k=2,nz-1
cmdir 3 2
            do 21 j=2,ny-1
               ioff = 1 - mod((j+k+2),2)
cmdir 3 3
               do 22 i=2+ioff,nx-1, 2
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
                  if (znerror .lt. zntol) goto 27
                  if (niters .gt. nitmax) goto 26
 25               continue
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
                     if (znerror .lt. zntol) goto 27
                     if (niters .gt. nitmax) goto 26
c*
c*                *** loop ***
                  goto 25
c*
c*                *** tolerance not reached ***
 26               continue
                  ifail_tol = ifail_tol + 1
c*
c*                *** tolerance reached ***
 27               continue
c*
c*                *** newton's method complete -- update solution value ***
                  x(i,j,k) = x(i,j,k)*fac + omega*UU
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
         print*,'% NSOR: 1d newton tolerance failures:     ',ifail_tol
      endif
c*
c*    *** return and end ***
      return
      end
      subroutine nsor27(nx,ny,nz,ipc,rpc,
     2   oC,cc,fc,
     3   oE,oN,uC,oNE,oNW,uE,uW,uN,uS,uNE,uNW,uSE,uSW,
     4   x,w1,w2,r,
     5   itmax,iters,errtol,omega,iresid,iadjoint)
c* *********************************************************************
c* purpose: 
c*
c*    27 diagonal nonlinear sor routine.
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
c* function [U,iterss,errors] = Zgs(U0,UTRUE,a,b,Llin,F,nx,tol,itmax,stp,plt);
c* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c* %%% nonlinear sor-newton iteration (ortega-rheinboldt, page 219)
c* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c*
c* %%% initial guess
c* iters = 0;
c* error = 1.0d0;
c* U     = U0;
c* Lnon  = c_vec(b,U,nx);
c* N     = Llin * U + Lnon - F ;
c* denom = 1.0d0; if (stp == 1) denom = norm(N); end;
c*
c* %%% loop
c* while ((error > tol) & (iters < itmax))
c*    iters = iters + 1;
c*
c*    %%% save old guys
c*    Lnon_OLD = Lnon;
c*    U_OLD    = U;
c*
c*    %%% relax each unknown
c*    for i = 1:nx
c*
c*       %%% define the 1d rhs
c*       FF = F(i) - Llin(i,1:nx) * U + Llin(i,i) * U(i);
c*       AA = Llin(i,i);
c*       UU = U(i);
c*       BB = b(i);
c*
c*       %%% perform 1d newton's method
c*       %%% note that in the linear case, this is equivalent to:  FF/AA
c*       nitmax = 10;
c*       ntol   = 1.0e-7;
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
c*    %%% new residual %%%
c*    Lnon  = c_vec(b,U,nx);
c*    N     = Llin * U + Lnon - F ;
c*    delta = zeros(U);
c*    delta_s = zeros(U);
c*
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
c*    *** do the sor iteration itmax times ***
      fac   = 1.0d0 - omega
      do 30 iters = 1, itmax
c*
c*       *** do all of the points ***
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
                  x(i,j,k) = x(i,j,k)*fac + omega*UU
 12            continue
 11         continue
 10      continue
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
         print*,'% NSOR: 1d newton tolerance failures:     ',ifail_tol
      endif
c*
c*    *** return and end ***
      return
      end
