c* ///////////////////////////////////////////////////////////////////////////
c* @file    cgd.f
c* @author  Michael Holst
c* @brief   Classical (Hestenes-Stiefel) Conjugate Gradient Method.
c* @version $Id: cgd.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

      subroutine cghs(nx,ny,nz,ipc,rpc,ac,cc,fc,x,p,ap,r,
     2   itmax,iters,errtol,omega,iresid,iadjoint,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    this routine solves the spd ax=b using conjugate gradients (cghs).
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),itmax,iters,iresid,iadjoint,nx,ny,nz
      integer          whichbc(3)
      double precision omega,errtol,rsnrm,pAp,denom,xnrm2
      double precision rpc(*),ac(nx,ny,nz,*),fc(nx,ny,nz),cc(nx,ny,nz)
      double precision x(nx,ny,nz),r(nx,ny,nz),p(nx,ny,nz),ap(nx,ny,nz)
      double precision rhok1,rhok2,alpha,beta,xdot,rp
c*
cmdir 0 0
c*
c*    *** setup for the looping ***
      iters = 0
      if ((iters .ge. itmax) .and. (iresid .eq. 0)) goto 99
      call mresid(nx,ny,nz,ipc,rpc,ac,cc,fc,x,r,whichbc)
      denom = xnrm2(nx,ny,nz,r)
      if (denom .eq. 0.0e0) goto 99
      if (iters .ge. itmax) goto 99
 30   continue
c*
c*       *** compute/check the current stopping test ***
         rhok2 = xdot(nx,ny,nz,r,r)
         rsnrm = dsqrt(rhok2)
c* ******print*,'% CGHS: iters, rsnrm = ',iters,rsnrm/denom
         !write(*,*) rsnrm/denom, denom, errtol
         if (rsnrm/denom .le. errtol) goto 99
         if (iters .ge. itmax) goto 99
c*
c*       *** form new direction vector from old one and residual ***
         if (iters .eq. 0) then
            call xcopy(nx,ny,nz,r,p)
         else
            beta = rhok2 / rhok1
            call xaxpy(nx,ny,nz,((1.0d0)/beta),r,p)
            call xscal(nx,ny,nz,beta,p)
         endif
c*
c*       *** linear case: alpha which minimizes energy norm of error ***
         call matvec(nx,ny,nz,ipc,rpc,ac,cc,p,ap,whichbc)
         pAp = xdot(nx,ny,nz,p,ap)
         rp = xdot(nx,ny,nz,r,p)
         alpha = rp / pAp
c*
c*       *** save rhok2 for next iteration ***
         rhok1 = rhok2
c*
c*       *** update solution in direction p of length alpha ***
         call xaxpy(nx,ny,nz,alpha,p,x)
c*
c*       *** update residual ***
         call xaxpy(nx,ny,nz,(-alpha),ap,r)
c*
c*       *** some bookkeeping ***
         iters = iters + 1
      goto 30
c*
c*    *** return and end ***
 99   continue
      return
      end

