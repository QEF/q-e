c* ///////////////////////////////////////////////////////////////////////////
c* @file    cgmgd.f
c* @author  Michael Holst
c* @brief   Multigrid-Preconditioned CG.
c* @version $Id: cgmgd.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

      subroutine fcgmg(nx,ny,nz,x,iz,w0,w1,w2,w3,
     2   istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3   iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4   w4,w5,w6,
     5   ipc,rpc,pc,ac,cc,fc,tru,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    nested iteration for multilevel preconditioned cg
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iz(50,*),iok,ilev,iinfo,nlev,itmax
      integer          iters,ierror,level,itmxd,nlevd,iterd,iokd,istop
      integer          nx,ny,nz,nxf,nyf,nzf,nxc,nyc,nzc,nlev_real,istpd
      integer          nu1,nu2,mgsmoo,mgsolv
      integer          whichbc(3)
      double precision epsiln,errd,errtol,omega
      double precision x(*),w0(*),w1(*),w2(*),w3(*)
      double precision rpc(*),pc(*),ac(*),cc(*),fc(*),tru(*)
      double precision w4(*),w5(*),w6(*)
c*
c*    *** recover gridsizes ***
      nxf = nx
      nyf = ny
      nzf = nz
      call mkcors(nlev-1,nxf,nyf,nzf,nxc,nyc,nzc)
c*
c*    *** move up grids: interpolate solution to finer, do cgmg ***
      if (iinfo.ne.0) then
         write(6,100)'% FCGMG: starting:  ',nxf,nyf,nzf,nxc,nyc,nzc
 100     format(a,2(2x,' [',i3,',',i3,',',i3,'] '))
      endif
      do 10 level = nlev_real, ilev+1, -1
c*
c*       *** call mv cycle ***
         errd  = 1.0e-5
         itmxd = 1
         nlevd = nlev_real - level + 1
         iterd = 0
         iokd  = 0
         istpd = 1
         if (iinfo .ge. 2) iokd = 2
         call cgmg(nxc,nyc,nzc,x,iz,w0,w1,w2,w3,
     2      istpd,itmxd,iterd,ierror,nlevd,level,nlev_real,mgsolv,
     3      iokd,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4      w4,w5,w6,
     5      ipc,rpc,pc,ac,cc,fc,tru,whichbc)
c*
c*       *** find new grid size ***
         call mkfine(1,nxc,nyc,nzc,nxf,nyf,nzf)
c*
c*       *** interpolate to next finer grid (use correct bc's) ***
         call interp(nxc,nyc,nzc,nxf,nyf,nzf,
     2      x(iz(1,level)),x(iz(1,level-1)),pc(iz(11,level-1)),
     3      whichbc)
c*
c*       *** new grid size ***
         nxc = nxf
         nyc = nyf
         nzc = nzf
 10   continue
c*
c*    *** call mv cycle ***
      level = ilev
      call cgmg(nxf,nyf,nzf,x,iz,w0,w1,w2,w3,
     2   istop,itmax,iters,ierror,nlev,level,nlev_real,mgsolv,
     3   iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4   w4,w5,w6,
     5   ipc,rpc,pc,ac,cc,fc,tru)
c*
c*    *** return and end ***
      return
      end
      subroutine cgmg(nx,ny,nz,x,iz,w0,w1,w2,w3,
     2   istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3   iok,iinfo,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     4   rr,zz,pp,
     5   ipc,rpc,pc,ac,cc,fc,tru,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    multilevel preconditioned cg
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** other declarations ***
      integer          ipc(*),iz(50,*),iok,ilev,iinfo,nlev,level,lev
      integer          itmax,iters,ierror,istop,nu1,nu2,mgsmoo
      integer          itmax_s,iters_s,ierror_s,iok_s,iinfo_s,istop_s
      integer          nx,ny,nz,mgsolv,nlev_real
      integer          whichbc(3)
      double precision errtol,epsiln,rhok1,rhok2,ztmp,omega
      double precision rsden,rsnrm,orsnrm,xnrm1,xnrm2,xdot
      double precision x(*),w0(*),w1(*),w2(*),w3(*)
      double precision rpc(*),pc(*),ac(*),cc(*),fc(*),tru(*)
      double precision rr(*),zz(*),pp(*)
c*
c*    *** recover level information ***
      level = 1
      lev   = (ilev-1)+level
c*
c*    *** do some i/o if requested ***
!      if (iinfo.ne.0) then
!         write(6,100)'% CGMG: starting:   ',nx,ny,nz
! 100     format(a,(2x,' [',i3,',',i3,',',i3,'] '))
!      endif
c*
c*    *** initial wall clock ***
!      if (iok.ne.0) then
!         call prtini(istop)
!         call prtstp(iok,-1,0.0d0,0.0d0,0.0d0)
!      endif
c*
c*    **************************************************************
c*    *** note: if (iok.ne.0) then:  use a stopping test.        ***
c*    ***       else:  use just the itmax to stop iteration.     ***
c*    **************************************************************
c*    *** istop=0 most efficient (whatever it is)                ***
c*    *** istop=1 relative residual                              ***
c*    *** istop=2 rms difference of successive iterates          ***
c*    *** istop=3 relative true error (provided for testing)     ***
c*    **************************************************************
c*
c*    *** compute denominator for stopping criterion ***
      if (istop .eq. 0) then
         rsden = 1.0d0
      elseif (istop .eq. 1) then
         rsden = xnrm1(nx,ny,nz,fc(iz(1,lev)))
      elseif (istop .eq. 2) then
         rsden = dsqrt(dble(nx*ny*nz))
      elseif (istop .eq. 3) then
         rsden = xnrm2(nx,ny,nz,tru(iz(1,lev)))
      elseif (istop .eq. 4) then
         rsden = xnrm2(nx,ny,nz,tru(iz(1,lev)))
      elseif (istop .eq. 5) then
         call matvec(nx,ny,nz,
     2      ipc(iz(5,lev)),rpc(iz(6,lev)),
     3      ac(iz(7,lev)),cc(iz(1,lev)),
     4      tru(iz(1,lev)),w1,whichbc)
         rsden = dsqrt(xdot(nx,ny,nz,tru(iz(1,lev)),w1))
      else
         print*,'% CGMG: bad istop value... '
      endif
      if (rsden.eq.0.0d0) then
         rsden = 1.0d0
         print*,'% CGMG: rhs is zero...'
      endif
      rsnrm = rsden
      orsnrm = rsnrm
!      if (iok.ne.0) call prtstp (iok,0,rsnrm,rsden,orsnrm)
c*
c* *********************************************************************
c* *** begin cg iteration 
c* *********************************************************************
c*
c*    *** compute the initial residual (always required) ***
      call mresid(nx,ny,nz,
     2   ipc(iz(5,lev)),rpc(iz(6,lev)),
     3   ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4   x(iz(1,lev)),rr(iz(1,lev)),whichbc)
c*
c* *********************************************************************
c* ** *** no preconditioning ***
c* ** call xcopy(nx,ny,nz,rr(iz(1,lev)),zz(iz(1,lev)))
c* *********************************************************************
c*    *** multilevel preconditioning ***
c*
c*    *** restrict residual to form rhs of coarse grid systems ***
      call getpre (nx,ny,nz,iz,ilev,nlev_real,rr,pc, whichbc)
c*
c*    *** do a linear multigrid solve of the precond equations ***
      call azeros(nx,ny,nz,zz(iz(1,lev)))
      istop_s   = 1
      itmax_s   = 1
      iters_s   = 0
      ierror_s  = 0
      iinfo_s   = 0
      iok_s     = 0
c* ***if ((iinfo .ge. 2) .and. (ilev .eq. 1)) iok_s = 2
      call mvcs(nx,ny,nz,zz,iz,w0,w1,w2,w3,
     2   istop_s,itmax_s,iters_s,ierror_s,
     3   nlev,ilev,nlev_real,mgsolv,
     4   iok_s,iinfo_s,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     5   ipc,rpc,pc,ac,cc,rr,tru,whichbc)
c* *********************************************************************
c*
c*    *** setup ***
      rhok1 = xdot(nx,ny,nz,zz(iz(1,lev)),rr(iz(1,lev)))
      if (rhok1 .eq. 0.0d0) goto 99
c*
c*    *** do the cg iteration ***
      iters = 0
 30   continue
         iters = iters + 1
c*
c*       *** save iterate if stop test will use it on next iter ***
         if (istop .eq. 2) call xcopy(nx,ny,nz,x(iz(1,lev)),
     2      tru(iz(1,lev)))
c*
c*       *** main computation ***
         if (iters .eq. 1) then
            call xcopy(nx,ny,nz,zz(iz(1,lev)),pp)
         else
            call xaxpy(nx,ny,nz,(rhok2/rhok1),zz(iz(1,lev)),pp)
            call xscal(nx,ny,nz,(rhok1/rhok2),pp)
         endif 
         call matvec(nx,ny,nz,
     2      ipc(iz(5,lev)),rpc(iz(6,lev)),
     3      ac(iz(7,lev)),cc(iz(1,lev)),
     4      pp,w1,whichbc)
         ztmp = rhok1 / xdot(nx,ny,nz,pp,w1)
         call xaxpy(nx,ny,nz,ztmp,pp,x(iz(1,lev)))
         call xaxpy(nx,ny,nz,(-ztmp),w1,rr(iz(1,lev)))
c*
c* *********************************************************************
c* ***** *** no preconditioning ***
c* ***** call xcopy(nx,ny,nz,rr(iz(1,lev)),zz(iz(1,lev)))
c* *********************************************************************
c*       *** multilevel preconditioning ***
c*
c*       *** restrict residual to form rhs of coarse grid systems ***
         call getpre (nx,ny,nz,iz,ilev,nlev_real,rr,pc,whichbc)
c*
c*       *** do a linear multigrid solve of the precond equations ***
         call azeros(nx,ny,nz,zz(iz(1,lev)))
         istop_s   = 1
         itmax_s   = 1
         iters_s   = 0
         ierror_s  = 0
         iinfo_s   = 0
         iok_s     = 0
c* ******if ((iinfo .ge. 2) .and. (ilev .eq. 1)) iok_s = 2
         call mvcs(nx,ny,nz,zz,iz,w0,w1,w2,w3,
     2      istop_s,itmax_s,iters_s,ierror_s,
     3      nlev,ilev,nlev_real,mgsolv,
     4      iok_s,iinfo_s,epsiln,errtol,omega,nu1,nu2,mgsmoo,
     5      ipc,rpc,pc,ac,cc,rr,tru)
c* *********************************************************************
c*
c*       *** new residual ***
         rhok2 = rhok1
         rhok1 = xdot(nx,ny,nz,zz(iz(1,lev)),rr(iz(1,lev)))
         if (rhok1 .eq. 0.0d0) goto 99
c*
c*       *** compute/check the current stopping test ***
         orsnrm = rsnrm
         if (istop .eq. 0) then
            rsnrm = xnrm1(nx,ny,nz,rr(iz(1,lev)))
         elseif (istop .eq. 1) then
            rsnrm = xnrm1(nx,ny,nz,rr(iz(1,lev)))
         elseif (istop .eq. 2) then
            call xcopy(nx,ny,nz,tru(iz(1,lev)),w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x(iz(1,lev)),w1)
            rsnrm = xnrm1(nx,ny,nz,w1)
         elseif (istop .eq. 3) then
            call xcopy(nx,ny,nz,tru(iz(1,lev)),w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x(iz(1,lev)),w1)
            rsnrm = xnrm2(nx,ny,nz,w1)
         elseif (istop .eq. 4) then
            call xcopy(nx,ny,nz,tru(iz(1,lev)),w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x(iz(1,lev)),w1)
            rsnrm = xnrm2(nx,ny,nz,w1)
         elseif (istop .eq. 5) then
            call xcopy(nx,ny,nz,tru(iz(1,lev)),w1)
            call xaxpy(nx,ny,nz,(-1.0d0),x(iz(1,lev)),w1)
            call matvec(nx,ny,nz,
     2         ipc(iz(5,lev)),rpc(iz(6,lev)),
     3         ac(iz(7,lev)),cc(iz(1,lev)),
     4         w1,w2,whichbc)
            rsnrm = dsqrt(xdot(nx,ny,nz,w1,w2))
         else
            print*,'% MVCS: bad istop value... '
         endif
!         if (iok.ne.0) then
!            call prtstp (iok,iters,rsnrm,rsden,orsnrm)
!         endif
c*
c*       *** check iteration count and tolerance ***
         if ((rsnrm/rsden) .le. errtol) goto 99
         if (iters .ge. itmax) goto 99
c*
c*    *** main loop ***
      goto 30
c*
c*    *** return and end ***
 99   continue
      return
      end
      subroutine getpre(nx,ny,nz,iz,lev,nlev_real,r,pc,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    form the residual on all levels to prepare for multilevel prec.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iz(50,*),nx,ny,nz,nlev_real,lev,level
      integer          nxx,nyy,nzz,nxold,nyold,nzold
      integer          whichbc(3)
      double precision pc(*),r(*)
c*
c*    *** setup ***
      nxx    = nx
      nyy    = ny
      nzz    = nz
c*
c*    *** build the (nlev-1) level operators ***
      do 10 level = lev+1, nlev_real
         nxold = nxx
         nyold = nyy
         nzold = nzz
         call mkcors(1,nxold,nyold,nzold,nxx,nyy,nzz)
c*
c*       *** make the coarse grid rhs functions ***
         call restrc(nxold,nyold,nzold,nxx,nyy,nzz,
     2      r(iz(1,level-1)),r(iz(1,level)),pc(iz(11,level-1)),whichbc)
 10   continue
c*
c*    *** end it ***
      return
      end
