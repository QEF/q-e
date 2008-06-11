c* ///////////////////////////////////////////////////////////////////////////
c* @file    cgmgdrvd.f
c* @author  Michael Holst
c* @brief   Driver for CG preconditioned with MG.
c* @version $Id: cgmgdrvd.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

      subroutine cgmgdriv(iparm,rparm,iwork,rwork,u,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    multilevel preconditioned cg driver
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** input parameters ***
      integer          iparm(*),iwork(*)
      double precision rparm(*),rwork(*),u(*)
      double precision xf(*),yf(*),zf(*),gxcf(*),gycf(*),gzcf(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*),tcf(*)
      integer          whichbc(3)
c*
c*    *** variables returned from mgsz ***
      integer          nxc,nyc,nzc,nf,nc,narr,narrc,n_rpc
      integer          n_iz,n_ipc,iretot,iintot
c*
c*    *** misc variables ***
      integer          nrwk,niwk,nx,ny,nz,nlev,ierror,maxlev,mxlv
      integer          mgcoar,mgdisc,mgsolv
      integer          k_iz,k_w1,k_w2
      integer          k_ipc,k_rpc,k_ac,k_cc,k_fc,k_pc
c*
c*    *** decode some parameters ***
      nrwk   = iparm(1)
      niwk   = iparm(2)
      nx     = iparm(3)
      ny     = iparm(4)
      nz     = iparm(5)
      nlev   = iparm(6)
c*
c*    *** some checks on input ***
      if ((nlev.le.0).or.(nx.le.0).or.(ny.le.0).or.(nz.le.0)) then
!         print*,'% CGMGDRIV:  nx,ny,nz, and nlev must be positive...'
         ierror = -1
         iparm(51) = ierror 
         return
      endif
      mxlv = maxlev(nx,ny,nz)
      if (nlev.gt.mxlv) then
!         print*,'% CGMGDRIV:  max levels for your grid size is: ',mxlv
         ierror = -2
         iparm(51) = ierror 
         return
      endif
c*
c*    *** basic grid sizes, etc. ***
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsolv = iparm(21)
      call mgsz(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev,nxc,nyc,nzc,
     2   nf,nc,narr,narrc,n_rpc,n_iz,n_ipc,iretot,iintot)
c*
c*    *** allocate space for two additional work vectors ***
!      print*,'% CGMGDRIV: real    work space must be: ', iretot + 2*nf
c*
c*    *** some more checks on input ***
      if ((nrwk.lt.iretot) .or. (niwk.lt.iintot)) then
!         print*,'% CGMGDRIV: real    work space must be: ',iretot
!         print*,'% CGMGDRIV: real    work space is: ',nrwk
!         print*,'% CGMGDRIV: integer work space must be: ',iintot
!         print*,'% CGMGDRIV: integer work space is: ',niwk
         ierror = -3
         iparm(51) = ierror 
         return
      endif
c*
c*    *** split up the integer work array ***
      k_iz   = 1
      k_ipc  = k_iz   + n_iz
c*
c*    *** split up the real work array ***
      k_rpc  = 1
      k_cc   = k_rpc  + n_rpc
      k_fc   = k_cc   + narr
      k_w1   = k_fc   + narr
      k_w2   = k_w1   + nf
      k_pc   = k_w2   + nf
      k_ac   = k_pc   + 27*narrc
c* ***k_ac_after = 4*nf + 14*narrc
c*
c*    *** call the multigrid driver ***
      call cgmgdriv2(iparm,rparm,nx,ny,nz,u,iwork(k_iz),
     2   rwork(k_w1),rwork(k_w2),
     3   iwork(k_ipc),rwork(k_rpc),
     4   rwork(k_pc),rwork(k_ac),rwork(k_cc),rwork(k_fc),
     5   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,whichbc)
c*
c*    *** return and end ***
      return
      end
      subroutine cgmgdriv2(iparm,rparm,nx,ny,nz,u,iz,
     2   w1,w2,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    this routine uses a conjugate gradient iteration combined with a
c*    multigrid v-cycle preconditioner to solve the three-dimensional, 
c*    2nd order elliptic partial differential equation:
c*
c*         lu = f, u in omega
c*          u = g, u on boundary of omega
c*    where
c*
c*         omega = [xmin,xmax]x[ymin,ymax]x[zmin,zmax]
c*
c*    the code requires the operator in the form:
c*
c*         - \nabla \cdot (a \nabla u) + c u = f
c*
c*    with
c*
c*        a(x,y,z),f(x,y,z), scalar functions (possibly discontinuous)
c*        on omega.  (discontinuities must be along fine grid lines).
c*        boundary function g(x,y,z) is smooth on boundary of omega.
c*
c*        the function c(x,y,z) and varies (possibly discontinuously) 
c*        with the spatial position.
c*
c* user inputs:
c*
c*    the user must provide the coefficients of the differential
c*    operator, some initial parameter settings in an integer and a
c*    real parameter array, and various work arrays.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
c*
c*    *** input parameters ***
      integer          iparm(*),ipc(*),iz(*)
      double precision rparm(*),rpc(*),pc(*),ac(*),cc(*),fc(*)
      double precision u(*),w1(*),w2(*)
      double precision xf(*),yf(*),zf(*),gxcf(*),gycf(*),gzcf(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*),tcf(*)
      integer          whichbc(3)
c*
c*    *** misc variables ***
      integer          mgkey,nlev,itmax,iok,iinfo,istop,ipkey,nu1,nu2
      integer          nx,ny,nz,ilev,ido,iters,ierror,nlev_real,ibound
      integer          mgprol,mgcoar,mgsolv,mgdisc,mgsmoo,mode
      double precision epsiln,epsmac,errtol,omegal,omegan
      double precision bf,oh,tsetupf,tsetupc,tsolve
c*
c*    *** decode the iparm array ***
      nlev   = iparm(6)
      nu1    = iparm(7)
      nu2    = iparm(8)
      mgkey  = iparm(9)
      itmax  = iparm(10)
      istop  = iparm(11)
      iinfo  = iparm(12)
      ipkey  = iparm(14)
      mgprol = iparm(17)
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsmoo = iparm(20)
      mgsolv = iparm(21)
      errtol = rparm(1)
      omegal = rparm(9)
      omegan = rparm(10)
c*
c*    *** intitialize the iteration timer ***
      call prtstp(0,-99,0.0d0,0.0d0,0.0d0)
c*
c*    *** build the multigrid data structure in iz ***
      call buildstr (nx,ny,nz,nlev,iz)
c*
c*    *** start timer ***
      call tstart(bf,oh)
c*
c*    *** build op and rhs on fine grid ***
      ido = 0
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call tstop(bf,oh,tsetupf)
!      print*,'% CGMGDRIV2: fine problem setup time: ',tsetupf
c*
c*    *** start timer ***
      call tstart(bf,oh)
c*
c*    *** build op and rhs on all coarse grids ***
      ido = 1
      call buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    *** stop timer ***
      call tstop(bf,oh,tsetupc)
!      print*,'% CGMGDRIV2: coarse problem setup time: ',tsetupc
c*
c* ******************************************************************
c* *** this overwrites the rhs array provided by pde specification
c* ****** compute an algebraically produced rhs for the given tcf ***
      mode = 0
      if ((istop .eq. 4) .or. (istop .eq. 5)) then
         call buildALG (nx,ny,nz,mode,nlev,iz,
     2      ipc,rpc,ac,cc,ccf,tcf,fc,fcf)
      endif
c* ******************************************************************
c*
c*    *** determine machine epsilon ***
      epsiln = epsmac(0)
c*
c*    *** impose zero dirichlet boundary conditions (now in source fcn) ***
      call fbound00(nx,ny,nz,u,whichbc)
c      call fbound(1,nx,ny,nz,u,gxcf,gycf,gzcf)
c*
c*    *** MATLAB ***
!      print*,' cgmg = [ '
c*
c*    *** start timer ***
      call tstart(bf,oh)
c*
c*    *** call specified multigrid method ***
      nlev_real = nlev
      iok  = 1
      ilev = 1
      if (mgkey .eq. 0) then
         call cgmg(nx,ny,nz,u,iz,ccf,fcf,w1,w2,
     2      istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3      iok,iinfo,epsiln,errtol,omegal,nu1,nu2,mgsmoo,
     4      a1cf,a2cf,a3cf,
     5      ipc,rpc,pc,ac,cc,fc,tcf,whichbc)
      else if (mgkey .eq. 1) then
         call fcgmg(nx,ny,nz,u,iz,ccf,fcf,w1,w2,
     2      istop,itmax,iters,ierror,nlev,ilev,nlev_real,mgsolv,
     3      iok,iinfo,epsiln,errtol,omegal,nu1,nu2,mgsmoo,
     4      a1cf,a2cf,a3cf,
     5      ipc,rpc,pc,ac,cc,fc,tcf)
      else
!         print*,'% CGMGDRIV2: bad mgkey given '
      endif
c*
c*    *** stop timer ***
      call tstop(bf,oh,tsolve)
!      print*,'% CGMGDRIV2: solve time: ',tsolve
c*
c*    *** MATLAB ***
!      write(*,100) 'cgmg_sf',tsetupf,'cgmg_sc',tsetupc,
!     2   'cgmg_st',(tsetupf+tsetupc),'cgmg_so',tsolve
! 100  format(' ];',4(' ',a7,'=',1pe9.3,';'))
c*
c*    *** restore boundary conditions ***
      !ibound = 1
      !call fbound(ibound,nx,ny,nz,u,gxcf,gycf,gzcf)
c*
c*    *** return and end ***
      return
      end
