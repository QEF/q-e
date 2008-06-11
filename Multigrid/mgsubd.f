c* ///////////////////////////////////////////////////////////////////////////
c* @file    mgsubd.f
c* @author  Michael Holst
c* @brief   Supporting routines for CS and FAS multigrid algorithms.
c* @version $Id: mgsubd.f,v 1.1 2008-06-11 10:47:38 degironc Exp $
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

      integer function maxlev(n1,n2,n3)
c* *********************************************************************
c* purpose:
c*
c*    find maximum multigrid possible coarsenning common to three grid
c*    sizes:  n1,n2,n3.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          n1,n2,n3,n1c,n2c,n3c,lev,iden,idone
c*
c*    *** fine the common level ***
      idone = 0
      lev = 0
 10   continue
         lev = lev + 1
         iden = 2**(lev-1)
         n1c = (n1-1)/iden + 1
         n2c = (n2-1)/iden + 1
         n3c = (n3-1)/iden + 1
         if ( ((n1c-1)*iden .ne. (n1-1)) .or. (n1c .le. 2) ) idone = 1
         if ( ((n2c-1)*iden .ne. (n2-1)) .or. (n2c .le. 2) ) idone = 1
         if ( ((n3c-1)*iden .ne. (n3-1)) .or. (n3c .le. 2) ) idone = 1 
      if (idone .ne. 1) goto 10
      maxlev = lev-1
c*
c*    *** return and end ***
      return
      end
      subroutine mkcors(numlev,nxold,nyold,nzold,nxnew,nynew,nznew)
c* *********************************************************************
c* purpose:
c*
c*    compute the number of grid points in the coarser grid, given the 
c*    number of grid points in a finger grid in each direction.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          numlev,nxold,nyold,nzold,nxnew,nynew,nznew
      integer          nxtmp,nytmp,nztmp,i
c*
c*    *** determine the coarser grid ***
      nxnew = nxold
      nynew = nyold
      nznew = nzold
      do 10 i = 1, numlev
         nxtmp = nxnew
         nytmp = nynew
         nztmp = nznew
         call corsr(nxtmp,nxnew)
         call corsr(nytmp,nynew)
         call corsr(nztmp,nznew)
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine corsr(nold,nnew)
c* *********************************************************************
c* purpose:
c*
c*    compute the number of grid points in the coarser grid, given the 
c*    number of grid points in a finer grid.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nold,nnew
c*
c*    *** find the coarser grid size ***
      nnew = (nold - 1) / 2 + 1
c*
c*    *** check a few things ***
      if (((nnew-1)*2).ne.(nold-1)) 
     2   print*,'% CORSR:  may not corsen grid this far... '
      if (nnew .lt. 1)
     2   print*,'% CORSR:  have corsenned grid below zero... '
c*
c*    *** return and end ***
      return
      end
      subroutine mkfine(numlev,nxold,nyold,nzold,nxnew,nynew,nznew)
c* *********************************************************************
c* purpose:
c*
c*    compute the number of grid points in the finer grid, given the 
c*    number of grid points in a coarser grid in each direction.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          numlev,nxold,nyold,nzold,nxnew,nynew,nznew
      integer          nxtmp,nytmp,nztmp,i
c*
c*    *** determine the finer grid ***
      nxnew = nxold
      nynew = nyold
      nznew = nzold
      do 10 i = 1, numlev
         nxtmp = nxnew
         nytmp = nynew
         nztmp = nznew
         call finer(nxtmp,nxnew)
         call finer(nytmp,nynew)
         call finer(nztmp,nznew)
 10   continue
c*
c*    *** return and end ***
      return
      end
      subroutine finer(nold,nnew)
c* *********************************************************************
c* purpose:
c*
c*    compute the number of grid points in the finer grid, given the 
c*    number of grid points in a coarser grid.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nold,nnew
c*
c*    *** find the finer grid size ***
      nnew = (nold - 1) * 2 + 1
c*
c*    *** return and end ***
      return
      end
      function ivariv (nu,level)
c* *********************************************************************
c* purpose:
c*
c*    this routine defines the number of smoothings for a particular
c*    level in a variable-v-cycle method.
c*
c*    possible definitions:
c*       ivariv = nu * 2**(level - 1)
c*       ivariv = nu * level
c*       ivariv = nu + (level - 1)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ivariv,nu,level
c*
c* ** *** variable V-cycle ***
c* ** ivariv = nu * 2**(level - 1)
c*
c* ** *** standard V-cycle ***
      ivariv = nu
c*
c*    *** return and end ***
      return
      end
      subroutine prtini(istop)
c* *********************************************************************
c* purpose:
c*
c*    this routine prints out some info and such from inside multigrid.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          istop
      character*20     str0a,str1a,str2a,str3a
      character*20     str0b,str1b,str2b,str3b
      character*20     str0c,str1c,str2c,str3c
c*
c*    *** do some i/o ***
      str0a = 'iteration  '
      str0b = 'count      '
      str0c = '---------  '
      if (istop .eq. 0) then
         str1a = 'absol resid'
         str1b = 'discr(1-nm)'
         str1c = '-----------'
      elseif (istop .eq. 1) then
         str1a = 'relat resid'
         str1b = 'discr(1-nm)'
         str1c = '-----------'
      elseif (istop .eq. 2) then
         str1a = 'rms change '
         str1b = 'discr(1-nm)'
         str1c = '---------- '
      elseif (istop .eq. 3) then
         str1a = 'relat error'
         str1b = 'conti(2-nm)'
         str1c = '-----------'
      elseif (istop .eq. 4) then
         str1a = 'relat error'
         str1b = 'discr(2-nm)'
         str1c = '-----------'
      elseif (istop .eq. 5) then
         str1a = 'relat error'
         str1b = 'discr(A-nm)'
         str1c = '-----------'
      else
         print*,'% PRTINI: bad istop value... '
      endif
      str2a = 'contraction'
      str2b = 'number     '
      str2c = '-----------'
      str3a = 'wall '
      str3b = 'clock'
      str3c = '-----'
      write(6,100) str0c,str1c,str2c,str3c
      write(6,100) str0a,str1a,str2a,str3a
      write(6,100) str0b,str1b,str2b,str3b
      write(6,100) str0c,str1c,str2c,str3c
c*
c*    *** format statements ***
 100  format('% ',a12,1x,a12,3x,a12,3x,a12)
c*
c*    *** return and end ***
      return
      end
      subroutine prtstp(iok,iters,rsnrm,rsden,orsnrm)
c* *********************************************************************
c* purpose:
c*
c*    this routine prints out some info and such from inside multigrid.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iok,iters
      double precision rsnrm,rsden,orsnrm
      double precision relres,contrac
      double precision bf,oh,cputme
      save             bf,oh,cputme
c*
c*    *** initializing timer ***
      if (iters .eq. -99) then
!         print*,'% PRTSTP: initializing timer...'
         call tstart(bf,oh)
         cputme = 0.0d0
         goto 99
c*
c*    *** setup for the iteration ***
      elseif (iters .eq. -1) then
         call tstop(bf,oh,cputme)
         if (iok .eq. 1) then
            write(6,100) -1,0.0d0,0.0d0,cputme
         elseif (iok .eq. 2) then
            write(6,110) -1,0.0d0,0.0d0,cputme
         endif
         goto 99
c*
c*    *** during the iteration ***
      else
c*
c*       *** stop the timer ***
         call tstop(bf,oh,cputme)
c*
c*       *** relative residual ***
         if (rsden .eq. 0.0d0) then
            relres = 1.0e6
            print*,'% PRTSTP: avoided division by zero'
         else
            relres = rsnrm/rsden
         endif
c*
c*       *** contraction number ***
         if (orsnrm .eq. 0.0d0) then
            contrac = 1.0e6
            print*,'% PRTSTP: avoided division by zero'
         else
            contrac = rsnrm/orsnrm
         endif
c*
c*       *** the i/o ***
         if (iok .eq. 1) then
            write(6,100) iters,relres,contrac,cputme
         elseif (iok .eq. 2) then
            write(6,110) iters,relres,contrac,cputme
         endif
      endif
c*
c*    *** format statements ***
 100  format(2x,  i5,8x,1pe11.5,4x,1pe11.5,4x,1pe8.2,20x,'%%%')
 110  format('% ',i5,8x,1pe11.5,4x,1pe11.5,4x,1pe8.2)
c*
c*    *** return and end ***
 99   continue
      return
      end
      subroutine buildstr (nx,ny,nz,nlev,iz)
c* *********************************************************************
c* purpose:  
c*
c*    build the nested operator framework in the array iz
c*
c*    note: iz(50,i) indexes into the gridfcn arrays 
c*       for each level i=(1,...,nlev+1) as follows:
c*
c*          fun(i)    = fun (iz(1,i))
c*          bndx(i)   = bndx(iz(2,i))
c*          bndy(i)   = bndy(iz(3,i))
c*          bndz(i)   = bndz(iz(4,i))
c*          ipc(i)    = ipc(iz(5,i))
c*          rpc(i)    = rpc(iz(6,i))
c*          oper(i)   = oper(iz(7,i))
c*          grdx(i)   = brdx(iz(8,i))
c*          grdy(i)   = brdy(iz(9,i))
c*          grdz(i)   = brdz(iz(10,i))
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iz(50,*),nx,ny,nz,nlev,lev,n
      integer          nxold,nyold,nzold,nxnew,nynew,nznew
c*
c*    *** setup ***
      nxnew  = nx
      nynew  = ny
      nznew  = nz
      n      = nxnew*nynew*nznew
c*
c*    *** start with level 1 ***
      lev = 1
c*
c*    *** mark beginning of everything at level 1 ***
      iz(1,lev)  = 1
      iz(2,lev)  = 1
      iz(3,lev)  = 1
      iz(4,lev)  = 1
      iz(5,lev)  = 1
      iz(6,lev)  = 1
      iz(7,lev)  = 1
      iz(8,lev)  = 1
      iz(9,lev)  = 1
      iz(10,lev) = 1
      iz(11,lev) = 1
c*
c*    *** mark beginning of everything at level 2 ***
      iz(1,lev+1)  = iz(1,lev)  + n
      iz(2,lev+1)  = iz(2,lev)  + 4*nynew*nznew
      iz(3,lev+1)  = iz(3,lev)  + 4*nxnew*nznew
      iz(4,lev+1)  = iz(4,lev)  + 4*nxnew*nynew
      iz(5,lev+1)  = iz(5,lev)  + 100
      iz(6,lev+1)  = iz(6,lev)  + 100
      iz(8,lev+1)  = iz(8,lev)  + nxnew
      iz(9,lev+1)  = iz(9,lev)  + nynew
      iz(10,lev+1) = iz(10,lev) + nznew
c* *********************************************************************
c* ***NOTE: we mark operator offsets as we build the operators ***
c* ***iz(7,lev+1)  = iz(7,lev)  + 4*n
c* *********************************************************************
c* ***NOTE: we mark prolongation operator offsets lagging a level ***
c* ***iz(11,lev)   = iz(11,lev-1) + 27*nsmall
c* *********************************************************************
c*
c*    *** mark the beginning of everything at (nlev-1) more ***
      do 10 lev = 2, nlev
         nxold = nxnew
         nyold = nynew
         nzold = nznew
         call mkcors(1,nxold,nyold,nzold,nxnew,nynew,nznew)
         n = nxnew*nynew*nznew
c*
c*       *** mark the beginning of everything at level (lev+1) ***
         iz(1,lev+1)  = iz(1,lev)  + n
         iz(2,lev+1)  = iz(2,lev)  + 4*nynew*nznew
         iz(3,lev+1)  = iz(3,lev)  + 4*nxnew*nznew
         iz(4,lev+1)  = iz(4,lev)  + 4*nxnew*nynew
         iz(5,lev+1)  = iz(5,lev)  + 100
         iz(6,lev+1)  = iz(6,lev)  + 100
         iz(7,lev+1)  = iz(7,lev)  + 4*n
         iz(8,lev+1)  = iz(8,lev)  + nxnew
         iz(9,lev+1)  = iz(9,lev)  + nynew
         iz(10,lev+1) = iz(10,lev) + nznew
c*       *** mark prolongation operator storage for previous level ***
         iz(11,lev)   = iz(11,lev-1) + 27*n
c*       ***************************************************************
c*       ***NOTE: we mark operator offsets as we build the operators ***
c*       *** iz(7,lev+1)  = iz(7,lev)  + 4*n
c*       ***************************************************************
 10   continue
c*
c*    *** end it ***
      return
      end
      subroutine buildops (nx,ny,nz,nlev,ipkey,iinfo,ido,iz,
     2   mgprol,mgcoar,mgsolv,mgdisc,
     3   ipc,rpc,pc,ac,cc,fc,
     4   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:  
c*
c*    build operators, boundary arrays, modify affine vectors.
c*    if (ido=0) do only fine level
c*    if (ido=1) do only coarse levels (including a second op at coarsest)
c*    if (ido=2) do all levels
c*    if (ido=3) rebuild the second operator at the coarsest level
c* 
c*    note:  the fine level must be build before any coarse levels.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),iz(50,*),nx,ny,nz,nlev,iinfo,lev,ido
      integer          nxx,nyy,nzz,nxold,nyold,nzold,numdia,key,ipkey
      integer          mgprol,mgcoar,mgsolv,mgdisc
      double precision rpc(*),pc(*),ac(*),cc(*),fc(*)
      double precision a1cf(*),a2cf(*),a3cf(*),ccf(*),fcf(*),tcf(*)
      double precision gxcf(*),gycf(*),gzcf(*)
      double precision xf(*),yf(*),zf(*)
c*
c*    *** setup ***
      nxx    = nx
      nyy    = ny
      nzz    = nz
c*
c*    *** build the operator a on the finest level ***
      if ((ido .eq. 0) .or. (ido .eq. 2)) then
         lev = 1
c*
c*       *** some i/o ***
!         if (iinfo .ne. 0) then
!            write(6,100)'% BUILDOPS: (FINE) :',nxx,nyy,nzz
! 100        format(a,(2x,' [',i3,',',i3,',',i3,'] '))
!         endif
c*
c*       *** finest level discretization ***
         call buildA (nxx,nyy,nzz,ipkey,mgdisc,numdia,
     2      ipc(iz(5,lev)),rpc(iz(6,lev)),
     3      ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4      xf(iz(8,lev)),yf(iz(9,lev)),zf(iz(10,lev)),
     5      gxcf(iz(2,lev)),gycf(iz(3,lev)),gzcf(iz(4,lev)),
     6      a1cf(iz(1,lev)),a2cf(iz(1,lev)),a3cf(iz(1,lev)),
     7      ccf(iz(1,lev)),fcf(iz(1,lev)))
c*
c*       *** now initialize the differential operator offset ***
!#ifdef _MULTIGRID_VERBOSE
!         print*,'% BUILDOPS: operator stencil (lev,numdia) = ',
!     2      lev,numdia
!#endif
         iz(7,lev+1)  = iz(7,lev) + numdia * nxx * nyy * nzz
c*
c*       *** debug ***
!         if (iinfo.gt.7) call prtmatd(nxx,nyy,nzz,
!     2       ipc(iz(5,lev)),rpc(iz(6,lev)),ac(iz(7,lev)))
      endif
c*
c*    *** build the (nlev-1) level operators ***
      if ((ido .eq. 1) .or. (ido .eq. 2) .or. (ido .eq. 3)) then
         do 10 lev = 2, nlev
            nxold = nxx
            nyold = nyy
            nzold = nzz
            call mkcors(1,nxold,nyold,nzold,nxx,nyy,nzz)
            if (ido .ne. 3) then
c*
c*             *** build the interpolation operator on this level ***
               call buildP (nxold,nyold,nzold,nxx,nyy,nzz,mgprol,
     2            ipc(iz(5,lev-1)),rpc(iz(6,lev-1)),
     3            pc(iz(11,lev-1)),ac(iz(7,lev-1)),
     4            xf(iz(8,lev-1)),yf(iz(9,lev-1)),zf(iz(10,lev-1)))
c*
c*             *** differential operator this level with standard disc. ***
               if (mgcoar .eq. 0) then
!                  if (iinfo .ne. 0) then
!                     write(6,100)'% BUILDOPS: (STAND) ',nxx,nyy,nzz
!                  endif
                  call buildcopy0 (nxx,nyy,nzz,nxold,nyold,nzold,
     2               xf(iz(8,lev)),yf(iz(9,lev)),zf(iz(10,lev)),
     3               gxcf(iz(2,lev)),gycf(iz(3,lev)),gzcf(iz(4,lev)),
     4               a1cf(iz(1,lev)),a2cf(iz(1,lev)),a3cf(iz(1,lev)),
     5               ccf(iz(1,lev)),fcf(iz(1,lev)),tcf(iz(1,lev)),
     6               xf(iz(8,lev-1)),yf(iz(9,lev-1)),zf(iz(10,lev-1)),
     7               gxcf(iz(2,lev-1)),gycf(iz(3,lev-1)),
     7               gzcf(iz(4,lev-1)),
     8               a1cf(iz(1,lev-1)),a2cf(iz(1,lev-1)),
     8               a3cf(iz(1,lev-1)),
     9               ccf(iz(1,lev-1)),fcf(iz(1,lev-1)),
     9               tcf(iz(1,lev-1)))
                  call buildA (nxx,nyy,nzz,ipkey,mgdisc,numdia,
     2               ipc(iz(5,lev)),rpc(iz(6,lev)),
     3               ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4               xf(iz(8,lev)),yf(iz(9,lev)),zf(iz(10,lev)),
     5               gxcf(iz(2,lev)),gycf(iz(3,lev)),gzcf(iz(4,lev)),
     6               a1cf(iz(1,lev)),a2cf(iz(1,lev)),a3cf(iz(1,lev)),
     7               ccf(iz(1,lev)),fcf(iz(1,lev)))
c*
c*             *** differential operator this level with harmonic disc. ***
               elseif (mgcoar .eq. 1) then
!                  if (iinfo .ne. 0) then
!                     write(6,100)'% BUILDOPS: (HARMO) ',nxx,nyy,nzz
!                  endif
                  call buildharm0 (nxx,nyy,nzz,nxold,nyold,nzold,
     2               xf(iz(8,lev)),yf(iz(9,lev)),zf(iz(10,lev)),
     3               gxcf(iz(2,lev)),gycf(iz(3,lev)),gzcf(iz(4,lev)),
     4               a1cf(iz(1,lev)),a2cf(iz(1,lev)),a3cf(iz(1,lev)),
     5               ccf(iz(1,lev)),fcf(iz(1,lev)),tcf(iz(1,lev)),
     6               xf(iz(8,lev-1)),yf(iz(9,lev-1)),zf(iz(10,lev-1)),
     7               gxcf(iz(2,lev-1)),gycf(iz(3,lev-1)),
     7               gzcf(iz(4,lev-1)),
     8               a1cf(iz(1,lev-1)),a2cf(iz(1,lev-1)),
     8               a3cf(iz(1,lev-1)),
     9               ccf(iz(1,lev-1)),fcf(iz(1,lev-1)),
     9               tcf(iz(1,lev-1)))
                  call buildA (nxx,nyy,nzz,ipkey,mgdisc,numdia,
     2               ipc(iz(5,lev)),rpc(iz(6,lev)),
     3               ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4               xf(iz(8,lev)),yf(iz(9,lev)),zf(iz(10,lev)),
     5               gxcf(iz(2,lev)),gycf(iz(3,lev)),gzcf(iz(4,lev)),
     6               a1cf(iz(1,lev)),a2cf(iz(1,lev)),a3cf(iz(1,lev)),
     7               ccf(iz(1,lev)),fcf(iz(1,lev)))
c*
c*             *** differential operator with galerkin formulation ***
               elseif (mgcoar .eq. 2) then
                  if (iinfo .ne. 0) then
!                     write(6,100)'% BUILDOPS: (GALER) ',nxx,nyy,nzz
                  endif
                  call buildgaler0 (nxold,nyold,nzold,nxx,nyy,nzz,
     2               ipkey,numdia,pc(iz(11,lev-1)),
     3               ipc(iz(5,lev-1)),rpc(iz(6,lev-1)),
     4               ac(iz(7,lev-1)),cc(iz(1,lev-1)),fc(iz(1,lev-1)),
     5               ipc(iz(5,lev)),rpc(iz(6,lev)),
     6               ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)))
                  call extrac(nxold,nyold,nzold,nxx,nyy,nzz,
     2               tcf(iz(1,lev-1)),tcf(iz(1,lev)))
               else
                  print*,'% BUILDOPS: bad mgcoar key given...'
               endif
c*
c*             *** now initialize the differential operator offset ***
!               print*,'% BUILDOPS: operator stencil (lev,numdia) = ',
!     2            lev,numdia
               iz(7,lev+1)  = iz(7,lev) + numdia * nxx * nyy * nzz
c*
c*             *** debug ***
               if (iinfo.gt.8) call prtmatd(nxx,nyy,nzz,
     2             ipc(iz(5,lev)),rpc(iz(6,lev)),ac(iz(7,lev)))
            endif
 10      continue
c*
c*       *** build a sparse format coarse grid operator ***
         if (mgsolv .eq. 1) then
            lev = nlev
            call buildband (key,nxx,nyy,nzz,
     2         ipc(iz(5,lev)),rpc(iz(6,lev)),ac(iz(7,lev)),
     3         ipc(iz(5,lev+1)),rpc(iz(6,lev+1)),ac(iz(7,lev+1)))
            if (key .eq. 1) then
!               print*,'% BUILDOPS: changing your MGSOLV to iterative...'
               mgsolv = 0
            endif
         endif
      endif
c*
c*    *** end it ***
      return
      end
      subroutine buildcopy0 (nx,ny,nz,nxf,nyf,nzf,
     2   xc,yc,zc,gxc,gyc,gzc,a1c,a2c,a3c,cc,fc,tc,
     3   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:   
c*
c*    produce information for a coarser grid.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nx,ny,nz,i,j,k,ii,jj,kk
      integer          iadd,jadd,kadd
      double precision xc(nx),yc(ny),zc(nz)
      double precision gxc(ny,nz,*),gyc(nx,nz,*),gzc(nx,ny,*)
      double precision a1c(nx,ny,nz),a2c(nx,ny,nz),a3c(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz),tc(nx,ny,nz)
      double precision xf(nxf),yf(nyf),zf(nzf)
      double precision gxcf(nyf,nzf,*),gycf(nxf,nzf,*),gzcf(nxf,nyf,*)
      double precision a1cf(nxf,nyf,nzf),a2cf(nxf,nyf,nzf)
      double precision a3cf(nxf,nyf,nzf),tcf(nxf,nyf,nzf)
      double precision ccf(nxf,nyf,nzf),fcf(nxf,nyf,nzf)
c*
cmdir 0 0
c*
c*    *** how far to step into the coefficient arrays ***
      iadd   = (nxf-1)/(nx-1)
      jadd   = (nyf-1)/(ny-1)
      kadd   = (nzf-1)/(nz-1)
      if ((iadd.ne.2).or.(jadd.ne.2).or.(kadd.ne.2)) then
         print*,'% BUILDCOPY0: problem with grid dimensions...'
      endif
c*
c*    *** compute the coarse grid pde coefficients ***
cmdir 3 1
      do 30 k = 1, nz
         kk = 2 * k - 1
         zc(k) = zf(kk)
cmdir 3 2
         do 31 j = 1, ny
            jj = 2 * j - 1
            yc(j) = yf(jj)
cmdir 3 3
            do 32 i = 1, nx
               ii = 2 * i - 1
               xc(i) = xf(ii)
c*
c*             *** true solution ***
               tc(i,j,k) = tcf(ii,jj,kk)
c*
c*             *** helmholtz coefficient ***
               cc(i,j,k) = ccf(ii,jj,kk)
c*
c*             *** source function ***
               fc(i,j,k) = fcf(ii,jj,kk)
c*
c*             *** east/west neighbor ***
               a1c(i,j,k) = a1cf(ii,jj,kk)
c*
c*             *** north/south neighbor ***
               a2c(i,j,k) = a2cf(ii,jj,kk)
c*
c*             *** up/down neighbor ***
               a3c(i,j,k) = a3cf(ii,jj,kk)
 32         continue
 31      continue
 30   continue
c*
c*    *** the (i=1) and (i=nx) boundaries ***
cmdir 2 1
      do 50 k = 1, nz
         kk = 2 * k - 1
cmdir 2 2
         do 51 j = 1, ny
            jj = 2 * j - 1
            gxc(j,k,1) = gxcf(jj,kk,1)
            gxc(j,k,2) = gxcf(jj,kk,2)
            gxc(j,k,3) = gxcf(jj,kk,3)
            gxc(j,k,4) = gxcf(jj,kk,4)
 51      continue
 50   continue
c*
c*    *** the (j=1) and (j=ny) boundaries ***
cmdir 2 1
      do 60 k = 1, nz
         kk = 2 * k - 1
cmdir 2 2
         do 61 i = 1, nx
            ii = 2 * i - 1
            gyc(i,k,1) = gycf(ii,kk,1)
            gyc(i,k,2) = gycf(ii,kk,2)
            gyc(i,k,3) = gycf(ii,kk,3)
            gyc(i,k,4) = gycf(ii,kk,4)
 61      continue
 60   continue
c*
c*    *** the (k=1) and (k=nz) boundaries ***
cmdir 2 1
      do 70 j = 1, ny
         jj = 2 * j - 1
cmdir 2 2
         do 71 i = 1, nx
            ii = 2 * i - 1
            gzc(i,j,1) = gzcf(ii,jj,1)
            gzc(i,j,2) = gzcf(ii,jj,2)
            gzc(i,j,3) = gzcf(ii,jj,3)
            gzc(i,j,4) = gzcf(ii,jj,4)
 71      continue
 70   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildharm0 (nx,ny,nz,nxf,nyf,nzf,
     2   xc,yc,zc,gxc,gyc,gzc,a1c,a2c,a3c,cc,fc,tc,
     3   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c* *********************************************************************
c* purpose:   
c*
c*    produce information for a coarser grid.
c*    but also harmonically average the problem coefficients.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          nxf,nyf,nzf,nx,ny,nz,i,j,k,ii,jj,kk
      integer          iadd,jadd,kadd
      double precision xc(nx),yc(ny),zc(nz)
      double precision gxc(ny,nz,*),gyc(nx,nz,*),gzc(nx,ny,*)
      double precision a1c(nx,ny,nz),a2c(nx,ny,nz),a3c(nx,ny,nz)
      double precision cc(nx,ny,nz),fc(nx,ny,nz),tc(nx,ny,nz)
      double precision xf(nxf),yf(nyf),zf(nzf)
      double precision gxcf(nyf,nzf,*),gycf(nxf,nzf,*),gzcf(nxf,nyf,*)
      double precision a1cf(nxf,nyf,nzf),a2cf(nxf,nyf,nzf)
      double precision a3cf(nxf,nyf,nzf),tcf(nxf,nyf,nzf)
      double precision ccf(nxf,nyf,nzf),fcf(nxf,nyf,nzf)
c*
c*    *** statement functions ***
      double precision harmo2 
c     ,harmo4,arith2,arith4,arith6,arith8
      double precision a,b,c,d,e,f,g,h
      harmo2(a,b)      = 2. * a * b / (a + b)
c      harmo4(a,b,c,d)  = 1. / ( 0.25 * ( 1./a + 1./b + 1./c + 1./d ) )
c      arith2(a,b)      = (a+b) / 2.
c      arith4(a,b,c,d)  = (a+b+c+d) / 4.
c      arith6(a,b,c,d,e,f) = (a+b+c+d+e+f) / 6.
c      arith8(a,b,c,d,e,f,g,h) = (a+b+c+d+e+f+g+h) / 8.
c*
cmdir 0 0
c*
c*    *** how far to step into the coefficient arrays ***
      iadd   = (nxf-1)/(nx-1)
      jadd   = (nyf-1)/(ny-1)
      kadd   = (nzf-1)/(nz-1)
      if ((iadd.ne.2).or.(jadd.ne.2).or.(kadd.ne.2)) then
         print*,'% BUILDHARM0: problem with grid dimensions...'
      endif
c*
c*    *** compute the coarse grid pde coefficients ***
cmdir 3 1
      do 30 k = 1, nz
         kk = 2 * k - 1
         zc(k) = zf(kk)
cmdir 3 2
         do 31 j = 1, ny
            jj = 2 * j - 1
            yc(j) = yf(jj)
cmdir 3 3
            do 32 i = 1, nx
               ii = 2 * i - 1
               xc(i) = xf(ii)
c*
c*             *** true solution ***
               tc(i,j,k) = tcf(ii,jj,kk)
c*
c*             *** helmholtz coefficient ***
               cc(i,j,k) = ccf(ii,jj,kk)
CZZZ           cc(i,j,k) = (
CZZZ 2            +0.5e0 * ccf(ii,jj,kk)
CZZZ 3            +0.5e0 * arith6( ccf(max0(1,ii-1),jj,kk),
CZZZ 4                             ccf(min0(nxf,ii+1),jj,kk),
CZZZ 5                             ccf(ii,max0(1,jj-1),kk),
CZZZ 6                             ccf(ii,min0(nyf,jj+1),kk),
CZZZ 7                             ccf(ii,jj,max0(nzf,kk-1)),
CZZZ 8                             ccf(ii,jj,min0(nzf,kk+1)) )
CZZZ 9            )
c*
c*             *** source function ***
               fc(i,j,k) = fcf(ii,jj,kk)
CZZZ           fc(i,j,k) = (
CZZZ 2            +0.5e0 * fcf(ii,jj,kk)
CZZZ 3            +0.5e0 * arith6( fcf(max0(1,ii-1),jj,kk),
CZZZ 4                             fcf(min0(nxf,ii+1),jj,kk),
CZZZ 5                             fcf(ii,max0(1,jj-1),kk),
CZZZ 6                             fcf(ii,min0(nyf,jj+1),kk),
CZZZ 7                             fcf(ii,jj,max0(nzf,kk-1)),
CZZZ 8                             fcf(ii,jj,min0(nzf,kk+1)) )
CZZZ 9            )
c*
c*             *** east/west neighbor ***
               a1c(i,j,k) = (
     2           +0.5e0  *harmo2( a1cf(ii,jj,kk),
     2                 a1cf(min0(nxf,ii+1),jj,kk) )
     3           +0.125e0*harmo2( a1cf(ii,jj,max0(1,kk-1)),
     3                 a1cf(min0(nxf,ii+1),jj,max0(1,kk-1)) )
     4           +0.125e0*harmo2( a1cf(ii,jj,min0(nzf,kk+1)),
     4                 a1cf(min0(nxf,ii+1),jj,min0(nzf,kk+1)) )
     5           +0.125e0*harmo2( a1cf(ii,max0(1,jj-1),kk),
     5                 a1cf(min0(nxf,ii+1),max0(1,jj-1),kk) )
     6           +0.125e0*harmo2( a1cf(ii,min0(nyf,jj+1),kk),
     6                 a1cf(min0(nxf,ii+1),min0(nyf,jj+1),kk) )
     7            )
c*
c*             *** north/south neighbor ***
               a2c(i,j,k) = (
     2           +0.5e0  *harmo2( a2cf(ii,jj,kk),
     2                 a2cf(ii,min0(nyf,jj+1),kk) )
     3           +0.125e0*harmo2( a2cf(ii,jj,max0(1,kk-1)),
     3                 a2cf(ii,min0(nyf,jj+1),max0(1,kk-1)) )
     4           +0.125e0*harmo2( a2cf(ii,jj,min0(nzf,kk+1)),
     4                 a2cf(ii,min0(nyf,jj+1),min0(nzf,kk+1)) )
     5           +0.125e0*harmo2( a2cf(max0(1,ii-1),jj,kk),
     5                 a2cf(max0(1,ii-1),min0(nyf,jj+1),kk) )
     6           +0.125e0*harmo2( a2cf(min0(nxf,ii+1),jj,kk),
     6                 a2cf(min0(nxf,ii+1),min0(nyf,jj+1),kk) )
     7            )
c*
c*             *** up/down neighbor ***
               a3c(i,j,k) = (
     2           +0.5e0  *harmo2( a3cf(ii,jj,kk),
     2                 a3cf(ii,jj,min0(nzf,kk+1)) )
     3           +0.125e0*harmo2( a3cf(ii,max0(1,jj-1),kk),
     3                 a3cf(ii,max0(1,jj-1),min0(nzf,kk+1)) )
     4           +0.125e0*harmo2( a3cf(ii,min0(nyf,jj+1),kk),
     4                 a3cf(ii,min0(nyf,jj+1),min0(nzf,kk+1)) )
     5           +0.125e0*harmo2( a3cf(max0(1,ii-1),jj,kk),
     5                 a3cf(max0(1,ii-1),jj,min0(nzf,kk+1)) )
     6           +0.125e0*harmo2( a3cf(min0(nxf,ii+1),jj,kk),
     6                 a3cf(min0(nxf,ii+1),jj,min0(nzf,kk+1)) )
     7            )
 32         continue
 31      continue
 30   continue
c*
c*    *** the (i=1) and (i=nx) boundaries ***
cmdir 2 1
      do 50 k = 1, nz
         kk = 2 * k - 1
cmdir 2 2
         do 51 j = 1, ny
            jj = 2 * j - 1
            gxc(j,k,1) = gxcf(jj,kk,1)
            gxc(j,k,2) = gxcf(jj,kk,2)
            gxc(j,k,3) = gxcf(jj,kk,3)
            gxc(j,k,4) = gxcf(jj,kk,4)
 51      continue
 50   continue
c*
c*    *** the (j=1) and (j=ny) boundaries ***
cmdir 2 1
      do 60 k = 1, nz
         kk = 2 * k - 1
cmdir 2 2
         do 61 i = 1, nx
            ii = 2 * i - 1
            gyc(i,k,1) = gycf(ii,kk,1)
            gyc(i,k,2) = gycf(ii,kk,2)
            gyc(i,k,3) = gycf(ii,kk,3)
            gyc(i,k,4) = gycf(ii,kk,4)
 61      continue
 60   continue
c*
c*    *** the (k=1) and (k=nz) boundaries ***
cmdir 2 1
      do 70 j = 1, ny
         jj = 2 * j - 1
cmdir 2 2
         do 71 i = 1, nx
            ii = 2 * i - 1
            gzc(i,j,1) = gzcf(ii,jj,1)
            gzc(i,j,2) = gzcf(ii,jj,2)
            gzc(i,j,3) = gzcf(ii,jj,3)
            gzc(i,j,4) = gzcf(ii,jj,4)
 71      continue
 70   continue
c*
c*    *** return and end ***
      return
      end
      subroutine buildgaler0 (nxf,nyf,nzf,nxc,nyc,nzc,ipkey,numdia,
     2   pcFF,
     3   ipcFF,rpcFF,acFF,ccFF,fcFF,
     4   ipc,rpc,ac,cc,fc)
c* *********************************************************************
c* purpose: 
c*
c*    form the galerkin coarse grid system.
c*
c*    notes: although the fine grid matrix may be 7 or 27 diagonal,
c*           the coarse grid matrix is always 27 diagonal.
c*           (only 14 stored due to symmetry.)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipcFF(*),ipc(*),nxf,nyf,nzf,nxc,nyc,nzc,ipkey
      integer          numdia,numdia_loc
      double precision pcFF(*),rpcFF(*),acFF(*),ccFF(*),fcFF(*)
      double precision rpc(*),ac(*),cc(*),fc(*)
c*
cmdir 0 0
c*
c*    *** call the algebraic galerkin routine ***
      write(*,*) "pb: in buildgaler0"
      numdia_loc = ipcFF(11)
      call buildG (nxf,nyf,nzf,nxc,nyc,nzc,numdia_loc,pcFF,acFF,ac)
c*
c*    *** note how many nonzeros in this new discretization stencil ***
      ipc(11) = 27
      numdia  = 14
c*
c*    *** save the problem key with this new operator ***
      ipc(10) = ipkey
c*
c*    *** restrict the helmholtz term and source function ***
      call restrc(nxf,nyf,nzf,nxc,nyc,nzc,ccFF,cc,pcFF)
      call restrc(nxf,nyf,nzf,nxc,nyc,nzc,fcFF,fc,pcFF)
c*
c*    *** return and end ***
      return
      end
      subroutine buildgaler1 (nxf,nyf,nzf,nxc,nyc,nzc,numdia,
     2   pcFF,
     3   ipcFF,rpcFF,ccFF,
     4   ipc,rpc,cc)
c* *********************************************************************
c* purpose: 
c*
c*    form the helmholtz term of a galerkin coarse grid system.
c*
c*    notes: although the fine grid matrix may be 1 or 27 diagonal,
c*           the coarse grid matrix is always 27 diagonal.
c*           (only 14 stored due to symmetry.)
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipcFF(*),ipc(*),nxf,nyf,nzf,nxc,nyc,nzc,numdia
      integer          numdia_loc
      double precision pcFF(*),rpcFF(*),ccFF(*),rpc(*),cc(*)
c*
cmdir 0 0
c*
c*    *** call the algebraic galerkin routine ***
      numdia_loc = ipcFF(12)
      call buildG (nxf,nyf,nzf,nxc,nyc,nzc,numdia_loc,pcFF,ccFF,cc)
c*
c*    *** note how many nonzeros in this new discretization stencil ***
      ipc(12) = 27
      numdia  = 14
c*
c*    *** return and end ***
      return
      end
      subroutine buildALG (nx,ny,nz,mode,nlev,iz,
     3   ipc,rpc,ac,cc,fc,x,y,tmp)
c* *********************************************************************
c* purpose:  
c*
c*    build rhs algebraically for analysis purposes.
c* 
c*    note:  the fine level must be build before any coarse levels.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          ipc(*),iz(50,*),nx,ny,nz,mode,nlev,lev
      integer          nxx,nyy,nzz,nxold,nyold,nzold
      double precision rpc(*),ac(*),cc(*),fc(*),x(*),y(*),tmp(*)
c*
c*    *** setup ***
      nxx    = nx
      nyy    = ny
      nzz    = nz
c*
c*    *** build the rhs the finest level ***
      lev = 1
      if ((mode .eq. 1) .or. (mode .eq. 2)) then
         call nmatvec(nxx,nyy,nzz,
     2      ipc(iz(5,lev)),rpc(iz(6,lev)),
     3      ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4      x(iz(1,lev)),y(iz(1,lev)),tmp)
      else
         call matvec(nxx,nyy,nzz,
     2      ipc(iz(5,lev)),rpc(iz(6,lev)),
     3      ac(iz(7,lev)),cc(iz(1,lev)),
     4      x(iz(1,lev)),y(iz(1,lev)))
      endif
c*
c*    *** build the (nlev-1) level rhs function ***
      do 10 lev = 2, nlev
         nxold = nxx
         nyold = nyy
         nzold = nzz
         call mkcors(1,nxold,nyold,nzold,nxx,nyy,nzz)
c*
c*       *** build the rhs on this level ***
         if ((mode .eq. 1) .or. (mode .eq. 2)) then
            call nmatvec(nxx,nyy,nzz,
     2         ipc(iz(5,lev)),rpc(iz(6,lev)),
     3         ac(iz(7,lev)),cc(iz(1,lev)),fc(iz(1,lev)),
     4         x(iz(1,lev)),y(iz(1,lev)),tmp)
         else
            call matvec(nxx,nyy,nzz,
     2         ipc(iz(5,lev)),rpc(iz(6,lev)),
     3         ac(iz(7,lev)),cc(iz(1,lev)),
     4         x(iz(1,lev)),y(iz(1,lev)))
         endif
 10   continue
c*
c*    *** end it ***
      return
      end
