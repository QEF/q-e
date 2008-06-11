c* ///////////////////////////////////////////////////////////////////////////
c* @file    mainsubd.f
c* @author  Michael Holst
c* @brief   Supporting I/O and other routines for main.c and maind.f.
c* @version $Id: mainsubd.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

      subroutine readit (iparm,rparm,nx,ny,nz,nlev,nrwk,niwk,key,meth,
     2                   errtol,itmax,whichbc)
c* *********************************************************************
c* purpose:
c*
c*    this routine reads in some initial values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iparm(*),nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey
      integer          istop,iinfo,key,meth,itmax,ipcon
      integer          ipkey,iperf, nlevmax
      integer          nonlin,mgprol,mgcoar,mgdisc,mgsmoo,mgsolv, maxlev
      double precision rparm(*),errtol,omegal,omegan
      integer          whichbc(3)
c*
c*    *** parameters ***
      integer          iread,irite
      parameter        (iread=7, irite=8)
c*
c*    *** if not interactive mode then open i/o files ***
!      open(unit=iread,  file='in',   status='unknown')
!      open(unit=irite,  file='outt', status='unknown')
!      rewind(iread)
!      rewind(irite)
c*
c*    *** input the controling parameters ***
!      read (iread,10)
!      read (iread,10)
!      read (iread,20) errtol
!      read (iread,30) itmax
!      read (iread,30) istop
!      read (iread,30) iinfo
!      read (iread,30) ipkey
!      read (iread,30) key
!      read (iread,30) iperf
!      read (iread,10)
!      read (iread,10)
!      read (iread,10)
!      read (iread,30) meth
!      read (iread,30) nonlin
!      read (iread,30) mgkey
!      read (iread,30) nlev
!      read (iread,30) nu1
!      read (iread,30) nu2
!      read (iread,30) mgsmoo
!      read (iread,30) mgprol
!      read (iread,30) mgcoar
!      read (iread,30) mgsolv
!      read (iread,30) mgdisc
!      read (iread,20) omegal
!      read (iread,20) omegan
!      read (iread,30) ipcon
!      read (iread,30) whichbc(1)
!      read (iread,30) whichbc(2)
!      read (iread,30) whichbc(3)

!   1.0e-22 errtol    error tolerance desired
!            errtol = 1.0e-22
!            print *,"errtol vale",errtol
!      1000 itmax     maximum number of iterations
!            itmax = 1000
!            print *,"itmax vale ",itmax
!         0 istop     stop crit (0=res, 1=relr, 2=dif, 3=errc, 4=errd, 5=aerrd)
            istop = 0
!         2 iinfo     runtime status messages (0=no, 1=some, 2=lots, 3=more)
            iinfo = 2
!         0 ipkey     problem setup option (specific to problem)
            ipkey = 0
!         0 key       print solution to ascii file (0=no, 1=yes)
            key = 0
!         0 iperf     analysis (0=no, 1=kappa(A), 2=rho(M), 3=1+2)
            iperf = 0
!%%%%%%%%%%%%%% method options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!         0 meth      (0=cgmg, 1=new, 2=mg, 3=cg, 4=sor, 5=rbgs, 6=wj, 7=rich)
            meth = 0
!         0 nonlin    (0=linear, 1=nonlinear, 2=linear-then-nonlinear)
            nonlin = 0
!         0 mgkey     method (0=variable v-cycle, 1=nested iteration)
            mgkey = 0
!         5 nlev      number of mesh levels to use
!            nlev = 5
!            print *,"nlev vale :",nlev
!        75 nu1       num of pre-smoothings
            nu1 = 75
!        75 nu2       num of post-smoothings
            nu2 = 75
!         0 mgsmoo    smoothing (0=wjac, 1=gsrb, 2=sor, 3=rich, 4=cghs)
            mgsmoo = 0
!
! mgsmoo 4 also should be available is it is called from the program
! independently of the value set here ...
!
!         0 mgprol    prolongation (0=trilinear, 1=op-based, 2=mod-op-based)
            mgprol = 0
!         0 mgcoar    coarse equations (0=stan, 1=harm, 2=galerkin)
            mgcoar = 0
!         0 mgsolv    coarse equation solver (0=cghs, 1=banded linpack)
            mgsolv = 0
!         0 mgdisc    discretization (0=fv, 1=fe)
            mgdisc = 0
!        .8 omegal    linear relaxation parameter
            omegal = .8
!  9.00e-01 omegan    nonlinear relaxation parameter
            omegan = 9.00e-01
!         3 ipcon     prec flag (0=diag, 1=iccg, 2=iccgdw, 3=miccgdw, 4=none)
            ipcon = 3
!         0 whichbcx  x boundary condition  (0=zero, 1=periodic)
!            print *,"whichbc vale",whichbc
!            whichbc(1) = 0
!         0 whichbcy  y boundary condition  (0=zero, 1=periodic)
!            whichbc(2) = 0
!         0 whichbcz  z boundary condition  (0=zero, 1=periodic)
!            whichbc(3) = 0

      nlevmax = maxlev(nx,ny,nz)
      IF( nlev .gt. nlevmax ) THEN
         nlev = nlevmax
         print*,"From multigrid. Warning: nlev = ", nlevmax
      END IF

c*
c*    *** pack iparm/rparm correctly for desired method ***
      call packmg (iparm,rparm,
     2   nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey,itmax,istop,ipcon,
     3   nonlin,mgsmoo,mgprol,mgcoar,mgsolv,mgdisc,iinfo,errtol,ipkey,
     4   omegal,omegan,irite,iperf)
c*
c*    *** do a little output now ***
!      write(*,40) '% READIT: done reading input file... '
c*
c*    *** format statements ***
! 10   format()
! 20   format(e10.1)
! 30   format(i10)
! 40   format (a)
! 41   format ('Warning: nlev = ',i10,';')
c*
c*    *** return and end ***
      return
      end
      subroutine packmg (iparm,rparm,
     2   nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey,itmax,istop,ipcon,
     3   nonlin,mgsmoo,mgprol,mgcoar,mgsolv,mgdisc,iinfo,errtol,ipkey,
     4   omegal,omegan,irite,iperf)
c* *********************************************************************
c* purpose:
c*
c*    this routine reads in some initial values.
c*
c* author:  michael holst
c* *********************************************************************
      implicit         none
      integer          iparm(*),nrwk,niwk,nx,ny,nz,nlev,nu1,nu2,mgkey
      integer          itmax,istop,ipcon,nonlin,iinfo,irite,ipkey,iperf
      integer          mgsmoo,mgprol,mgcoar,mgsolv,mgdisc
      double precision rparm(*),errtol,omegal,omegan
c*
c*    *** encode iparm parameters ***
      iparm(1)  = nrwk
      iparm(2)  = niwk
      iparm(3)  = nx
      iparm(4)  = ny
      iparm(5)  = nz
      iparm(6)  = nlev
      iparm(7)  = nu1
      iparm(8)  = nu2
      iparm(9)  = mgkey
      iparm(10) = itmax
      iparm(11) = istop
      iparm(12) = iinfo
      iparm(13) = irite
      iparm(14) = ipkey
      iparm(15) = ipcon
      iparm(16) = nonlin
      iparm(17) = mgprol
      iparm(18) = mgcoar
      iparm(19) = mgdisc
      iparm(20) = mgsmoo
      iparm(21) = mgsolv
      iparm(22) = iperf
c*
c*    *** encode rparm parameters ***
      rparm(1)  = errtol
      rparm(9)  = omegal
      rparm(10) = omegan
c*
c*    *** return and end ***
     
      return
      end
