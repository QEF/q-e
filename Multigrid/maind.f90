!* ///////////////////////////////////////////////////////////////////////////
!* @file    maind.f
!* @author  Michael Holst
!* @brief   Main FORTRAN driver for PMG.
!* @version $Id: maind.f90,v 1.2 2008-09-23 09:19:31 giannozz Exp $
!* @attention
!* @verbatim
!*
!* PMG -- Parallel algebraic MultiGrid
!* Copyright (c) 1994-2006.  Michael Holst.
!*
!* Michael Holst <mholst@math.ucsd.edu>
!* University of California, San Diego
!* Department of Mathematics, 5739 AP&M
!* 9500 Gilman Drive, Dept. 0112
!* La Jolla, CA 92093-0112 USA                                                  
!* http://math.ucsd.edu/~mholst
!*
!* This file is part of PMG.
!*
!* PMG is free software; you can redistribute it and/or modify
!* it under the terms of the GNU General Public License as published by
!* the Free Software Foundation; either version 2 of the License, or
!* (at your option) any later version.
!*
!* PMG is distributed in the hope that it will be useful,
!* but WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!* GNU General Public License for more details.
!*
!* You should have received a copy of the GNU General Public License
!* along with PMG; if not, write to the Free Software
!* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
!*
!* Linking PMG statically or dynamically with other modules is making a
!* combined work based on PMG. Thus, the terms and conditions of the GNU
!* General Public License cover the whole combination.
!* 
!* SPECIAL GPL EXCEPTION
!* In addition, as a special exception, the copyright holders of PMG
!* give you permission to combine the PMG program with free software
!* programs and libraries that are released under the GNU LGPL or with
!* code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision.
!* Such combined software may be linked with PMG and redistributed together 
!* in original or modified form as mere aggregation without requirement that 
!* the entire work be under the scope of the GNU General Public License.
!* This special exception permission is also extended to any software listed
!* in the SPECIAL GPL EXCEPTION clauses by the FEtk and APBS libraries.
!* 
!* Note that people who make modified versions of PMG are not obligated
!* to grant this special exception for their modified versions; it is
!* their choice whether to do so. The GNU General Public License gives
!* permission to release a modified version without this exception; this
!* exception also makes it possible to release a modified version which
!* carries forward this exception.
!*
!* @endverbatim
!* ///////////////////////////////////////////////////////////////////////////

!* ///////////////////////////////////////////////////////////////////////////
!*                           Welcome to PMG (MG/XMG).
!* ///////////////////////////////////////////////////////////////////////////
!*
!* ======================
!* The Author Information
!* ======================
!*
!* The computer codes PMG (MG/XMG, "the Code"), including FORTRAN/C/C++
!* language multilevel subroutines called "MG" and the C-language X-Window
!* interface called "XMG", were developed by:
!*
!*    Michael Holst                TELE:  (858) 534-4899
!*    Department of Mathematics    FAX:   (858) 534-5273
!*    UC San Diego, AP&M 5739      EMAIL: mholst@math.ucsd.edu
!*    La Jolla, CA 92093 USA       WEB:   http://www.scicomp.ucsd.edu/~mholst
!*
!* ====================================================
!* Licensing, permission to use, modify, and distribute
!* ====================================================
!*
!* The code has been placed under the GNU General Public License (GPL).  
!*
!* This means that essentially you can use and modify the code as you like, 
!* under a few conditions, such as acknowledging the authors of the original 
!* code in any derived works, such as new programs or research results using 
!* the code.
!*
!* Please have a look at the GPL which should have accompanied the code.
!* If it is missing, please write to the Free Software Foundation, Inc., 
!* 675 Mass Ave, Cambridge, MA 02139, USA.
!*
!* It would have been easier for me to just keep the code completely private, 
!* or to hand out very restricted use copies.  However, in the interests of 
!* freedom of information, the advancement of science, the future of the 
!* planet, and so on, I decided to distribute the code in this very free way.
!* The GNU world functions very much on the honor system.  I want you to have 
!* all of my source code because I think it will do more toward helping you 
!* understand the ideas in there than anything else.  The GNU idea is that in 
!* return for all of my work on the source code and the ideas behind it, all I 
!* ask for is some acknowledgment in your derived works for using my original 
!* code and the original ideas that went into the code.  For GNU to be viable, 
!* we all need to abide by this very minimal requirement.
!*
!* I'm typing this on a nearly all-GNU system right now, namely Linux 1.2.8, 
!* a GNU version of UN*X available under the same conditions as PMG.
!* I.e., you get all of the source code for the entire UN*X operating system.
!* For more information about Linux, see: http://sunsite.unc.edu/mdw/linux.html
!*
!* ====================================================================
!* Pointers to Papers and other Reference Materials Related to the Code
!* ====================================================================
!*
!* The methods and techniques  implemented  in  the  Code  were  developed  by 
!* Michael Holst as a part of his PhD research, and are described in detail in 
!* the PhD thesis:
!*
!*    @phdthesis{Hols93a,
!*       author    = "M. Holst",
!*       title     = "Multilevel Methods for the {Poisson-Boltzmann} Equation",
!*       note      = "Also published as Tech. Rep. UIUCDCS-R-03-1821",
!*       school    = "Numerical Computing Group,
!*                    Department of Computer Science,
!*                    University of Illinois at Urbana-Champaign",
!*       year      = 1993 }
!*
!* Also refer to the supporting papers:
!*
!*    @manual{Hols94d,
!*       author    = "M. Holst",
!*       title     = "The {Poisson-Boltzmann} Equation:
!*                    Analysis and Multilevel Numerical Solution",
!*       note      = "Monograph (updated and extended form of the
!*                    Ph.D. thesis~\cite{Hols93a})." }
!*    
!*    @article{SHHN94,
!*       author    = "R. Sampogna and J. Hecht and M. Holst and A. Nicholls
!*                    and B. Honig",
!*       title     = "Nonlinear {Poisson-Boltzmann} Calculation of {pKa} 
!*                    Values",
!*       journal   = "Biophysics",
!*       note      = "(In Progress?)",
!*       year      = 1994 }
!*    
!*    @article{HoVa94b,
!*       author    = "M. Holst and S. Vandewalle",
!*       title     = "Schwarz Methods: to Symmetrize or not to Symmetrize",
!*       journal   = SINUM,
!*       note      = "(Accepted)",
!*       year      = 1994 }
!*    
!*    @article{Hols94e,
!*       author    = "M. Holst and F. Saied",
!*       title     = "Numerical Solution of the Nonlinear {Poisson-Boltzmann}
!*                    Equation: Developing More Robust and Efficient Methods",
!*       journal   = JCC,
!*       volume    = "16",
!*       number    = "3",
!*       pages     = "337--364",
!*       year      = 1995 }
!*    
!*    @article{HoSa93a,
!*       author    = "M. Holst and F. Saied",
!*       title     = "Multigrid Solution of the {Poisson-Boltzmann} Equation",
!*       journal   = JCC,
!*       volume    = "14",
!*       number    = "1",
!*       pages     = "105--113",
!*       year      = 1993 }
!*    
!*    @article{HKSS94,
!*       author    = "M. Holst and R. Kozack and F. Saied and S. Subramaniam",
!*       title     = "Protein Electrostatics: Rapid Multigrid-based {Newton}
!*                    Algorithm for Solution of the Full Nonlinear
!*                    {Poisson-Boltzmann} Equation",
!*       journal   = "J. Biomol. Struct. Dyn.",
!*       volume    = "11",
!*       pages     = "1437--1445",
!*       year      = 1994 }
!*    
!*    @article{HKSS93b,
!*       author    = "M. Holst and R. Kozack and F. Saied and S. Subramaniam",
!*       title     = "Treatment of Electrostatic Effects in Proteins:
!*                    Multigrid-based-{Newton} Iterative Method for Solution 
!*                    of the Full Nonlinear {Poisson-Boltzmann} Equation",
!*       journal   = "Proteins: Structure, Function, and Genetics",
!*       volume    = "18",
!*       number    = "3",
!*       pages     = "231--245",
!*       year      = 1994 }
!*    
!*    @techreport{Hols94h,
!*       author    = "M. Holst",
!*       title     = "A Theoretical Analysis of the {Poisson-Boltzmann} 
!*                    Equation: A priori Estimates and Well-posedness",
!*       institution = "Applied Mathematics and CRPC,
!*                      California Institute of Technology",
!*       year      = 1994 }
!*    
!*    @inproceedings{HoVa95a,
!*       author    = "M. Holst and S. Vandewalle",
!*       title     = "Schwarz Methods: To Symmetrize or not to Symmetrize",
!*       booktitle = "Proceedings of the Seventh Copper Mountain Conference on
!*                    Multigrid Methods, April 2-7, 1995, Copper Mountain, 
!*                    Colorado",
!*       editor    = "J. Mandel and S. McCormick",
!*       publisher = "NASA Langley Research Center",
!*       year      = "1995" }
!*    
!*    @inproceedings{HoSa93b,
!*       author    = "M. Holst and F. Saied",
!*       title     = "Multigrid and Domain Decomposition Methods for 
!*                    Electrostatics Problems",
!*       booktitle = "Domain Decomposition Methods in Science and Engineering
!*                    (Proceedings of the Seventh International Conference on
!*                    Domain Decomposition, October 27-30, 1993, 
!*                    The Pennsylvania State University)",
!*       editor    = "D. E. Keyes and J. Xu",
!*       publisher = "American Mathematical Society, Providence",
!*       year      = "1995" }
!*    
!*    @techreport{HoVa94a,
!*       author    = "M. Holst and S. Vandewalle",
!*       title     = "Schwarz Methods: to Symmetrize or not to Symmetrize",
!*       institution = "Applied Mathematics and CRPC,
!*                      California Institute of Technology",
!*       number    = "CRPC-94-13",
!*       year      = 1994 }
!*    
!*    @techreport{Hols94c,
!*       author    = "M. Holst",
!*       title     = "{An} {Algebraic} {Schwarz} {Theory}",
!*       institution = "Applied Mathematics and CRPC,
!*                      California Institute of Technology",
!*       number    = "CRPC-94-12",
!*       year      = 1994 }
!*    
!*    @techreport{Hols94f,
!*       author    = "M. Holst",
!*       title     = "A Robust and Efficient Numerical Method for Nonlinear
!*                    Protein Modeling Equations",
!*       institution = "Applied Mathematics and CRPC,
!*                      California Institute of Technology",
!*       number    = "CRPC-94-9",
!*       year      = 1994 }
!*    
!*    @techreport{HoSa93c,
!*       author    = "M. Holst and F. Saied",
!*       title     = "A Short Note Comparing Multigrid and Domain Decomposition
!*                    for Protein Modeling Equations",
!*       institution = "Applied Mathematics and CRPC,
!*                      California Institute of Technology",
!*       number    = "CRPC-94-10",
!*       year      = 1994 }
!*
!* =====================
!* PMG Bugs and Notes
!* =====================
!*
!* The  timing  routine,  found  in  "secd.c", is machine dependent.  We  have 
!* implemented  a  number  of  routines  for  different  machines, including a
!* generic  routine  for  a  standard  unix  machine with the getrusage system
!* routine.
!*
!* =============================
!* PMG Last Modification Date
!* =============================
!*
!* MG:  10-01-95
!* XMG: 10-01-95
!*
!* ///////////////////////////////////////////////////////////////////////////



      subroutine mgmain(source,v,eps,kap,mr1,mr2,mr3,d1,d2,d3,       &
                                            whichbc,nlev,errtol,itmax)
!***********************************************************************
!* pupose:
!* -------
!*
!*    this program uses a numerical method to solve the following
!*    three-dimensional, 2nd order linear elliptic partial 
!*    differential equation:
!*
!*            lu = f, u in omega
!*             u = g, u on boundary of omega
!*
!*    where
!*
!*            omega = [xmin,xmax] x [ymin,ymax] x [zmin,zmax]
!*
!*    the multigrid code requires the operator in the form:
!*
!*           - \nabla \cdot (a \nabla u) + c u = f
!*
!*    with;
!*
!*        a=(a1(x,y,z),a2(x,y,z),a3(x,y,z)), c(x,y,z), f(x,y,z), 
!*        all scalar functions (possibly discontinuous) on omega.
!*        (the discontinuities must be along grid lines on fine grid.)
!*        the boundary function g(x,y,z) is smooth on boundary of omega.
!*
!* user inputs:
!* ------------
!*
!*    the user selects one of several possible methods:
!*       (1) multigrid
!*       (2) ...
!*
!*    the calls to each of the routines are as follows:
!*
!*       call mgdriv (iparm,rparm,iwork,rwork,u,
!*   .      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
!*
!*    where the parameters are as follows:
!*
!*       iparm     = integer parameters for the solver
!*       rparm     = real parameters for the solver
!*       iwork     = some integer work space
!*       rwork     = some real work space
!*       u         = solution (perhaps some initial approximation)
!*       xf        = the x-coordinates of the tensor-product mesh
!*       yf        = the y-coordinates of the tensor-product mesh
!*       zf        = the z-coordinates of the tensor-product mesh
!*       gxcf      = the x=xmin and x=xmax boundary types and conditions
!*       gycf      = the y=ymin and y=ymax boundary types and conditions
!*       gzcf      = the z=zmin and z=zmax boundary types and conditions
!*       a1cf      = entry (1,1) of the diagonal tensor "a" above
!*       a2cf      = entry (2,2) of the diagonal tensor "a" above
!*       a3cf      = entry (3,3) of the diagonal tensor "a" above
!*       ccf       = the coefficient "c" above
!*       fcf       = the sourse function "f" above
!*       tcf       = the true solution of the pde if available (for testing)
!*                   (must always be provided as may be used for storage)
!*
!*    the sizes of the arrays must be as follows:
!*
!*       nxm    = number of x-points in the tensor-product mesh
!*       nym    = number of y-points in the tensor-product mesh
!*       nzm    = number of z-points in the tensor-product mesh
!*       nf     = nxm * nym * nzm (total number of nodes)
!*
!*       ifudge = 1000
!*       nfarr  = nf + nf/6 + ifudge
!*       niwk   = calculated by multigrid code based on input parameters
!*       nrwk   = calculated by multigrid code based on input parameters
!*
!*       iparm  = iparm(100)
!*       rparm  = rparm(100)
!*       iwork  = iwork(niwk)
!*       rwork  = rwork(nrwk)
!*       a1cf   = a1cf(nfarr)
!*       a2cf   = a2cf(nfarr)
!*       a3cf   = a3cf(nfarr)
!*       ccf    = ccf(nfarr)
!*       fcf    = fcf(nfarr)
!*       tcf    = tcf(nfarr)
!*       u      = u(nfarr)
!*       xf     = xf(nxm*5)
!*       yf     = yf(nym*5)
!*       zf     = zf(nzm*5)
!*       gxcf   = gxcf(nym*nzm*10)
!*       gycf   = gycf(nxm*nzm*10)
!*       gzcf   = gzcf(nxm*nym*10)
!*
!*   note that all of the above arrays are written over by the multigrid
!*   code during processing; there is no extra storage here.
!*
!*   the parameters that must be initialized in the iparm and rparm
!*   arrays in different ways for the three methods.  for example,
!*   the multigrid code requires the following:
!*
!*       *** integer parameters ***
!*       iparm(1)  = nrwk
!*       iparm(2)  = niwk
!*       iparm(3)  = nx
!*       iparm(4)  = ny
!*       iparm(5)  = nz
!*       iparm(6)  = nlev
!*       iparm(7)  = nu1
!*       iparm(8)  = nu2
!*       iparm(9)  = mgkey
!*       iparm(10) = itmax
!*       iparm(11) = istop
!*       iparm(12) = iinfo
!*       iparm(13) = irite
!*       iparm(14) = ipkey
!*       iparm(15) = ipcon
!*       iparm(16) = nonlin
!*       iparm(17) = mgprol
!*       iparm(18) = mgcoar
!*       iparm(19) = mgdisc
!*       iparm(20) = mgsmoo
!*       iparm(21) = mgsolv
!*       iparm(22) = iperf
!*
!*       *** real parameters ***
!*       rparm(1)  = errtol
!*       rparm(3)  = xmin
!*       rparm(4)  = xmax
!*       rparm(5)  = ymax
!*       rparm(6)  = ymax
!*       rparm(7)  = zmax
!*       rparm(8)  = zmax
!*       rparm(9)  = omegal
!*       rparm(10) = omegan
!*
!*    the following are returned from the multigrid solver:
!*
!*       u         = the computed solution
!*       iparm(50) = iters  (number of iters taken to reach tolerance)
!*       iparm(51) = ierror (=0 if tolerance reached, <>0 if problems)
!*
!* author:  michael holst
!***********************************************************************
      implicit         none
!*
!*    *** input parameter ***
      integer          mr1, mr2, mr3, itmax
      double precision d1, d2, d3 
      double precision source(mr1*mr2*mr3), v(mr1*mr2*mr3), &
         eps(mr1, mr2, mr3, 3), kap(mr1, mr2, mr3), &
         errtol

!*
!*    *** the three main storage allocation parameters ***
      integer          nrwk,niwk
!***********************************************************************
!* ***parameter        (ngrid=5,   nrwk=2465,      niwk=2000)
!* ***parameter        (ngrid=9,   nrwk=11885,     niwk=2000)
!* ***parameter        (ngrid=17,  nrwk=130149,    niwk=2000)
!* ***parameter        (ngrid=25,  nrwk=363620,    niwk=2000)
!* ***parameter        (ngrid=33,  nrwk=837892,    niwk=2000)
!      parameter        (ngrid=45,  nrwk=1187935,   niwk=600)
!*    parameter        (ngrid=49,  nrwk=2827484,   niwk=2000)
!***********************************************************************
!* ***parameter        (ngrid=65,  nrwk=3025573,   niwk=2000)
!* ***parameter        (ngrid=65,  nrwk=3456041,   niwk=2000)
!* ***parameter        (ngrid=65,  nrwk=6202291,   niwk=2000)
!***********************************************************************
!* ***parameter        (ngrid=97,  nrwk=9947081,   niwk=750)
!* ***parameter        (ngrid=97,  nrwk=11320640,  niwk=2000)
!* ***parameter        (ngrid=97,  nrwk=20447370,  niwk=2000)
!***********************************************************************
!* ***parameter        (ngrid=129, nrwk=23315466,  niwk=2000)
!***********************************************************************
!* ***parameter        (ngrid=193, nrwk=77735113,  niwk=2000)
!***********************************************************************
!* ***parameter        (ngrid=257, nrwk=183123727, niwk=2000)
!***********************************************************************
!* ***parameter        (ngrid=65,  nrwk=6202291,   niwk=2000)
!* ***parameter        (ngrid=65,  nrwk=4005291,   niwk=2000)
!* ***parameter        (ngrid=97,  nrwk=11320640,  niwk=2000)
!*
!*    *** other parameters ***
      integer          nxm,nym,nzm,nf,nfarr,ifudge
!*
!*    *** main storage *** 
      integer          nx,ny,nz,key,meth,nlev,mgcoar,mgdisc,mgsolv
      integer          nrwk_tot,niwk_tot
      integer          iparm(100)
      double precision rparm(100)
      integer, allocatable :: iwork(:)
      integer              :: whichbc(3)
      real(kind=8), allocatable :: rwork(:)
      real(kind=8), allocatable :: a1cf (:)
      real(kind=8), allocatable :: a2cf (:)
      real(kind=8), allocatable :: a3cf (:)
      real(kind=8), allocatable :: ccf  (:)
      real(kind=8), allocatable :: fcf  (:)
      real(kind=8), allocatable :: tcf  (:)
      real(kind=8), allocatable :: u    (:)
      real(kind=8), allocatable :: xf   (:,:)
      real(kind=8), allocatable :: yf   (:,:)
      real(kind=8), allocatable :: zf   (:,:)
      real(kind=8), allocatable :: gxcf (:,:,:)
      real(kind=8), allocatable :: gycf (:,:,:)
      real(kind=8), allocatable :: gzcf (:,:,:)
      real(kind=8) :: bf,oh,cputme
!*    *** start timer ***
      call tstart(bf,oh)
!*
!*    *** read input options ***
      nx=mr1
      ny=mr2
      nz=mr3
      nxm=mr1
      nym=mr2
      nzm=mr3
      call readit (iparm,rparm,nx,ny,nz,nlev,nrwk,niwk,key,meth,  &
                         errtol, itmax,whichbc)
!*
!*    *** print out the required sizes ***
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsolv = iparm(21)
      call mgsize(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev,nrwk,niwk)
!*
      nf=nxm*nym*nzm
      ifudge=100
      nfarr=nf+nf/6+ifudge
      allocate( iwork(niwk), rwork(nrwk) )
      allocate( a1cf (nfarr), a2cf (nfarr), a3cf (nfarr) )
      allocate(  ccf  (nfarr), fcf  (nfarr), tcf  (nfarr), u    (nfarr))
      allocate( xf   (nxm,5), yf   (nym,5), zf   (nzm,5) )
      allocate( gxcf (nym,nzm,10), gycf (nxm,nzm,10), gzcf (nxm,nym,10))
!*    *** some i/o ***
      niwk_tot = 100 + niwk
      nrwk_tot = 100 + nrwk + 7*nfarr &
         + nxm*5 + nym*5 + nzm*5      &
         + nym*nzm*10 + nxm*nzm*10 + nxm*nym*10


!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      print*,'% MGMAIN: TOTAL FP  STORAGE DECLARED:        ',nrwk_tot
!      print*,'% MGMAIN: TOTAL INT STORAGE DECLARED:        ',niwk_tot
!      print*,'% MGMAIN: max number of unknowns:            ',nf
!      print*,'% MGMAIN: max vector length for all levels:  ',nfarr
!      print*,'% MGMAIN: DECLARED fp  work array size:      ',nrwk
!      print*,'% MGMAIN: DECLARED int work array size:      ',niwk
!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

!*    *** this is to reinitialize iparm, rparm -- to be modified
      call readit (iparm,rparm,nx,ny,nz,nlev,nrwk,niwk,key,meth,errtol, &
                                                 itmax,whichbc)
!*    *** fill the coefficient arrays ***
      call fillco (iparm,rparm,nx,ny,nz,                                &
         xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,source,v,   &
         eps,kap,mr1,mr2,mr3,d1,d2,d3)
!*
!*    *** stop timer ***
      call tstop(bf,oh,cputme)
!      print*,'% MGMAIN: time to read in problem description: ',cputme
!*
!*    *** start timer ***
      call tstart(bf,oh)
!*
!*    *** call linear solvers ***
!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      if (meth .eq. 0) then
!         print*,'% MGMAIN: calling cgmg (linear) ... '
         call cgmgdriv (iparm,rparm,iwork,rwork,v,                      &
            xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,         &
            whichbc)
!      elseif (meth .eq. 1) then
!         print*,'% MGMAIN: calling newton (nonlinear) ... '
!         call newdriv (iparm,rparm,iwork,rwork,u,
!     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
!      elseif (meth .eq. 2) then
!         print*,'% MGMAIN: calling nmg (linear/nonlinear) ... '
!         call mgdriv (iparm,rparm,iwork,rwork,u,
!     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
!      elseif (meth .eq. 3) then
!         print*,'% MGMAIN: calling ncg (linear/nonlinear) ... '
!         call ncghsdriv (iparm,rparm,iwork,rwork,v,
!     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,
!     3      whichbc)
!      elseif (meth .eq. 4) then
!         print*,'% MGMAIN: calling nsor (linear/nonlinear) ... '
!         call nsordriv (iparm,rparm,iwork,rwork,u,
!     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
!      elseif (meth .eq. 5) then
!         print*,'% MGMAIN: calling ngsrb (linear/nonlinear) ... '
!         call ngsrbdriv (iparm,rparm,iwork,rwork,u,
!     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
!      elseif (meth .eq. 6) then
!         print*,'% MGMAIN: calling nwjac (linear/nonlinear) ... '
!         call nwjacdriv (iparm,rparm,iwork,rwork,u,
!     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
!      elseif (meth .eq. 7) then
!         print*,'% MGMAIN: calling nrich (linear/nonlinear) ... '
!         call nrichdriv (iparm,rparm,iwork,rwork,u,
!     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
      else
         print*,'From multigrid.  % MGMAIN: bad method key ... '
      endif
!      print*,'% MGMAIN: returned from iterative solver.'
!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!*
!*    *** stop timer ***
      call tstop(bf,oh,cputme)
!      print*,'% MGMAIN: total time spent in solve routine: ',cputme
!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!*
!*    *** write soln ***
!      call writit (iparm,rparm,nx,ny,nz,u,
!     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,key)

!      v( 1:nf ) = u( 1:nf )
    
!*
      deallocate( iwork )
      deallocate( rwork )
      deallocate( a1cf )
      deallocate( a2cf )
      deallocate( a3cf )
      deallocate( ccf  )
      deallocate( fcf  )
      deallocate( tcf  )
      deallocate( u    )
      deallocate( xf   )
      deallocate( yf   )
      deallocate( zf   )
      deallocate( gxcf )
      deallocate( gycf )
      deallocate( gzcf )

!*    *** end it ***
 99   continue
      return
      end subroutine

