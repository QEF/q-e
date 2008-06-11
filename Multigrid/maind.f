c* ///////////////////////////////////////////////////////////////////////////
c* @file    maind.f
c* @author  Michael Holst
c* @brief   Main FORTRAN driver for PMG.
c* @version $Id: maind.f,v 1.1 2008-06-11 10:47:37 degironc Exp $
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

c* ///////////////////////////////////////////////////////////////////////////
c*                           Welcome to PMG (MG/XMG).
c* ///////////////////////////////////////////////////////////////////////////
c*
c* ======================
c* The Author Information
c* ======================
c*
c* The computer codes PMG (MG/XMG, "the Code"), including FORTRAN/C/C++
c* language multilevel subroutines called "MG" and the C-language X-Window
c* interface called "XMG", were developed by:
c*
c*    Michael Holst                TELE:  (858) 534-4899
c*    Department of Mathematics    FAX:   (858) 534-5273
c*    UC San Diego, AP&M 5739      EMAIL: mholst@math.ucsd.edu
c*    La Jolla, CA 92093 USA       WEB:   http://www.scicomp.ucsd.edu/~mholst
c*
c* ====================================================
c* Licensing, permission to use, modify, and distribute
c* ====================================================
c*
c* The code has been placed under the GNU General Public License (GPL).  
c*
c* This means that essentially you can use and modify the code as you like, 
c* under a few conditions, such as acknowledging the authors of the original 
c* code in any derived works, such as new programs or research results using 
c* the code.
c*
c* Please have a look at the GPL which should have accompanied the code.
c* If it is missing, please write to the Free Software Foundation, Inc., 
c* 675 Mass Ave, Cambridge, MA 02139, USA.
c*
c* It would have been easier for me to just keep the code completely private, 
c* or to hand out very restricted use copies.  However, in the interests of 
c* freedom of information, the advancement of science, the future of the 
c* planet, and so on, I decided to distribute the code in this very free way.
c* The GNU world functions very much on the honor system.  I want you to have 
c* all of my source code because I think it will do more toward helping you 
c* understand the ideas in there than anything else.  The GNU idea is that in 
c* return for all of my work on the source code and the ideas behind it, all I 
c* ask for is some acknowledgment in your derived works for using my original 
c* code and the original ideas that went into the code.  For GNU to be viable, 
c* we all need to abide by this very minimal requirement.
c*
c* I'm typing this on a nearly all-GNU system right now, namely Linux 1.2.8, 
c* a GNU version of UN*X available under the same conditions as PMG.
c* I.e., you get all of the source code for the entire UN*X operating system.
c* For more information about Linux, see: http://sunsite.unc.edu/mdw/linux.html
c*
c* ====================================================================
c* Pointers to Papers and other Reference Materials Related to the Code
c* ====================================================================
c*
c* The methods and techniques  implemented  in  the  Code  were  developed  by 
c* Michael Holst as a part of his PhD research, and are described in detail in 
c* the PhD thesis:
c*
c*    @phdthesis{Hols93a,
c*       author    = "M. Holst",
c*       title     = "Multilevel Methods for the {Poisson-Boltzmann} Equation",
c*       note      = "Also published as Tech. Rep. UIUCDCS-R-03-1821",
c*       school    = "Numerical Computing Group,
c*                    Department of Computer Science,
c*                    University of Illinois at Urbana-Champaign",
c*       year      = 1993 }
c*
c* Also refer to the supporting papers:
c*
c*    @manual{Hols94d,
c*       author    = "M. Holst",
c*       title     = "The {Poisson-Boltzmann} Equation:
c*                    Analysis and Multilevel Numerical Solution",
c*       note      = "Monograph (updated and extended form of the
c*                    Ph.D. thesis~\cite{Hols93a})." }
c*    
c*    @article{SHHN94,
c*       author    = "R. Sampogna and J. Hecht and M. Holst and A. Nicholls
c*                    and B. Honig",
c*       title     = "Nonlinear {Poisson-Boltzmann} Calculation of {pKa} 
c*                    Values",
c*       journal   = "Biophysics",
c*       note      = "(In Progress?)",
c*       year      = 1994 }
c*    
c*    @article{HoVa94b,
c*       author    = "M. Holst and S. Vandewalle",
c*       title     = "Schwarz Methods: to Symmetrize or not to Symmetrize",
c*       journal   = SINUM,
c*       note      = "(Accepted)",
c*       year      = 1994 }
c*    
c*    @article{Hols94e,
c*       author    = "M. Holst and F. Saied",
c*       title     = "Numerical Solution of the Nonlinear {Poisson-Boltzmann}
c*                    Equation: Developing More Robust and Efficient Methods",
c*       journal   = JCC,
c*       volume    = "16",
c*       number    = "3",
c*       pages     = "337--364",
c*       year      = 1995 }
c*    
c*    @article{HoSa93a,
c*       author    = "M. Holst and F. Saied",
c*       title     = "Multigrid Solution of the {Poisson-Boltzmann} Equation",
c*       journal   = JCC,
c*       volume    = "14",
c*       number    = "1",
c*       pages     = "105--113",
c*       year      = 1993 }
c*    
c*    @article{HKSS94,
c*       author    = "M. Holst and R. Kozack and F. Saied and S. Subramaniam",
c*       title     = "Protein Electrostatics: Rapid Multigrid-based {Newton}
c*                    Algorithm for Solution of the Full Nonlinear
c*                    {Poisson-Boltzmann} Equation",
c*       journal   = "J. Biomol. Struct. Dyn.",
c*       volume    = "11",
c*       pages     = "1437--1445",
c*       year      = 1994 }
c*    
c*    @article{HKSS93b,
c*       author    = "M. Holst and R. Kozack and F. Saied and S. Subramaniam",
c*       title     = "Treatment of Electrostatic Effects in Proteins:
c*                    Multigrid-based-{Newton} Iterative Method for Solution 
c*                    of the Full Nonlinear {Poisson-Boltzmann} Equation",
c*       journal   = "Proteins: Structure, Function, and Genetics",
c*       volume    = "18",
c*       number    = "3",
c*       pages     = "231--245",
c*       year      = 1994 }
c*    
c*    @techreport{Hols94h,
c*       author    = "M. Holst",
c*       title     = "A Theoretical Analysis of the {Poisson-Boltzmann} 
c*                    Equation: A priori Estimates and Well-posedness",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       year      = 1994 }
c*    
c*    @inproceedings{HoVa95a,
c*       author    = "M. Holst and S. Vandewalle",
c*       title     = "Schwarz Methods: To Symmetrize or not to Symmetrize",
c*       booktitle = "Proceedings of the Seventh Copper Mountain Conference on
c*                    Multigrid Methods, April 2-7, 1995, Copper Mountain, 
c*                    Colorado",
c*       editor    = "J. Mandel and S. McCormick",
c*       publisher = "NASA Langley Research Center",
c*       year      = "1995" }
c*    
c*    @inproceedings{HoSa93b,
c*       author    = "M. Holst and F. Saied",
c*       title     = "Multigrid and Domain Decomposition Methods for 
c*                    Electrostatics Problems",
c*       booktitle = "Domain Decomposition Methods in Science and Engineering
c*                    (Proceedings of the Seventh International Conference on
c*                    Domain Decomposition, October 27-30, 1993, 
c*                    The Pennsylvania State University)",
c*       editor    = "D. E. Keyes and J. Xu",
c*       publisher = "American Mathematical Society, Providence",
c*       year      = "1995" }
c*    
c*    @techreport{HoVa94a,
c*       author    = "M. Holst and S. Vandewalle",
c*       title     = "Schwarz Methods: to Symmetrize or not to Symmetrize",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-13",
c*       year      = 1994 }
c*    
c*    @techreport{Hols94c,
c*       author    = "M. Holst",
c*       title     = "{An} {Algebraic} {Schwarz} {Theory}",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-12",
c*       year      = 1994 }
c*    
c*    @techreport{Hols94f,
c*       author    = "M. Holst",
c*       title     = "A Robust and Efficient Numerical Method for Nonlinear
c*                    Protein Modeling Equations",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-9",
c*       year      = 1994 }
c*    
c*    @techreport{HoSa93c,
c*       author    = "M. Holst and F. Saied",
c*       title     = "A Short Note Comparing Multigrid and Domain Decomposition
c*                    for Protein Modeling Equations",
c*       institution = "Applied Mathematics and CRPC,
c*                      California Institute of Technology",
c*       number    = "CRPC-94-10",
c*       year      = 1994 }
c*
c* =====================
c* PMG Bugs and Notes
c* =====================
c*
c* The  timing  routine,  found  in  "secd.c", is machine dependent.  We  have 
c* implemented  a  number  of  routines  for  different  machines, including a
c* generic  routine  for  a  standard  unix  machine with the getrusage system
c* routine.
c*
c* =============================
c* PMG Last Modification Date
c* =============================
c*
c* MG:  10-01-95
c* XMG: 10-01-95
c*
c* ///////////////////////////////////////////////////////////////////////////



      subroutine mgmain(source,v,eps,kap,mr1,mr2,mr3,d1,d2,d3,        
     2                                       whichbc,nlev,errtol,itmax)
c***********************************************************************
c* pupose:
c* -------
c*
c*    this program uses a numerical method to solve the following
c*    three-dimensional, 2nd order linear elliptic partial 
c*    differential equation:
c*
c*            lu = f, u in omega
c*             u = g, u on boundary of omega
c*
c*    where
c*
c*            omega = [xmin,xmax] x [ymin,ymax] x [zmin,zmax]
c*
c*    the multigrid code requires the operator in the form:
c*
c*           - \nabla \cdot (a \nabla u) + c u = f
c*
c*    with;
c*
c*        a=(a1(x,y,z),a2(x,y,z),a3(x,y,z)), c(x,y,z), f(x,y,z), 
c*        all scalar functions (possibly discontinuous) on omega.
c*        (the discontinuities must be along grid lines on fine grid.)
c*        the boundary function g(x,y,z) is smooth on boundary of omega.
c*
c* user inputs:
c* ------------
c*
c*    the user selects one of several possible methods:
c*       (1) multigrid
c*       (2) ...
c*
c*    the calls to each of the routines are as follows:
c*
c*       call mgdriv (iparm,rparm,iwork,rwork,u,
c*   .      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf)
c*
c*    where the parameters are as follows:
c*
c*       iparm     = integer parameters for the solver
c*       rparm     = real parameters for the solver
c*       iwork     = some integer work space
c*       rwork     = some real work space
c*       u         = solution (perhaps some initial approximation)
c*       xf        = the x-coordinates of the tensor-product mesh
c*       yf        = the y-coordinates of the tensor-product mesh
c*       zf        = the z-coordinates of the tensor-product mesh
c*       gxcf      = the x=xmin and x=xmax boundary types and conditions
c*       gycf      = the y=ymin and y=ymax boundary types and conditions
c*       gzcf      = the z=zmin and z=zmax boundary types and conditions
c*       a1cf      = entry (1,1) of the diagonal tensor "a" above
c*       a2cf      = entry (2,2) of the diagonal tensor "a" above
c*       a3cf      = entry (3,3) of the diagonal tensor "a" above
c*       ccf       = the coefficient "c" above
c*       fcf       = the sourse function "f" above
c*       tcf       = the true solution of the pde if available (for testing)
c*                   (must always be provided as may be used for storage)
c*
c*    the sizes of the arrays must be as follows:
c*
c*       nxm    = number of x-points in the tensor-product mesh
c*       nym    = number of y-points in the tensor-product mesh
c*       nzm    = number of z-points in the tensor-product mesh
c*       nf     = nxm * nym * nzm (total number of nodes)
c*
c*       ifudge = 1000
c*       nfarr  = nf + nf/6 + ifudge
c*       niwk   = calculated by multigrid code based on input parameters
c*       nrwk   = calculated by multigrid code based on input parameters
c*
c*       iparm  = iparm(100)
c*       rparm  = rparm(100)
c*       iwork  = iwork(niwk)
c*       rwork  = rwork(nrwk)
c*       a1cf   = a1cf(nfarr)
c*       a2cf   = a2cf(nfarr)
c*       a3cf   = a3cf(nfarr)
c*       ccf    = ccf(nfarr)
c*       fcf    = fcf(nfarr)
c*       tcf    = tcf(nfarr)
c*       u      = u(nfarr)
c*       xf     = xf(nxm*5)
c*       yf     = yf(nym*5)
c*       zf     = zf(nzm*5)
c*       gxcf   = gxcf(nym*nzm*10)
c*       gycf   = gycf(nxm*nzm*10)
c*       gzcf   = gzcf(nxm*nym*10)
c*
c*   note that all of the above arrays are written over by the multigrid
c*   code during processing; there is no extra storage here.
c*
c*   the parameters that must be initialized in the iparm and rparm
c*   arrays in different ways for the three methods.  for example,
c*   the multigrid code requires the following:
c*
c*       *** integer parameters ***
c*       iparm(1)  = nrwk
c*       iparm(2)  = niwk
c*       iparm(3)  = nx
c*       iparm(4)  = ny
c*       iparm(5)  = nz
c*       iparm(6)  = nlev
c*       iparm(7)  = nu1
c*       iparm(8)  = nu2
c*       iparm(9)  = mgkey
c*       iparm(10) = itmax
c*       iparm(11) = istop
c*       iparm(12) = iinfo
c*       iparm(13) = irite
c*       iparm(14) = ipkey
c*       iparm(15) = ipcon
c*       iparm(16) = nonlin
c*       iparm(17) = mgprol
c*       iparm(18) = mgcoar
c*       iparm(19) = mgdisc
c*       iparm(20) = mgsmoo
c*       iparm(21) = mgsolv
c*       iparm(22) = iperf
c*
c*       *** real parameters ***
c*       rparm(1)  = errtol
c*       rparm(3)  = xmin
c*       rparm(4)  = xmax
c*       rparm(5)  = ymax
c*       rparm(6)  = ymax
c*       rparm(7)  = zmax
c*       rparm(8)  = zmax
c*       rparm(9)  = omegal
c*       rparm(10) = omegan
c*
c*    the following are returned from the multigrid solver:
c*
c*       u         = the computed solution
c*       iparm(50) = iters  (number of iters taken to reach tolerance)
c*       iparm(51) = ierror (=0 if tolerance reached, <>0 if problems)
c*
c* author:  michael holst
c***********************************************************************
      implicit         none
c*
c*    *** input parameter ***
      integer          mr1, mr2, mr3, itmax
      double precision d1, d2, d3 
      double precision source(mr1*mr2*mr3), v(mr1*mr2*mr3), 
     2   eps(mr1, mr2, mr3, 3), kap(mr1, mr2, mr3),
     3   errtol

c*
c*    *** the three main storage allocation parameters ***
      integer          nrwk,niwk
c***********************************************************************
c* ***parameter        (ngrid=5,   nrwk=2465,      niwk=2000)
c* ***parameter        (ngrid=9,   nrwk=11885,     niwk=2000)
c* ***parameter        (ngrid=17,  nrwk=130149,    niwk=2000)
c* ***parameter        (ngrid=25,  nrwk=363620,    niwk=2000)
c* ***parameter        (ngrid=33,  nrwk=837892,    niwk=2000)
c      parameter        (ngrid=45,  nrwk=1187935,   niwk=600)
c*    parameter        (ngrid=49,  nrwk=2827484,   niwk=2000)
c***********************************************************************
c* ***parameter        (ngrid=65,  nrwk=3025573,   niwk=2000)
c* ***parameter        (ngrid=65,  nrwk=3456041,   niwk=2000)
c* ***parameter        (ngrid=65,  nrwk=6202291,   niwk=2000)
c***********************************************************************
c* ***parameter        (ngrid=97,  nrwk=9947081,   niwk=750)
c* ***parameter        (ngrid=97,  nrwk=11320640,  niwk=2000)
c* ***parameter        (ngrid=97,  nrwk=20447370,  niwk=2000)
c***********************************************************************
c* ***parameter        (ngrid=129, nrwk=23315466,  niwk=2000)
c***********************************************************************
c* ***parameter        (ngrid=193, nrwk=77735113,  niwk=2000)
c***********************************************************************
c* ***parameter        (ngrid=257, nrwk=183123727, niwk=2000)
c***********************************************************************
c* ***parameter        (ngrid=65,  nrwk=6202291,   niwk=2000)
c* ***parameter        (ngrid=65,  nrwk=4005291,   niwk=2000)
c* ***parameter        (ngrid=97,  nrwk=11320640,  niwk=2000)
c*
c*    *** other parameters ***
      integer          nxm,nym,nzm,nf,nfarr,ifudge
c*
c*    *** main storage *** 
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
c*    *** start timer ***
      call tstart(bf,oh)
c*
c*    *** read input options ***
      nx=mr1
      ny=mr2
      nz=mr3
      nxm=mr1
      nym=mr2
      nzm=mr3
      call readit (iparm,rparm,nx,ny,nz,nlev,nrwk,niwk,key,meth, 
     2                   errtol, itmax,whichbc)
c*
c*    *** print out the required sizes ***
      mgcoar = iparm(18)
      mgdisc = iparm(19)
      mgsolv = iparm(21)
      call mgsize(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev,nrwk,niwk)
c*
      nf=nxm*nym*nzm
      ifudge=100
      nfarr=nf+nf/6+ifudge
      allocate( iwork(niwk), rwork(nrwk) )
      allocate( a1cf (nfarr), a2cf (nfarr), a3cf (nfarr) )
      allocate(  ccf  (nfarr), fcf  (nfarr), tcf  (nfarr), u    (nfarr))
      allocate( xf   (nxm,5), yf   (nym,5), zf   (nzm,5) )
      allocate( gxcf (nym,nzm,10), gycf (nxm,nzm,10), gzcf (nxm,nym,10))
c*    *** some i/o ***
      niwk_tot = 100 + niwk
      nrwk_tot = 100 + nrwk + 7*nfarr
     2   + nxm*5 + nym*5 + nzm*5
     3   + nym*nzm*10 + nxm*nzm*10 + nxm*nym*10


!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!      print*,'% MGMAIN: TOTAL FP  STORAGE DECLARED:        ',nrwk_tot
!      print*,'% MGMAIN: TOTAL INT STORAGE DECLARED:        ',niwk_tot
!      print*,'% MGMAIN: max number of unknowns:            ',nf
!      print*,'% MGMAIN: max vector length for all levels:  ',nfarr
!      print*,'% MGMAIN: DECLARED fp  work array size:      ',nrwk
!      print*,'% MGMAIN: DECLARED int work array size:      ',niwk
!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

c*    *** this is to reinitialize iparm, rparm -- to be modified
      call readit (iparm,rparm,nx,ny,nz,nlev,nrwk,niwk,key,meth,errtol, 
     2                                           itmax,whichbc)
c*    *** fill the coefficient arrays ***
      call fillco (iparm,rparm,nx,ny,nz,
     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,source,v,
     3   eps,kap,mr1,mr2,mr3,d1,d2,d3)
c*
c*    *** stop timer ***
      call tstop(bf,oh,cputme)
!      print*,'% MGMAIN: time to read in problem description: ',cputme
c*
c*    *** start timer ***
      call tstart(bf,oh)
c*
c*    *** call linear solvers ***
!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      if (meth .eq. 0) then
!         print*,'% MGMAIN: calling cgmg (linear) ... '
         call cgmgdriv (iparm,rparm,iwork,rwork,v,
     2      xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,
     3      whichbc)
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
c*
c*    *** stop timer ***
      call tstop(bf,oh,cputme)
!      print*,'% MGMAIN: total time spent in solve routine: ',cputme
!      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
c*
c*    *** write soln ***
!      call writit (iparm,rparm,nx,ny,nz,u,
!     2   xf,yf,zf,gxcf,gycf,gzcf,a1cf,a2cf,a3cf,ccf,fcf,tcf,key)

!      v( 1:nf ) = u( 1:nf )
    
c*
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

c*    *** end it ***
 99   continue
      return
      end subroutine

