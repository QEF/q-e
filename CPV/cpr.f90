!
! Copyright (C) 2002-2004 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "../include/machine.h"

!
!=======================================================================
!
   subroutine cprmain( tau, fion_out, etot_out )
!
!=======================================================================
!***  Molecular Dynamics using Density-Functional Theory   ****
!***  this is a Car-Parrinello program using Vanderbilt pseudopotentials
!***********************************************************************
!***  based on version 11 of cpv code including ggapw 07/04/99
!***  copyright Alfredo Pasquarello 10/04/1996
!***  parallelized and converted to f90 by Paolo Giannozzi (2000),
!***  using parallel FFT written for PWSCF by Stefano de Gironcoli
!***  PBE added by Michele Lazzeri (2000)
!***  variable-cell dynamics by Andrea Trave (1998-2000)
!***********************************************************************
!***  appropriate citation for use of this code:
!***  Car-Parrinello method    R. Car and M. Parrinello, PRL 55, 2471 (1985) 
!***  current implementation   A. Pasquarello, K. Laasonen, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRL 69, 1982 (1992);
!***                           K. Laasonen, A. Pasquarello, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRB 47, 10142 (1993).
!***  implementation gga       A. Dal Corso, A. Pasquarello, A. Baldereschi,
!***                           and R. Car, PRB 53, 1180 (1996).
!***********************************************************************
!***  
!***  f90 version, with dynamical allocation of memory
!***  Variables that do not change during the dynamics are in modules
!***  (with some exceptions) All other variables are passed as arguments
!***********************************************************************
!***
!*** fft : uses machine's own complex fft routines, two real fft at the time
!*** ggen: g's only for positive halfspace (g>)
!*** all routines : keep equal c(g) and [c(-g)]*
!***
!***********************************************************************
!    general variables:
!     delt           = delta t
!     emass          = electron mass (fictitious)
!     dt2bye         = 2*delt/emass
!***********************************************************************
!
      use control_flags, only: iprint, thdyn, tpre, tbuff, iprsta, trhor, &
            tfor, tvlocw, trhow, taurdr, tprnfor
      use control_flags, only: ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
            tnosep, trane, tranp, tsdp, tcp, tcap, ampre, amprp, tnoseh

      use core, only: nlcc
      use core, only: deallocate_core
      use cvan, only: nvb, nhx, nhsa
      use cvan, only: deallocate_cvan
      use energies, only: eht, epseu, exc, etot, eself, enl, ekin
      use elct, only: nx, n, ispin, f, nspin, nel, iupdwn, nupdwn
      use elct, only: deallocate_elct
      use gvec, only: tpiba2, ng
      use gvec, only: deallocate_gvec
      use gvecs, only: ngs
      use gvecb, only: ngb
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use ions_base, only: na, nat, pmass, nas => nax, nsp, ipp, rcmax
      use ions_base, only: ind_srt
      use grid_dimensions, only: nnr => nnrx, nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3
      use cell_base, only: omega, alat
      use cell_base, only: h, hold, deth, wmass
      use smooth_grid_dimensions, only: nnrsx, nr1s, nr2s, nr3s
      use smallbox_grid_dimensions, only: nnrb => nnrbx, nr1b, nr2b, nr3b
      use pseu, only: vps, rhops
      use pseu, only: deallocate_pseu
      use work
      use work_box, only: qv, deallocate_work_box
      use io_global, ONLY: io_global_start, stdout, ionode
      use mp_global, ONLY: mp_global_start
      use mp, ONLY: mp_sum, mp_barrier
      use para_mod
      use dener
      use derho
      use dpseu
      use cdvan
      use stre
      use gvecw, only: ggp, agg => ecutz, sgg => ecsig, e0gg => ecfix
      use restart
      use parameters, only: nacx, natx, nsx, nbndxx
      use constants, only: pi, factem
      use io_files, only: psfile, pseudo_dir
      use input_cp, only: iosys
      use qgb_mod, only: deallocate_qgb_mod
      use dqgb_mod, only: deallocate_dqgb_mod
      use qradb_mod, only: deallocate_qradb_mod
      use dqrad_mod, only: deallocate_dqrad_mod
      use betax, only: deallocate_betax
      use input_parameters, only: outdir
      use wave_base, only: wave_steepest, wave_verlet
      use wave_base, only: wave_speed2
      USE control_flags, ONLY : conv_elec, tconvthrs
      USE check_stop, ONLY : check_stop_now

! wavefunctions
!
      use wavefunctions_module, only: c0, cm, phi => cp
      use wavefunctions_module, only: deallocate_wavefunctions
!
!
      implicit none
!
! input variables
!
      real(kind=8) :: tau(3,*)
      real(kind=8) :: fion_out(3,*)
      real(kind=8) :: etot_out

!
!
! control variables
!
      logical twall, tbump
      logical thdiag
      logical tfirst, tlast
      logical tstop
      logical tconv
      real(kind=8) :: delta_etot
!
! structure factors e^{-ig*R}
!
      complex(kind=8), allocatable:: ei1(:,:,:),  ei2(:,:,:),  ei3(:,:,:)
      complex(kind=8), allocatable:: eigr(:,:,:)
!
! structure factors (summed over atoms of the same kind)
!
      complex(kind=8), allocatable:: sfac(:,:)
!
! indexes, positions, and structure factors for the box grid
!
      integer irb(3,natx,nsx)
      real(kind=8) taub(3,natx,nsx)
      complex(kind=8), allocatable:: eigrb(:,:,:)
! 
! charge densities and potentials
!     rhog  = charge density in g space
!     rhor  = charge density in r space (dense grid)
!     rhos  = charge density in r space (smooth grid)
!     rhoc  = core charge density in real space (dense grid)
!
      complex(kind=8), allocatable:: rhog(:,:)
      real(kind=8), allocatable:: rhor(:,:), rhos(:,:), rhoc(:)
!
! nonlocal projectors:
!     bec   = scalar product of projectors and wave functions
!     betae = nonlocal projectors in g space = beta x e^(-ig.R) 
!     becdr = <betae|g|psi> used in force calculation
!     rhovan= \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
!     deeq  = \int V_eff(r) q_lm(r) dr
!
      complex(kind=8), allocatable::   betae(:,:)
      real(kind=8), allocatable:: bec(:,:), becdr(:,:,:)
      real(kind=8), allocatable:: bephi(:,:), becp(:,:)
      real(kind=8), allocatable:: rhovan(:,:,:), deeq(:,:,:,:)
!
!  mass preconditioning
!
      real(kind=8), allocatable:: ema0bg(:)
      real(kind=8), allocatable:: emadt2(:)
      real(kind=8), allocatable:: emaver(:)
      real(kind=8), allocatable:: emainv(:)
      real(kind=8)  emaec
!
!  constraints (lambda at t, lambdam at t-dt, lambdap at t+dt)
!
      real(kind=8), allocatable:: lambda(:,:), lambdam(:,:), lambdap(:,:)
!
!  ionic positions, center of mass position
!
      real(kind=8) tau0(3,natx,nsx), taum(3,natx,nsx), taup(3,natx,nsx)
      real(kind=8) cdm0(3),cdmvel(3)
!
!  forces on ions
!
      real(kind=8) fion(3,natx,nsx), fionm(3,natx,nsx)
      integer iforce(3,natx,nsx)
!
! for variable cell dynamics: scaled tau
!
      real(kind=8) taus(3,natx,nsx)
      real(kind=8) f_(nbndxx)
      real(kind=8) ispin_(nbndxx)
      integer iforceh(3,3)
!
      integer maxit
!
! work variables
!
      real(kind=8) celldm(6), ecut, ecutw
      real(kind=8) acc(nacx)
      complex(kind=8), allocatable:: c2(:), c3(:)
      complex(kind=8)  speed
      real(kind=8)                                                      & 
     &       tempp, xnhe0, vnhp, xnhp0, xnhpm, verl1, verl2, verl3,     &
     &       fccc, xnhem, vnhe, savee, saveh, savep, press,             &
     &       enthal, epot, xnhpp, xnhep, epre, enow, tps, econs, econt, &
     &       fricp, greasp, eps, qnp, tempw, qne,                       &
     &       frice,  grease, emass, delt, ccc, bigr, dt2,               &
     &       dt2by2, twodel, gausp, dt2bye, gkbt, dt2hbe
      real(kind=8) ekinc0, ekinp, ekinpr, ekincm, ekinc, ekincw
      integer nnn, is, nacc, ia, j, iter, nfi, i, isa, ipos
!
! work variables, 2
!
      real(kind=8) tausm(3,natx,nsx),tausp(3,natx,nsx)
      real(kind=8) vels(3,natx,nsx),velsm(3,natx,nsx),velsp(3,natx,nsx)
      real(kind=8) hnew(3,3),velh(3,3),hgamma(3,3)
      real(kind=8) cdm(3)
      real(kind=8) qr(3)
      real(kind=8) xnhh0(3,3),xnhhm(3,3),xnhhp(3,3),vnhh(3,3),temphh(3,3)
!
      integer ibrav, k, ii, l, m
      real(kind=8) ekinh, alfar, temphc, alfap, frich, tolp,    &
     &     factp, temp1, temp2, temph, greash, qnh, randy
      real(kind=8) ftmp

      logical :: twmass
      character(len=256) :: filename
      character(len=256) :: dirname
      integer :: strlen, dirlen
!
!     CP loop starts here
!
      call start_clock( 'initialize' ) 
      etot_out = 0.0d0

!
!     ==================================================================
!     ====  units and constants                                     ====
!     ====  1 hartree           = 1 a.u.                            ====
!     ====  1 bohr radius       = 1 a.u. = 0.529167 Angstrom        ====
!     ====  1 rydberg           = 1/2 a.u.                          ====
!     ====  1 electron volt     = 1/27.212 a.u.                     ====
!     ====  1 kelvin *k-boltzm. = 1/(27.212*11606) a.u.'='3.2e-6 a.u====
!     ====  1 second            = 1/(2.4189)*1.e+17 a.u.            ====
!     ====  1 proton mass       = 1822.89 a.u.                      ====
!     ====  1 tera              = 1.e+12                            ====
!     ====  1 pico              = 1.e-12                            ====
!     ==================================================================

      factp   = 3.3989 * 0.00001
!
!     ==================================================================
!     read input from standard input (unit 5)
!     ==================================================================

      call iosys( nbeg , ndr , ndw , nomore , iprint                       &
     & , delt , emass , emaec  , tsde , frice , grease , twall              &
     & , tortho , eps , maxit , trane , ampre , tranp , amprp                &
     & , tfor , tsdp , fricp , greasp , tcp , tcap , tolp , trhor , trhow , tvlocw &
     & , tnosep , qnp , tempw , tnosee , qne , ekincw                 &
     & , tpre , thdyn , thdiag , twmass , wmass , frich , greash , press   &
     & , tnoseh , qnh , temph , celldm , ibrav , tau0 , ecutw , ecut , iforce &
     & , nat , nsp , na , pmass , rcmax , ipp , f_ , nel , nspin , nupdwn  &
     & , iupdwn , n , nx , nr1 , nr2 , nr3 , omega , alat , a1 , a2 , a3  &
     & , nr1b , nr2b , nr3b , nr1s , nr2s , nr3s , agg , sgg , e0gg &
     & , psfile , pseudo_dir, iprsta, ispin_ )
      allocate( f( nx ) )
      f( :   ) = 0.0d0
      f( 1:n ) = f_( 1:n )
      allocate( ispin( nx ) )
      ispin( :   ) = 0.0d0
      ispin( 1:n ) = ispin_( 1:n )

!     ==================================================================
!
!     general variables
!
      tfirst = .true.
      tlast  = .false.
      nacc = 5
!
      gausp = delt * sqrt(tempw/factem)
      twodel = 2.d0 * delt
      dt2 = delt * delt
      dt2by2 = .5d0 * dt2
      dt2bye = dt2/emass
      dt2hbe = dt2by2/emass

      if (ionode) then

         dirlen = index(outdir,' ') - 1
         filename = 'fort.8'
         if( dirlen >= 1 ) then
           filename = outdir(1:dirlen) // '/' // filename
         end if
         strlen  = index(filename,' ') - 1
         WRITE( stdout, * ) ' UNIT8 = ', filename
         OPEN(unit=8, file=filename(1:strlen), status='unknown')

         filename = 'fort.77'
         if( dirlen >= 1 ) then
           filename = outdir(1:dirlen) // '/' // filename
         end if
         strlen  = index(filename,' ') - 1
         OPEN(unit=77, file=filename(1:strlen), status='unknown')

         filename = 'fort.78'
         if( dirlen >= 1 ) then
           filename = outdir(1:dirlen) // '/' // filename
         end if
         strlen  = index(filename,' ') - 1
         OPEN(unit=78, file=filename(1:strlen), status='unknown')

         filename = 'fort.79'
         if( dirlen >= 1 ) then
           filename = outdir(1:dirlen) // '/' // filename
         end if
         strlen  = index(filename,' ') - 1
         OPEN(unit=79, file=filename(1:strlen), status='unknown')

      end if

!
!     ==================================================================
!     initialize g-vectors, fft grids
!     ==================================================================

      call init( ibrav, celldm, ecut, ecutw, tranp, amprp, ndr, nbeg, tfirst,  &
           twmass, thdiag, iforceh, tau0, taus, delt )

      WRITE( stdout,*) ' out from init'
!
!     more initialization requiring atomic positions
!
      nas = MAXVAL( na( 1 : nsp ) )
      if( iprsta > 1 ) then
         WRITE( stdout,*) ' tau0 '
         WRITE( stdout,'(3f14.8)') (((tau0(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
      endif

!
!     ==================================================================
!     allocate and initialize nonlocal potentials
!     ==================================================================

      call nlinit

      WRITE( stdout,*) ' out from nlinit'

!
!     ==================================================================
!     allocation of all arrays not already allocated in init and nlinit
!     ==================================================================
!
      allocate(c0(ngw,nx,1,1))
      allocate(cm(ngw,nx,1,1))
      allocate(phi(ngw,nx,1,1))
      allocate(wrk2(ngw,max(nas,n)))
      allocate(eigr(ngw,nas,nsp))
      allocate(eigrb(ngb,nas,nsp))
      allocate(sfac(ngs,nsp))
      allocate(rhops(ngs,nsp))
      allocate(vps(ngs,nsp))
      allocate(rhor(nnr,nspin))
      allocate(rhos(nnrsx,nspin))
      allocate(rhog(ng,nspin))
      if (nlcc > 0) allocate(rhoc(nnr))
      allocate(wrk1(nnr))
      allocate(qv(nnrb))
      allocate(c2(ngw))
      allocate(c3(ngw))
      allocate(ema0bg(ngw))
      allocate(lambda(nx,nx))
      allocate(lambdam(nx,nx))
      allocate(lambdap(nx,nx))
      allocate(ei1(-nr1:nr1,nas,nsp))
      allocate(ei2(-nr2:nr2,nas,nsp))
      allocate(ei3(-nr3:nr3,nas,nsp))
      allocate(betae(ngw,nhsa))
      allocate(becdr(nhsa,n,3))
      allocate(bec  (nhsa,n))
      allocate(bephi(nhsa,n))
      allocate(becp (nhsa,n))
      allocate(deeq(nat,nhx,nhx,nspin))
      allocate(rhovan(nat,nhx*(nhx+1)/2,nspin))
      allocate(dbec (nhsa,n,3,3))
      allocate(dvps(ngs,nsp))
      allocate(drhops(ngs,nsp))
      allocate(drhog(ng,nspin,3,3))
      allocate(drhor(nnr,nspin,3,3))
      allocate(drhovan(nat,nhx*(nhx+1)/2,nspin,3,3))
#ifdef __PARA
      allocate(aux(nnr))
#endif
      deeq(:,:,:,:) = 0.d0
!
 666  continue

!
!
      temp1=tempw+tolp
      temp2=tempw-tolp
      gkbt = 3.*nat*tempw/factem
      press = press*factp
!     ==========================================================
      do is=1,nsp
         do ia=1,na(is)
            do i=1,3
               tausm(i,ia,is)=taus(i,ia,is)
               tausp(i,ia,is)=0.
               taum(i,ia,is)=tau0(i,ia,is)
               taup(i,ia,is)=0.
               vels(i,ia,is)=0.
               velsm(i,ia,is)=0.
               velsp(i,ia,is)=0.
            end do
         end do
      enddo
!
      hnew=h
!
      lambda(:,:)=0.d0
      cm(:,:,1,1) = (0.d0, 0.d0)
      c0(:,:,1,1) = (0.d0, 0.d0)
!
!     mass preconditioning: ema0bg(i) = ratio of emass(g=0) to emass(g)
!     for g**2>emaec the electron mass ema0bg(g) rises quadratically
!
      do i=1,ngw
         ema0bg(i)=1./max(1.d0,tpiba2*ggp(i)/emaec)
         if(iprsta.ge.10)print *,i,' ema0bg(i) ',ema0bg(i)
      end do

      !WRITE( stdout, * ) 'NBEG = ', nbeg
      !fion_out(1:3,1:nat) = 0.0d0
      !etot_out = 0.0d0
      !RETURN
!
      if ( nbeg < 0 ) then

         !======================================================================
         !    nbeg = -1 or nbeg = -2 or nbeg = -3
         !======================================================================

         if( nbeg == -1 ) then
           call readfile_new                                            &
     &     ( 0, ndr,h,hold,nfi,cm(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc,         &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion)
         endif
!     
         call phfac( tau0, ei1, ei2, ei3, eigr )
!
         call initbox ( tau0, taub, irb )
!
         call phbox( taub, eigrb )
!
         if( iprsta > 2 ) then
            do is=1,nvb
               WRITE( stdout,'(/,2x,''species= '',i2)') is 
               do ia=1,na(is)
                  WRITE( stdout,2000) ia, (irb(i,ia,is),i=1,3)
 2000             format(2x,'atom= ',i3,' irb1= ',i3,' irb2= ',i3,      &
     &                 ' irb3= ',i3) 
               end do
            end do
         endif
!
         if(trane) then
!       
!     random initialization
!
            call randin(1,n,gstart,ngw,ampre,cm)

         else if(nbeg.eq.-3) then
!       
!     gaussian initialization
!
            call gausin(eigr,cm)
         end if
!
!     prefor calculates betae (used by graham)
!
         call prefor(eigr,betae)
         call graham(betae,bec,cm)
         if(iprsta.ge.3) call dotcsc(eigr,cm)
!     
         nfi=0
!
!     strucf calculates the structure factor sfac
!
         call strucf(ei1,ei2,ei3,sfac)
         call formf(tfirst,eself)
         call calbec (1,nsp,eigr,cm,bec)
         if (tpre) call caldbec(1,nsp,eigr,cm)
         call rhoofr (nfi,cm,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
         if(iprsta.gt.0) WRITE( stdout,*) ' out from rhoofr'
!
!     put core charge (if present) in rhoc(r)
!
         if (nlcc.gt.0) call set_cc(irb,eigrb,rhoc)
!
         call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
     &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
         do i=1,3
            do j=1,3
               stress(i,j)=-1.d0/omega*(detot(i,1)*h(j,1)+              &
     &                      detot(i,2)*h(j,2)+detot(i,3)*h(j,3))
            enddo
         enddo
         if(iprsta.gt.0) WRITE( stdout,*) ' out from vofrho'
         if(iprsta.gt.2) then
            WRITE( stdout,*) ' fion '
            WRITE( stdout,'(3f14.8)')                                         &
     &                   (((fion(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
         end if
! 
!     forces for eigenfunctions
!
!     newd calculates deeq and a contribution to fion
!
         call newd(rhor,irb,eigrb,rhovan,deeq,fion)
         WRITE( stdout,*) ' out from newd'
         call prefor(eigr,betae)
!
!     if n is odd => c(*,n+1)=0
!
         ALLOCATE( emadt2( ngw ) )
         ccc = dt2hbe
         if(tsde) ccc = dt2bye
         emadt2 = ccc * ema0bg

         do i=1,n,2
            call dforce(bec,deeq,betae,i,cm(1,i,1,1),cm(1,i+1,1,1),c2,c3,rhos)
            call wave_steepest( c0(:,i,1,1), cm(:,i,1,1), emadt2, c2 )
            call wave_steepest( c0(:,i+1,1,1), cm(:,i+1,1,1), emadt2, c3 )
         end do

         DEALLOCATE( emadt2 )
         WRITE( stdout,*) ' out from dforce'
!
!     buffer for wavefunctions is unit 21
!
         if(tbuff) rewind 21
!
!     nlfq needs deeq calculated in newd
!
         if ( tfor .or. tprnfor ) call nlfq(cm,deeq,eigr,bec,becdr,fion)
         WRITE( stdout,*) ' out from nlfq'
! 
!     imposing the orthogonality
!     ==========================================================
!
         call calphi(cm,ema0bg,bec,betae,phi)
         WRITE( stdout,*) ' out from calphi'
!     ==========================================================
!
         if(tortho) then
            call ortho  (eigr,c0,phi,lambda,                            &
     &                   bigr,iter,ccc,eps,maxit,delt,bephi,becp)
         else
            call graham(betae,bec,c0)
            WRITE( stdout,*) ' graham  c0 '
         endif
!
!     nlfl needs lambda becdr and bec
!
         if ( tfor .or. tprnfor ) call nlfl(bec,becdr,lambda,fion)
         WRITE( stdout,*) ' out from nlfl'
!
         if(iprsta.ge.3) then
            nnn=min(12,n)
            WRITE( stdout,*) 'from main:'
            do i=1,nnn
               WRITE( stdout,'(12f8.5)') (lambda(i,j)*ccc,j=1,nnn)
            end do
            WRITE( stdout,*)
         endif
!
         if(tpre) then
            call nlfh(bec,dbec,lambda)
            WRITE( stdout,*) ' out from nlfh'
            WRITE( stdout,*) 
            WRITE( stdout,*) ' internal stress tensor:'
            WRITE( stdout,5555) ((stress(i,j),j=1,3),i=1,3)
         endif
 5555    format(1x,f12.5,1x,f12.5,1x,f12.5/                             &
     &          1x,f12.5,1x,f12.5,1x,f12.5/                             &
     &          1x,f12.5,1x,f12.5,1x,f12.5//)
!
         if(tortho) then
            call updatc(ccc,lambda,phi,bephi,becp,bec,c0)
            WRITE( stdout,*) ' out from updatc'
         endif
         call calbec (nvb+1,nsp,eigr,c0,bec)
         if (tpre) call caldbec(1,nsp,eigr,cm)
         WRITE( stdout,*) ' out from calbec'
!     ==============================================================
!     cm now orthogonalized
!     ==============================================================
         if(iprsta.ge.3) call dotcsc(eigr,c0)
!     
         if(thdyn) then
            do i=1,3
               do j=1,3
                  h(i,j)=hold(i,j)+dt2by2/wmass*omega*iforceh(i,j)*     &
     &                  (stress(i,1)*ainv(j,1)+stress(i,2)*ainv(j,2)+   &
     &                   stress(i,3)*ainv(j,3)-press*ainv(j,i))
               end do
            end do
            call invmat( 3, h, ainv, deth )
         endif
!
         if( tfor ) then
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     taus(i,ia,is)=tausm(i,ia,is) +                     &
     &                    iforce(i,ia,is)*dt2by2/pmass(is) *            &
     &                    (fion(1,ia,is)*ainv(i,1) +                    &
     &                     fion(2,ia,is)*ainv(i,2) +                    &
     &                     fion(3,ia,is)*ainv(i,3) )
                  end do
               end do
            end do
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     tau0(i,ia,is)=h(i,1)*taus(1,ia,is)                 &
     &                           +h(i,2)*taus(2,ia,is)                  &
     &                           +h(i,3)*taus(3,ia,is)
                  end do
               end do
            end do
            call phfac(tau0,ei1,ei2,ei3,eigr)
            call calbec (1,nsp,eigr,c0,bec)
            if (tpre) call caldbec(1,nsp,eigr,c0)
         endif
!     
         xnhp0=0.
         xnhpm=0.
         vnhp =0.
         fionm(:,:,:)=0.
         do is=1,nsp
            do ia=1,na(is)
               do i=1,3
                  vels (i,ia,is)=(taus(i,ia,is)-tausm(i,ia,is))/delt
               end do
            end do
         end do
         xnhh0(:,:)=0.
         xnhhm(:,:)=0.
         vnhh (:,:) =0.
         velh (:,:)=(h(:,:)-hold(:,:))/delt
!     
!     ======================================================
!     kinetic energy of the electrons
!     ======================================================
         ALLOCATE( emainv( ngw ) )

         emainv = 1.0d0 / ema0bg
         ftmp = 1.0d0
         if( gstart == 2 ) ftmp = 0.5d0 
         
         ekincm=0.0
         do i=1,n
           ekincm = ekincm + 2.0d0 * &
                    wave_speed2( c0(:,i,1,1), cm(:,i,1,1), emainv, ftmp )
         end do
         ekincm = ekincm * emass / dt2
         
         CALL mp_sum( ekincm )

         DEALLOCATE( emainv )

         xnhe0=0.
         xnhem=0.
         vnhe =0.
!     
         lambdam(:,:)=lambda(:,:)
!     
      else

!======================================================================
!       nbeg = 0, nbeg = 1 or nbeg = 2
!======================================================================

            call readfile_new                                           &
     &     ( 1, ndr,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc,         &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion)
!
         do is=1,nsp
            do ia=1,na(is)
               do i=1,3
                  tau0(i,ia,is)=h(i,1)*taus(1,ia,is)                    &
     &                         +h(i,2)*taus(2,ia,is)                    &
     &                         +h(i,3)*taus(3,ia,is)
               end do
            end do
         end do
!
         if(trane.and.trhor) then
            call prefor(eigr,betae)
            call graham(betae,bec,c0)
            cm(:, 1:n,1,1)=c0(:, 1:n,1,1)
         endif
!
         if(iprsta.gt.2) then
            WRITE( stdout,*) ' read: taus '
            WRITE( stdout,'(3f14.8)')  (((taus(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
            WRITE( stdout,*) ' read: cell parameters h '
            WRITE( stdout,*)  (h(1,j),j=1,3)
            WRITE( stdout,*)  (h(2,j),j=1,3)
            WRITE( stdout,*)  (h(3,j),j=1,3)
         endif
!
         call phfac(tau0,ei1,ei2,ei3,eigr)
         call strucf(ei1,ei2,ei3,sfac)
         call formf(tfirst,eself)
         call calbec (1,nsp,eigr,c0,bec)
         if (tpre) call caldbec(1,nsp,eigr,c0)
!
      end if
!==============================================end of if(nbeg.lt.0)====
!
!     =================================================================
!     restart with new averages and nfi=0
!     =================================================================
      if( nbeg <= 0 ) then
         acc = 0.0d0
         nfi=0
      end if
!
      if( ( .not. tfor ) .and. ( .not. tprnfor ) ) then
         fion (:,:,:) = 0.d0
      end if
!
      if( .not. tpre ) then
         stress (:,:) = 0.d0
      endif
!         
      fccc = 1.0d0
      !
      nomore = nomore + nfi
!
      call cofmass( taus, cdm0 )
!
!======================================================================
!
!           basic loop for molecular dynamics starts here
!
!======================================================================
!
    call stop_clock( 'initialize' ) 

    MAIN_LOOP: DO

      call start_clock( 'total_time' )
!
!     calculation of velocity of nose-hoover variables
!
      if(.not.tsde) fccc=1./(1.+frice)
      if(tnosep)then
         vnhp=2.*(xnhp0-xnhpm)/delt-vnhp
      endif
      if(tnosee)then
         vnhe=2.*(xnhe0-xnhem)/delt-vnhe
         fccc=1./(1.+0.5*delt*vnhe)
      endif
      if(tnoseh) then
         vnhh(:,:)=2.*(xnhh0(:,:)-xnhhm(:,:))/delt-vnhh(:,:)
         velh(:,:)=2.*(h(:,:)-hold(:,:))/delt-velh(:,:)
      endif
! 
      if ( tfor .or. thdyn .or. tfirst ) then 
         call initbox ( tau0, taub, irb )
         call phbox(taub,eigrb)
      endif
!
      if( tfor .or. thdyn ) call phfac(tau0,ei1,ei2,ei3,eigr) 
!
!     strucf calculates the structure factor sfac
!
      call strucf(ei1,ei2,ei3,sfac)
      if (thdyn) call formf(tfirst,eself)
!
      nfi=nfi+1
      tlast=(nfi.eq.nomore)
!
      call rhoofr (nfi,c0,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
!
#ifdef __PARA     
      if(trhow .and. tlast) call write_rho(47,nspin,rhor)
#else
      if(trhow .and. tlast) write(47) ((rhor(i,is),i=1,nnr),is=1,nspin)
#endif
!
!     put core charge (if present) in rhoc(r)
!
      if (nlcc.gt.0) call set_cc(irb,eigrb,rhoc)
!
      call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                 &
     &            ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
      do i=1,3
         do j=1,3
            stress(i,j)=-1.d0/omega*(detot(i,1)*h(j,1)+                 &
     &              detot(i,2)*h(j,2)+detot(i,3)*h(j,3))
         enddo
      enddo
!
      enthal=etot+press*omega
!
!=======================================================================
!
!              verlet algorithm
!
!     loop which updates electronic degrees of freedom
!     cm=c(t+dt) is obtained from cm=c(t-dt) and c0=c(t)
!     the electron mass rises with g**2
!
!=======================================================================
!
      call newd(rhor,irb,eigrb,rhovan,deeq,fion)
      call prefor(eigr,betae)
!
!==== set friction ====
!
      if( tnosee ) then
        verl1 = 2.0d0 * fccc
        verl2 = 1.0d0 - verl1
        verl3 = 1.0d0 * fccc
      else
        verl1=2./(1.+frice) 
        verl2=1.-verl1
        verl3=1./(1.+frice)
      end if
!
!==== start loop ====
!

      ALLOCATE( emadt2( ngw ) )
      ALLOCATE( emaver( ngw ) )
      emadt2 = dt2bye * ema0bg
      emaver = emadt2 * verl3

      do i=1,n,2
         call dforce(bec,deeq,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
         if(tsde) then
            CALL wave_steepest( cm(:, i  , 1, 1), c0(:, i  , 1, 1 ), emadt2, c2 )
            CALL wave_steepest( cm(:, i+1, 1, 1), c0(:, i+1, 1, 1 ), emadt2, c3 )
         else 
            CALL wave_verlet( cm(:, i  , 1, 1), c0(:, i  , 1, 1 ), &
                 verl1, verl2, emaver, c2 )
            CALL wave_verlet( cm(:, i+1, 1, 1), c0(:, i+1, 1, 1 ), &
                 verl1, verl2, emaver, c3 )
         endif
         if ( gstart == 2 ) then
            cm(1,  i,1,1)=cmplx(real(cm(1,  i,1,1)),0.0)
            cm(1,i+1,1,1)=cmplx(real(cm(1,i+1,1,1)),0.0)
         end if
      end do

      ccc = fccc * dt2bye
      DEALLOCATE( emadt2 )
      DEALLOCATE( emaver )
!
!==== end of loop which updates electronic degrees of freedom
!
!     buffer for wavefunctions is unit 21
!
      if(tbuff) rewind 21
!
!----------------------------------------------------------------------
!                 contribution to fion due to lambda
!----------------------------------------------------------------------
!
!     nlfq needs deeq bec
!
      if ( tfor .or. tprnfor ) call nlfq(c0,deeq,eigr,bec,becdr,fion)
!
      if( tfor .or. thdyn ) then
!
! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
!
         lambdap(:,:) = 2.d0*lambda(:,:)-lambdam(:,:)
         lambdam(:,:)=lambda (:,:)
         lambda (:,:)=lambdap(:,:)
      endif
!
!     calphi calculates phi
!     the electron mass rises with g**2
!
      call calphi(c0,ema0bg,bec,betae,phi)
!
!     begin try and error loop (only one step!)
!
!       nlfl and nlfh need: lambda (guessed) becdr
!
      if ( tfor .or. tprnfor ) call nlfl(bec,becdr,lambda,fion)
      if(tpre) then
         if(iprsta.ge.4) then
            if((nfi.eq.0).or.tfirst.or.tlast                            &
     &                   .or.(mod(nfi-1,iprint).eq.0)) then
               WRITE( stdout,*) 
               WRITE( stdout,*) ' internal stress tensor (before nlfh):'
               WRITE( stdout,5555) ((stress(i,j),j=1,3),i=1,3)
            endif
         endif
         call nlfh(bec,dbec,lambda)
         if(iprsta.ge.4) then
            if((nfi.eq.0).or.tfirst.or.tlast                            &
     &                   .or.(mod(nfi-1,iprint).eq.0)) then
               WRITE( stdout,*) 
               WRITE( stdout,*) ' internal stress tensor (after nlfh):'
               WRITE( stdout,5555) ((stress(i,j),j=1,3),i=1,3)
            endif
         endif
         do i=1,3
            do j=1,3
               do is=1,nsp
                  do ia=1,na(is)
                     stress(i,j)=stress(i,j)+pmass(is)/omega*           &
     &                  ((h(i,1)*vels(1,ia,is)+h(i,2)*vels(2,ia,is)+    &
     &                    h(i,3)*vels(3,ia,is))*(h(j,1)*vels(1,ia,is)+  &
     &                    h(j,2)*vels(2,ia,is)+h(j,3)*vels(3,ia,is)))
                  enddo
               enddo
            enddo
         enddo
         if((nfi.eq.0).or.tfirst.or.tlast.or.(mod(nfi-1,iprint).eq.0))  &
     &        then
            WRITE( stdout,*) 
            WRITE( stdout,*) ' internal stress tensor:'
            WRITE( stdout,5555) ((stress(i,j),j=1,3),i=1,3)
         endif
      endif
!
!=======================================================================
!
!              verlet algorithm
!
!     loop which updates cell parameters and ionic degrees of freedom
!     hnew=h(t+dt) is obtained from hold=h(t-dt) and h=h(t)
!     tausp=pos(t+dt) from tausm=pos(t-dt) taus=pos(t) h=h(t)
!
!           guessed displacement of ions
!=======================================================================
!
      hgamma(:,:) = 0.d0
      if(thdyn) then
         verl1=2./(1.+frich)
         verl2=1.-verl1
         verl3=dt2/(1.+frich)
!     
         if (tnoseh) then
            do j=1,3
               do i=1,3
                  hnew(i,j) = h(i,j) +                                  &
     &                 (h(i,j) - hold(i,j) + dt2/wmass*omega*           &
     &           (ainv(j,1)*stress(i,1) + ainv(j,2)*stress(i,2) +       &
     &            ainv(j,3)*stress(i,3) - ainv(j,i)*press) -            &
     &            dt2*vnhh(i,j)*velh(i,j))*iforceh(i,j)
               enddo
            enddo
         else
            do j=1,3
               do i=1,3
                  hnew(i,j) = h(i,j) + ((verl1-1.)*h(i,j)               &
     &                + verl2*hold(i,j)                                 &
     &                + verl3/wmass*omega                               &
     &                *(ainv(j,1)*stress(i,1)+ainv(j,2)*stress(i,2)     &
     &                + ainv(j,3)*stress(i,3)-ainv(j,i)*press))         &
     &                * iforceh(i,j)
               enddo
            enddo
         endif
         !
         velh(:,:) = (hnew(:,:)-hold(:,:))/twodel
         !
         do i=1,3
            do j=1,3
               do k=1,3
                  do l=1,3
                     do m=1,3
                        hgamma(i,j)=hgamma(i,j)+ainv(i,l)*ainv(k,l)*    &
     &                       (velh(m,k)*h(m,j)+h(m,k)*velh(m,j))
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
!
!======================================================================
      if( tfor ) then
!
!==== set friction ====
!
         verl1=2./(1.+fricp)
         verl2=1.-verl1
         verl3=dt2/(1.+fricp)
!
         if(tsdp) then
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     tausp(i,ia,is) = taus(i,ia,is) +                   &
     &                    iforce(i,ia,is)*dt2by2/pmass(is)*             &
     &        (ainv(i,1)*fion(1,ia,is)+ainv(i,2)*fion(2,ia,is)+         &
     &         ainv(i,3)*fion(3,ia,is) ) -                              &
     &                    pmass(is)*(hgamma(i,1)*vels(1,ia,is)+         &
     &         hgamma(i,2)*vels(2,ia,is)+hgamma(i,3)*vels(3,ia,is))
                  end do
               end do
            end do
         else if (tnosep) then
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     fionm(i,ia,is) = (ainv(i,1)*fion(1,ia,is)          &
     &                                +ainv(i,2)*fion(2,ia,is)          &
     &                                +ainv(i,3)*fion(3,ia,is))         &
     &                              - vnhp*vels(i,ia,is)*pmass(is)      &
     &                    - pmass(is)*(hgamma(i,1)*vels(1,ia,is)        &
     &                                +hgamma(i,2)*vels(2,ia,is)        &
     &                                +hgamma(i,3)*vels(3,ia,is))
                     tausp(i,ia,is)=-tausm(i,ia,is)+2.*taus(i,ia,is)+   &
     &                   iforce(i,ia,is)*dt2*fionm(i,ia,is)/pmass(is)
                     velsp(i,ia,is) = velsm(i,ia,is) +                  &
     &                    twodel*fionm(i,ia,is)/pmass(is)
                  end do
               end do
            end do
         else 
            do is=1,nsp
               do ia=1,na(is)
                  do i=1,3
                     tausp(i,ia,is) = verl1*taus(i,ia,is)               &
     &                    + verl2*tausm(i,ia,is)                        &
     &        + verl3/pmass(is)*iforce(i,ia,is) * (ainv(i,1)*fion(1,ia,is)&
     &        + ainv(i,2)*fion(2,ia,is) + ainv(i,3)*fion(3,ia,is))      &
     &        - verl3*iforce(i,ia,is) * (hgamma(i,1)*vels(1,ia,is)      &
     &        + hgamma(i,2)*vels(2,ia,is) + hgamma(i,3)*vels(3,ia,is))
                     velsp(i,ia,is)=velsm(i,ia,is)                      &
     &        - 4.*fricp*vels(i,ia,is)                                  &
     &        + twodel/pmass(is)*iforce(i,ia,is)*(ainv(i,1)*fion(1,ia,is) &
     &        + ainv(i,2)*fion(2,ia,is) + ainv(i,3)*fion(3,ia,is))      &
     &        - twodel*iforce(i,ia,is) * (hgamma(i,1)*vels(1,ia,is)     &
     &        + hgamma(i,2)*vels(2,ia,is) + hgamma(i,3)*vels(3,ia,is))
                  end do
               end do
            end do
         endif
!cc   call cofmass(velsp,cdmvel)
         call cofmass(tausp,cdm)
         do is=1,nsp
            do ia=1,na(is)
               do i=1,3
!cc   velsp(i,ia,is)=velsp(i,ia,is)-cdmvel(i)
                  tausp(i,ia,is)=tausp(i,ia,is)+cdm0(i)-cdm(i)
               enddo
            enddo
         enddo
         do is=1,nsp
            do ia=1,na(is)
               do i=1,3
                  taup(i,ia,is) = hnew(i,1)*tausp(1,ia,is)+             &
     &                 hnew(i,2)*tausp(2,ia,is)+hnew(i,3)*tausp(3,ia,is)
               enddo
            enddo
         enddo
      endif
!     
!---------------------------------------------------------------------------
!              initialization with guessed positions of ions
!---------------------------------------------------------------------------
!
!  if thdyn=true g vectors and pseudopotentials are recalculated for 
!  the new cell parameters
!
      if ( tfor .or. thdyn ) then
         if( thdyn ) then
            hold = h
            h = hnew
            call newinit(ibrav)
            call newnlinit
         else
            hold = h
         endif
!
!       phfac calculates eigr
!
         call phfac(taup,ei1,ei2,ei3,eigr)
!
!       prefor calculates betae
!
         call prefor(eigr,betae)
      end if
!
!---------------------------------------------------------------------------
!                    imposing the orthogonality
!---------------------------------------------------------------------------
!
      if(tortho) then
         call ortho                                                     &
     &         (eigr,cm,phi,lambda,bigr,iter,ccc,eps,maxit,delt,bephi,becp)
      else
         call graham(betae,bec,cm)
         if(iprsta.gt.4) call dotcsc(eigr,cm)
      endif
!
!---------------------------------------------------------------------------
!                   correction to displacement of ions
!---------------------------------------------------------------------------
!
      if(iprsta.ge.3) then
         nnn=min(12,n)
         do i=1,nnn
            WRITE( stdout,*)' main lambda  = ',(lambda(i,k),k=1,nnn)
         end do
         WRITE( stdout,*)
      endif
!
      if(tortho) call updatc(ccc,lambda,phi,bephi,becp,bec,cm)
      call calbec (nvb+1,nsp,eigr,cm,bec)
      if (tpre) call caldbec(1,nsp,eigr,cm)
!
      if(iprsta.ge.3)  call dotcsc(eigr,cm)
!
!---------------------------------------------------------------------------
!                  temperature monitored and controlled
!---------------------------------------------------------------------------
!
      ekinp=0.0
      ekinpr=0.0
!
!     ionic kinetic energy 
!
      if( tfor ) then
         do is=1,nsp
            do ia=1,na(is)
               do i=1,3
                  vels(i,ia,is)=(tausp(i,ia,is)-tausm(i,ia,is))/twodel
               enddo
               do i=1,3
                  do j=1,3
                     do ii=1,3
                        ekinp=ekinp+pmass(is)*                          &
     &                       hold(j,i)*vels(i,ia,is)*                   &
     &                       hold(j,ii)*vels(ii,ia,is)
                     end do
                  end do
               end do
            end do
         end do
      endif
      ekinp=0.5*ekinp
!
!     ionic temperature
!
      call cofmass(vels,cdmvel)
      if( tfor ) then
         do i=1,3
            do j=1,3
               do ii=1,3
                  do is=1,nsp
                     do ia=1,na(is)
                        ekinpr=ekinpr+pmass(is)*hold(j,i)*              &
     &                       (vels(i,ia,is)-cdmvel(i))*                 &
     &                       hold(j,ii)*(vels(ii,ia,is)-cdmvel(ii))
                     end do
                  end do
               end do
            end do
         end do
      endif
      ekinpr=0.5*ekinpr
      tempp=ekinpr*factem/(1.5d0*nat)
!
!     fake electronic kinetic energy
!
      ekinc0 = 0.0d0

      ALLOCATE( emainv( ngw ) )

      emainv = 1.0d0 / ema0bg
      ftmp = 1.0d0
      if( gstart == 2 ) ftmp = 0.5d0

      do i=1,n
         ekinc0 = ekinc0 + 2.0d0 * &
                  wave_speed2( cm(:,i,1,1), c0(:,i,1,1), emainv, ftmp )
      end do

      CALL mp_sum( ekinc0 )

      ekinc0 = ekinc0 * emass / dt2
      ekinc = 0.5 * ( ekinc0 + ekincm )

      DEALLOCATE( emainv )
!
!     fake cell-parameters kinetic energy
!
      ekinh=0.
      if(thdyn) then
         do j=1,3 
            do i=1,3
               ekinh=ekinh+0.5*wmass*velh(i,j)*velh(i,j)
               temphh(i,j)=factem*wmass*velh(i,j)*velh(i,j)
            end do
         end do
      endif
      if(thdiag) then
         temphc=2.*factem*ekinh/3.
      else
         temphc=2.*factem*ekinh/9.
      endif
!
!     udating nose-hoover friction variables
!
      if(tnosep)then
         xnhpp=2.*xnhp0-xnhpm+2.*(dt2/qnp)*(ekinpr-gkbt/2.)
         vnhp =(xnhpp-xnhpm)/twodel
      endif
      if(tnosee)then
         xnhep=2.*xnhe0-xnhem+2.*(dt2/qne)*(ekinc-ekincw)
         vnhe =(xnhep-xnhem)/twodel
      endif
      if(tnoseh)then
         do j=1,3
            do i=1,3
               xnhhp(i,j)=2.*xnhh0(i,j)-xnhhm(i,j)+                     &
     &              (dt2/qnh)/factem*(temphh(i,j)-temph)
               vnhh(i,j) =(xnhhp(i,j)-xnhhm(i,j))/twodel
            end do
         end do
      endif
!
! warning! thdyn and tcp/tcap are not compatible yet!!!
!
      if(tcp.or.tcap.and.tfor.and.(.not.thdyn)) then
         if(tempp.gt.temp1.or.tempp.lt.temp2.and.tempp.ne.0.d0) then
            if(.not.tcap) then
               alfap=.5d0*sqrt(tempw/tempp)
               do is=1,nsp
                  do ia=1,na(is)
                     do i=1,3
                        taup(i,ia,is) = tau0(i,ia,is) +                 &
     &                       alfap*(taup(i,ia,is)-taum(i,ia,is)) +      &
     &                      dt2by2/pmass(is)*fion(i,ia,is)*iforce(i,ia,is)
                     end do
                  end do
               end do
            else
               do i=1,3
                  qr(i)=0.d0
                  do is=1,nsp
                     do ia=1,na(is)
                        alfar=gausp/sqrt(pmass(is))*cos(2.d0*pi*randy())&
     &                       *sqrt(-2.d0*log(randy()))
                        taup(i,ia,is)=alfar
                        qr(i)=qr(i)+alfar
                     end do
                  end do
                  qr(i)=qr(i)/nat
               end do
               do is=1,nsp
                  do ia=1,na(is)
                     do i=1,3
                        alfar=taup(i,ia,is)-qr(i)
                        taup(i,ia,is)=tau0(i,ia,is)+iforce(i,ia,is)*     &
     &                             (alfar+dt2by2/pmass(is)*fion(i,ia,is))
                     end do
                  end do
               end do
            end if
         end if
      end if
!
      if(mod(nfi-1,iprint).eq.0 .or. (nfi.eq.(nomore))) then
         call eigs(nspin,nx,nupdwn,iupdwn,f,lambda)
         WRITE( stdout,*)
      endif
!
      epot=eht+epseu+exc
!
      acc(1)=acc(1)+ekinc
      acc(2)=acc(2)+ekin
      acc(3)=acc(3)+epot
      acc(4)=acc(4)+etot
      acc(5)=acc(5)+tempp
!
      econs=ekinp+ekinh+enthal
      econt=econs+ekinc
      if(tnosep)then
         econt=econt+0.5*qnp*vnhp*vnhp+     gkbt*xnhp0
      endif
      if(tnosee)then
         econt=econt+0.5*qne*vnhe*vnhe+2.*ekincw*xnhe0
      endif
      if(tnoseh)then
         do i=1,3
            if(thdiag) then
               econt=econt+0.5*qnh*vnhh(i,i)*vnhh(i,i)+                 &
     &                temph/factem*xnhh0(i,i)
            else
               do j=1,3
                  econt=econt+0.5*qnh*vnhh(i,j)*vnhh(i,j)+              &
     &                 temph/factem*xnhh0(i,j)
               enddo
            endif
         enddo
      endif
!
      if(mod(nfi-1,iprint).eq.0.or.tfirst)  then
         WRITE( stdout,*)
         WRITE( stdout,1947)
      end if
!
      tps=nfi*delt*2.4189d-5
      WRITE( stdout,1948) nfi,ekinc,int(temphc),int(tempp),enthal,econs,      &
     &              econt,                                              &
     &              vnhh(3,3),xnhh0(3,3),vnhp,xnhp0
      write(8,2948) tps,ekinc,int(temphc),int(tempp),enthal,econs,      &
     &              econt,                                              &
     &              vnhh(3,3),xnhh0(3,3),vnhp,xnhp0
!     c              frice,frich,fricp
! 
 1947 format(2x,'nfi',4x,'ekinc',2x,'temph',1x,'tempp',6x,'enthal',     &
     &       7x,'econs',7x,'econt',4x,'vnhh',3x,'xnhh0',4x,'vnhp',      &
     &       3x,'xnhp0')
!cc     f       7x,'econs',7x,'econt',3x,'frice',2x,'frich',2x,'fricp')
 1948 format(i5,1x,f8.5,1x,i6,1x,i5,3(1x,f11.5),4(1x,f7.4))
 2948 format(f8.5,1x,f8.5,1x,i6,1x,i5,3(1x,f11.5),4(1x,f7.4))
!
      if( tfor ) then
        if ( ionode ) then
          write(77,3340) ((h(i,j),i=1,3),j=1,3)
          write(77,'(3f12.8)') (((taus(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
          write(78,'(3f12.8)') (((fion(i,ia,is),i=1,3),ia=1,na(is)),is=1,nsp)
          write(79,3340) ((stress(i,j),i=1,3),j=1,3)
#if defined (__ORIGIN) || defined (__T3E)
          call flush(77)
          call flush(78)
#elif defined (__AIX) || defined (__ABSOFT)
          call flush_(77)
          call flush_(78)
#endif
 3340     format(9(1x,f9.5))
        endif
!
!     new variables for next step
!
         do is=1,nsp
            do ia=1,na(is)
               do i=1,3
                  tausm(i,ia,is)=taus(i,ia,is)
                  taus(i,ia,is)=tausp(i,ia,is)
                  taum(i,ia,is)=tau0(i,ia,is)
                  tau0(i,ia,is)=taup(i,ia,is)
                  velsm(i,ia,is)=vels(i,ia,is)
                  vels(i,ia,is)=velsp(i,ia,is)
               end do
            end do
         end do
         if(tnosep) then
            xnhpm = xnhp0
            xnhp0 = xnhpp
         endif
         if(tnosee) then
            xnhem = xnhe0
            xnhe0 = xnhep
         endif
         if(tnoseh) then
            xnhhm(:,:) = xnhh0(:,:)
            xnhh0(:,:) = xnhhp(:,:)
         endif
      end if
!
      if(thdyn)then
         do i=1,ngw
            ema0bg(i)=1./max(1.d0,tpiba2*ggp(i)/emaec) 
         enddo
      endif
!
      ekincm=ekinc0
!  
!     cm=c(t+dt) c0=c(t)
!
      call DSWAP(2*ngw*n,c0,1,cm,1)
!
!     now:  cm=c(t) c0=c(t+dt)
!
      if (tfirst) then
         epre = etot
         enow = etot
      endif
!
      tfirst=.false.
!
!     write on file ndw each iprint
!
      if( ( mod( nfi, iprint ) == 0 ) .and. ( nfi < nomore ) ) then

         call writefile_new                                         &
     &     ( ndw,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc,               &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion)

      endif
!
!     =====================================================================
!     automatic adapt of friction using grease and twall
!     =====================================================================

      epre = enow
      enow = etot
      tbump=.false.
      if (enow .gt. (epre+0.00002)) then
         if (tsde) then
            WRITE( stdout,'(''etot rising with tsde - program stopped'')')
            delt = delt - 1.
            WRITE( stdout,'(''new delt = delt - 1. = '',f12.6)') delt       
            WRITE( stdout,*)
            if (delt .le. 0.) stop
            if(nbeg.lt.0) goto 666
         endif
         if (frice .lt. grease) then
            savee = frice
            savep = fricp
         endif
         if (twall) then
            tbump = .true.
            frice = 1./grease
            fricp = 1./greasp
            frich = 1./greash
         endif
      else
         if (tbump) then
            tbump = .false.
            frice = savee
            fricp = savep
            frich = saveh
         endif
      endif
      frice = frice * grease
      fricp = fricp * greasp
      frich = frich * greash
!     =====================================================
      call stop_clock( 'total_time' )
!
      delta_etot = ABS( epre - enow )
      tstop = check_stop_now()
      tconv = .FALSE.
      IF( tconvthrs%active ) THEN
        tconv = ( delta_etot < tconvthrs%derho ) .AND. ( ekinc < tconvthrs%ekin )
      END IF
      IF( tconv ) THEN
        IF( ionode ) THEN
          WRITE( stdout,fmt= &
            "(/,3X,'MAIN:',10X,'EKINC   (thr)',10X,'DETOT   (thr)',7X,'MAXFORCE   (thr)')" )
          WRITE( stdout,fmt="(3X,'MAIN: ',3(D14.6,1X,D8.1))" ) &
            ekinc, tconvthrs%ekin, delta_etot, tconvthrs%derho, 0.0d0, tconvthrs%force
          WRITE( stdout,fmt="(3X,'MAIN: convergence achieved for system relaxation')")
        END IF
      END IF

      tstop = tstop .OR. tconv

      if( (nfi >= nomore) .OR. tstop ) EXIT MAIN_LOOP

    END DO MAIN_LOOP

!
!=============================end of main loop of molecular dynamics====
!

      ! 
      !  Here copy relevant physical quantities into the output arrays/variables
      !

      etot_out = etot
      isa = 0
      do is = 1, nsp
        do ia = 1, na(is)
          isa = isa + 1
          ipos = ind_srt( isa )
          tau( 1, ipos ) = tau0( 1, ia, is )
          tau( 2, ipos ) = tau0( 2, ia, is )
          tau( 3, ipos ) = tau0( 3, ia, is )
          ! ftmp = ainv(1,1)*fion(1,ia,is)+ainv(1,2)*fion(2,ia,is)+ainv(1,3)*fion(3,ia,is)
          fion_out( 1, ipos ) = fion( 1, ia, is )
          ! ftmp = ainv(2,1)*fion(1,ia,is)+ainv(2,2)*fion(2,ia,is)+ainv(2,3)*fion(3,ia,is)
          fion_out( 2, ipos ) = fion( 2, ia, is )
          ! ftmp = ainv(3,1)*fion(1,ia,is)+ainv(3,2)*fion(2,ia,is)+ainv(3,3)*fion(3,ia,is)
          fion_out( 3, ipos ) = fion( 3, ia, is )
        end do
      end do

      !  Calculate statistics

      do i=1,nacc
         acc(i)=acc(i)/dble(nfi)
      end do
!
      WRITE( stdout,1949)
 1949 format(//'              averaged quantities :',/,                 &
     &       9x,'ekinc',10x,'ekin',10x,'epot',10x,'etot',5x,'tempp')
      WRITE( stdout,1950) (acc(i),i=1,nacc)
 1950 format(4f14.5,f10.1)
!
      call print_clock( 'initialize' )
      call print_clock( 'total_time' )
      call print_clock( 'formf' )
      call print_clock( 'rhoofr' )
      call print_clock( 'vofrho' )
      call print_clock( 'dforce' )
      call print_clock( 'calphi' )
      call print_clock( 'ortho' )
      call print_clock( 'updatc' )
      call print_clock( 'graham' )
      call print_clock( 'newd' )
      call print_clock( 'calbec' )
      call print_clock( 'prefor' )
      call print_clock( 'strucf' )
      call print_clock( 'nlfl' )
      call print_clock( 'nlfq' )
      call print_clock( 'set_cc' )
      call print_clock( 'rhov' )
      call print_clock( 'nlsm1' )
      call print_clock( 'nlsm2' )
      call print_clock( 'forcecc' )
      call print_clock( 'fft' )
      call print_clock( 'ffts' )
      call print_clock( 'fftw' )
      call print_clock( 'fftb' )
      call print_clock( 'rsg' )
      call print_clock( 'setfftpara' )
      call print_clock( 'reduce' )

!
         call writefile_new                                         &
     &     ( ndw,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc,               &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion)

!
      if(iprsta.gt.1) then
         WRITE( stdout,*)
         WRITE( stdout,3370)'    lambda   n = ',n
         do i=1,n
            WRITE( stdout,3380) (lambda(i,j),j=1,n)
         end do
3370     format(26x,a,i4)
3380     format(9f8.4)      
      endif
!
      if( tfor .or. tprnfor ) then
         WRITE( stdout,1970) ibrav, alat
         WRITE( stdout,1971)
         do i=1,3
            WRITE( stdout,1972) (h(i,j),j=1,3)
         enddo
         WRITE( stdout,1973)
         do is=1,nsp
            do ia=1,na(is)
               WRITE( stdout,1974) is,ia,(tau0(i,ia,is),i=1,3),               &
     &            ((ainv(j,1)*fion(1,ia,is)+ainv(j,2)*fion(2,ia,is)+    &
     &              ainv(j,3)*fion(3,ia,is)),j=1,3)
            end do
         end do
         WRITE( stdout,1975)
         do is=1,nsp
            do ia=1,na(is)
               WRITE( stdout,1976) is,ia,(taus(i,ia,is),i=1,3)
            end do
         end do
      endif
      conv_elec = .TRUE.



 1970 format(1x,'ibrav :',i4,'  alat : ',f10.4,/)
 1971 format(1x,'lattice vectors',/)
 1972 format(1x,3f10.4)
 1973 format(/1x,'Cartesian coordinates (a.u.)              forces' &
     &       /1x,'species',' atom #', &
     &           '   x         y         z      ', &
     &           '   fx        fy        fz'/)
 1974 format(1x,2i5,3f10.4,2x,3f10.4)
 1975 format(/1x,'Scaled coordinates '/1x,'species',' atom #')
 1976 format(1x,2i5,3f10.4)
      WRITE( stdout,1977) 

!
 600  continue
!
      call memory
!      
 1977 format(5x,//'====================== end cprvan ',                 &
     &            '======================',//)

      IF( ALLOCATED( ei1 ) ) DEALLOCATE( ei1 )
      IF( ALLOCATED( ei2 ) ) DEALLOCATE( ei2 )
      IF( ALLOCATED( ei3 ) ) DEALLOCATE( ei3 )
      IF( ALLOCATED( eigr ) ) DEALLOCATE( eigr )
      IF( ALLOCATED( sfac ) ) DEALLOCATE( sfac )
      IF( ALLOCATED( eigrb ) ) DEALLOCATE( eigrb )
      IF( ALLOCATED( rhor ) ) DEALLOCATE( rhor )
      IF( ALLOCATED( rhos ) ) DEALLOCATE( rhos )
      IF( ALLOCATED( rhog ) ) DEALLOCATE( rhog )
      IF( ALLOCATED( rhoc ) ) DEALLOCATE( rhoc )
      IF( ALLOCATED( betae ) ) DEALLOCATE( betae )
      IF( ALLOCATED( bec ) ) DEALLOCATE( bec )
      IF( ALLOCATED( becdr ) ) DEALLOCATE( becdr )
      IF( ALLOCATED( bephi ) ) DEALLOCATE( bephi )
      IF( ALLOCATED( becp ) ) DEALLOCATE( becp )
      IF( ALLOCATED( rhovan ) ) DEALLOCATE( rhovan )
      IF( ALLOCATED( deeq ) ) DEALLOCATE( deeq )
      IF( ALLOCATED( ema0bg ) ) DEALLOCATE( ema0bg )
      IF( ALLOCATED( lambda ) ) DEALLOCATE( lambda )
      IF( ALLOCATED( lambdam ) ) DEALLOCATE( lambdam )
      IF( ALLOCATED( lambdap ) ) DEALLOCATE( lambdap )
      IF( ALLOCATED( c2 ) ) DEALLOCATE( c2 )
      IF( ALLOCATED( c3 ) ) DEALLOCATE( c3 )

      CALL deallocate_elct()
      CALL deallocate_core()
      CALL deallocate_cvan()
      CALL deallocate_gvec()
      CALL deallocate_pseu()
      CALL deallocate_qgb_mod()
      CALL deallocate_qradb_mod()
      CALL deallocate_work()
      CALL deallocate_work_box()
      CALL deallocate_derho()
      CALL deallocate_dqgb_mod()
      CALL deallocate_dpseu()
      CALL deallocate_cdvan()
      CALL deallocate_dqrad_mod()
      CALL deallocate_betax()
      CALL deallocate_para_mod()
      CALL deallocate_wavefunctions()

      if( ionode ) then
        CLOSE( 8 )
        CLOSE( 77 )
        CLOSE( 78 )
        CLOSE( 79 )
      end if

!
      return
      end subroutine
