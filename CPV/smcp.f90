!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Version OSCP.01_050404
!
! main subroutine for SMD by Yosuke Kanai
!

#include "../include/f_defs.h"

!
!=======================================================================
!
subroutine smdmain( tau, fion_out, etot_out )
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
  use control_flags, only: iprint, isave, thdyn, tpre, tbuff, iprsta, trhor, &
       tfor, tvlocw, trhow, taurdr, tprnfor, tsdc
  use control_flags, only: ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
       tnosep, trane, tranp, tsdp, tcp, tcap, ampre, amprp, tnoseh
  !
  use atom, only: nlcc
  use core, only: deallocate_core
  use uspp_param, only: nhm
  use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq
  use uspp, only: deallocate_uspp
  use cvan, only: nvb
  use energies, only: eht, epseu, exc, etot, eself, enl, ekin, esr
  use elct, only: nx, n, ispin, f, nspin, nel, iupdwn, nupdwn
  use elct, only: deallocate_elct
  use gvec, only: tpiba2, ng
  use gvec, only: deallocate_gvec
  use gvecs, only: ngs
  use gvecb, only: ngb
  use gvecw, only: ngw
  use reciprocal_vectors, only: ng0 => gstart
  use ions_base, only: na, nat, pmass, nas => nax, nsp, rcmax
  use ions_base, only: ind_srt, ions_vel, ions_cofmass, ions_kinene, ions_temp
  use ions_base, only: ions_thermal_stress, ions_vrescal
  use grid_dimensions, only: nnr => nnrx, nr1, nr2, nr3
  use cell_base, only: ainv, a1, a2, a3, r_to_s
  use cell_base, only: omega, alat, frich, greash, press
  use cell_base, only: h, hold, deth, wmass, s_to_r, iforceh, cell_force
  use smooth_grid_dimensions, only: nnrsx, nr1s, nr2s, nr3s
  use smallbox_grid_dimensions, only: nnrb => nnrbx, nr1b, nr2b, nr3b
  use pseu, only: vps, rhops
  use pseu, only: deallocate_pseu
  use work
  use work_box, only: qv, deallocate_work_box
  use io_global, ONLY: io_global_start, stdout, ionode
  use mp_global, ONLY: mp_global_start
  use mp, ONLY: mp_sum
  use para_mod
  use dener
  use derho
  use dpseu
  use cdvan
  use stre
  use gvecw, only: ggp, agg => ecutz, sgg => ecsig, e0gg => ecfix
  !use restartsm
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
  use input_parameters, only: outdir, ene_ini => smd_ene_ini, &
       ene_fin => smd_ene_fin, num_of_images
  use wave_base, only: wave_steepest, wave_verlet
  use wave_base, only: wave_speed2
  USE control_flags, ONLY : conv_elec, tconvthrs
  USE check_stop, ONLY : check_stop_now
  USE smd_variables
  USE smd_rep
  USE smd_ene
  USE cpr_subroutines
  use ions_positions, only: tau0, velsp
  use ions_positions, only: ions_hmove, ions_move
  use ions_nose, only: gkbt, qnp, ions_nosevel, ions_noseupd, tempw
  USE cell_base, ONLY: cell_kinene, cell_move, cell_gamma, cell_hmove
  USE cell_nose, ONLY: cell_nosevel, cell_noseupd


  !
  !
  implicit none
  !
  ! output variables
  !
  real(kind=8) :: tau(3,nat,0:num_of_images-1)
  real(kind=8) :: fion_out(3,nat,0:num_of_images-1)
  real(kind=8) :: etot_out(0:num_of_images-1)
  !
  !
  ! control variables
  !
  logical tbump
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
  real(kind=8) taub(3,natx)
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
  real(kind=8), allocatable:: becdr(:,:,:)
  real(kind=8), allocatable:: bephi(:,:), becp(:,:)
  !
  !  mass preconditioning
  !
  real(kind=8), allocatable:: ema0bg(:)
  real(kind=8), allocatable:: emadt2(:)
  real(kind=8), allocatable:: emaver(:)
  real(kind=8), allocatable:: emainv(:)
  real(kind=8)  emaec

  !
  !  iforce, 
  ! 
  integer :: iforce(3,natx)
  !
  ! for variable cell dynamics: scaled tau
  !
  real(kind=8) f_(nbndxx)
  real(kind=8) ispin_(nbndxx)
  !
  integer maxit
  !
  !
  ! work variables
  !
  real(kind=8) celldm(6), ecut, ecutw
  complex(kind=8), allocatable:: c2(:), c3(:)
  complex(kind=8)  speed
  real(kind=8)                                                      &
       &       tempp(smx), xnhe0(smx), vnhp(smx), xnhp0(smx), xnhpm(smx), &
       &       fccc(smx), xnhem(smx), vnhe(smx),  &
       &       epot(smx), xnhpp(smx), xnhep(smx), epre(smx), enow(smx),   &
       &       econs(smx), econt(smx),                                    &
       &       ccc(smx)
  real(kind=8) temps(nsx) 
  real(kind=8) verl1, verl2, verl3, anor, saveh, tps,        &
       &       eps, qne, emass, delt, bigr, dt2,              &
       &       dt2by2, twodel, gausp, dt2bye, dt2hbe, fricp, greasp,&
       &        frice, grease, savee, savep
  real(kind=8) ekinc0(smx), ekinp(smx), ekinpr(smx), ekincm(smx),   &
       &       ekinc(smx), pre_ekinc(smx), enthal(smx),ekincw
  integer is, nacc, ia, j, iter, nfi, i, isa, ipos
  !
  !
  ! work variables, 2
  !
  real(kind=8) fcell(3,3), hnew(3,3),velh(3,3),hgamma(3,3)
  real(kind=8) cdm(3)
  real(kind=8) qr(3)
  real(kind=8) xnhh0(3,3),xnhhm(3,3),xnhhp(3,3),vnhh(3,3),temphh(3,3)
  real(kind=8) thstress(3,3) 
  !
  integer ibrav, k, ii, l, m
  real(kind=8) ekinh, alfar, temphc, alfap, tolp,    &
       &     factp, temp1, temp2, temph, qnh, randy
  real(kind=8) ftmp

  character(len=256) :: filename
  character(len=256) :: dirname
  integer :: strlen, dirlen
  !
  character :: unitsmw*2, unitnum*2  
  !
  !
  !     SMD 
  !
  real(kind=8) :: t_arc_pre,t_arc_now,t_arc_tot
  integer :: sm_k,sm_file,sm_ndr,sm_ndw,unico,unifo,unist
  integer :: smpm,con_ite 

  real(kind=8) :: workvec(3,natx,nsx) 
  real(kind=8), allocatable :: deviation(:)
  real(kind=8), allocatable :: maxforce(:)
  real(kind=8), allocatable :: arc_now(:)
  real(kind=8), allocatable :: arc_tot(:)
  real(kind=8), allocatable :: arc_pre(:)
  real(kind=8), allocatable :: paraforce(:)
  real(kind=8), allocatable :: err_const(:)
  integer, allocatable :: pvvcheck(:)
  !
  type(ptr), allocatable :: p_tau0(:)
  type(ptr), allocatable :: p_taup(:)
  type(ptr), allocatable :: p_tan(:)
  !
  !
  !     CP loop starts here
  !
  call start_clock( 'initialize' )
  !
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
  tps     = 0.0d0

  !
  !     ==================================================================
  !     read input from standard input (unit 5)
  !     ==================================================================

  call iosys( nbeg , ndr , ndw , nomore , iprint, isave                  &
       & , delt , emass , emaec  , tsde , frice , grease          &
       & , tortho , eps , maxit , trane , ampre , tranp , amprp          &
       & , tfor , tsdp , fricp , greasp , tcp , tcap , tolp , trhor , trhow , tvlocw &
       & , tnosep , qnp , tempw , tnosee , qne , ekincw                 &
       & , tpre , thdyn , thdiag , iforceh , wmass , frich , greash , press   &
       & , tnoseh , qnh , temph , celldm , ibrav , tau0, ecutw , ecut , iforce &
       & , nat , nsp , na , pmass , rcmax , f_ , nel , nspin , nupdwn  &
       & , iupdwn , n , nx , nr1 , nr2 , nr3 , omega , alat , a1 , a2 , a3  &
       & , nr1b , nr2b , nr3b , nr1s , nr2s , nr3s , agg , sgg , e0gg &
       & , psfile , pseudo_dir, iprsta, ispin_ &
       & , sm_p, smcp, smlm, smopt, linr, polm, kwnp, codfreq, forfreq, smwfreq &
       & , tol, lmfreq, maxlm )

  !
  !
  allocate( f( nx ) )
  f( :   ) = 0.0d0
  f( 1:n ) = f_( 1:n )
  allocate( ispin( nx ) )
  ispin( :   ) = 0.0d0
  ispin( 1:n ) = ispin_( 1:n )
  !
  !
  !     ==================================================================
  !
  write(stdout,*) '                                                 '
  write(stdout,*) ' ================================================'
  write(stdout,*) '                                                 '
  write(stdout,*) ' CPSMD : replica = 0 .. ', sm_p
  write(stdout,*) '                                                 '
  write(stdout,*) '                                                 '
  IF(smopt) THEN
     write(stdout,*) ' CP optimizations for replicas 1, and 2  '
     write(stdout,*) '                                                 '
  ENDIF
  write(stdout,*) ' String Method : '
  write(stdout,*) '                                                 '
  write(stdout,*) '                                                 '
  write(stdout,*) '      SMLM    = ', smlm
  write(stdout,*) '      CPMD    = ', smcp
  write(stdout,*) '      LINR    = ', linr 
  write(stdout,*) '      POLM    = ', polm, 'Pt = ', kwnp
  write(stdout,*) '      STCD    = ', stcd
  write(stdout,*) '      SMfile  = ', smwout
  write(stdout,*) '      TFOR  = ', tfor
  write(stdout,*) '                                                 '
  IF(smlm) THEN
     write(stdout,*) '      CONST. TOL  = ', tol
     write(stdout,*) '      Frequency   = ', lmfreq
  ENDIF
  write(stdout,*) '                                                 '
  write(stdout,*) ' ================================================'
  !
  !
  !
  !     ===================================================================
  !
  !     general variables
  !
  tfirst = .true.
  tlast  = .false.
  nacc = 5
  !
  !
  !
  gausp = delt * sqrt(tempw/factem)
  twodel = 2.d0 * delt
  dt2 = delt * delt
  dt2by2 = .5d0 * dt2
  dt2bye = dt2/emass
  dt2hbe = dt2by2/emass
  smpm = sm_p -1
  !
  !
  IF(tbuff) THEN
     write(stdout,*) "TBUFF set to .FALSE., Option not implemented with SMD"
     tbuff = .FALSE.
  ENDIF
  !
  !
  !     Opening files
  !
  !
  unico = 41
  unifo = 42
  unist = 43
  !
  if (ionode) then

     dirlen = index(outdir,' ') - 1
     filename = 'fort.8'
     if( dirlen >= 1 ) then
        filename = outdir(1:dirlen) // '/' // filename
     end if
     strlen  = index(filename,' ') - 1
     WRITE( stdout, * ) ' UNIT8 = ', filename
     OPEN(unit=8, file=filename(1:strlen), status='unknown')

     filename = 'cod.out'
     if( dirlen >= 1 ) then
        filename = outdir(1:dirlen) // '/' // filename
     end if
     strlen  = index(filename,' ') - 1
     OPEN(unit=unico, file=filename(1:strlen), status='unknown')

     filename = 'for.out'
     if( dirlen >= 1 ) then
        filename = outdir(1:dirlen) // '/' // filename
     end if
     strlen  = index(filename,' ') - 1
     OPEN(unit=unifo, file=filename(1:strlen), status='unknown')

     filename = 'str.out'
     if( dirlen >= 1 ) then
        filename = outdir(1:dirlen) // '/' // filename
     end if
     strlen  = index(filename,' ') - 1
     OPEN(unit=unist, file=filename(1:strlen), status='unknown')


     WRITE(unitsmw, '(i2)') smwout

     DO sm_k=1,smpm
        sm_file = smwout + sm_k

        unitnum = "00"

        IF(sm_k < 10) THEN
           write(unitnum(2:2), '(i1)') sm_k
        ELSE
           write(unitnum,'(i2)') sm_k
        ENDIF

        filename = 'rep.'
        if( dirlen >= 1 ) then
           filename = outdir(1:dirlen) // '/' // filename
        end if
        strlen  = index(filename,' ') - 1

        OPEN(unit=sm_file,file=filename(1:strlen)//unitnum, status='unknown')
     ENDDO

  end if

  !
  !     ==================================================================
  !     initialize g-vectors, fft grids
  !     ==================================================================
  !
  !
  !
  !      ... taus is calculated here.
  !

  call sminit( ibrav, celldm, ecut, ecutw, ndr, nbeg, tfirst, delt, tps, iforce )

  !
  !

  call r_to_s(rep(0)%tau0,rep(0)%taus, na, nsp, ainv)
  call r_to_s(rep(sm_p)%tau0,rep(sm_p)%taus, na, nsp, ainv)

  !
  !
  WRITE( stdout,*) ' out from init'
  !
  !     more initialization requiring atomic positions
  !
  nas = MAXVAL( na( 1 : nsp ) )

  DO sm_k=1,smpm
     sm_file = smwout + sm_k
     if( iprsta > 1 ) then
        if(ionode) WRITE( sm_file,*) ' tau0 '
        if(ionode) WRITE( sm_file,'(3f14.8)') ((rep(sm_k)%tau0(i,ia),i=1,3),ia=1,SUM(na(1:nsp)))
     endif
  ENDDO

  !
  !     ==================================================================
  !     allocate and initialize nonlocal potentials
  !     ==================================================================

  call nlinit
  !

  WRITE( stdout,*) ' out from nlinit'

  !     ==================================================================
  !     allocation of all arrays not already allocated in init and nlinit
  !     ==================================================================
  !
  !
  WRITE( stdout,*) ' Allocation begun, smpm = ', smpm
  !
#if defined (__ORIGIN) || defined (__T3E)
  call flush(stdout)
#elif defined (__AIX) || defined (__ABSOFT)
  call flush_(stdout)
#endif
  !
  !
  DO sm_k=1,smpm
     allocate(rep_el(sm_k)%c0(ngw,nx))
     allocate(rep_el(sm_k)%cm(ngw,nx))
     allocate(rep_el(sm_k)%phi(ngw,nx))
     allocate(rep_el(sm_k)%lambda(nx,nx))
     allocate(rep_el(sm_k)%lambdam(nx,nx))
     allocate(rep_el(sm_k)%lambdap(nx,nx))
     allocate(rep_el(sm_k)%bec  (nhsa,n))
     allocate(rep_el(sm_k)%rhovan(nhm*(nhm+1)/2,nat,nspin))
  ENDDO
  !
  WRITE( stdout,*) " Allocation for W.F. s : successful " 

  !
  !
  allocate(wrk2(ngw,max(nas,n)))
  allocate(eigr(ngw,nas,nsp))
  allocate(eigrb(ngb,nas,nsp))
  allocate(sfac(ngs,nsp))
  allocate(rhops(ngs,nsp))
  allocate(vps(ngs,nsp))
  allocate(rhor(nnr,nspin))
  allocate(rhos(nnrsx,nspin))
  allocate(rhog(ng,nspin))
  if ( ANY( nlcc ) ) allocate(rhoc(nnr))
  allocate(wrk1(nnr))
  allocate(qv(nnrb))
  allocate(c2(ngw))
  allocate(c3(ngw))
  allocate(ema0bg(ngw))
  allocate(ei1(-nr1:nr1,nas,nsp))
  allocate(ei2(-nr2:nr2,nas,nsp))
  allocate(ei3(-nr3:nr3,nas,nsp))
  allocate(betae(ngw,nhsa))
  allocate(becdr(nhsa,n,3))
  allocate(bephi(nhsa,n))
  allocate(becp (nhsa,n))
  allocate(deeq(nhm,nhm,nat,nspin))
  allocate(dbec (nhsa,n,3,3))
  allocate(dvps(ngs,nsp))
  allocate(drhops(ngs,nsp))
  allocate(drhog(ng,nspin,3,3))
  allocate(drhor(nnr,nspin,3,3))
  allocate(drhovan(nhm*(nhm+1)/2,nat,nspin,3,3))
  !
  WRITE( stdout,*) ' Allocation for CP core : successful '
  !
  allocate(etot_ar(0:sm_p))
  allocate(ekin_ar(0:sm_p))
  allocate(eht_ar(0:sm_p))
  allocate(epseu_ar(0:sm_p))
  allocate(exc_ar(0:sm_p))
  allocate(esr_ar(0:sm_p))
  !
  allocate(deviation(smpm))
  allocate(maxforce(smpm))
  allocate(arc_now(0:sm_p))
  allocate(arc_tot(0:sm_p))
  allocate(arc_pre(0:sm_p))
  allocate(paraforce(smpm))
  allocate(err_const(sm_p))
  allocate(pvvcheck(smpm))
  !
  WRITE( stdout,*) ' Allocation for SM variables : successful '
  !
  allocate(p_tau0(0:sm_p))
  allocate(p_taup(0:sm_p))
  allocate(p_tan(0:sm_p))
  !
  WRITE( stdout,*) ' Allocation for pointers : successful '
  !
#ifdef __PARA
  allocate(aux(nnr))
#endif
  !
  !
  con_ite = 0
  deeq(:,:,:,:) = 0.d0
  deviation(:) = 0.d0
  maxforce(:) = 0.d0
  pre_ekinc(:) = 0.d0
  paraforce(:) = 0.d0
  err_const(:) = 0.d0
  ekinc(:) = 0.d0
  arc_now(:) = 0.d0
  arc_tot(:) = 0.d0
  arc_pre(:) = 0.d0
  !
  !
  WRITE( stdout,*) ' Allocation ended '
  !
#if defined (__ORIGIN) || defined (__T3E)
  call flush(stdout)
#elif defined (__AIX) || defined (__ABSOFT)
  call flush_(stdout)
#endif

  !
666 continue
  !
  !
  DO sm_k=0,sm_p
     p_tau0(sm_k)%d3 => rep(sm_k)%tau0
     IF(tfor) THEN
        p_taup(sm_k)%d3 => rep(sm_k)%taup
     ELSE
        p_taup(sm_k)%d3 => rep(sm_k)%tau0
     ENDIF
     p_tan(sm_k)%d3 => rep(sm_k)%tan
  ENDDO
  !
  !
  temp1=tempw+tolp
  temp2=tempw-tolp
  gkbt = 3.*nat*tempw/factem
  press = press*factp
  !     ==========================================================

  etot_ar(0) = ene_ini
  etot_ar(sm_p) = ene_fin

  DO sm_k=0,sm_p
     rep(sm_k)%tausm=rep(sm_k)%taus
     rep(sm_k)%tausp=0.
     rep(sm_k)%taum=rep(sm_k)%tau0
     rep(sm_k)%taup=0.
     rep(sm_k)%vels  = 0.0d0
     rep(sm_k)%velsm = 0.0d0
  ENDDO
  velsp = 0.
  !
  hnew=h
  !

  DO sm_k=1,smpm
     rep_el(sm_k)%lambda = 0.d0
     rep_el(sm_k)%cm = (0.d0, 0.d0)
     rep_el(sm_k)%c0 = (0.d0, 0.d0)
  ENDDO
  !
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
  !
  ! ... calculating tangent for force transformation.
  !

  IF(smlm) call TANGENT(p_tau0,p_tan)

  !
  !
  !
  INI_REP_LOOP : DO sm_k=1,smpm  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
     !

     sm_file = smwout + sm_k
     sm_ndr = ndr + sm_k
     !
     !
     if ( nbeg < 0 ) then      ! >>>>>>>>>>>>>>>>>>>>>>>>>> !

        !======================================================================
        !    nbeg = -1 or nbeg = -2 or nbeg = -3
        !======================================================================

        if( nbeg == -1 ) then
           call readfile_new                                            &
                &     ( 0, sm_ndr,h,hold,nfi,rep_el(sm_k)%cm,rep_el(sm_k)%cm,rep(sm_k)%taus,  &
                &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,         &
                &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,                             &
                &       xnhe0(sm_k),xnhem(sm_k),vnhe(sm_k),xnhp0(sm_k),xnhpm(sm_k),vnhp(sm_k),&
                &       ekincm(sm_k),                           &
                &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,rep(sm_k)%fion, &
                &       tps)
        endif


        !     
        call phfac( rep(sm_k)%tau0, ei1, ei2, ei3, eigr )
        !
        call initbox ( rep(sm_k)%tau0, taub, irb )
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
        if(trane) then            ! >>>>>>>>>>>>>>>>>> ! 
           if(sm_k == 1) then
              !       
              !     random initialization
              !
              call randin(1,n,ng0,ngw,ampre,rep_el(sm_k)%cm)
           else
              rep_el(sm_k)%cm = rep_el(1)%cm
           endif
        else if(nbeg.eq.-3) then
           if(sm_k == 1) then
              !       
              !     gaussian initialization
              !
              call gausin(eigr,rep_el(sm_k)%cm)
           else
              rep_el(sm_k)%cm = rep_el(1)%cm 
           endif
        end if            ! <<<<<<<<<<<<<<<<<<<< !

        !
        !
        !     prefor calculates betae (used by graham)
        !
        call prefor(eigr,betae)
        call graham(betae,rep_el(sm_k)%bec,rep_el(sm_k)%cm)
        if(iprsta.ge.3) call dotcsc(eigr,rep_el(sm_k)%cm)
        !     
        nfi=0

        !
        !     strucf calculates the structure factor sfac
        !
        call strucf(ei1,ei2,ei3,sfac)

        IF(sm_k == 1) call formf(tfirst,eself) 

        call calbec (1,nsp,eigr,rep_el(sm_k)%cm,rep_el(sm_k)%bec)
        if (tpre) call caldbec(1,nsp,eigr,rep_el(sm_k)%cm)
        !
        !
        WRITE(stdout,*) " "
        WRITE(stdout,*) " --------------------- Replica : ", sm_k
        WRITE(stdout,*) " "

        !

        call rhoofr (nfi,rep_el(sm_k)%cm,irb,eigrb,rep_el(sm_k)%bec, & 
             & rep_el(sm_k)%rhovan,rhor,rhog,rhos,enl,ekin_ar(sm_k))

        ekin = ekin_ar(sm_k)

        !
        if(iprsta.gt.0) WRITE( stdout,*) ' out from rhoofr'
        !
        !     put core charge (if present) in rhoc(r)
        !
        if ( ANY( nlcc ) ) call set_cc(irb,eigrb,rhoc)
        !
        !


        call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
             &        ei1,ei2,ei3,irb,eigrb,sfac,rep(sm_k)%tau0,rep(sm_k)%fion)

        !WRITE( stdout, * ) 'debug 1, rcmax: ', rcmax(1:nsp) 
        !WRITE( stdout, * ) 'debug 1, esr  : ', esr, eself

        !
        etot_ar(sm_k) = etot
        eht_ar(sm_k) = eht
        epseu_ar(sm_k) = epseu
        exc_ar(sm_k) = exc
        !
        !
        call compute_stress( stress, detot, h, omega )
        !
        !
        if(iprsta.gt.0 .AND. ionode) WRITE( sm_file,*) ' out from vofrho'
        if(iprsta.gt.2) call print_atomic_var( rep(sm_k)%fion, na, nsp, ' fion ', iunit = sm_file )
        ! 
        !     forces for eigenfunctions
        !
        !     newd calculates deeq and a contribution to fion
        !


        !call newd(rhor,irb,eigrb,rep_el(sm_k)%rhovan,deeq,rep(sm_k)%fion)
        call newd(rhor,irb,eigrb,rep_el(sm_k)%rhovan,rep(sm_k)%fion)

        if(ionode) WRITE( sm_file,*) ' out from newd'
        call prefor(eigr,betae)
        !
        !     if n is odd => c(*,n+1)=0
        !
        ALLOCATE( emadt2( ngw ) )
        ccc(sm_k) = dt2hbe
        if(tsde) ccc(sm_k) = dt2bye
        emadt2 = ccc(sm_k) * ema0bg


        do i=1,n,2

           call dforce(rep_el(sm_k)%bec,betae,i,rep_el(sm_k)%cm(1,i), &
                & rep_el(sm_k)%cm(1,i+1),c2,c3,rhos)
           call wave_steepest( rep_el(sm_k)%c0(:,i), rep_el(sm_k)%cm(:,i), emadt2, c2 )
           call wave_steepest( rep_el(sm_k)%c0(:,i+1), rep_el(sm_k)%cm(:,i+1), emadt2, c3 )
        end do

        DEALLOCATE( emadt2 )
        !
        !
        if(ionode) WRITE( sm_file,*) ' out from dforce'
        !
        !     buffer for wavefunctions is unit 21
        !
        if(tbuff) rewind 21

        !call hangup
        !stop 'qui'
        !
        !     nlfq needs deeq calculated in newd
        !
        if ( tfor .or. tprnfor ) call nlfq(rep_el(sm_k)%cm,eigr, &
             & rep_el(sm_k)%bec,becdr,rep(sm_k)%fion)
        !
        if(ionode) WRITE( sm_file,*) ' out from nlfq'
        ! 
        !     imposing the orthogonality
        !     ==========================================================
        !
        call calphi(rep_el(sm_k)%cm,ema0bg,rep_el(sm_k)%bec, &
             & betae,rep_el(sm_k)%phi)
        !
        !
        if(ionode) WRITE( sm_file,*) ' out from calphi'
        !     ==========================================================
        !
        if(tortho) then
           call ortho  (eigr,rep_el(sm_k)%c0,rep_el(sm_k)%phi,rep_el(sm_k)%lambda, &
                &                   bigr,iter,ccc(sm_k),eps,maxit,delt,bephi,becp)
        else
           call graham(betae,rep_el(sm_k)%bec,rep_el(sm_k)%c0)
           !
           if(ionode) WRITE( sm_file,*) ' graham  c0 '
        endif
        !
        !     nlfl needs lambda becdr and bec
        !
        if ( tfor .or. tprnfor ) call nlfl(rep_el(sm_k)%bec,becdr, &
             & rep_el(sm_k)%lambda,rep(sm_k)%fion)
        !
        if((tfor .or. tprnfor) .AND. ionode) WRITE( sm_file,*) ' out from nlfl'
        !
        if(iprsta.ge.3) CALL print_lambda( rep_el(sm_k)%lambda, n, 9, ccc(sm_k), iunit = sm_file )
        !
        if(tpre) then
           call nlfh(rep_el(sm_k)%bec,dbec,rep_el(sm_k)%lambda)
           WRITE( stdout,*) ' out from nlfh'
        endif
        !
        if(tortho) then
           call updatc(ccc(sm_k),rep_el(sm_k)%lambda,rep_el(sm_k)%phi, &
                & bephi,becp,rep_el(sm_k)%bec,rep_el(sm_k)%c0)
           !
           if(ionode) WRITE( sm_file,*) ' out from updatc'
        endif
        call calbec (nvb+1,nsp,eigr,rep_el(sm_k)%c0,rep_el(sm_k)%bec)
        if (tpre) call caldbec(1,nsp,eigr,rep_el(sm_k)%cm)
        !
        if(ionode) WRITE( sm_file,*) ' out from calbec'
        !
        !     ==============================================================
        !     cm now orthogonalized
        !     ==============================================================
        if(iprsta.ge.3) call dotcsc(eigr,rep_el(sm_k)%c0)
        !     
        if(thdyn) then
           call cell_force( fcell, ainv, stress, omega, press, wmass )
           call cell_hmove( h, hold, delt, iforceh, fcell )
           call invmat( 3, h, ainv, deth )
        endif
        !
        if( tfor ) then

           if(smlm) call PERP(rep(sm_k)%fion,rep(sm_k)%tan,paraforce(sm_k)) 

           CALL ions_hmove( rep(sm_k)%taus, rep(sm_k)%tausm, iforce, pmass, rep(sm_k)%fion, ainv, delt, na, nsp )
           CALL s_to_r( rep(sm_k)%taus, rep(sm_k)%tau0, na, nsp, h )

        endif
        !
        !     
     else         ! <<<<<<<<<<< |||||||||||| >>>>>>>>>>>>!

        !
        !======================================================================
        !       nbeg = 0, nbeg = 1 or nbeg = 2
        !======================================================================

        call readfile_new                                           &
             &     ( 1, sm_ndr,h,hold,nfi,rep_el(sm_k)%c0,rep_el(sm_k)%cm,rep(sm_k)%taus, &
             &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,         &
             &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,                   &
             &       xnhe0(sm_k),xnhem(sm_k),vnhe(sm_k),xnhp0(sm_k),xnhpm(sm_k),vnhp(sm_k),&
             &       ekincm(sm_k),                                                         &
             &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,rep(sm_k)%fion, &
             &       tps)


        !
        CALL s_to_r(  rep(sm_k)%taus,  rep(sm_k)%tau0, na, nsp, h )

        !
        if(trane.and.trhor) then
           call prefor(eigr,betae)
           call graham(betae,rep_el(sm_k)%bec,rep_el(sm_k)%c0)
           rep_el(sm_k)%cm(:, 1:n)=rep_el(sm_k)%c0(:, 1:n)
        endif
        !
        if(iprsta.gt.2) then
           call print_atomic_var( rep(sm_k)%taus, na, nsp, ' read: taus ' )
           WRITE( stdout,*) ' read: cell parameters h '
           WRITE( stdout,*)  (h(1,j),j=1,3)
           WRITE( stdout,*)  (h(2,j),j=1,3)
           WRITE( stdout,*)  (h(3,j),j=1,3)
        endif
        !
        call phfac(rep(sm_k)%tau0,ei1,ei2,ei3,eigr)
        call strucf(ei1,ei2,ei3,sfac)
        call formf(tfirst,eself)
        call calbec (1,nsp,eigr,rep_el(sm_k)%c0,rep_el(sm_k)%bec)
        if (tpre) call caldbec(1,nsp,eigr,rep_el(sm_k)%c0)
        !
     end if               ! <<<<<<<<<<<<<<<<<<<<<<<< !


     !==============================================end of if(nbeg.lt.0)====
     !=============== END of NBEG selection ================================

     !
     !
     !     =================================================================
     !     restart with new averages and nfi=0
     !     =================================================================
     if( nbeg <= 0 ) then
        rep(sm_k)%acc = 0.0d0
        nfi=0
     end if
     !
     if( ( .not. tfor ) .and. ( .not. tprnfor ) ) then
        rep(sm_k)%fion = 0.d0
     end if
     !
     if( .not. tpre ) then
        stress (:,:) = 0.d0
     endif
     !         
     fccc = 1.0d0
     !
     IF(sm_k==1) nomore = nomore + nfi
     !
     !
     !      call cofmass( rep(sm_k)%taus, rep(sm_k)%cdm0 )
     !
     !
  ENDDO  INI_REP_LOOP      ! <<<<<<<<<<<<<<<<<<<<<<< ! 
  !
  !
  IF(smlm .and. nbeg < 0) THEN

     !  .... temp assignment ...

     DO sm_k=0,sm_p
        Nullify(p_taup(sm_k)%d3)
        p_taup(sm_k)%d3 => rep(sm_k)%taum
     ENDDO

     !
     !    ... calculate LM for smd .
     !

     call SMLAMBDA(p_tau0,p_taup,p_tan,con_ite,err_const)

     IF(maxlm <= con_ite ) write(stdout, *) "Warning ! : ", maxlm, con_ite 
     !
     !    ... to reduced (crystal) coordinates.
     ! 

     DO sm_k =1,smpm
        call r_to_s(rep(sm_k)%tau0,rep(sm_k)%taus, na, nsp, ainv)
     ENDDO

     !
     !  .... back to regular assignemnt
     !
     DO sm_k=0,sm_p
        Nullify(p_taup(sm_k)%d3)
        p_taup(sm_k)%d3 => rep(sm_k)%taup
     ENDDO

     !
  ENDIF
  !
  !
  INI2_REP_LOOP : DO sm_k=1,smpm   ! >>>>>>>>>>>>>>>>>>>>>> !

     IF(tfor .OR. smlm) THEN  
        call phfac(rep(sm_k)%tau0,ei1,ei2,ei3,eigr)
        call calbec (1,nsp,eigr,rep_el(sm_k)%c0,rep_el(sm_k)%bec)
        if (tpre) call caldbec(1,nsp,eigr,rep_el(sm_k)%c0)
     ENDIF
     !
     xnhp0(sm_k)=0.
     xnhpm(sm_k)=0.
     vnhp(sm_k) =0.
     rep(sm_k)%fionm=0.d0
     CALL ions_vel(  rep(sm_k)%vels,  rep(sm_k)%taus,  rep(sm_k)%tausm, na, nsp, delt )
     xnhh0(:,:)=0.
     xnhhm(:,:)=0.
     vnhh (:,:) =0.
     velh (:,:)=(h(:,:)-hold(:,:))/delt
     !
     !     ======================================================
     !     kinetic energy of the electrons
     !     ======================================================
     !
     ekincm(sm_k)=0.0
     CALL elec_fakekine2( ekincm(sm_k), ema0bg, emass, rep_el(sm_k)%c0, rep_el(sm_k)%cm, ngw, n, delt )

     xnhe0(sm_k)=0.
     xnhem(sm_k)=0.
     vnhe(sm_k) =0.
     !
     rep_el(sm_k)%lambdam(:,:)=rep_el(sm_k)%lambda(:,:)
     !
  ENDDO INI2_REP_LOOP  ! <<<<<<<<<<<<<<<<<<<  !

  !
  !
  ! ... Copying the init and final state.
  !
  rep(0)%taup = rep(0)%tau0
  rep(sm_p)%taup = rep(sm_p)%tau0

  ! ... Center of mass
  !
  WRITE(stdout,*) " "
  WRITE(stdout,*) " Center of mass "
  WRITE(stdout,*) " "
  DO sm_k=0,sm_p
     call ions_cofmass( rep(sm_k)%tau0, pmass, na, nsp, rep(sm_k)%cdm0)
     WRITE(stdout,'(i4,1x,3f8.5)') sm_k, (rep(sm_k)%cdm0(i),i=1,3)
  ENDDO
  WRITE(stdout,*) " "
  !
  ! ... Initial geometry ..
  !
  IF(ionode) THEN
     IF(nbeg < 0) THEN
        WRITE(unico,'(9(1x,f9.5))') &
             & (((rep(j)%taum(i,ia),i=1,3),ia=1,SUM(na(1:nsp))),j=0,sm_p)
     ELSE
        WRITE(unico,'(9(1x,f9.5))') &
             & (((rep(j)%tau0(i,ia),i=1,3),ia=1,SUM(na(1:nsp))),j=0,sm_p)
     ENDIF
     !
#if defined (__ORIGIN) || defined (__T3E)
     call flush(stdout)
#elif defined (__AIX) || defined (__ABSOFT)
     call flush_(stdout)
#endif
     !
  ENDIF

  !
  !
  !======================================================================
  !
  !           basic loop for molecular dynamics starts here
  !
  !======================================================================
  !

  !
  !
  !
  WRITE(stdout,*) " _____________________________________________"
  WRITE(stdout,*) " "
  WRITE(stdout,*) " *****    Entering SMD LOOP   ***** "
  WRITE(stdout,*) " _____________________________________________"
  !
  !
  !
  call stop_clock( 'initialize' )

  MAIN_LOOP: DO

     !
     !
     EL_REP_LOOP: DO sm_k=1,smpm        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
        !
        sm_file = smwout + sm_k
        sm_ndr = ndr + sm_k
        !
        !
        !     calculation of velocity of nose-hoover variables
        !
        if(.not.tsde) fccc(sm_k)=1./(1.+frice)
        if(tnosep)then
           CALL ions_nosevel( vnhp(sm_k), xnhp0(sm_k), xnhpm(sm_k), delt )
        endif
        if(tnosee)then
           CALL elec_nosevel( vnhe(sm_k), xnhe0(sm_k), xnhem(sm_k), delt, fccc(sm_k) )
        endif
        if(tnoseh) then
           CALL cell_nosevel( vnhh, xnhh0, xnhhm, delt, velh, h, hold )
        endif
        ! 
        !
        call initbox ( rep(sm_k)%tau0, taub, irb )
        call phbox(taub,eigrb)
        call phfac( rep(sm_k)%tau0,ei1,ei2,ei3,eigr) 
        !
        !     strucf calculates the structure factor sfac
        !
        call strucf(ei1,ei2,ei3,sfac)
        if (thdyn) call formf(tfirst,eself)
        !
        IF(sm_k==1) THEN
           nfi=nfi+1
           tlast=(nfi.eq.nomore)
        ENDIF
        !
        !
        IF((nfi.eq.0) .OR. (mod(nfi-1,iprint).eq.0)) THEN
           WRITE(stdout,*) " "
           WRITE(stdout,*) " -----------  REPLICA : ", sm_k
           WRITE(stdout,*) " "
        ENDIF
        !
        !

        call rhoofr (nfi,rep_el(sm_k)%c0,irb,eigrb,rep_el(sm_k)%bec, &
             & rep_el(sm_k)%rhovan,rhor,rhog,rhos,enl,ekin_ar(sm_k))

        ekin = ekin_ar(sm_k)


        IF(trhow) THEN
           WRITE(stdout,*) " TRHOW set to .FASLE., it is not implemented with SMD" 
           trhow = .false.   
        ENDIF
        !
        ! Y.K.
        !#ifdef __PARA     
        !      if(trhow .and. tlast) call write_rho(47,nspin,rhor)
        !#else
        !      if(trhow .and. tlast) write(47) ((rhor(i,is),i=1,nnr),is=1,nspin)
        !#endif
        !
        !     put core charge (if present) in rhoc(r)
        !
        if ( ANY( nlcc ) ) call set_cc(irb,eigrb,rhoc)
        !
        !
        IF(tlast) THEN
           WRITE(stdout,*) " "
           WRITE(stdout,*) " -----------  REPLICA : ", sm_k
           WRITE(stdout,*) " "
        ENDIF
        !
        IF(.NOT. tfirst) esr = esr_ar(sm_k)
        !

        call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                 &
             &            ei1,ei2,ei3,irb,eigrb,sfac,rep(sm_k)%tau0,rep(sm_k)%fion)

        !WRITE( stdout, * ) 'debug 2, rcmax: ', rcmax(1:nsp) 
        !WRITE( stdout, * ) 'debug 2, esr  : ', esr, eself

        !
        etot_ar(sm_k)  = etot
        eht_ar(sm_k)   = eht
        epseu_ar(sm_k) = epseu   
        exc_ar(sm_k)   = exc   
        !
        IF(tfirst) esr_ar(sm_k) = esr
        !
        !
        call compute_stress( stress, detot, h, omega )

        !
        enthal(sm_k)=etot+press*omega
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

        !call newd(rhor,irb,eigrb,rep_el(sm_k)%rhovan,deeq,rep(sm_k)%fion)
        call newd(rhor,irb,eigrb,rep_el(sm_k)%rhovan,rep(sm_k)%fion)
        call prefor(eigr,betae)


        !
        !==== set friction ====
        !
        if( tnosee ) then
           verl1 = 2.0d0 * fccc(sm_k)
           verl2 = 1.0d0 - verl1
           verl3 = 1.0d0 * fccc(sm_k)
        else
           verl1=2./(1.+frice) 
           verl2=1.-verl1
           verl3=1./(1.+frice)
        end if
        !
        !==== start loop over electronic degrees of freedom ====
        !

        ALLOCATE( emadt2( ngw ) )
        ALLOCATE( emaver( ngw ) )
        emadt2 = dt2bye * ema0bg
        emaver = emadt2 * verl3
        !
        do i=1,n,2
           call dforce(rep_el(sm_k)%bec,betae, &
                & i,rep_el(sm_k)%c0(1,i),rep_el(sm_k)%c0(1,i+1),c2,c3,rhos)
           !
           if(tsde) then
              CALL wave_steepest( rep_el(sm_k)%cm(:, i), rep_el(sm_k)%c0(:, i), emadt2, c2 )
              CALL wave_steepest( rep_el(sm_k)%cm(:, i+1), rep_el(sm_k)%c0(:, i+1), emadt2, c3 )
           else 
              CALL wave_verlet( rep_el(sm_k)%cm(:, i), rep_el(sm_k)%c0(:, i), &
                   verl1, verl2, emaver, c2 )
              CALL wave_verlet( rep_el(sm_k)%cm(:, i+1), rep_el(sm_k)%c0(:, i+1), &
                   verl1, verl2, emaver, c3 )
           endif
           if (ng0.eq.2) then
              rep_el(sm_k)%cm(1,  i)=cmplx(real(rep_el(sm_k)%cm(1,  i)),0.0)
              rep_el(sm_k)%cm(1,i+1)=cmplx(real(rep_el(sm_k)%cm(1,i+1)),0.0)
           end if
        end do

        ccc(sm_k) = fccc(sm_k) * dt2bye
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
        if ( tfor .or. tprnfor ) call nlfq(rep_el(sm_k)%c0,eigr, &
             & rep_el(sm_k)%bec,becdr,rep(sm_k)%fion)
        !
        if( tfor .or. thdyn ) then
           !
           ! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
           !
           rep_el(sm_k)%lambdap(:,:) = 2.d0*rep_el(sm_k)%lambda(:,:)-rep_el(sm_k)%lambdam(:,:)
           rep_el(sm_k)%lambdam(:,:)= rep_el(sm_k)%lambda (:,:)
           rep_el(sm_k)%lambda (:,:)= rep_el(sm_k)%lambdap(:,:)
        endif
        !
        !     calphi calculates phi
        !     the electron mass rises with g**2
        !
        call calphi(rep_el(sm_k)%c0,ema0bg,rep_el(sm_k)%bec,betae,rep_el(sm_k)%phi)
        !
        !     begin try and error loop (only one step!)
        !
        !       nlfl and nlfh need: lambda (guessed) becdr
        !
        if ( tfor .or. tprnfor ) call nlfl(rep_el(sm_k)%bec,becdr,rep_el(sm_k)%lambda, &
             & rep(sm_k)%fion)
        !
        !
        !     This part is not compatible with SMD.
        !
        if(tpre) then
           call nlfh(rep_el(sm_k)%bec,dbec,rep_el(sm_k)%lambda)
           call ions_thermal_stress( thstress, pmass, omega, h, rep(sm_k)%vels, nsp, na )
           stress = stress + thstress
        endif

        !
     ENDDO EL_REP_LOOP   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

     !_________________________________________________________________________!
     !
     !
     IF(smlm) THEN ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
        !
        !
        ! ... Transforming ionic forces to perp force. 
        !

        call TANGENT(p_tau0,p_tan)

        DO sm_k=1,smpm
           call PERP(rep(sm_k)%fion,rep(sm_k)%tan,paraforce(sm_k))
        ENDDO


        DO sm_k=1,smpm 

           deviation(sm_k) = 0.d0
           isa = 0
           DO is=1,nsp
              DO ia=1,na(is)
                 isa = isa + 1
                 DO i=1,3

                    deviation(sm_k) = deviation(sm_k) + (rep(sm_k)%fion(i,isa)* &
                         &    iforce(i,isa))**2.d0   

                    workvec(i,ia,is) = rep(sm_k)%fion(i,isa)*iforce(i,isa)

                 ENDDO
              ENDDO
           ENDDO

           deviation(sm_k) = DSQRT(deviation(sm_k))
           maxforce(sm_k) = MAX(ABS(MAXVAL(workvec)),ABS(MINVAL(workvec)))
        ENDDO
        !
     ENDIF ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
     !
     !
     !  ... Writing main output ...
     !
     !
     IF(tfirst) THEN
        write(stdout,*) " " 
        write(stdout,1997)  
        write(stdout,*) " " 
     ENDIF
     !
     IF(smlm) THEN
        write(stdout,1998) nfi,SUM(etot_ar(0:sm_p)),con_ite,MAXVAL(err_const), &
             & MAX(MAXLOC(err_const),0),MAXVAL(ekinc(1:smpm)), &
             & MAXVAL(ekinc(1:smpm)-pre_ekinc(1:smpm)), & 
             & SUM(deviation)/(sm_p-1),MAXVAL(maxforce)
     ENDIF
     !
1997 format(/2x,'nfi',3x,'SUM(E)',3x,'con_ite',5x,'con_usat',5x, &
          &     'cons_pls',5x,'max_ekinc',3x,   &
          &     'deviation',3x,'maxforce')
     !
1998 format(i5,1x,F14.5,1x,i5,1x,E12.5, &
          & 1x,i3,1x,f8.5, &
          & 1x,f8.5, &
          & 1x,E12.5,1x,E12.5)
     !
     !
     !________________________________________________________________________!
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
     !     CELL_dynamics
     !
     hgamma(:,:) = 0.d0
     if(thdyn) then

         call cell_force( fcell, ainv, stress, omega, press, wmass )

         call cell_move( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, vnhh, velh, tsdc )

         !
         velh(:,:) = (hnew(:,:)-hold(:,:))/twodel
         !
         call cell_gamma( hgamma, ainv, h, velh )

     endif
     !
     !======================================================================
     !
     TFOR_IF : if( tfor ) then

        ION_REP_LOOP : DO sm_k=1,smpm  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !


           CALL ions_move( rep(sm_k)%tausp, rep(sm_k)%taus, rep(sm_k)%tausm, iforce, pmass, rep(sm_k)%fion, &
             ainv, delt, na, nsp, fricp, hgamma, rep(sm_k)%vels, tsdp, tnosep, rep(sm_k)%fionm, vnhp(sm_k), velsp, &
             rep(sm_k)%velsm )

           !
           !cc   call cofmass(velsp,rep(sm_k)%cdmvel)
           !         call cofmass(rep(sm_k)%tausp,cdm)
           !         do is=1,nsp
           !            do ia=1,na(is)
           !               do i=1,3
           !cc   velsp(i,ia,is)=velsp(i,ia,is)-cdmvel(i)
           !                  tausp(i,ia,is)=tausp(i,ia,is) ! +cdm0(i)-cdm(i)
           !               enddo
           !            enddo
           !         enddo
           !

           !
           !  ... taup is obtained from tausp ...
           !
           CALL  s_to_r( rep(sm_k)%tausp, rep(sm_k)%taup, na, nsp, hnew )

        ENDDO ION_REP_LOOP             ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !

     endif TFOR_IF
     !
     !---------------------------------------------------------------------------
     !              String method const  .. done in real coordiantes ...
     !---------------------------------------------------------------------------
     !
     !

     call ARC(p_taup,arc_pre,t_arc_pre,1)

     !
     !
     IF(mod(nfi,lmfreq) == 0 ) THEN
        IF(smlm) THEN
           !

           call SMLAMBDA(p_taup,p_tau0,p_tan,con_ite,err_const)

           !
        ENDIF
     ENDIF
     !

     call ARC(p_taup,arc_now,t_arc_now,1)

     call ARC(p_taup,arc_tot,t_arc_tot,0)
     !
     !
     !     ... move back to reduced coordiinates
     !     

     DO sm_k=1,smpm 
        call r_to_s(rep(sm_k)%taup,rep(sm_k)%tausp, na, nsp, ainv)
     ENDDO

     !
     ! 
     !    
     POST_REP_LOOP : DO sm_k = 1,smpm          ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
        !
        sm_file  =  smwout + sm_k
        sm_ndw = ndw + sm_k
        !
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
           !        ... phfac calculates eigr
           !
           call phfac(rep(sm_k)%taup,ei1,ei2,ei3,eigr)

        else 
           !
           call phfac(rep(sm_k)%tau0,ei1,ei2,ei3,eigr)
           !
        end if
        !
        !
        !        ... prefor calculates betae
        !
        call prefor(eigr,betae)
        !
        !
        !---------------------------------------------------------------------------
        !                    imposing the orthogonality
        !---------------------------------------------------------------------------
        !
        if(tortho) then
           call ortho                                                     &
                &         (eigr,rep_el(sm_k)%cm,rep_el(sm_k)%phi,rep_el(sm_k)%lambda, &
                & bigr,iter,ccc(sm_k),eps,maxit,delt,bephi,becp)
        else
           call graham(betae,rep_el(sm_k)%bec,rep_el(sm_k)%cm)
           if(iprsta.gt.4) call dotcsc(eigr,rep_el(sm_k)%cm)
        endif
        !
        !---------------------------------------------------------------------------
        !                   correction to displacement of ions
        !---------------------------------------------------------------------------
        !
        if(iprsta.ge.3) CALL print_lambda( rep_el(sm_k)%lambda, n, 9, 1.0d0 )
        !
        if(tortho) call updatc(ccc(sm_k),rep_el(sm_k)%lambda,rep_el(sm_k)%phi,bephi, &
             & becp,rep_el(sm_k)%bec,rep_el(sm_k)%cm)
        !
        call calbec (nvb+1,nsp,eigr,rep_el(sm_k)%cm,rep_el(sm_k)%bec)
        if (tpre) call caldbec(1,nsp,eigr,rep_el(sm_k)%cm)
        !
        if(iprsta.ge.3)  call dotcsc(eigr,rep_el(sm_k)%cm)
        !
        !---------------------------------------------------------------------------
        !                  temperature monitored and controlled
        !---------------------------------------------------------------------------
        !
        ekinp(sm_k)  = 0.d0
        ekinpr(sm_k) = 0.d0
        !
        !     ionic kinetic energy 
        !
        if( tfor ) then
         CALL ions_vel( rep(sm_k)%vels, rep(sm_k)%tausp, rep(sm_k)%tausm, na, nsp, delt )
         CALL ions_kinene( ekinp(sm_k), rep(sm_k)%vels, na, nsp, hold, pmass )
        endif
        !
        !     ionic temperature
        !
        if( tfor ) then
          CALL ions_temp( tempp(sm_k), temps, ekinpr(sm_k), rep(sm_k)%vels, na, nsp, hold, pmass )
        endif
        !
        !     fake electronic kinetic energy
        !
        call elec_fakekine2( ekinc0(sm_k), ema0bg, emass, rep_el(sm_k)%c0, rep_el(sm_k)%cm, ngw, n, delt )

        !     ... previous ekinc
        !
        pre_ekinc(sm_k) = ekinc(sm_k)

        ekinc(sm_k) = 0.5 * ( ekinc0(sm_k) + ekincm(sm_k) )

        !
        !
        !     fake cell-parameters kinetic energy
        !
        ekinh=0.
        if(thdyn) then
          call cell_kinene( ekinh, temphh, velh )
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
           CALL ions_noseupd( xnhpp(sm_k), xnhp0(sm_k), xnhpm(sm_k), delt, qnp, ekinpr(sm_k), gkbt, vnhp(sm_k) )
        endif
        if(tnosee)then
           call elec_noseupd( xnhep(sm_k), xnhe0(sm_k), xnhem(sm_k), delt, qne, ekinc(sm_k), ekincw, vnhe(sm_k) )
        endif
        if(tnoseh)then
           call cell_noseupd( xnhhp, xnhh0, xnhhm, delt, qnh, temphh, temph, vnhh )
        endif
        !
        ! warning! thdyn and tcp/tcap are not compatible yet!!!
        !
        if(tcp.or.tcap.and.tfor.and.(.not.thdyn)) then
           if(tempp(sm_k).gt.temp1.or.tempp(sm_k).lt.temp2.and.tempp(sm_k).ne.0.d0) then
             call  ions_vrescal( tcap, tempw, tempp(sm_k), rep(sm_k)%taup, rep(sm_k)%tau0, rep(sm_k)%taum, &
               na, nsp, rep(sm_k)%fion, iforce, pmass, delt )
           end if
        end if
        !
        !
        ! ------------------------------
        !
        !
        IF(mod(nfi-1,iprint).eq.0 .or. tlast) then
           write(stdout,*) " "
           write(stdout,*) " ---------------------- Replica : ", sm_k
           write(stdout,*) " "
        ENDIF
        !
        !
        if(mod(nfi-1,iprint).eq.0 .or. (nfi.eq.(nomore))) then
           call eigs(nspin,nx,nupdwn,iupdwn,f,rep_el(sm_k)%lambda)
           WRITE( stdout,*)
        endif
        !
        epot(sm_k)=eht_ar(sm_k)+epseu_ar(sm_k)+exc_ar(sm_k)
        !
        rep(sm_k)%acc(1)=rep(sm_k)%acc(1)+ekinc(sm_k)
        rep(sm_k)%acc(2)=rep(sm_k)%acc(2)+ekin_ar(sm_k)
        rep(sm_k)%acc(3)=rep(sm_k)%acc(3)+epot(sm_k)
        rep(sm_k)%acc(4)=rep(sm_k)%acc(4)+etot_ar(sm_k)
        rep(sm_k)%acc(5)=rep(sm_k)%acc(5)+tempp(sm_k)
        !
        econs(sm_k)=ekinp(sm_k)+ekinh+enthal(sm_k)
        econt(sm_k)=econs(sm_k)+ekinc(sm_k)
        !
        if(tnosep)then
           econt(sm_k)=econt(sm_k)+0.5*qnp*vnhp(sm_k)*vnhp(sm_k)+     gkbt*xnhp0(sm_k)
        endif
        if(tnosee)then
           econt(sm_k)=econt(sm_k)+0.5*qne*vnhe(sm_k)*vnhe(sm_k)+2.*ekincw*xnhe0(sm_k)
        endif
        if(tnoseh)then
           do i=1,3
              if(thdiag) then
                 econt(sm_k)=econt(sm_k)+0.5*qnh*vnhh(i,i)*vnhh(i,i)+                 &
                      &                temph/factem*xnhh0(i,i)
              else
                 do j=1,3
                    econt(sm_k)=econt(sm_k)+0.5*qnh*vnhh(i,j)*vnhh(i,j)+              &
                         &                 temph/factem*xnhh0(i,j)
                 enddo
              endif
           enddo
        endif
        !
        !     ... Writing the smfiles ...
        !
        if(mod(nfi-1,iprint).eq.0.or.tfirst)  then
           if(ionode) WRITE( sm_file,*)
           if(ionode) WRITE( sm_file,1949)
        end if
        !
        tps=nfi*delt*2.4189d-5
        !
        if(ionode) WRITE( sm_file,1950) nfi, ekinc(sm_k), int(tempp(sm_k)), &
             &              etot_ar(sm_k), econs(sm_k), econt(sm_k),              &
             &              arc_now(sm_k),t_arc_now,arc_pre(sm_k),arc_tot(sm_k),  &
             &              deviation(sm_k),maxforce(sm_k),paraforce(sm_k)
        !
        !

        ! Y.K.
        !      write(8,2948) tps,ekinc,int(temphc),int(tempp),enthal,econs,      &
        !     &              econt,                                              &
        !     &              vnhh(3,3),xnhh0(3,3),vnhp,xnhp0

        !
        !     c              frice,frich,fricp
        ! 
1949    format(/2x,'nfi',4x,'ekinc',1x,'tempp',8x,'etot',7x,'econs',7x,'econt',   &
             &     3x,'arc_diff',2x,'real',4x,'ori_arc',4x,'tot_arc',3x,'dev',4x,'maxF',4x,'paraF')
        !
        !cc     f       7x,'econs',7x,'econt',3x,'frice',2x,'frich',2x,'fricp')
        !
1950    format(i5,1x,f8.5,1x,i5,1x,f11.5,1x,f11.5,1x,f11.5, &
             & 1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5)
        !
#if defined (__ORIGIN) || defined (__T3E)
        call flush(sm_file)
#elif defined (__AIX) || defined (__ABSOFT)
        call flush_(sm_file)
#endif
        !
        !
        ! 
2948    format(f8.5,1x,f8.5,1x,i6,1x,i5,3(1x,f11.5),4(1x,f7.4))
        !
        If( tfor ) then
           if ( ionode ) then
              IF(tlast .OR. (mod(nfi,codfreq) == 0)) THEN
                 write(unico,*) "== COORD ==  rep : ", sm_k
                 write(unico,3340) ((h(i,j),i=1,3),j=1,3)
                 write(unico,'(3f12.8)') ((rep(sm_k)%tau0(i,ia),i=1,3),ia=1,SUM(na(1:nsp)))
              ENDIF
              !
              IF(tlast .OR. (mod(nfi,forfreq) == 0)) THEN
                 write(unifo,*) "== FORCE ==  rep : ", sm_k
                 write(unifo,'(3f12.8)') ((rep(sm_k)%fion(i,ia),i=1,3),ia=1,SUM(na(1:nsp)))
              ENDIF
              !
              IF(tlast) THEN
                 write(unist,3340) ((stress(i,j),i=1,3),j=1,3)
              ENDIF
              !
#if defined (__ORIGIN) || defined (__T3E)
              call flush(unico)
              call flush(unifo)
#elif defined (__AIX) || defined (__ABSOFT)
              call flush_(unico)
              call flush_(unifo)
#endif
3340          format(9(1x,f9.5))
           endif
           !
           !     new variables for next step
           !
           rep(sm_k)%tausm=rep(sm_k)%taus
           rep(sm_k)%taus=rep(sm_k)%tausp
           rep(sm_k)%taum=rep(sm_k)%tau0
           rep(sm_k)%tau0=rep(sm_k)%taup
           rep(sm_k)%velsm = rep(sm_k)%vels
           rep(sm_k)%vels  = velsp
           if(tnosep) then
              xnhpm(sm_k) = xnhp0(sm_k)
              xnhp0(sm_k) = xnhpp(sm_k)
           endif
           if(tnosee) then
              xnhem(sm_k) = xnhe0(sm_k)
              xnhe0(sm_k) = xnhep(sm_k)
           endif
           if(tnoseh) then
              xnhhm(:,:) = xnhh0(:,:)
              xnhh0(:,:) = xnhhp(:,:)
           endif
        End if
        !
        if(thdyn)then
           do i=1,ngw
              ema0bg(i)=1./max(1.d0,tpiba2*ggp(i)/emaec) 
           enddo
        endif
        !
        ekincm(sm_k)=ekinc0(sm_k)
        !  
        !     cm=c(t+dt) c0=c(t)
        !
        call DSWAP(2*ngw*n,rep_el(sm_k)%c0,1,rep_el(sm_k)%cm,1)
        !
        !     now:  cm=c(t) c0=c(t+dt)
        !
        if (tfirst) then
           epre(sm_k) = etot_ar(sm_k)
           enow(sm_k) = etot_ar(sm_k)
        endif
        !
        IF(sm_k == smpm) tfirst=.false.
        !
        !     write on file ndw each isave
        !
        if( ( mod( nfi, isave ) == 0 ) .and. ( nfi < nomore ) ) then

           call writefile_new                                                    &
                &     ( sm_ndw,h,hold,nfi,rep_el(sm_k)%c0,rep_el(sm_k)%cm,rep(sm_k)%taus, &
                &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,     &
                &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,xnhe0(sm_k),xnhem(sm_k), &
                &       vnhe(sm_k),xnhp0(sm_k),xnhpm(sm_k),vnhp(sm_k),ekincm(sm_k),       &
                &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,         &
                &       rep(sm_k)%fion, tps)

        endif
        !
        !
        epre(sm_k) = enow(sm_k)
        enow(sm_k) = etot_ar(sm_k)
        !
        frice = frice * grease
        fricp = fricp * greasp
        frich = frich * greash

        !     =====================================================
        !
        delta_etot = ABS( epre(sm_k) - enow(sm_k) )
        tstop = check_stop_now()
        tconv = .FALSE.

        ! Y.K.
        !      IF( tconvthrs%active ) THEN
        !       tconv = ( delta_etot < tconvthrs%derho ) .AND. ( ekinc(sm_k) < tconvthrs%ekin )
        !      END IF

        !
        !
        IF(tconv) THEN
           WRITE(stdout,*) "TCONV set to .F. : not implemented with calculation = smd "
           tconv = .false. 
        ENDIF
        !
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

     ENDDO POST_REP_LOOP ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !


     if( (nfi >= nomore) .OR. tstop ) EXIT MAIN_LOOP

  END DO MAIN_LOOP

  !
  !============================= END of main LOOP ======================
  !

  !
  !  Here copy relevant physical quantities into the output arrays/variables
  !  for 0 and sm_p replicas.
  !

  etot_out(0) = etot_ar(0)
  etot_out(sm_p) = etot_ar(sm_p)

  isa = 0
  do is = 1, nsp
     do ia = 1, na(is)
        isa = isa + 1
        ipos = ind_srt( isa )
        tau( 1, ipos , 0 ) = rep(0)%tau0( 1, isa )
        tau( 2, ipos , 0 ) = rep(0)%tau0( 2, isa )
        tau( 3, ipos , 0 ) = rep(0)%tau0( 3, isa )
        tau( 1, ipos , sm_p) = rep(sm_p)%tau0( 1, isa )
        tau( 2, ipos , sm_p) = rep(sm_p)%tau0( 2, isa )
        tau( 3, ipos , sm_p) = rep(sm_p)%tau0( 3, isa )
        fion_out( 1, ipos, 0 ) = rep(0)%fion( 1, isa )
        fion_out( 2, ipos, 0 ) = rep(0)%fion( 2, isa )
        fion_out( 3, ipos, 0 ) = rep(0)%fion( 3, isa )
        fion_out( 1, ipos, sm_p ) = rep(sm_p)%fion( 1, isa )
        fion_out( 2, ipos, sm_p ) = rep(sm_p)%fion( 2, isa )
        fion_out( 3, ipos, sm_p ) = rep(sm_p)%fion( 3, isa )
     end do
  end do


  FIN_REP_LOOP : DO sm_k=1,smpm     ! >>>>>>>>>>>>>>>>>>>>>>>>>> !
     !
     sm_file = smwout + sm_k
     sm_ndw = ndw + sm_k
     !
     ! 
     !  Here copy relevant physical quantities into the output arrays/variables
     !
     !
     etot_out(sm_k) = etot_ar(sm_k)
     !
     isa = 0
     do is = 1, nsp
        do ia = 1, na(is)
           isa = isa + 1
           ipos = ind_srt( isa )
           tau( 1, ipos , sm_k) = rep(sm_k)%tau0( 1, isa )
           tau( 2, ipos , sm_k) = rep(sm_k)%tau0( 2, isa )
           tau( 3, ipos , sm_k) = rep(sm_k)%tau0( 3, isa )
           fion_out( 1, ipos, sm_k ) = rep(sm_k)%fion( 1, isa )
           fion_out( 2, ipos, sm_k ) = rep(sm_k)%fion( 2, isa )
           fion_out( 3, ipos, sm_k ) = rep(sm_k)%fion( 3, isa )
        end do
     end do

     !  Calculate statistics

     anor=1.d0/dfloat(nfi)
     do i=1,nacc
        rep(sm_k)%acc(i)=rep(sm_k)%acc(i)*anor
     end do
     !
     !
     if(ionode) WRITE( sm_file,1951)
1951 format(//'              averaged quantities :',/,                 &
          &       9x,'ekinc',10x,'ekin',10x,'epot',10x,'etot',5x,'tempp')
     if(ionode) WRITE( sm_file,1952) (rep(sm_k)%acc(i),i=1,nacc)
1952 format(4f14.5,f10.1)
     !
     !
#if defined (__ORIGIN) || defined (__T3E)
     call flush(sm_file)
#elif defined (__AIX) || defined (__ABSOFT)
     call flush_(sm_file)
#endif
     !
     !
     IF( sm_k == smpm ) THEN
        call print_clock( 'initialize' )
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
     END IF
     !
     !

     call writefile_new                                                    &
          &     ( sm_ndw,h,hold,nfi,rep_el(sm_k)%c0,rep_el(sm_k)%cm,rep(sm_k)%taus, &
          &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,     &
          &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,xnhe0(sm_k),xnhem(sm_k), &
          &       vnhe(sm_k),xnhp0(sm_k),xnhpm(sm_k),vnhp(sm_k),ekincm(sm_k),       &
          &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,         &
          &       rep(sm_k)%fion, tps)

     !
     !
     !
     if(iprsta.gt.1) CALL print_lambda( rep_el(sm_k)%lambda, n, n, 1.0d0, iunit = sm_file )
     !
     if( tfor .or. tprnfor ) then
        if(ionode) WRITE( sm_file,1970) ibrav, alat
        if(ionode) WRITE( sm_file,1971)
        do i=1,3
           if(ionode) WRITE( sm_file,1972) (h(i,j),j=1,3)
        enddo
        if(ionode) WRITE( sm_file,1973)
        isa = 0
        do is=1,nsp
           do ia=1,na(is)
              isa = isa + 1
              if(ionode) WRITE( sm_file,1974) is,ia,(rep(sm_k)%tau0(i,isa),i=1,3),       &
                   &            ((ainv(j,1)*rep(sm_k)%fion(1,isa)+ainv(j,2)*rep(sm_k)%fion(2,isa)+    &
                   &              ainv(j,3)*rep(sm_k)%fion(3,isa)),j=1,3)
           end do
        end do
        if(ionode) WRITE( sm_file,1975)
        isa = 0
        do is=1,nsp
           do ia=1,na(is)
              isa = isa + 1
              if(ionode) WRITE( sm_file,1976) is,ia,(rep(sm_k)%taus(i,isa),i=1,3)
           end do
        end do
     endif
     conv_elec = .TRUE.



1970 format(1x,'ibrav :',i4,'  alat : ',f10.4,/)
1971 format(1x,'lattice vectors',/)
1972 format(1x,3f10.4)
1973 format(/1x,'Cartesian coordinates (a.u.)              forces (redu. units)' &
          &       /1x,'species',' atom #', &
          &           '   x         y         z      ', &
          &           '   fx        fy        fz'/)
1974 format(1x,2i5,3f10.4,2x,3f10.4)
1975 format(/1x,'Scaled coordinates '/1x,'species',' atom #')
1976 format(1x,2i5,3f10.4)
     if(ionode) WRITE( sm_file,1977) 

     !
600  continue
     !
  ENDDO FIN_REP_LOOP               ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
  !
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
  IF( ALLOCATED( becdr ) ) DEALLOCATE( becdr )
  IF( ALLOCATED( bephi ) ) DEALLOCATE( bephi )
  IF( ALLOCATED( becp ) ) DEALLOCATE( becp )
  IF( ALLOCATED( deeq ) ) DEALLOCATE( deeq )
  IF( ALLOCATED( ema0bg ) ) DEALLOCATE( ema0bg )
  IF( ALLOCATED( c2 ) ) DEALLOCATE( c2 )
  IF( ALLOCATED( c3 ) ) DEALLOCATE( c3 )

  IF( ALLOCATED( deviation )) DEALLOCATE( deviation )
  IF( ALLOCATED( maxforce )) DEALLOCATE( maxforce )
  IF( ALLOCATED( arc_now )) DEALLOCATE( arc_now )
  IF( ALLOCATED( arc_tot )) DEALLOCATE( arc_tot )
  IF( ALLOCATED( arc_pre )) DEALLOCATE( arc_pre )
  IF( ALLOCATED( paraforce )) DEALLOCATE( paraforce )
  IF( ALLOCATED( err_const )) DEALLOCATE( err_const )
  IF( ALLOCATED( pvvcheck )) DEALLOCATE( pvvcheck  )

  DO sm_k=0,sm_p  
     IF( ASSOCIATED( p_tau0(sm_k)%d3 )) NULLIFY( p_tau0(sm_k)%d3 ) 
     IF( ASSOCIATED( p_taup(sm_k)%d3 )) NULLIFY( p_taup(sm_k)%d3 ) 
     IF( ASSOCIATED( p_tan(sm_k)%d3 )) NULLIFY( p_tan(sm_k)%d3 ) 
  ENDDO

  IF( ALLOCATED( p_tau0 )) DEALLOCATE( p_tau0  )
  IF( ALLOCATED( p_taup )) DEALLOCATE( p_taup )
  IF( ALLOCATED( p_tan )) DEALLOCATE( p_tan )

  CALL deallocate_elct()
  CALL deallocate_core()
  CALL deallocate_uspp()
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

  ! CALL deallocate_wavefunctions()


  CALL deallocate_smd_rep()
  CALL deallocate_smd_ene()


  if( ionode ) then
     CLOSE( 8 )
     CLOSE( unico )
     CLOSE( unifo )
     CLOSE( unist )
     DO sm_k=1,smpm
        sm_file =  smwout + sm_k
        CLOSE(sm_file)
     ENDDO
  end if

  !
  return
end subroutine smdmain
