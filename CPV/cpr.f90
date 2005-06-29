!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
!***  Makov Payne Correction for charged systems by Filippo De Angelis
!***********************************************************************
!***  appropriate citation for use of this code:
!***  Car-Parrinello method    R. Car and M. Parrinello, PRL 55, 2471 (1985) 
!***  current implementation   A. Pasquarello, K. Laasonen, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRL 69, 1982 (1992);
!***                           K. Laasonen, A. Pasquarello, R. Car, 
!***                           C. Lee, and D. Vanderbilt, PRB 47, 10142 (1993).
!***  implementation gga       A. Dal Corso, A. Pasquarello, A. Baldereschi,
!***                           and R. Car, PRB 53, 1180 (1996).
!***  implementation Wannier   M. Sharma, Y. Wu and R. Car, Int.J.Quantum.Chem.
!***  function dynamics        95, 821, (2003).
!***
!***  implementation           M. Sharma and R.Car, ???
!***  Electric Field
!***  ensemble-DFT
! cf. "Ensemble Density-Functional Theory for Ab Initio Molecular Dynamics
!      of Metals and Finite-Temperature Insulators"  PRL v.79,nbsp.7 (1997)
!      nbsp. Marzari, D. Vanderbilt and M.C. Payne
!***
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
!----------------------------------------------------------------------------
SUBROUTINE cprmain( tau, fion_out, etot_out )
  !----------------------------------------------------------------------------
  !
  USE kinds,                    ONLY : dbl
  USE control_flags,            ONLY : iprint, isave, thdyn, tpre, tbuff, &
                                       iprsta, trhor, tfor, tvlocw, trhow, &
                                       taurdr, tprnfor, tsdc, lconstrain
  USE control_flags,            ONLY : ndr, ndw, nbeg, nomore, tsde, tortho, &
                                       tnosee, tnosep, trane, tranp, tsdp,   &
                                       tcp, tcap, ampre, amprp, tnoseh, tolp
  USE control_flags,            ONLY : lwf, ortho_eps, ortho_max, printwfc
  USE core,                     ONLY : nlcc_any
  USE uspp_param,               ONLY : nhm, nh
  USE cvan,                     ONLY : nvb, ish
  USE uspp,                     ONLY : nkb, vkb, becsum, deeq
  USE energies,                 ONLY : eht, epseu, exc, etot, eself, enl, &
                                       ekin, atot, entropy, egrand
  USE electrons_base,           ONLY : nbspx, nbsp, ispin => fspin, f, nspin
  USE electrons_base,           ONLY : nel, iupdwn, nupdwn, nudx, nelt
  USE efield_module,            ONLY : efield, epol, tefield, allocate_efield, &
                                       efield_update, ipolp, qmat, gqq,        &
                                       evalue, berry_energy
  USE ensemble_dft,             ONLY : tens, tgrand, ninner, ismear, etemp,   &
                                       ef, tdynz, tdynf, zmass, fmass, fricz, &
                                       fricf, allocate_ensemble_dft,          &
                                       id_matrix_init, z0, c0diag, becdiag,   &
                                       bec0, v0s, vhxcs, becdrdiag, gibbsfe
  USE cg_module,                ONLY : tcg, maxiter, etresh, passop, &
                                       allocate_cg, cg_update, &
                                       itercg, c0old
  USE gvecp,                    ONLY : ngm
  USE gvecs,                    ONLY : ngs
  USE gvecb,                    ONLY : ngb
  USE gvecw,                    ONLY : ngw
  USE reciprocal_vectors,       ONLY : gstart, mill_l
  USE ions_base,                ONLY : na, nat, pmass, nax, nsp, rcmax
  USE ions_base,                ONLY : ind_srt, ions_cofmass, ions_kinene, &
                                       ions_temp, ions_thermal_stress
  USE ions_base,                ONLY : ions_vrescal, fricp, greasp, &
                                       iforce, ions_shiftvar, ityp, &
                                       atm, ind_bck
  USE cell_base,                ONLY : ainv, a1, a2, a3, frich, greash, tpiba2
  USE cell_base,                ONLY : omega, alat, ibrav, celldm
  USE cell_base,                ONLY : h, hold, deth, wmass, press
  USE cell_base,                ONLY : iforceh, cell_force, thdiag
  USE grid_dimensions,          ONLY : nnrx, nr1, nr2, nr3
  USE smooth_grid_dimensions,   ONLY : nnrsx, nr1s, nr2s, nr3s
  USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE io_global,                ONLY : io_global_start, stdout, ionode
  USE mp_global,                ONLY : mp_global_start
  USE mp,                       ONLY : mp_sum, mp_barrier
  USE dener,                    ONLY : detot
  USE derho,                    ONLY : drhor, drhog
  USE cdvan,                    ONLY : dbec, drhovan
  USE stre,                     ONLY : stress
  USE gvecw,                    ONLY : ggp
  USE parameters,               ONLY : nacx, natx, nsx, nbndxx
  USE constants,                ONLY : pi, factem, au_gpa, au_ps, gpa_au
  USE io_files,                 ONLY : psfile, pseudo_dir
  USE input,                    ONLY : iosys
  USE wave_base,                ONLY : wave_steepest, wave_verlet
  USE wave_base,                ONLY : wave_speed2, frice, grease
  USE control_flags,            ONLY : conv_elec, tconvthrs
  USE check_stop,               ONLY : check_stop_now
  USE efcalc,                   ONLY : clear_nbeg, ef_force  !Electric Field (M.S)
  USE ions_base,                ONLY : zv, ions_vel
  USE cp_electronic_mass,       ONLY : emass, emass_cutoff, emass_precond
  USE ions_positions,           ONLY : tau0, taum, taup, taus, tausm, tausp, &
                                       vels, velsm, velsp, ions_hmove, ions_move
  USE ions_nose,                ONLY : gkbt, kbt, ndega, nhpcl, qnp, vnhp, &
                                       xnhp0, xnhpm, xnhpp, ions_nosevel,  &
                                       ions_noseupd, tempw, ions_nose_nrg, &
                                       ions_nose_shiftvar
  USE electrons_nose,           ONLY : qne, ekincw, xnhe0, xnhep, xnhem,  &
                                       vnhe, electrons_nose_nrg,          &
                                       electrons_nose_shiftvar,           &
                                       electrons_nosevel, electrons_noseupd
  USE from_scratch_module,      ONLY : from_scratch
  USE from_restart_module,      ONLY : from_restart
  USE wavefunctions_module,     ONLY : c0, cm, phi => cp
  USE wannier_module,           ONLY : allocate_wannier
  USE printout_base,            ONLY : printout_base_open, &
                                       printout_base_close, &
                                       printout_pos, printout_cell, &
                                       printout_stress, print_pos_in
  USE cell_nose,                ONLY : xnhh0, xnhhm, xnhhp, vnhh, temph, &
                                       qnh, cell_nosevel, cell_noseupd,  &
                                       cell_nose_nrg, cell_nose_shiftvar
  USE cell_base,                ONLY : cell_kinene, cell_gamma, &
                                       cell_move, cell_hmove
  USE gvecw,                    ONLY : ecutw
  USE gvecp,                    ONLY : ecutp
  USE time_step,                ONLY : delt, tps, dt2, dt2by2, twodelt
  USE electrons_module,         ONLY : cp_eigs
  USE print_out_module,         ONLY : cp_print_rho
  USE cp_main_variables,        ONLY : allocate_mainvar, &
                                       acc, bec, lambda, lambdam, lambdap, &
                                       ema0bg, sfac, eigr, ei1, ei2, ei3,  &
                                       irb, becdr, taub, eigrb, rhog, rhos, &
                                       rhor, rhoc, bephi, becp
  !
  USE cell_base,                ONLY : s_to_r, r_to_s
  USE phase_factors_module,     ONLY : strucf
  USE cpr_subroutines,          ONLY : print_lambda, print_atomic_var, &
                                       ions_cofmsub, elec_fakekine
  USE wannier_subroutines,      ONLY : wannier_init, wf_closing_options, &
                                       read_efwan_param, ef_enthalpy
  USE restart_file,             ONLY : readfile, writefile
  USE constraints_module,       ONLY : check_constraint
  !
  IMPLICIT NONE
  !
  ! ... input/output variables
  !
  REAL(KIND=dbl), INTENT(INOUT) :: tau(3,nat)
  REAL(KIND=dbl), INTENT(OUT)   :: fion_out(3,nat)
  REAL(KIND=dbl), INTENT(OUT)   :: etot_out
  !
  ! ... control variables
  !
  LOGICAL :: tfirst, tlast, tstop, tconv
  LOGICAL :: ttprint    
    !  logical variable used to control printout
  !
  ! ... ionic positions, center of mass position
  !
  REAL(KIND=dbl) :: cdm0(3)
  REAL(KIND=dbl) :: cdm(3)
  !
  ! ... forces on ions
  !
  REAL(KIND=dbl) :: fion(3,natx), fionm(3,natx)
  REAL(KIND=dbl) :: maxfion
  !
  ! ... work variables
  !
  REAL(KIND=dbl) :: tempp,  fccc, savee, saveh, savep, enthal, epot, epre, &
                    enow, econs, econt, ettt, ccc, bigr, dt2bye
  REAL(KIND=dbl) :: ekinc0, ekinp, ekinpr, ekincm, ekinc
  REAL(KIND=dbl) :: temps(nsx)
  REAL(KIND=dbl) :: ekinh, temphc, temp1, temp2, randy
  REAL(KIND=dbl) :: delta_etot
  REAL(KIND=dbl) :: ftmp, enb, enbi, one, toang
  INTEGER        :: is, nacc, ia, j, iter, nfi, i, isa, ipos
  INTEGER        :: k, ii, l, m, ibeg
  !
  REAL(KIND=dbl) :: hnew(3,3), velh(3,3), hgamma(3,3), temphh(3,3)
  REAL(KIND=dbl) :: fcell(3,3)
  !
  REAL(KIND=dbl) :: b1(3), b2(3), b3(3)
  REAL(KIND=dbl) :: stress_gpa(3,3), thstress(3,3)
  !
  REAL(KIND=dbl), ALLOCATABLE :: tauw(:,:)  
    ! temporary array used to printout positions
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
  !     ====  1 Volt / meter      = 1/(5.1412*1.e+11) a.u.            ====
  !     ==================================================================
  !
  ! ... CP starts here
  !
  CALL start_clock( 'initialize' )
  !
  ! ... general variables
  !
  !=======================================================================
  !     copy-in input parameters from input_parameter module
  !=======================================================================
  !
  CALL iosys()
  !
  IF ( lwf ) CALL read_efwan_param( nbeg )
  !
  !=======================================================================
  !     initialize g-vectors, fft grids
  !=======================================================================
  !
  CALL init_dimensions()
  !
  dt2bye   = dt2 / emass
  etot_out = 0.D0
  one = 1.0d0
  toang = 0.5291772083

  tfirst   = .TRUE.
  tlast    = .FALSE.
  nacc     = 5
  !
  CALL init_geometry()
  !
  IF ( lwf ) CALL clear_nbeg( nbeg )
  !
  ! ... more initialization requiring atomic positions
  !
  nax = MAXVAL( na(1:nsp) )
  !
  IF ( iprsta > 1 ) CALL print_atomic_var( tau0, na, nsp, ' tau0 ' )
  !
  !=======================================================================
  !     allocate and initialize nonlocal potentials
  !=======================================================================
  !
  CALL nlinit()
  !
  !=======================================================================
  !     allocation of all arrays not already allocated in init and nlinit
  !=======================================================================
  !
  CALL allocate_mainvar( ngw, ngb, ngs, ngm, nr1, nr2, nr3, nnrx, nnrsx, &
                         nat, nax, nsp, nspin, nbsp, nbspx, nkb, nlcc_any )
  !
  ALLOCATE( c0(  ngw, nbspx, 1, 1 ) )
  ALLOCATE( cm(  ngw, nbspx, 1, 1 ) )
  ALLOCATE( phi( ngw, nbspx, 1, 1 ) )
  !
  CALL allocate_local_pseudo( ngs, nsp )
  !
  ALLOCATE( vkb( ngw, nkb ) )
  ALLOCATE( deeq( nhm, nhm, nat, nspin ) )
  ALLOCATE( becsum( nhm*(nhm+1)/2, nat, nspin ) )
  ALLOCATE( dbec( nkb, nbsp, 3, 3 ) )
  ALLOCATE( drhog( ngm, nspin, 3, 3 ) )
  ALLOCATE( drhor( nnrx, nspin, 3, 3 ) )
  ALLOCATE( drhovan( nhm*(nhm+1)/2, nat, nspin, 3, 3 ) )
  !
  IF ( lwf ) CALL allocate_wannier(  nbsp, nnrsx, nspin, ngm )
  !
  IF ( tens .or. tcg) CALL allocate_ensemble_dft( nkb, nbsp, ngw, nudx, &
                                          nspin, nbspx, nnrsx, natx )
  !
  IF( tcg ) CALL allocate_cg( ngw, nbspx )
  !
  IF( tefield ) CALL allocate_efield( ngw, nbspx, nhm, nax, nsp )
  !
  deeq(:,:,:,:) = 0.D0
  !
  temp1 = tempw + tolp
  temp2 = tempw - tolp
  !
  !=======================================================================
  !
  taum  = tau0
  taup  = 0.D0
  tausm = taus
  tausp = 0.D0
  vels  = 0.D0
  velsm = 0.D0
  velsp = 0.D0
  !
  hnew=h
  !
  lambda(:,:) = 0.D0
  cm(:,:,1,1) = ( 0.D0, 0.D0 )
  c0(:,:,1,1) = ( 0.D0, 0.D0 )
  !
  IF ( tens ) CALL id_matrix_init( nupdwn, nspin )
  !
  CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
  !
  IF ( lwf ) CALL wannier_init( ibrav, alat, a1, a2, a3, b1, b2, b3 )
  !
  IF ( nbeg < -1 ) THEN
     !
     !======================================================================
     !     Initialize from scratch nbeg = -2 or nbeg = -3
     !======================================================================
     !
     nfi = 0
     !
     CALL from_scratch ( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst,    &
                         eself, fion, taub, irb, eigrb, b1, b2, b3, nfi,   &
                         rhog, rhor, rhos, rhoc, enl, ekin, stress, detot, &
                         enthal, etot, lambda, lambdam, lambdap, ema0bg,   &
                         dbec, delt, bephi, becp, velh, dt2bye, iforce,    &
                         fionm, nbeg, xnhe0, xnhem, vnhe, ekincm )
     !
  ELSE
     !
     !======================================================================
     !     nbeg = -1, nbeg = 0, nbeg = 1 or nbeg = 2
     !======================================================================
     !
     IF( nbeg == -1 ) THEN
        !
        ! ... read only wavefunction cm from restart
        !
        ! ibeg = 0  
        ibeg = 1  !  DEBUG
        !
     ELSE
        !
        ! ... read all dynamic variable from restart
        !
        ibeg = 1
        !
     END IF
     !
     CALL readfile( ibeg, ndr, h, hold, nfi, c0(:,:,1,1), cm(:,:,1,1), taus, &
                    tausm, vels, velsm, acc, lambda, lambdam, xnhe0, xnhem,  &
                    vnhe, xnhp0, xnhpm, vnhp, nhpcl, ekincm, xnhh0, xnhhm,   &
                    vnhh, velh, ecutp, ecutw, delt, pmass, ibrav, celldm,    &
                    fion, tps, z0, f )
     !
     CALL from_restart( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, fion, &
                        taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, &
                        rhoc, stress, detot, enthal, lambda, lambdam, lambdap, &
                        ema0bg, dbec, bephi, becp, velh, dt2bye, fionm, ekincm )
     !
  END IF
  !
  !=======================================================================
  !     restart with new averages and nfi=0
  !=======================================================================
  !
  ! ... Fix. Center of Mass - M.S
  !
  IF ( lwf ) CALL ions_cofmass( tau0, pmass, na, nsp, cdm0 )
  !
  IF ( nbeg <= 0 ) THEN
     !
     acc = 0.D0
     nfi = 0
     !
  END IF
  !
  IF ( .NOT. tfor .AND. .NOT. tprnfor ) fion = 0.D0
  !
  IF ( .NOT. tpre ) stress = 0.D0
  !         
  fccc = 1.D0
  !
  nomore = nomore + nfi
  !
  CALL ions_cofmass( taus, pmass, na, nsp, cdm0 )
  !
  CALL stop_clock( 'initialize' )
  !
  !======================================================================
  !
  !           basic loop for molecular dynamics starts here
  !
  !======================================================================
  !
  main_loop: DO
     !
     CALL start_clock( 'total_time' )
     !
     nfi     = nfi + 1
     tlast   = ( nfi == nomore )
     ttprint = ( MOD( nfi, iprint ) == 0 )
     !
     IF( ionode .AND. ttprint ) THEN
       WRITE( stdout, fmt = '( /, " * Step ",  I6 )' ) nfi
     END IF
     !
     ! ... calculation of velocity of nose-hoover variables
     !
     IF ( .NOT. tsde ) fccc = 1.D0 / ( 1.D0 + frice )
     !
     IF ( tnosep ) CALL ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl )
     !
     IF ( tnosee ) THEN
        !
        CALL electrons_nosevel( vnhe, xnhe0, xnhem, delt )
        !
        fccc = 1.D0 / ( 1.D0 + 0.5D0 * delt * vnhe )
        !
     END IF
     !
     IF ( tnoseh ) THEN
        !
        CALL cell_nosevel( vnhh, xnhh0, xnhhm, delt )
        !
        velh(:,:) = 2.D0 * ( h(:,:) - hold(:,:) ) / delt - velh(:,:)
        !
     END IF
     ! 
     IF ( tfor .OR. thdyn .OR. tfirst ) THEN
        !
        CALL initbox( tau0, taub, irb )
        CALL phbox( taub, eigrb )
        !
     END IF
     !
     IF ( tfor .OR. thdyn ) CALL phfac( tau0, ei1, ei2, ei3, eigr ) 
     !
     ! ... strucf calculates the structure factor sfac
     !
     CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
     !
     IF ( thdyn ) CALL formf( tfirst, eself )
     !
     ! ... why this call ??? from Paolo Umary
     !
     IF ( tefield ) CALL calbec( 1, nsp, eigr, c0, bec ) ! ATTENZIONE  
     !
     IF ( ( tfor .OR. tfirst ) .AND. tefield ) CALL efield_update( eigr )
     !
     !=======================================================================
     !
     !    electronic degrees of freedom are updated here
     !
     !=======================================================================
     !
     CALL move_electrons( nfi, tfirst, tlast, b1, b2, b3, fion, &
                          enthal, enb, enbi, fccc, ccc, dt2bye )
     !
     IF ( tpre ) THEN
        !
        CALL nlfh( bec, dbec, lambda )
        !
        CALL ions_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
        !
     END IF
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
     hgamma(:,:) = 0.D0
     !
     IF ( thdyn ) THEN
        !
        CALL cell_force( fcell, ainv, stress, omega, press, wmass )
        !
        CALL cell_move( hnew, h, hold, delt, iforceh, &
                        fcell, frich, tnoseh, vnhh, velh, tsdc )
        !
        velh(:,:) = ( hnew(:,:) - hold(:,:) ) / twodelt
        !
        CALL cell_gamma( hgamma, ainv, h, velh )
        !
     END IF
     !
     !======================================================================
     !
     IF ( tfor ) THEN
        !
        IF ( lwf ) CALL ef_force( fion, na, nsp, zv )
        !
        CALL ions_move( tausp, taus, tausm, iforce, pmass, fion, &
                        ainv, delt, na, nsp, fricp, hgamma, vels, &
                        tsdp, tnosep, fionm, vnhp, velsp, velsm )
        !
        IF ( lconstrain ) THEN
           !
           ! ... constraints are imposed here
           !
           CALL s_to_r( tausp, taup, na, nsp, hnew )
           !
           CALL check_constraint( nat, taup, tau0, &
                                  fion, iforce, ityp, 1.D0, delt )
           !
           CALL r_to_s( taup, tausp, na, nsp, ainv )
           !
        END IF
        !
        CALL ions_cofmass( tausp, pmass, na, nsp, cdm )
        !
        CALL ions_cofmsub( tausp, na, nsp, cdm, cdm0 )
        !
        CALL s_to_r( tausp, taup, na, nsp, hnew )
        !
     END IF
     !     
     !---------------------------------------------------------------------------
     !              initialization with guessed positions of ions
     !---------------------------------------------------------------------------
     !
     ! ... if thdyn=true g vectors and pseudopotentials are recalculated for 
     ! ... the new cell parameters
     !
     IF ( tfor .OR. thdyn ) THEN
        !
        IF ( thdyn ) THEN
           !
           hold = h
           h    = hnew
           !
           CALL newinit( h )
           CALL newnlinit()
           !
        ELSE
           !
           hold = h
           !
        END IF
        !
        ! ... phfac calculates eigr
        !
        CALL phfac( taup, ei1, ei2, ei3, eigr )
        !
        ! ... prefor calculates vkb
        !
        CALL prefor( eigr, vkb )
        !
     END IF
     !
     !---------------------------------------------------------------------------
     !                    imposing the orthogonality
     !---------------------------------------------------------------------------
     !
     IF ( .NOT. tcg ) THEN
        !
        IF ( tortho ) THEN
           !
           CALL ortho( eigr, cm, phi, lambda, bigr, iter, ccc, &
                       ortho_eps, ortho_max, delt, bephi, becp )
           !
        ELSE
           !
           CALL gram( vkb, bec,cm )
           !
           IF ( iprsta > 4 ) CALL dotcsc( eigr, cm )
           !
        END IF
        !
     END IF
     !
     !---------------------------------------------------------------------------
     !                   correction to displacement of ions
     !---------------------------------------------------------------------------
     !
     IF ( .NOT. tcg ) THEN
        !
        IF ( iprsta >= 3 ) CALL print_lambda( lambda, nbsp, 9, 1.D0 )
        !
        IF ( tortho ) CALL updatc( ccc, lambda, phi, bephi, becp, bec, cm )
        !
        CALL calbec( nvb+1, nsp, eigr, cm, bec )
        !
        IF ( tpre ) CALL caldbec( 1, nsp, eigr, cm )
        !
        IF ( iprsta >= 3 ) CALL dotcsc( eigr, cm )
        !
     END IF
     !
     !---------------------------------------------------------------------------
     !                  temperature monitored and controlled
     !---------------------------------------------------------------------------
     !
     ekinp  = 0.D0
     ekinpr = 0.D0
     tempp  = 0.D0
     !
     ! ... ionic kinetic energy 
     !
     IF ( tfor ) THEN
        !
        CALL ions_vel( vels, tausp, tausm, na, nsp, delt )
        !
        CALL ions_kinene( ekinp, vels, na, nsp, hold, pmass )
        !
     END IF
     !
     ! ... ionic temperature
     !
     IF ( tfor ) CALL ions_temp( tempp, temps, ekinpr, vels, &
                                 na, nsp, hold, pmass, ndega )
     !
     ! ... fake electronic kinetic energy
     !
     IF ( .NOT. tcg ) THEN
        !
        CALL elec_fakekine( ekinc0, ema0bg, emass, c0, cm, ngw, nbsp, delt )
        !
        ekinc = ekinc0
        !
     END IF
     !
     ! ... fake cell-parameters kinetic energy
     !
     ekinh = 0.D0
     !
     IF ( thdyn ) CALL cell_kinene( ekinh, temphh, velh )
     !
     IF ( COUNT( iforceh == 1 ) > 0 ) THEN
        !
        temphc = 2.D0 * factem * ekinh / DBLE( COUNT( iforceh == 1 ) )
        !
     ELSE
        !
        temphc = 0.D0
        !
     END IF
     !
     ! ... udating nose-hoover friction variables
     !
     IF ( tnosep ) CALL ions_noseupd( xnhpp, xnhp0, xnhpm, delt, &
                                      qnp, ekinpr, gkbt, vnhp, kbt, nhpcl )
     !
     IF ( tnosee ) CALL electrons_noseupd( xnhep, xnhe0, xnhem, &
                                           delt, qne, ekinc, ekincw, vnhe )
     !
     IF ( tnoseh ) CALL cell_noseupd( xnhhp, xnhh0, xnhhm, &
                                      delt, qnh, temphh, temph, vnhh )
     !
     ! ... warning:  thdyn and tcp/tcap are not compatible yet!!!
     !
     IF ( tcp .OR. tcap .AND. tfor .AND. .NOT.thdyn ) THEN
        !
        IF ( tempp > temp1 .OR. tempp < temp2 .AND. tempp /= 0.D0 ) THEN
           !
           CALL  ions_vrescal( tcap, tempw, tempp, taup, &
                               tau0, taum, na, nsp, fion, iforce, pmass, delt )
           !
        END IF
        !
     END IF
     !
     IF( ( MOD(nfi-1,iprint) == 0 ) .OR. ( nfi == nomore ) ) THEN
        !
        CALL cp_eigs( nfi, bec, c0, irb, eigrb, rhor, &
                      rhog, rhos, lambdap, lambda, tau0, h )
        !
        IF( printwfc >= 0 ) CALL cp_print_rho( nfi, bec, c0, eigr, irb, &
                                               eigrb, rhor, rhog, rhos, &
                                               lambdap, lambda, tau0, h )
        !
     END IF
     !
     IF ( lwf ) CALL ef_enthalpy( enthal, tau0 )
     !
     IF ( tens ) THEN
        !
        IF ( MOD( nfi-1, iprint ) == 0 .OR. ( nfi == nomore ) ) THEN
           !
           WRITE( stdout, '("Occupations  :")' )
           WRITE( stdout, '(10F9.6)' ) ( f(i), i = 1, nbsp )
           !
        END IF
        !
     END IF
     !
     epot = eht + epseu + exc
     !
     acc(1) = acc(1) + ekinc
     acc(2) = acc(2) + ekin
     acc(3) = acc(3) + epot
     acc(4) = acc(4) + etot
     acc(5) = acc(5) + tempp
     !
     IF ( .NOT. tcg ) THEN
        !
        econs = ekinp + ekinh + enthal
        econt = econs + ekinc
        !
     ELSE
        !
        IF ( .NOT. tens ) THEN
           !
           econs = ekinp + etot
           atot  = etot
           econt = econs
           !
        ELSE
           !
           gibbsfe = atot
           econs   = ekinp + atot
           econt   = econs
           !
        END IF
        !
     END IF
     !
     ! ... add energies of thermostats
     !
     IF ( tnosep ) econt = econt + &
                           ions_nose_nrg( xnhp0, vnhp, qnp, gkbt, kbt, nhpcl )
     IF ( tnosee ) econt = econt + &
                           electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
     IF ( tnoseh ) econt = econt + &
                           cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
     !
     IF( ( MOD( nfi-1, iprint ) == 0 ) .OR. tfirst )  THEN
        !
        WRITE( stdout, * )
        WRITE( stdout, 1947 )
        !
     END IF
     !
     tps = tps + delt * au_ps
     !
     WRITE( stdout, 1948 ) nfi, ekinc, temphc, tempp, etot, enthal, econs, &
                           econt, vnhh(3,3), xnhh0(3,3), vnhp(1),  xnhp0(1)
     !
     IF( tcg ) THEN
        !
        IF ( MOD( nfi-1, iprint ) == 0 .OR. tfirst ) THEN
           !
           WRITE( stdout, * )
           WRITE( stdout, 255 ) 'nfi','tempp','E','-T.S-mu.nbsp','+K_p'
           !
        END IF
        !
        WRITE( stdout, 256 ) nfi, INT( tempp ), etot, atot, econs, itercg
        !
     END IF
     !
1947 FORMAT( 2X,'nfi',4X,'ekinc',2X,'temph',2X,'tempp',8X,'etot',6X,'enthal', &
          &  7X,'econs',7X,'econt',4X,'vnhh',3X,'xnhh0',4X,'vnhp',3X,'xnhp0' )
1948 FORMAT( I5,1X,F8.5,1X,F6.1,1X,F6.1,4(1X,F11.5),4(1X,F7.4) )
     !
     IF ( ionode .AND. ttprint ) THEN
        !
        ! ... Open units 30, 31, ... 40 for simulation output
        !
        CALL printout_base_open()
        !
        WRITE( stdout, 10 )
        !
        CALL printout_cell( stdout, nfi, hold, tps )
        CALL printout_cell( 36, nfi, hold, tps )
        !
        WRITE( stdout, 17 )
        !
        stress_gpa = stress * au_gpa
        !
        CALL printout_stress( stdout, nfi, stress_gpa, tps )
        CALL printout_stress( 38, nfi, stress_gpa, tps )
        !
        WRITE( stdout,11)
        !
        CALL printout_pos( stdout, nfi, tau0, nat, tps )
        CALL printout_pos( 35    , nfi, tau0, nat, tps )
        !
        ! ... write out a standard XYZ file in angstroms
        !
!        CALL print_pos_in( stdout, nfi, tau0 , nat, tps, ityp, atm,ind_bck,one)
!        WRITE( 35, '(I5)') nat
!        CALL print_pos_in( 35    , nfi, tau0 , nat, tps, ityp, atm,ind_bck,toang)

        !
        ALLOCATE( tauw( 3, natx ) )
        !
        isa = 0 
        !
        DO is = 1, nsp
           !
           DO ia = 1, na(is)
              !
              isa = isa + 1
              !
              CALL s_to_r( vels(:,isa), tauw(:,isa), hold )
              !
           END DO
           !
        END DO
        !
        WRITE( stdout, 12 )
        !
        CALL printout_pos( stdout, nfi, tauw, nat, tps )
        CALL printout_pos( 34    , nfi, tauw, nat, tps )
!        CALL print_pos_in( stdout, nfi, tauw , nat, tps, ityp, atm,ind_bck,one)
!        CALL print_pos_in( 34    , nfi, tauw , nat, tps, ityp, atm,ind_bck,one)
        !
        WRITE( stdout, 13 )
        !
        CALL printout_pos( stdout, nfi, fion, nat, tps )
        CALL printout_pos( 37    , nfi, fion, nat, tps )
!        CALL print_pos_in( stdout, nfi, fion , nat, tps, ityp, atm,ind_bck,one)
!        CALL print_pos_in( 37    , nfi, fion , nat, tps, ityp, atm,ind_bck,one)
        !
        DEALLOCATE( tauw )
        !
        WRITE( 33, 2948 ) tps, ekinc, temphc, tempp, etot, enthal, econs, econt
        WRITE( 39, 2949 ) tps, vnhh(3,3), xnhh0(3,3), vnhp(1), xnhp0(1)
        !
        ! ... Close and flush unit 30, ... 40
        !
        CALL printout_base_close()
        !
     END IF
     !
10   FORMAT( /,3X,'Cell Variables (AU)',/ )
11   FORMAT( /,3X,'Atomic Positions (AU)',/ )
12   FORMAT( /,3X,'Atomic Velocities (AU)',/ )
13   FORMAT( /,3X,'Atomic Forces (AU)',/ )
17   FORMAT( /,3X,'Total Stress (GPa)',/ )
255  FORMAT( '     ',5(1X,A12) )
256  FORMAT( 'Step ',I5,1X,I7,1X,F12.5,1X,F12.5,1X,F12.5,1X,I5 )
2948 FORMAT( F8.5,1X,F8.5,1X,F6.1,1X,F6.1,3(1X,F11.5) )
2949 FORMAT( F8.5,1X,4(1X,F7.4) )
     !
     IF( tfor ) THEN
        !
        ! ... new variables for next step
        !
        CALL ions_shiftvar( taup,  tau0, taum  )   !  real positions
        CALL ions_shiftvar( tausp, taus, tausm )   !  scaled positions         
        CALL ions_shiftvar( velsp, vels, velsm )   !  scaled velocities
        !
        IF ( tnosep ) CALL ions_nose_shiftvar( xnhpp, xnhp0, xnhpm )
        IF ( tnosee ) CALL electrons_nose_shiftvar( xnhep, xnhe0, xnhem )
        IF ( tnoseh ) CALL cell_nose_shiftvar( xnhhp, xnhh0, xnhhm )
        !
     END IF
     !
     IF ( thdyn ) CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
     !
     ekincm = ekinc0
     !  
     ! ... cm=c(t+dt) c0=c(t)
     !
     IF( .NOT. tcg ) THEN
        !
        CALL dswap( 2*ngw*nbsp, c0, 1, cm, 1 )
        !
     ELSE
        !
        CALL cg_update( tfirst, nfi, c0 )
        !
     END IF
     !
     ! ... now:  cm=c(t) c0=c(t+dt)
     !
     IF ( tfirst ) THEN
        !
        epre = etot
        enow = etot
        !
     END IF
     !
     tfirst = .FALSE.
     !
     ! ... write on file ndw each isave
     !
     IF ( ( MOD( nfi, isave ) == 0 ) .AND. ( nfi < nomore ) ) THEN
        !
        IF ( tcg ) THEN
           !
           CALL writefile( ndw, h, hold ,nfi, c0(:,:,1,1), c0old, taus, tausm, &
                           vels, velsm, acc, lambda, lambdam, xnhe0, xnhem,    &
                           vnhe, xnhp0, xnhpm, vnhp, nhpcl, ekincm, xnhh0,     &
                           xnhhm, vnhh, velh, ecutp, ecutw, delt, pmass, ibrav,&
                           celldm, fion, tps, z0, f )
           !
        ELSE
           !
           CALL writefile( ndw, h, hold, nfi, c0(:,:,1,1), cm(:,:,1,1), taus,  &
                           tausm, vels, velsm, acc,  lambda, lambdam, xnhe0,   &
                           xnhem, vnhe, xnhp0, xnhpm, vnhp, nhpcl, ekincm,     &
                           xnhh0, xnhhm, vnhh, velh, ecutp, ecutw, delt, pmass,&
                           ibrav, celldm, fion, tps, z0, f )
           !
        END IF
        !
     END IF
     !
     epre = enow
     enow = etot
     !
     frice = frice * grease
     fricp = fricp * greasp
     frich = frich * greash
     !
     !======================================================================
     !
     CALL stop_clock( 'total_time' )
     !
     delta_etot = ABS( epre - enow )
     !
     tstop = check_stop_now()
     tconv = .FALSE.
     !
     IF ( tconvthrs%active ) THEN
        !
        ! ... electrons
        !
        tconv = ( ekinc < tconvthrs%ekin .AND. delta_etot < tconvthrs%derho )
        !
        IF ( tfor ) THEN
           !
           ! ... ions
           !
           maxfion = MAXVAL( ABS( fion(:,1:nat) ) )
           !
           PRINT *, maxfion, tconvthrs%force
           !
           tconv = tconv .AND. ( maxfion < tconvthrs%force )
           !
        END IF
        !
     END IF
     !
     ! ... in the case cp-wf the check on convergence is done starting
     ! ... from the second step 
     !
     IF ( lwf .AND. tfirst ) tconv = .FALSE.
     !
     IF ( tconv ) THEN
        !
        IF ( ionode ) THEN
           !
           WRITE( stdout, &
                & "(/,3X,'MAIN:',10X,'EKINC   (thr)', &
                & 10X,'DETOT   (thr)',7X,'MAXFORCE   (thr)')" )
           WRITE( stdout, "(3X,'MAIN: ',3(D14.6,1X,D8.1))" ) &
               ekinc, tconvthrs%ekin, delta_etot,                  &
               tconvthrs%derho, 0.D0, tconvthrs%force
           WRITE( stdout, &
                  "(3X,'MAIN: convergence achieved for system relaxation')" )
           !
        END IF
        !
     END IF
     !
     tstop = tstop .OR. tconv
     !
     IF ( lwf ) &
        CALL wf_closing_options( nfi, c0, cm, bec, becdr, eigr, eigrb, taub, &
                                 irb, ibrav, b1, b2, b3, taus, tausm, vels,  &
                                 velsm, acc, lambda, lambdam, xnhe0, xnhem,  &
                                 vnhe, xnhp0, xnhpm, vnhp, nhpcl, ekincm,    &
                                 xnhh0, xnhhm, vnhh, velh, ecutp, ecutw,     &
                                 delt, celldm, fion, tps, z0, f )
     !
     IF ( ( nfi >= nomore ) .OR. tstop ) EXIT main_loop
     !
  END DO main_loop
  !
  !===================== end of main loop of molecular dynamics ===============
  ! 
  ! ... Here copy relevant physical quantities into the output arrays/variables
  !
  etot_out = etot
  !
  isa = 0
  !
  DO is = 1, nsp
     !
     DO ia = 1, na(is)
        !
        isa = isa + 1
        !
        ipos = ind_srt( isa )
        !
        tau(:,ipos) = tau0(:,isa)
        !
        fion_out(:,ipos) = fion(:,isa)
        !
     END DO
     !
  END DO
  !
  ! ...  Calculate statistics
  !
  acc = acc / DBLE( nfi )
  !
  IF ( ionode ) THEN
     !
     WRITE( stdout,1949)
     WRITE( stdout,1950) ( acc(i), i = 1, nacc )
     !
  END IF
  !
1949 FORMAT( //'              averaged quantities :',/,9X,&
           & 'ekinc',10X,'ekin',10X,'epot',10X,'etot',5X,'tempp' )
1950 FORMAT( 4F14.5,F10.1 )
  !
  CALL print_clock( 'initialize' )
  CALL print_clock( 'total_time' )
  CALL print_clock( 'formf' )
  CALL print_clock( 'rhoofr' )
  CALL print_clock( 'vofrho' )
  CALL print_clock( 'dforce' )
  CALL print_clock( 'calphi' )
  CALL print_clock( 'ortho' )
  CALL print_clock( 'updatc' )
  CALL print_clock( 'gram' )
  CALL print_clock( 'newd' )
  CALL print_clock( 'calbec' )
  CALL print_clock( 'prefor' )
  CALL print_clock( 'strucf' )
  CALL print_clock( 'nlfl' )
  CALL print_clock( 'nlfq' )
  CALL print_clock( 'set_cc' )
  CALL print_clock( 'rhov' )
  CALL print_clock( 'nlsm1' )
  CALL print_clock( 'nlsm2' )
  CALL print_clock( 'forcecc' )
  CALL print_clock( 'fft' )
  CALL print_clock( 'ffts' )
  CALL print_clock( 'fftw' )
  CALL print_clock( 'fftb' )
  CALL print_clock( 'rsg' )
  CALL print_clock( 'reduce' )
  !
  IF ( tcg ) THEN
     !
     CALL writefile( ndw, h, hold, nfi, c0(:,:,1,1), c0old, taus, tausm, vels, &
                     velsm, acc, lambda, lambdam, xnhe0, xnhem, vnhe, xnhp0,   &
                     xnhpm, vnhp, nhpcl, ekincm, xnhh0, xnhhm, vnhh, velh,     &
                     ecutp, ecutw, delt, pmass, ibrav, celldm, fion, tps,      &
                     z0, f )
     !
  ELSE
     !
     CALL writefile( ndw, h, hold, nfi, c0(:,:,1,1), cm(:,:,1,1), taus, tausm, &
                     vels, velsm, acc, lambda, lambdam, xnhe0, xnhem, vnhe,    &
                     xnhp0, xnhpm, vnhp, nhpcl, ekincm, xnhh0, xnhhm, vnhh,    &
                     velh, ecutp, ecutw, delt, pmass, ibrav, celldm, fion, tps,&
                     z0, f )
     !
  END IF
  !
  IF( iprsta > 1 ) CALL print_lambda( lambda, nbsp, nbsp, 1.D0 )
  !
  conv_elec = .TRUE.
  !
1974 FORMAT( 1X,2I5,3F10.4,2X,3F10.4 )
1975 FORMAT( /1X,'Scaled coordinates '/1X,'species',' atom #' )
1976 FORMAT( 1X,2I5,3F10.4 )
  !
  IF( ionode ) WRITE( stdout, 1977 ) 
  !
  CALL memory()
  !      
1977 FORMAT(5X,//'====================== end cprvan ======================',//)
  !
  CALL deallocate_modules_var()
  !
  RETURN
  !
END SUBROUTINE cprmain
