!
! Copyright (C) 2002-2004 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

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
!***ensemble-DFT
! cf. "Ensemble Density-Functional Theory for Ab Initio Molecular Dynamics
!      of Metals and Finite-Temperature Insulators"  PRL v.79,n.7 (1997)
!      N. Marzari, D. Vanderbilt and M.C. Payne
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
      use control_flags, only: iprint, isave, thdyn, tpre, tbuff, iprsta, trhor, &
            tfor, tvlocw, trhow, taurdr, tprnfor, tsdc
      use control_flags, only: ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
            tnosep, trane, tranp, tsdp, tcp, tcap, ampre, amprp, tnoseh, tolp
      use control_flags, only: lwf, ortho_eps, ortho_max, printwfc

      use core, only: nlcc_any
      use core, only: deallocate_core
      use uspp_param, only: nhm, nh
      use cvan, only: nvb, ish, deallocate_cvan
      use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, &
           deeq
      use uspp, only: deallocate_uspp
      use energies, only: eht, epseu, exc, etot, eself, enl, ekin, &
                          atot, entropy, egrand

      use electrons_base, only: nx => nbspx, n => nbsp, ispin => fspin, f, nspin
      use electrons_base, only: deallocate_elct, nel, iupdwn, nupdwn, nudx, nelt

      use efield_module, only: efield, epol, tefield, allocate_efield, deallocate_efield, &
            efield_update, ipolp, qmat, gqq, evalue, berry_energy

      use ensemble_dft, only: tens, tgrand, ninner, ismear, etemp, ef,              &
            tdynz, tdynf, zmass, fmass, fricz, fricf, allocate_ensemble_dft,        &
            deallocate_ensemble_dft, id_matrix_init, z0, c0diag, becdiag, bec0, &
            v0s, vhxcs, becdrdiag, compute_entropy2, gibbsfe

      use cg_module, only: tcg, maxiter, etresh, passop, allocate_cg, deallocate_cg, &
            cg_update, itercg, c0old


      use gvec, only: tpiba2, ng
      use gvec, only: deallocate_gvec
      use gvecs, only: ngs
      use gvecb, only: ngb
      use gvecw, only: ngw
      use reciprocal_vectors, only: gstart
      use ions_base, only: na, nat, pmass, nax, nsp, rcmax
      use ions_base, only: ind_srt, ions_cofmass, ions_kinene, ions_temp, ions_thermal_stress
      use ions_base, only: ions_vrescal, fricp, greasp, iforce, ions_shiftvar
      use grid_dimensions, only: nnr => nnrx, nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3, frich, greash
      use cell_base, only: omega, alat, ibrav, celldm
      use cell_base, only: h, hold, deth, wmass, press
      use cell_base, only: s_to_r, r_to_s
      use cell_base, only: iforceh, cell_force, thdiag
      use smooth_grid_dimensions, only: nnrsx, nr1s, nr2s, nr3s
      use smallbox_grid_dimensions, only: nnrb => nnrbx, nr1b, nr2b, nr3b
      use pseu, only: vps, rhops
      use pseu, only: deallocate_pseu
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
      use restart_file
      use parameters, only: nacx, natx, nsx, nbndxx
      use constants, only: pi, factem, au_gpa, au_ps, gpa_au
      use io_files, only: psfile, pseudo_dir
      use input, only: iosys
      use qgb_mod, only: deallocate_qgb_mod
      use dqgb_mod, only: deallocate_dqgb_mod
      use qradb_mod, only: deallocate_qradb_mod
      use dqrad_mod, only: deallocate_dqrad_mod
      use betax, only: deallocate_betax
      use wave_base, only: wave_steepest, wave_verlet
      use wave_base, only: wave_speed2, frice, grease
      USE control_flags, ONLY : conv_elec, tconvthrs
      USE check_stop, ONLY : check_stop_now
      use efcalc, ONLY: clear_nbeg, ef_force    !Electric Field (M.S)
      use ions_base, only: zv, ions_vel           !
      use cp_electronic_mass, only: emass, emaec => emass_cutoff
      use cp_electronic_mass, only: emass_precond
      use cpr_subroutines
      use ions_positions, only: tau0, taum, taup, taus, tausm, tausp, vels, velsm, velsp
      use ions_positions, only: ions_hmove, ions_move
      use ions_nose, only: gkbt, kbt, ndega, nhpcl, qnp, vnhp, xnhp0, xnhpm, xnhpp, &
            ions_nosevel, ions_noseupd, tempw, ions_nose_nrg, ions_nose_shiftvar
      use electrons_nose, only: qne, ekincw, xnhe0, xnhep, xnhem, vnhe, &
            electrons_nose_nrg, electrons_nose_shiftvar, electrons_nosevel, electrons_noseupd
      use from_scratch_module, only: from_scratch
      use from_restart_module, only: from_restart

! wavefunctions
!
      use wavefunctions_module, only: c0, cm, phi => cp
      use wavefunctions_module, only: deallocate_wavefunctions

      use wannier_module, only: allocate_wannier, deallocate_wannier
      use wannier_subroutines
      USE printout_base, ONLY: printout_base_open, printout_base_close, &
            printout_pos, printout_cell, printout_stress
      USE cell_nose, ONLY: xnhh0, xnhhm, xnhhp, vnhh, temph, qnh, &
            cell_nosevel, cell_noseupd, cell_nose_nrg, cell_nose_shiftvar
      USE cell_base, ONLY: cell_kinene, cell_gamma, cell_move, cell_hmove
      USE gvecw, ONLY: ecutw
      USE gvecp, ONLY: ecutp

      USE time_step, ONLY: delt
      USE runcp_module, ONLY: runcp_uspp
      USE wave_constrains, ONLY: interpolate_lambda
      USE electrons_module, ONLY: cp_eigs
      USE print_out_module, ONLY: cp_print_rho

      USE cp_main_variables

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
! control variables
!
      logical :: tfirst, tlast, tstop, tconv
      logical :: ttprint    !  logical variable used to control printout
!
!  ionic positions, center of mass position
!
      real(kind=8) :: cdm0(3)
      real(kind=8) :: cdm(3)
!
!  forces on ions
!
      real(kind=8) :: fion(3,natx), fionm(3,natx)
!
! work variables
!
      real(kind=8)                                                      & 
     &       tempp,  fccc, savee, saveh, savep,             &
     &       enthal, epot, epre, enow, tps, econs, econt, &
     &       ettt, ccc, bigr, dt2, dt2by2, twodel, dt2bye, dt2hbe
      real(kind=8) ekinc0, ekinp, ekinpr, ekincm, ekinc
      real(kind=8) temps(nsx)
      real(kind=8) ekinh, temphc, temp1, temp2, randy
      real(kind=8) :: delta_etot
      real(kind=8) ftmp, enb, enbi
      integer :: is, nacc, ia, j, iter, nfi, i, isa, ipos
      integer :: k, ii, l, m, ibeg
!
! work variables, 2
!
      real(kind=8) hnew(3,3), velh(3,3), hgamma(3,3), temphh(3,3)
      real(kind=8) fcell(3,3)
!

      real(kind=8) :: b1(3), b2(3), b3(3)
      real(kind=8) :: stress_gpa(3,3), thstress(3,3)

      !   temporary array used to printout positions
      real(kind=8), allocatable :: tauw( :, : )  

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
!     CP starts here
!

      call start_clock( 'initialize' )

      etot_out = 0.0d0

      tps     = 0.0d0

!     general variables
!
      tfirst = .true.
      tlast  = .false.
      nacc = 5
!
!     ==================================================================
!     read input from standard input (unit 5)
!     ==================================================================

      call iosys( )

      if( lwf ) then
        call read_efwan_param( nbeg )
      end if

!     ==================================================================
!
      twodel = 2.d0 * delt
      dt2    = delt * delt
      dt2by2 = .5d0 * dt2
      dt2bye = dt2/emass
      dt2hbe = dt2by2/emass

!
!     ==================================================================
!     initialize g-vectors, fft grids
!     ==================================================================

      call init_dimensions( )

      call init1 ( tau0 )

      call init( ibrav, celldm, ecutp, ecutw, ndr, nbeg, tfirst,  &
           tau0, taus, delt, tps, iforce )

      WRITE( stdout,*) ' out from init'

      if( lwf ) then
        call clear_nbeg( nbeg )
      end if
!
!     more initialization requiring atomic positions
!
      nax = MAXVAL( na( 1 : nsp ) )
      if( iprsta > 1 ) then
         call print_atomic_var( tau0, na, nsp, ' tau0 ' )
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
      CALL allocate_mainvar &
         ( ngw, ngb, ngs, ng, nr1, nr2, nr3, nnr, nnrsx, nax, nsp, nspin, n, nx, nhsa, nlcc_any )

      allocate(c0(ngw,nx,1,1))
      allocate(cm(ngw,nx,1,1))
      allocate(phi(ngw,nx,1,1))
      allocate(rhops(ngs,nsp))
      allocate(vps(ngs,nsp))
      allocate(betae(ngw,nhsa))
      allocate(deeq(nhm,nhm,nat,nspin))
      allocate(rhovan(nhm*(nhm+1)/2,nat,nspin))
      allocate(dbec (nhsa,n,3,3))
      allocate(dvps(ngs,nsp))
      allocate(drhops(ngs,nsp))
      allocate(drhog(ng,nspin,3,3))
      allocate(drhor(nnr,nspin,3,3))
      allocate(drhovan(nhm*(nhm+1)/2,nat,nspin,3,3))

      if( lwf ) then
        call allocate_wannier(  n, nnrsx, nspin, ng )
      end if

      IF( tens ) THEN
        CALL allocate_ensemble_dft( nhsa, n, ngw, nudx, nspin, nx, nnrsx, nat )
      END IF

      IF( tcg ) THEN
        CALL allocate_cg( ngw, nx)
      END IF

      IF( tefield ) THEN
        CALL allocate_efield( ngw, nx, nhm, nax, nsp )
      END IF

      deeq(:,:,:,:) = 0.d0
      !
      !
 666  continue
      !
      !
      temp1=tempw+tolp
      temp2=tempw-tolp

!     ==========================================================

      taum  = tau0
      taup  = 0.0d0
      tausm = taus
      tausp = 0.0d0
      vels  = 0.0d0
      velsm = 0.0d0
      velsp = 0.0d0
!
      hnew=h
!
      lambda(:,:)=0.d0
      cm(:,:,1,1) = (0.d0, 0.d0)
      c0(:,:,1,1) = (0.d0, 0.d0)

      IF( tens ) THEN
        CALL id_matrix_init( nupdwn, nspin )
      END IF
!
      CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emaec )

      if( lwf ) then
        call wannier_init( ibrav, alat, a1, a2, a3, b1, b2, b3 )
      end if
!
      if ( nbeg < -1 ) then

         !======================================================================
         !    Initialize from scratch nbeg = -2 or nbeg = -3
         !======================================================================
!
         nfi = 0

         CALL from_scratch &
            ( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, eself, fion, &
              taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, rhoc, enl, ekin, stress,  &
              detot, enthal, etot, lambda, lambdam, lambdap, ema0bg, dbec, delt,  &
              bephi, becp, velh, dt2bye, iforce, fionm, nbeg, xnhe0, xnhem, vnhe, ekincm )
!
      else

!======================================================================
!        nbeg = -1, nbeg = 0, nbeg = 1 or nbeg = 2
!======================================================================

         IF( nbeg == -1 ) THEN
           ! ibeg = 0  !  read only wavefunction cm from restart
           ibeg = 1  !  DEBUG
         ELSE
           ibeg = 1  !  read all dynamic variable from restart
         END IF

         call readfile                                           &
     &     ( ibeg, ndr,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc,         &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,fion, tps, z0, f )
!

         call from_restart     &
            ( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, fion, &
              taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, rhoc, stress,  &
              detot, enthal, lambda, lambdam, lambdap, ema0bg, dbec, &
              bephi, becp, velh, dt2bye, fionm, ekincm )
!
      end if
!==============================================end of if(nbeg.lt.0)====
!
!     =================================================================
!     restart with new averages and nfi=0
!     =================================================================
!       Fix. Center of Mass - M.S
!
      if( lwf ) then
        call ions_cofmass(tau0, pmass, na, nsp, cdm0)
      end if

      if( nbeg <= 0 ) then
         acc = 0.0d0
         nfi = 0
      end if
!
      if( ( .not. tfor ) .and. ( .not. tprnfor ) ) then
         fion = 0.d0
      end if
!
      if( .not. tpre ) then
         stress = 0.d0
      endif
!         
      fccc = 1.0d0
      !
      nomore = nomore + nfi
!
      call ions_cofmass(taus, pmass, na, nsp, cdm0)
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
      if( .not. tsde ) fccc = 1.0d0 / ( 1.0d0 + frice )

      if(tnosep)then
         call ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl )
      endif
      if(tnosee)then
         call electrons_nosevel( vnhe, xnhe0, xnhem, delt )
         fccc=1./(1.+0.5*delt*vnhe)
      endif
      if(tnoseh) then
         call cell_nosevel( vnhh, xnhh0, xnhhm, delt, velh, h, hold )
      endif
      ! 
      if ( tfor .or. thdyn .or. tfirst ) then 
         call initbox ( tau0, taub, irb )
         call phbox( taub, eigrb )
      endif
      !
      if( tfor .or. thdyn ) call phfac(tau0,ei1,ei2,ei3,eigr) 
      !
      !     strucf calculates the structure factor sfac
      !
      call strucf(ei1,ei2,ei3,sfac)
      if (thdyn) call formf(tfirst,eself)

      if( tefield ) then
        ! why this call ??? from Paolo Umary
        call calbec(1,nsp,eigr,c0,bec) ! ATTENZIONE  
      end if
!
      nfi = nfi + 1
      tlast   = nfi == nomore
      ttprint = MOD(nfi, iprint) == 0


      if( ( tfor .or. tfirst ) .and. tefield ) call efield_update( eigr )


      ELECTRON_DYNAMIC: if ( tcg ) then

        call runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
             rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, fion, ema0bg, becdr, &
             lambdap, lambda  )

      else

        if( lwf ) then
          call get_wannier_center( tfirst, cm, bec, becdr, eigr, eigrb, taub, irb, ibrav, b1, b2, b3 )
        end if
!
        if( .not. tens ) then
          call rhoofr (nfi,c0,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
        else
          !     calculation of the rotated quantities
          call rotate(z0,c0(:,:,1,1),bec,c0diag,becdiag)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
          bec0(:,:)=bec(:,:)
        end if
!
#ifdef __PARA     
        if(trhow .and. tlast) call write_rho(47,nspin,rhor)
#else
        if(trhow .and. tlast) write(47) ((rhor(i,is),i=1,nnr),is=1,nspin)
#endif
        !
        !     put core charge (if present) in rhoc(r)
        !
        if ( nlcc_any ) call set_cc(irb,eigrb,rhoc)

        if( lwf ) then
          call write_charge_and_exit( rhog )
          call ef_tune( rhog, tau0 )
        end if
  !
        if( .not. tens ) then
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                 &
                 ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        else
          if( tfirst ) then
            call compute_entropy2( entropy, f, n, nspin )
          end if
          call vofrho2(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                &
                 ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion,v0s,vhxcs)
        end if

        if( lwf ) then
          call wf_options( tfirst, nfi, cm, rhovan, bec, becdr, eigr, eigrb, taub, irb, &
               ibrav, b1, b2, b3, rhor, rhog, rhos, enl, ekin  )
        end if

        call compute_stress( stress, detot, h, omega )
!
        enthal = etot + press * omega

        if( tefield .and. ( ( mod(nfi-1,10) == 0 ) .or. tfirst ) ) then
          call berry_energy( enb, enbi, bec, c0( :, :, 1, 1), fion )
          etot = etot + enb + enbi
        end if
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
        call newd(rhor,irb,eigrb,rhovan,fion)
        call prefor(eigr,betae)

        CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0(:,:,1,1), cm(:,:,1,1) )

        !
        !----------------------------------------------------------------------
        !                 contribution to fion due to lambda
        !----------------------------------------------------------------------
        !
        !     nlfq needs deeq bec
        !
        if ( .not. tens ) then
          if ( tfor .or. tprnfor ) call nlfq(c0,eigr,bec,becdr,fion)
        else
          if ( tfor .or. tprnfor ) call nlfq(c0diag,eigr,becdiag,becdrdiag,fion)
        end if 

        if( tfor .and. tefield ) then
           call bforceion( fion, tfor, ipolp, qmat, bec, becdr, gqq, evalue )
        endif
!
        if( tfor .or. thdyn ) then
          CALL interpolate_lambda( lambdap, lambda, lambdam )
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
        if ( .not. tens ) then
          if ( tfor .or. tprnfor ) call nlfl( bec, becdr, lambda, fion )
        else
          if (tfor .or. tprnfor) then
            call nlsm2( eigr, c0(:,:,1,1), becdr )
            call nlfl( bec, becdr, lambda, fion )
          end if
        end if

      end if ELECTRON_DYNAMIC

!
!
      if(tpre) then
         call nlfh(bec,dbec,lambda)
         call ions_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
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
      !
      if(thdyn) then
         !
         call cell_force( fcell, ainv, stress, omega, press, wmass )
         !
         ! WRITE(6,*) 'DEBUG fcell = ', fcell
         ! WRITE(6,*) 'DEBUG stress = ', stress
         ! WRITE(6,*) 'DEBUG omega = ', omega
         ! WRITE(6,*) 'DEBUG press = ', press
         ! WRITE(6,*) 'DEBUG wmass = ', wmass
         !
         call cell_move( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, vnhh, velh, tsdc )
         !
         velh(:,:) = (hnew(:,:)-hold(:,:))/twodel
         !
         call cell_gamma( hgamma, ainv, h, velh )
         !
      endif
      !
      !======================================================================
      !
      if( tfor ) then

        if( lwf ) then
          call ef_force( fion, na, nsp, zv )
        end if

        ! WRITE(6,*) 'DEBUG fion = '
        ! WRITE(6,fmt='(3F20.14)') fion(:,1:nat)
        ! WRITE(6,*) 'DEBUG iforce = '
        ! WRITE(6,fmt='(3I4)') iforce(:,1:nat)

        call ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp, &
                        fricp, hgamma, vels, tsdp, tnosep, fionm, vnhp, velsp, velsm )

        ! WRITE(6,*) 'DEBUG tausp = '
        ! WRITE(6,fmt='(3F20.14)') tausp(:,1:nat)
        ! WRITE(6,*) 'DEBUG taus  = '
        ! WRITE(6,fmt='(3F20.14)') tausp(:,1:nat)
        ! WRITE(6,*) 'DEBUG tausm = '
        ! WRITE(6,fmt='(3F20.14)') tausm(:,1:nat)

        !cc   call cofmass(velsp,cdmvel)
        !cc   velsp(i,isa)=velsp(i,isa)-cdmvel(i)

        call ions_cofmass(tausp, pmass, na, nsp, cdm)

        call ions_cofmsub( tausp, na, nsp, cdm, cdm0 )

        CALL s_to_r( tausp, taup, na, nsp, hnew )

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
      if ( .not. tcg ) then
        if(tortho) then
          call ortho(eigr,cm,phi,lambda,bigr,iter,ccc,ortho_eps,ortho_max,delt,bephi,becp)
        else
          call graham(betae,bec,cm)
          if(iprsta.gt.4) call dotcsc(eigr,cm)
        endif
      end if
!
      !---------------------------------------------------------------------------
      !                   correction to displacement of ions
      !---------------------------------------------------------------------------
!
      if( .not. tcg ) then
        if( iprsta >= 3 ) CALL print_lambda( lambda, n, 9, 1.0d0 )
!
        if(tortho) call updatc(ccc,lambda,phi,bephi,becp,bec,cm)
        call calbec (nvb+1,nsp,eigr,cm,bec)
        if (tpre) call caldbec(1,nsp,eigr,cm)
!
        if(iprsta.ge.3)  call dotcsc(eigr,cm)
      end if
!
      !---------------------------------------------------------------------------
      !                  temperature monitored and controlled
      !---------------------------------------------------------------------------
!
      ekinp=0.0
      ekinpr=0.0
      tempp=0.0

      !
      !     ionic kinetic energy 
      !
      if( tfor ) then
         CALL ions_vel( vels, tausp, tausm, na, nsp, delt )
         CALL ions_kinene( ekinp, vels, na, nsp, hold, pmass )
      endif
      !
      !     ionic temperature
      !
      if( tfor ) then
         CALL ions_temp( tempp, temps, ekinpr, vels, na, nsp, hold, pmass, ndega )
      endif
      !
      !     fake electronic kinetic energy
      !
      if( .not. tcg ) then
        call elec_fakekine( ekinc0, ema0bg, emass, c0, cm, ngw, n, delt )
        ! ekinc = 0.5 * ( ekinc0 + ekincm )  !  what ???
        ekinc = ekinc0
      end if
      !
      !     fake cell-parameters kinetic energy
      !
      ekinh=0.
      if(thdyn) then
         call cell_kinene( ekinh, temphh, velh )
      endif
      if( COUNT( iforceh == 1 ) > 0 ) then
         temphc = 2.0d0 * factem * ekinh / DBLE( COUNT( iforceh == 1 ) )
      else
         temphc = 0.0d0 
      endif
      !
      !     udating nose-hoover friction variables
      !
      if(tnosep)then
        call ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, ekinpr, gkbt, vnhp, kbt, nhpcl )
      endif
      if(tnosee)then
        call electrons_noseupd( xnhep, xnhe0, xnhem, delt, qne, ekinc, ekincw, vnhe )
      endif
      if(tnoseh)then
        call cell_noseupd( xnhhp, xnhh0, xnhhm, delt, qnh, temphh, temph, vnhh )
      endif
      !
      ! warning! thdyn and tcp/tcap are not compatible yet!!!
      !
      if(tcp.or.tcap.and.tfor.and.(.not.thdyn)) then
         if(tempp.gt.temp1.or.tempp.lt.temp2.and.tempp.ne.0.d0) then
            call  ions_vrescal( tcap, tempw, tempp, taup, tau0, taum, na, nsp, fion, iforce, pmass, delt )
         end if
      end if
!
      if( ( mod(nfi-1,iprint) .eq. 0 ) .or. ( nfi .eq. nomore ) ) then
        CALL cp_eigs( nfi, bec, c0, irb, eigrb, rhor, rhog, rhos, lambdap, lambda, tau0, h )
        if( printwfc >= 0 ) then
          CALL cp_print_rho(nfi, bec, c0, eigr, irb, eigrb, rhor, rhog, rhos, lambdap, lambda, tau0, h )
        end if
      endif
!
      if( lwf ) then
        call ef_enthalpy( enthal, tau0 )
      end if

      if(tens) then
         if(mod(nfi-1,iprint).eq.0 .or. (nfi.eq.(nomore))) then
            write(stdout,'("Occupations  :")')
            write(stdout,'(10f9.6)') (f(i),i=1,n)
         endif
      endif

      epot = eht + epseu + exc
!
      acc(1)=acc(1)+ekinc
      acc(2)=acc(2)+ekin
      acc(3)=acc(3)+epot
      acc(4)=acc(4)+etot
      acc(5)=acc(5)+tempp
!

      if ( .not. tcg ) then
        econs = ekinp + ekinh + enthal
        econt = econs + ekinc
      else
        if (.not.tens) then
          econs=ekinp+etot
          atot=etot
          econt=econs
        else
          gibbsfe=atot
          econs=ekinp+atot
          econt=econs
        endif
      end if
      !
      !   add energies of thermostats
      !
      if ( tnosep ) econt = econt + ions_nose_nrg( xnhp0, vnhp, qnp, gkbt, kbt, nhpcl )
      if ( tnosee ) econt = econt + electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
      if ( tnoseh ) econt = econt + cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
!
      if( ( MOD( nfi-1, iprint ) == 0 ) .or. tfirst )  then
         WRITE( stdout,*)
         WRITE( stdout,1947)
      end if
!
      tps = tps + delt * AU_PS
      WRITE( stdout,1948) nfi, ekinc, temphc, tempp, etot, enthal, econs,      &
     &              econt, vnhh(3,3), xnhh0(3,3), vnhp(1),  xnhp0(1)

      if( tcg ) then
        if(mod(nfi-1,iprint).eq.0.or.tfirst)  then
          WRITE( stdout,*)
          WRITE( stdout,255) 'nfi','tempp','E','-T.S-mu.N','+K_p'
        end if
        WRITE(stdout,256)  nfi,int(tempp),etot,atot,econs, itercg
      end if

 1947 format(2x,'nfi',4x,'ekinc',2x,'temph',2x,'tempp',8x,'etot',6x,'enthal',     &
     &       7x,'econs',7x,'econt',4x,'vnhh',3x,'xnhh0',4x,'vnhp',      &
     &       3x,'xnhp0')
 1948 format(i5,1x,f8.5,1x,f6.1,1x,f6.1,4(1x,f11.5),4(1x,f7.4))
!

      IF ( ionode .AND. ttprint ) THEN
        ! ...  Open units 30, 31, ... 40 for simulation output
        CALL printout_base_open()
        WRITE( stdout, 10 )
        CALL printout_cell( stdout, nfi, hold, tps )
        CALL printout_cell( 36, nfi, hold, tps )
        WRITE( stdout, 17 )
        stress_gpa = stress * au_gpa
        CALL printout_stress( stdout, nfi, stress_gpa, tps )
        CALL printout_stress( 38, nfi, stress_gpa, tps )

        WRITE( stdout,11)
        CALL printout_pos( stdout, nfi, tau0     , nat, tps )
        CALL printout_pos( 35    , nfi, tau0     , nat, tps )

        ALLOCATE( tauw( 3, natx ) )
        isa = 0 
        DO is = 1, nsp
          DO ia = 1, na(is)
            isa = isa + 1
            CALL s_to_r( vels(:,isa), tauw(:,isa), hold )
          END DO
        END DO
        WRITE( stdout, 12 )
        CALL printout_pos( stdout, nfi, tauw     , nat, tps )
        CALL printout_pos( 34    , nfi, tauw     , nat, tps )

        WRITE( stdout, 13 )
        CALL printout_pos( stdout, nfi, fion     , nat, tps )
        CALL printout_pos( 37    , nfi, fion     , nat, tps )

        DEALLOCATE( tauw )

        WRITE( 33, 2948 ) tps, ekinc, temphc, tempp, etot, enthal, econs, econt
        WRITE( 39, 2949 ) tps, vnhh(3,3),xnhh0(3,3),vnhp(1),xnhp0(1)

        ! ...   Close and flush unit 30, ... 40
        CALL printout_base_close()
      END IF

 10  FORMAT(/,3X,'Cell Variables (AU)',/)
 11  FORMAT(/,3X,'Atomic Positions (AU)',/)
 12  FORMAT(/,3X,'Atomic Velocities (AU)',/)
 13  FORMAT(/,3X,'Atomic Forces (AU)',/)
 17  FORMAT(/,3X,'Total Stress (GPa)',/)
  255 format('     ',5(1x,a12  ))
  256 format('Step ',i5,1x,i7,1x,f12.5,1x,f12.5,1x,f12.5,1x,i5)
 2948 format(f8.5,1x,f8.5,1x,f6.1,1x,f6.1,3(1x,f11.5))
 2949 format(f8.5,1x,4(1x,f7.4))



      if( tfor ) then
         !
         !     new variables for next step
         !
         call ions_shiftvar( tausp, taus, tausm )   !  scaled positions 
         call ions_shiftvar( taup,  tau0,  taum )   !  real positions
         call ions_shiftvar( velsp, vels, velsm )   !  scaled velocities
         !
         if( tnosep ) call ions_nose_shiftvar( xnhpp, xnhp0, xnhpm )
         if( tnosee ) call electrons_nose_shiftvar( xnhep, xnhe0, xnhem )
         if( tnoseh ) call cell_nose_shiftvar( xnhhp, xnhh0, xnhhm )
         !
      end if
!
      if(thdyn)then
         CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emaec )
      endif
!
      ekincm=ekinc0
!  
!     cm=c(t+dt) c0=c(t)
!
      if( .not. tcg ) then
        call DSWAP(2*ngw*n,c0,1,cm,1)
      else
        call cg_update( tfirst, nfi, c0 )
      end if
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
!     write on file ndw each isave
!
      if( ( mod( nfi, isave ) == 0 ) .and. ( nfi < nomore ) ) then

         ! WRITE(6,*) 'DEBUG write taus  = ', taus(:,1)
         ! WRITE(6,*) 'DEBUG write tausm = ', tausm(:,1)

         if( tcg ) then
           call writefile                                         &
     &       ( ndw,h,hold,nfi,c0(:,:,1,1),c0old,taus,tausm,vels,velsm,acc,               &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,fion, tps, z0, f )
         else
           call writefile                                         &
     &       ( ndw,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc,               &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,fion, tps, z0, f )
         end if

      endif
!
      epre = enow
      enow = etot

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

      if( lwf ) then
        CALL wf_closing_options( nfi, c0, cm, bec, becdr, eigr, eigrb, taub, irb, &
             ibrav, b1, b2, b3, taus, tausm, vels, velsm, acc, lambda, lambdam, xnhe0, &
             xnhem, vnhe, xnhp0, xnhpm, vnhp, nhpcl, ekincm, xnhh0, xnhhm, vnhh, velh, &
             ecutp, ecutw, delt, celldm, fion, tps, z0, f )
      end if

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
        tau( 1, ipos ) = tau0( 1, isa )
        tau( 2, ipos ) = tau0( 2, isa )
        tau( 3, ipos ) = tau0( 3, isa )
        fion_out( 1, ipos ) = fion( 1, isa )
        fion_out( 2, ipos ) = fion( 2, isa )
        fion_out( 3, ipos ) = fion( 3, isa )
      end do
    end do

    !  Calculate statistics

    acc = acc / dble( nfi )
!
    if( ionode ) then
      WRITE( stdout,1949)
      WRITE( stdout,1950) (acc(i),i=1,nacc)
    end if

 1949 format(//'              averaged quantities :',/,                 &
     &       9x,'ekinc',10x,'ekin',10x,'epot',10x,'etot',5x,'tempp')
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
    call print_clock( 'reduce' )

!
    ! WRITE(6,*) 'DEBUG write taus  = ', taus(:,1)
    ! WRITE(6,*) 'DEBUG write tausm = ', tausm(:,1)

    if( tcg ) then
      call writefile &
          ( ndw,h,hold,nfi,c0(:,:,1,1),c0old,taus,tausm,vels,velsm,acc, &
            lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,ekincm,   &
            xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,fion, tps, z0, f )
    else
      call writefile &
          ( ndw,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc, &
            lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,ekincm,   &
            xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,fion, tps, z0, f )
    end if

!
    if( iprsta > 1 ) CALL print_lambda( lambda, n, n, 1.0d0 )
!
    conv_elec = .TRUE.


! 1970 format(1x,'ibrav :',i4,'  alat : ',f10.4,/)
! 1971 format(1x,'lattice vectors',/)
! 1972 format(1x,3f10.4)
! 1973 format(/1x,'Cartesian coordinates (a.u.)              forces' &
!     &       /1x,'species',' atom #', &
!     &           '   x         y         z      ', &
!     &           '   fx        fy        fz'/)
 1974 format(1x,2i5,3f10.4,2x,3f10.4)
 1975 format(/1x,'Scaled coordinates '/1x,'species',' atom #')
 1976 format(1x,2i5,3f10.4)

      if( ionode ) WRITE( stdout,1977) 
!
      call memory
!      
 1977 format(5x,//'====================== end cprvan ======================',//)

      CALL deallocate_mainvar()
      CALL deallocate_cvan()
      CALL deallocate_efield( )
      CALL deallocate_ensemble_dft()
      CALL deallocate_cg( )
      CALL deallocate_elct()
      CALL deallocate_core()
      CALL deallocate_uspp()
      CALL deallocate_gvec()
      CALL deallocate_pseu()
      CALL deallocate_qgb_mod()
      CALL deallocate_qradb_mod()
      CALL deallocate_derho()
      CALL deallocate_dqgb_mod()
      CALL deallocate_dpseu()
      CALL deallocate_cdvan()
      CALL deallocate_dqrad_mod()
      CALL deallocate_betax()
      CALL deallocate_para_mod()
      CALL deallocate_wavefunctions()
      if( lwf ) then
        CALL deallocate_wannier()
      end if

!
      return
      end subroutine
