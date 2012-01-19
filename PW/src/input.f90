
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE iosys()
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine reads input data from standard input ( unit 5 )
  ! ...  Use "-input filename" to read input from file "filename":
  ! ...  may be useful if you have trouble reading from standard input
  ! ...  ---------------------------------------------------------------
  !
  ! ...  access the modules renaming the variables that have the same name
  ! ...  as the input parameters, this is required in order to use a code
  ! ...  independent input parser
  !
  !
  USE kinds,         ONLY : DP
  USE funct,         ONLY : dft_has_finite_size_correction, &
                            set_finite_size_volume, get_inlc 
#if defined(EXX)
  USE funct,         ONLY: set_exx_fraction, set_screening_parameter
  USE control_flags, ONLY: adapt_thr, tr2_init, tr2_multi
#endif
  USE constants,     ONLY : autoev, eV_to_kelvin, pi, rytoev, &
                            ry_kbar, amconv, bohr_radius_angs, eps8
  USE mp_global,     ONLY : npool, nproc_pool
  !
  USE io_global,     ONLY : stdout, ionode, ionode_id
  !
  USE kernel_table,  ONLY : initialize_kernel_table
  !
  USE mp,            ONLY : mp_bcast
  !
  USE bp,            ONLY : nppstr_    => nppstr, &
                            gdir_      => gdir, &
                            lberry_    => lberry, &
                            lelfield_  => lelfield, &
                            efield_    => efield, &
                            nberrycyc_ => nberrycyc, &
                            efield_cart_ => efield_cart
  !
  USE cell_base,     ONLY : at, alat, omega, &
                            cell_base_init, init_dofree
  !
  USE ions_base,     ONLY : if_pos, ityp, tau, extfor, &
                            ntyp_ => nsp, &
                            nat_  => nat, &
                            amass, tau_format
  !
  USE basis,         ONLY : startingconfig, starting_wfc, starting_pot
  !
  USE run_info,      ONLY : title_ => title
  !
  USE cellmd,        ONLY : cmass, omega_old, at_old, ntcheck, &
                            cell_factor_ => cell_factor , &
                            press_       => press, &
                            calc, lmovecell
  !
  USE dynamics_module, ONLY : control_temp, temperature, thermostat, &
                              dt_         => dt, &
                              delta_t_    => delta_t, &
                              nraise_     => nraise, &
                              refold_pos_ => refold_pos
  !
  USE extfield,      ONLY : tefield_  => tefield, &
                            dipfield_ => dipfield, &
                            edir_     => edir, &
                            emaxpos_  => emaxpos, &
                            eopreg_   => eopreg, &
                            eamp_     => eamp, &
                            forcefield
  !
  USE io_files,      ONLY : input_drho, output_drho, &
                            psfile, tmp_dir, wfc_dir, &
                            prefix_     => prefix, &
                            pseudo_dir_ => pseudo_dir
  !
  USE force_mod,     ONLY : lforce, lstres, force
  !
  USE gvecs,         ONLY : dual
  USE gvect,         ONLY : ecutrho_ => ecutrho
  !
  USE fft_base, ONLY : dfftp
  USE fft_base, ONLY : dffts
  !
  USE klist,         ONLY : lgauss, ngauss, two_fermi_energies, &
                            smearing_          => smearing, &
                            degauss_           => degauss, &
                            tot_charge_        => tot_charge, &
                            tot_magnetization_ => tot_magnetization
  !
  USE ktetra,        ONLY : ltetra
  USE start_k,       ONLY : init_start_k
  !
  USE ldaU,          ONLY : Hubbard_U_     => hubbard_u, &
                            Hubbard_alpha_ => hubbard_alpha, &
                            lda_plus_u_    => lda_plus_u, &
                            niter_with_fixed_ns, starting_ns, U_projection
  !
  USE martyna_tuckerman, ONLY: do_comp_mt
#ifdef __SOLVENT
  USE control_flags, ONLY : save_vltot
  USE constants,     ONLY : rydberg_si, bohr_radius_si
  USE solvent_base,  ONLY : do_solvent_ => do_solvent,        &
                            verbose_ => verbose,              &
                            solvent_thr_ => solvent_thr,      &
                            stype_ => stype,                  &
                            rhomax_ => rhomax,                &
                            rhomin_ => rhomin,                &
                            tbeta_ => tbeta,                  &
                            epszero_ => epszero,              &
                            eps_mode_ => eps_mode,            &
                            solvationrad_ => solvationrad,    &
                            atomicspread_ => atomicspread,    &
                            ifdtype_ => ifdtype,              &
                            nfdpoint_ => nfdpoint,            &
                            mixrhopol_ => mixrhopol,          &
                            tolrhopol_ => tolrhopol,          &
                            gamma_ => gamma,                  &
                            delta_ => delta,                  &
                            extpressure_ => extpressure
#endif
  !
  USE esm,           ONLY: do_comp_esm, &
                           esm_bc_ => esm_bc, &
                           esm_nfit_ => esm_nfit, &
                           esm_efield_ => esm_efield, &
                           esm_w_ => esm_w
  !
  USE a2F,           ONLY : la2F_ => la2F
  !
  USE exx,           ONLY : x_gamma_extrapolation_ => x_gamma_extrapolation, &
                            nqx1_ => nq1, &
                            nqx2_ => nq2, &
                            nqx3_ => nq3, &
                            exxdiv_treatment_ => exxdiv_treatment, &
                            yukawa_           => yukawa, &
                            ecutvcut_         => ecutvcut
  !
  !
  USE lsda_mod,      ONLY : nspin_                  => nspin, &
                            starting_magnetization_ => starting_magnetization, &
                            lsda
  !
  USE kernel_table,  ONLY : vdw_table_name_ => vdw_table_name
  !
  USE relax,         ONLY : epse, epsf, epsp, starting_scf_threshold
  !
  USE control_flags, ONLY : isolve, max_cg_iter, david, tr2, imix, gamma_only,&
                            nmix, iverbosity, niter, pot_order, wfc_order, &
                            remove_rigid_rot_ => remove_rigid_rot, &
                            diago_full_acc_   => diago_full_acc, &
                            tolp_             => tolp, &
                            upscale_          => upscale, &
                            mixing_beta_      => mixing_beta, &
                            nstep_            => nstep, &
                            iprint_           => iprint, &
                            noinv_            => noinv, &
                            lkpoint_dir_      => lkpoint_dir, &
                            tqr_              => tqr, &
                            io_level, ethr, lscf, lbfgs, lmd, &
                            ldamped, lbands, llang,           &
                            lconstrain, restart, twfcollect, &
                            llondon, do_makov_payne, &
                            lecrpa_           => lecrpa
  !
  USE wvfct,         ONLY : nbnd_ => nbnd, &
                            ecutwfc_ => ecutwfc, &
                            ecfixed_ => ecfixed, &
                            qcutz_   => qcutz, &
                            q2sigma_ => q2sigma
  !
  USE fixed_occ,     ONLY : tfixed_occ, f_inp, &
                            one_atom_occupations_ => one_atom_occupations
  !
  USE noncollin_module, ONLY : i_cons, mcons, bfield, &
                               noncolin_  => noncolin, &
                               lambda_    => lambda, &
                               angle1_    => angle1, &
                               angle2_    => angle2, &
                               report_    => report
  !
  USE spin_orb, ONLY : lspinorb_ => lspinorb,  &
                       starting_spin_angle_ => starting_spin_angle

  !
  USE symm_base, ONLY : no_t_rev_ => no_t_rev, nofrac, allfrac, &
                        nosym_ => nosym, nosym_evc_=> nosym_evc
  !
  USE bfgs_module,   ONLY : bfgs_ndim_        => bfgs_ndim, &
                            trust_radius_max_ => trust_radius_max, &
                            trust_radius_min_ => trust_radius_min, &
                            trust_radius_ini_ => trust_radius_ini, &
                            w_1_              => w_1, &
                            w_2_              => w_2
  USE wannier_new, ONLY :   use_wannier_      => use_wannier, &
                            use_energy_int_   => use_energy_int, &
                            nwan_             => nwan, &
                            print_wannier_coeff_    => print_wannier_coeff

  USE realus,                ONLY : real_space_ => real_space

  USE read_pseudo_mod,       ONLY : readpp

#if defined __MS2
  USE MS2,                   ONLY : MS2_enabled_ => MS2_enabled, &
                                    MS2_handler_ => MS2_handler
#endif
  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : title, calculation, verbosity, restart_mode,    &
                               nstep, iprint, tstress, tprnfor, dt, outdir,    &
                               wfcdir, prefix, etot_conv_thr, forc_conv_thr,   &
                               pseudo_dir, disk_io, tefield, dipfield, lberry, &
                               gdir, nppstr, wf_collect,lelfield, efield,      &
                               nberrycyc, lkpoint_dir, efield_cart, lecrpa,    &
                               vdw_table_name

#if defined __MS2
  USE input_parameters, ONLY : MS2_enabled, MS2_handler
#endif
  !
  ! ... SYSTEM namelist
  !
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               nat, ntyp, nbnd,tot_charge,tot_magnetization,&
                               ecutwfc, ecutrho, nr1, nr2, nr3, nr1s, nr2s, &
                               nr3s, noinv, nosym, nosym_evc, no_t_rev,     &
                               use_all_frac, force_symmorphic,              &
                               starting_magnetization,                      &
                               occupations, degauss, smearing, nspin,       &
                               ecfixed, qcutz, q2sigma, lda_plus_U,         &
                               Hubbard_U, Hubbard_alpha, input_dft, la2F,   &
                               starting_ns_eigenvalue, U_projection_type,   &
#if defined (EXX)
                               x_gamma_extrapolation, nqx1, nqx2, nqx3,     &
                               exxdiv_treatment, yukawa, ecutvcut,          &
                               exx_fraction, screening_parameter,           &
#endif
#ifdef __SOLVENT
                               do_solvent,                                  &
#endif
                               edir, emaxpos, eopreg, eamp, noncolin, lambda, &
                               angle1, angle2, constrained_magnetization,     &
                               B_field, fixed_magnetization, report, lspinorb,&
                               starting_spin_angle,                           &
                               assume_isolated, spline_ps, london, london_s6, &
                               london_rcut, one_atom_occupations, &
                               esm_bc, esm_efield, esm_w, esm_nfit
#ifdef __SOLVENT
  !
  ! ... SOLVENT namelist
  !
  USE input_parameters, ONLY : verbose, solvent_thr,                          &
                               stype, rhomax, rhomin, tbeta,                  &
                               epszero, eps_mode, solvationrad, atomicspread, &
                               ifdtype, nfdpoint,                             &
                               mixrhopol, tolrhopol,                          &
                               gamma, delta,                                  &
                               extpressure 
#endif
  !
  ! ... ELECTRONS namelist
  !
  USE input_parameters, ONLY : electron_maxstep, mixing_mode, mixing_beta, &
                               mixing_ndim, mixing_fixed_ns, conv_thr,     &
                               tqr, diago_thr_init, diago_cg_maxiter,      &
                               diago_david_ndim, diagonalization,          &
                               diago_full_acc, startingwfc, startingpot,   &
                               real_space
#if defined (EXX)
  USE input_parameters, ONLY : adaptive_thr, conv_thr_init, conv_thr_multi
#endif
  !
  ! ... IONS namelist
  !
  USE input_parameters, ONLY : phase_space, ion_dynamics, ion_positions, tolp, &
                               tempw, delta_t, nraise, ion_temperature,        &
                               refold_pos, remove_rigid_rot, upscale,          &
                               pot_extrapolation,  wfc_extrapolation,          &
                               w_1, w_2, trust_radius_max, trust_radius_min,   &
                               trust_radius_ini, bfgs_ndim
  !
  ! ... CELL namelist
  !
  USE input_parameters, ONLY : cell_parameters, cell_dynamics, press, wmass,  &
                               cell_temperature, cell_factor, press_conv_thr, &
                               cell_dofree
  !
  ! ... WANNIER_NEW namelist
  !
  USE input_parameters, ONLY : use_wannier, nwan, constrain_pot, &
                               use_energy_int, print_wannier_coeff
  !
  ! ... CARDS
  !
  USE input_parameters,   ONLY : k_points, xk, wk, nk1, nk2, nk3,  &
                                 k1, k2, k3, nkstot
  USE input_parameters, ONLY : nconstr_inp, ncolvar_inp, trd_ht, rd_ht, &
                               cell_units
  !
  USE constraints_module,    ONLY : init_constraint
  USE read_namelists_module, ONLY : read_namelists, sm_not_set
  USE london_module,         ONLY : init_london, lon_rcut, scal6
  USE us, ONLY : spline_ps_ => spline_ps
  !
  USE input_parameters,       ONLY : deallocate_input_parameters
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER  :: ia, image, nt, inlc
  REAL(DP) :: theta, phi
  !
  !
  ! ... various initializations of control variables
  !
  lforce    = tprnfor
  !
  SELECT CASE( trim( calculation ) )
  CASE( 'scf' )
     !
     lscf  = .true.
     nstep = 1
     !
  CASE( 'nscf' )
     !
     lforce = .false.
     nstep  = 1
     !
  CASE( 'bands' )
     !
     lforce = .false.
     lbands = .true.
     nstep  = 1
     !
  CASE( 'relax' )
     !
     lscf   = .true.
     lforce = .true.
     !
     epse = etot_conv_thr
     epsf = forc_conv_thr
     !
     SELECT CASE( trim( ion_dynamics ) )
     CASE( 'bfgs' )
        !
        lbfgs = .true.
        !
     CASE ( 'damp' )
        !
        lmd     = .true.
        ldamped = .true.
        !
        ntcheck = nstep + 1
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // trim( calculation ) // &
                   & ': ion_dynamics=' // trim( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
  CASE( 'md' )
     !
     lscf   = .true.
     lmd    = .true.
     lforce = .true.
     !
     SELECT CASE( trim( ion_dynamics ) )
     CASE( 'verlet' )
        !
        CONTINUE
        !
     CASE( 'langevin' )
        !
        llang       = .true.
        temperature = tempw
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys ', 'calculation=' // trim( calculation ) // &
                   & ': ion_dynamics=' // trim( ion_dynamics ) // &
                   & ' not supported', 1 )
     END SELECT
     !
  CASE( 'vc-relax' )
     !
     lscf      = .true.
     lmd       = .true.
     lmovecell = .true.
     lforce    = .true.
     ldamped   = .true.
     !
     epse =  etot_conv_thr
     epsf =  forc_conv_thr
     epsp = press_conv_thr
     !
     SELECT CASE( trim( cell_dynamics ) )
     CASE( 'none' )
        !
        calc    = 'mm'
        ntcheck = nstep + 1
        !
     CASE( 'damp-pr' )
        !
        calc    = 'cm'
        ntcheck = nstep + 1
        !
     CASE( 'damp-w' )
        !
        calc    = 'nm'
        ntcheck = nstep + 1
        !
     CASE( 'bfgs' )
        !
        lbfgs = .true.
        lmd   = .false.
        ldamped = .false.
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // trim( calculation ) // &
                   & ': cell_dynamics=' // trim( cell_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( .not. ldamped .and. .not. lbfgs) &
        CALL errore( 'iosys', 'calculation='// trim( calculation ) // &
                   & ': incompatible ion (' // trim( ion_dynamics )// &
                   & ') and cell dynamics ('// trim(cell_dynamics )// ')', 1 )
     !
  CASE( 'vc-md' )
     !
     lscf      = .true.
     lmd       = .true.
     lmovecell = .true.
     lforce    = .true.
     !
     ntcheck = nstep + 1
     !
     SELECT CASE( trim( cell_dynamics ) )
     CASE( 'none' )
        !
        calc = 'md'
        !
     CASE( 'pr' )
        !
        calc = 'cd'
        !
     CASE( 'w' )
        !
        calc = 'nd'
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // trim( calculation ) // &
                   & ': ion_dynamics=' // trim( ion_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( trim( ion_dynamics ) /= 'beeman' ) &
        CALL errore( 'iosys', 'calculation=' // trim( calculation ) // &
                   & ': ion_dynamics=' // trim( ion_dynamics ) // &
                   & ' not supported', 1 )
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'calculation ' // &
                & trim( calculation ) // ' not implemented', 1 )
     !
  END SELECT
  !
  lstres = lmovecell .OR. ( tstress .and. lscf )
  !
  IF ( tefield .and. ( .not. nosym ) ) THEN
     nosym = .true.
     WRITE( stdout, &
            '(5x,"Presently no symmetry can be used with electric field",/)' )
  ENDIF
  IF ( tefield .and. tstress ) THEN
     tstress = .false.
     WRITE( stdout, &
            '(5x,"Presently stress not available with electric field",/)' )
  ENDIF
  IF ( tefield .and. ( nspin > 2 ) ) THEN
     CALL errore( 'iosys', 'LSDA not available with electric field' , 1 )
  ENDIF
  !
  twfcollect = wf_collect
  !
  ! ... Set Values for electron and bands
  !
  tfixed_occ = .false.
  ltetra     = .false.
  lgauss     = .false.
  !
  SELECT CASE( trim( occupations ) )
  CASE( 'fixed' )
     !
     ngauss = 0
     IF ( degauss /= 0.D0 ) THEN
        CALL errore( ' iosys ', &
                   & ' fixed occupations, gauss. broadening ignored', -1 )
        degauss = 0.D0
     ENDIF
     !
  CASE( 'smearing' )
     !
     lgauss = ( degauss > 0.0_dp ) 
     IF ( .NOT. lgauss ) &
        CALL errore( ' iosys ', &
                   & ' smearing requires gaussian broadening', 1 )
     !
     SELECT CASE ( trim( smearing ) )
     CASE ( 'gaussian', 'gauss', 'Gaussian', 'Gauss' )
        ngauss = 0
        smearing_ = 'gaussian'
     CASE ( 'methfessel-paxton', 'm-p', 'mp', 'Methfessel-Paxton', 'M-P', 'MP' )
        ngauss = 1
        smearing_ = 'Methfessel-Paxton'
     CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt', 'M-V', 'MV')
        ngauss = -1
        smearing_ = 'Marzari-Vanderbilt'
     CASE ( 'fermi-dirac', 'f-d', 'fd', 'Fermi-Dirac', 'F-D', 'FD')
        ngauss = -99
        smearing_ = 'Fermi-Dirac'
     CASE DEFAULT
        CALL errore( ' iosys ', ' smearing '//trim(smearing)//' unknown', 1 )
     END SELECT
     !
  CASE( 'tetrahedra' )
     !
     ! replace "errore" with "infomsg" in the next line if you really want
     ! to perform a calculation with forces using tetrahedra 
     !
     IF( lforce ) CALL errore( 'iosys', &
        'force calculation with tetrahedra not recommanded: use smearing',1)
     !
     ! as above, for stress
     !
     IF( lstres ) CALL errore( 'iosys', &
        'stress calculation with tetrahedra not recommanded: use smearing',1)
     ngauss = 0
     ltetra = .true.
     !
  CASE( 'from_input' )
     !
     ngauss     = 0
     tfixed_occ = .true.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','occupations ' // trim( occupations ) // &
                & 'not implemented', 1 )
     !
  END SELECT
  !
  IF( nbnd < 1 ) &
     CALL errore( 'iosys', 'nbnd less than 1', nbnd )
  !
  SELECT CASE( nspin )
  CASE( 1 )
     !
     lsda = .false.
     IF ( noncolin ) nspin = 4
     !
  CASE( 2 )
     !
     lsda = .true.
     IF ( noncolin ) &
        CALL errore( 'iosys', &
                     'noncolin .and. nspin==2 are conflicting flags', 1 )
     !
  CASE( 4 )
     !
     lsda = .false.
     noncolin = .true.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'wrong input value for nspin', 1 )
     !
  END SELECT
  !
  IF ( noncolin .and. lda_plus_u ) CALL errore('iosys', &
       'LDA+U not implemented with noncollinear magnetization', 1)
  !
  two_fermi_energies = ( tot_magnetization /= -1._DP)
  IF ( two_fermi_energies .and. tot_magnetization < 0._DP) &
     CALL errore( 'iosys', 'tot_magnetization only takes positive values', 1 )
  IF ( two_fermi_energies .and. .not. lsda ) &
     CALL errore( 'iosys', 'tot_magnetization requires nspin=2', 1 )
  !
  IF ( occupations == 'fixed' .and. lsda  .and. lscf ) THEN
     !
     IF ( two_fermi_energies ) THEN
        !
        IF ( abs( nint(tot_magnetization ) - tot_magnetization ) > eps8 ) &
           CALL errore( 'iosys', &
                 & 'fixed occupations requires integer tot_magnetization', 1 )
        IF ( abs( nint(tot_charge ) - tot_charge ) > eps8 ) &
           CALL errore( 'iosys', &
                      & 'fixed occupations requires integer charge', 1 )
        !
     ELSE
        !
        CALL errore( 'iosys', &
                   & 'fixed occupations and lsda need tot_magnetization', 1 )
        !
     ENDIF
     !
  ENDIF
  !
  IF (noncolin) THEN
     DO nt = 1, ntyp
        !
        angle1(nt) = pi * angle1(nt) / 180.D0
        angle2(nt) = pi * angle2(nt) / 180.D0
        !
     ENDDO
  ELSE
     angle1=0.d0
     angle2=0.d0
  ENDIF
  !
  SELECT CASE( trim( constrained_magnetization ) )
  CASE( 'none' )
     !
     ! ... starting_magnetization(nt) = sm_not_set means "not set"
     ! ... if no constraints are imposed on the magnetization, 
     ! ... starting_magnetization must be set for at least one atomic type
     !
     IF ( lscf .AND. lsda .AND. ( .NOT. tfixed_occ ) .AND. &
          ( .not. two_fermi_energies )  .AND. &
          ALL (starting_magnetization(1:ntyp) == sm_not_set) ) &
        CALL errore('iosys','some starting_magnetization MUST be set', 1 )
     !
     ! ... bring starting_magnetization between -1 and 1
     !
     DO nt = 1, ntyp
        !
        IF ( starting_magnetization(nt) == sm_not_set ) THEN
           starting_magnetization(nt) = 0.0_dp
        ELSEIF ( starting_magnetization(nt) > 1.0_dp ) THEN
          starting_magnetization(nt) = 1.0_dp
        ELSEIF ( starting_magnetization(nt) <-1.0_dp ) THEN
          starting_magnetization(nt) =-1.0_dp
        ENDIF
        !
     ENDDO
     !
     i_cons = 0
     !
  CASE( 'atomic' )
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys','constrained atomic magnetizations ' // &
                   & 'require nspin=2 or 4 ', 1 )
     IF ( ALL (starting_magnetization(1:ntyp) == sm_not_set) ) &
        CALL errore( 'iosys','constrained atomic magnetizations ' // &
                   & 'require that some starting_magnetization is set', 1 )
     !
     i_cons = 1
     !
     IF (nspin == 4) THEN
        ! non-collinear case
        DO nt = 1, ntyp
           !
           theta = angle1(nt)
           phi   = angle2(nt)
           !
           mcons(1,nt) = starting_magnetization(nt) * sin( theta ) * cos( phi )
           mcons(2,nt) = starting_magnetization(nt) * sin( theta ) * sin( phi )
           mcons(3,nt) = starting_magnetization(nt) * cos( theta )
           !
        ENDDO
     ELSE
        ! collinear case
        DO nt = 1, ntyp
           !
           mcons(1,nt) = starting_magnetization(nt)
           !
        ENDDO
     ENDIF
     !
  CASE( 'atomic direction' )
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys','constrained atomic magnetization ' // &
                   & 'directions require nspin=2 or 4 ', 1 )
     !
     i_cons = 2
     !
     DO nt = 1, ntyp
        !
        ! ... angle between the magnetic moments and the z-axis is
        ! ... constrained
        !
        theta = angle1(nt)
        mcons(3,nt) = cos(theta)
        !
     ENDDO
     !
  CASE( 'total' )
     !
     IF ( nspin == 4 ) THEN
        !
        i_cons = 3
        !
        mcons(1,1) = fixed_magnetization(1)
        mcons(2,1) = fixed_magnetization(2)
        mcons(3,1) = fixed_magnetization(3)
        !
     ELSE
        !
        CALL errore( 'iosys','constrained total magnetization ' // &
                   & 'requires nspin= 4 ', 1 )
        !
     ENDIF
     !
  CASE( 'total direction' )
     i_cons = 6
     mcons(3,1) = fixed_magnetization(3)
     IF ( mcons(3,1) < 0.D0 .or. mcons(3,1) > 180.D0 ) &
        CALL errore( 'iosys','constrained magnetization angle: ' // &
                   & 'theta must be within [0,180] degrees', 1 )
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','constrained magnetization ' // &
                & trim( constrained_magnetization ) // 'not implemented', 1 )
     !
  END SELECT
  !
  IF ( B_field(1) /= 0.D0 .or. &
       B_field(2) /= 0.D0 .or. &
       B_field(3) /= 0.D0 ) THEN
     !
     IF ( nspin == 1 ) CALL errore( 'iosys', &
          & 'non-zero external B_field requires nspin=2 or 4', 1 )
     IF ( TRIM( constrained_magnetization ) /= 'none' ) &
          CALL errore( 'iosys', 'constrained_magnetization and ' // &
                     & 'non-zero external B_field are conflicting flags', 1 )
     IF ( nspin == 2 .AND. ( B_field(1) /= 0.D0 .OR. B_field(2) /= 0.D0 ) ) &
        CALL errore('iosys','only B_field(3) can be specified with nspin=2', 1)
     IF ( i_cons /= 0 ) CALL errore( 'iosys', &
          & 'non-zero external B_field and constrained magnetization?', i_cons)
     !
     ! i_cons=4 signals the presence of an external B field
     ! this should be done in a cleaner way
     !
     i_cons = 4
     bfield(:)=B_field(:)
     !
  ENDIF
  !
  IF ( ecutrho <= 0.D0 ) THEN
     !
     dual = 4.D0
     ecutrho = dual*ecutwfc
     !
  ELSE
     !
     dual = ecutrho / ecutwfc
     IF ( dual <= 1.D0 ) &
        CALL errore( 'iosys', 'invalid dual?', 1 )
     !
  ENDIF
  !
  SELECT CASE( trim( restart_mode ) )
  CASE( 'from_scratch' )
     !
     restart        = .false.
     IF ( lscf ) THEN
        startingconfig = 'input'
     ELSE
        startingconfig = 'file'
     ENDIF
     !
  CASE( 'restart' )
     !
     restart = .true.
     !
     IF ( trim( ion_positions ) == 'from_input' ) THEN
        !
        startingconfig = 'input'
        !
     ELSE
        !
        startingconfig = 'file'
        !
     ENDIF
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', &
                & 'unknown restart_mode ' // trim( restart_mode ), 1 )
     !
  END SELECT
  !
  SELECT CASE( trim( disk_io ) )
  CASE( 'high' )
     !
     io_level = 2
     !
  CASE ( 'low' )
     !
     io_level = 0
     restart  = .false.
     !
  CASE ( 'none' )
     !
     io_level = -1
     restart  = .false.
     IF ( twfcollect ) THEN
        CALL infomsg('iosys', 'minimal I/O required, wf_collect reset to FALSE')
        twfcollect= .false.
     ENDIF
     !
  CASE DEFAULT
     !
     io_level = 1
     !
     IF ( lscf ) restart  = .false.
     !
  END SELECT
  !
  Hubbard_U(:)    = Hubbard_U(:) / rytoev
  Hubbard_alpha(:)= Hubbard_alpha(:) / rytoev
  !
  ethr = diago_thr_init
  !
  IF ( startingpot /= 'atomic' .and. startingpot /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingpot: use default (1)' )
     !
     IF ( lscf ) THEN
        startingpot = 'atomic'
     ELSE 
        startingpot = 'file'
     END IF
     !
  ENDIF
  !
  IF ( .not. lscf .and. startingpot /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingpot: use default (2)' )
     !
     startingpot = 'file'
     !
  ENDIF
  !
  IF (      startingwfc /= 'atomic' .and. &
            startingwfc /= 'random' .and. &
            startingwfc /= 'atomic+random' .and. &
            startingwfc /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingwfc: use default' )
     !
     startingwfc = 'atomic'
     !
  ENDIF
  !
  SELECT CASE( trim( diagonalization ) )
  CASE ( 'cg' )
     !
     isolve = 1
     max_cg_iter = diago_cg_maxiter
     !
  CASE ( 'david', 'davidson' )
     !
     isolve = 0
     david = diago_david_ndim
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'diagonalization ' // &
                & trim( diagonalization ) // ' not implemented', 1 )
     !
  END SELECT
  !
  tr2   = conv_thr
  niter = electron_maxstep
#if defined (EXX)
  adapt_thr = adaptive_thr
  tr2_init  = conv_thr_init
  tr2_multi = conv_thr_multi
#endif
  !
  pot_order = 1
  SELECT CASE( trim( pot_extrapolation ) )
  CASE( 'from_wfcs', 'from-wfcs' )
     ! not actually implemented
     pot_order =-1
     !
  CASE( 'none' )
     !
     pot_order = 0
     !
  CASE( 'first_order', 'first-order', 'first order' )
     !
     IF ( lmd  ) THEN
        pot_order = 2
     ELSE
        CALL infomsg('iosys', "pot_extrapolation='"//trim(pot_extrapolation)//&
                     "' not available, using 'atomic'")
     ENDIF
     !
  CASE( 'second_order', 'second-order', 'second order' )
     !
     IF ( lmd  ) THEN
        pot_order = 3
     ELSE
        CALL infomsg('iosys', "pot_extrapolation='"//trim(pot_extrapolation)//&
                     "' not available, using 'atomic'")
     ENDIF
     !
  CASE DEFAULT
     !
     pot_order = 1
     !
  END SELECT
  !
  wfc_order = 0
  SELECT CASE( trim( wfc_extrapolation ) )
     !
  CASE( 'first_order', 'first-order', 'first order' )
     !
     IF ( lmd  ) THEN
        wfc_order = 2
     ELSE
        CALL infomsg('iosys', "wfc_extrapolation='"//trim(pot_extrapolation)//&
                     "' not available, using 'atomic'")
     ENDIF
     !
  CASE( 'second_order', 'second-order', 'second order' )
     !
     IF ( lmd  ) THEN
        wfc_order = 3
     ELSE
        CALL infomsg('iosys', "wfc_extrapolation='"//trim(pot_extrapolation)//&
                     "' not available, using 'atomic'")
     ENDIF
     !
  END SELECT
  !
  SELECT CASE( trim( ion_temperature ) )
  CASE( 'not_controlled', 'not-controlled', 'not controlled' )
     !
     control_temp = .false.
     !
  CASE( 'initial' )
     !
     control_temp = .TRUE.
     thermostat   = TRIM( ion_temperature )
     temperature  = tempw
     !
  CASE( 'rescaling' )
     !
     control_temp = .true.
     thermostat   = trim( ion_temperature )
     temperature  = tempw
     tolp_        = tolp
     ntcheck      = nraise
     !
  CASE( 'rescale-v', 'rescale-V', 'rescale_v', 'rescale_V' )
     !
     control_temp = .true.
     thermostat   = trim( ion_temperature )
     temperature  = tempw
     nraise_      = nraise
     !
  CASE( 'reduce-T', 'reduce-t', 'reduce_T', 'reduce_t' )
     !
     control_temp = .true.
     thermostat   = trim( ion_temperature )
     temperature  = tempw
     delta_t_     = delta_t
     nraise_      = nraise
     !
  CASE( 'rescale-T', 'rescale-t', 'rescale_T', 'rescale_t' )
     !
     control_temp = .true.
     thermostat   = trim( ion_temperature )
     temperature  = tempw
     delta_t_     = delta_t
     !
  CASE( 'berendsen', ' Berendsen' )
     !
     control_temp = .true.
     thermostat   = trim( ion_temperature )
     temperature  = tempw
     nraise_      = nraise
     !
  CASE( 'andersen', 'Andersen' )
     !
     control_temp = .true.
     thermostat   = trim( ion_temperature )
     temperature  = tempw
     nraise_      = nraise
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', &
                & 'unknown ion_temperature ' // trim( ion_temperature ), 1 )
     !
  END SELECT
  !
  SELECT CASE( trim( mixing_mode ) )
  CASE( 'plain' )
     !
     imix = 0
     !
  CASE( 'TF' )
     !
     imix = 1
     !
  CASE( 'local-TF' )
     !
     imix = 2
     !
  CASE( 'potential' )
     !
     CALL errore( 'iosys', 'potential mixing no longer implemented', 1 )
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'unknown mixing ' // trim( mixing_mode ), 1 )
     !
  END SELECT
  !
  starting_scf_threshold = tr2
  nmix = mixing_ndim
  niter_with_fixed_ns = mixing_fixed_ns
  !
  IF ( ion_dynamics == ' bfgs' .and. epse <= 20.D0 * ( tr2 / upscale ) ) &
       CALL errore( 'iosys', 'required etot_conv_thr is too small:' // &
                     & ' conv_thr must be reduced', 1 )
  !
  SELECT CASE( trim( verbosity ) )
  CASE( 'debug', 'high', 'medium' )
     !
     iverbosity = 1
     !
  CASE( 'low', 'default', 'minimal' )
     !
     iverbosity = 0 
     !
  CASE DEFAULT
     !
     iverbosity = 0
     !
  END SELECT
  !
  tmp_dir = trimcheck ( outdir )
  !
  IF ( lberry ) THEN
     IF ( npool > 1 ) &
        CALL errore( 'iosys', 'Berry Phase not implemented with pools', 1 )
     IF ( noncolin )  &
        CALL errore( 'iosys', 'Noncolinear Berry Phase not implemented', 1 )
     IF ( lgauss .OR. ltetra ) &
        CALL errore( 'iosys', 'Berry Phase only for insulators!', 1 )
  END IF
  !
  IF ( lelfield ) THEN
     IF ( lgauss .OR. ltetra ) &
        CALL errore( 'iosys', 'electric fields only for insulators!', 1 )
  END IF
  !
  ! ... Copy values from input module to PW internals
  !
  nppstr_     = nppstr
  gdir_       = gdir
  lberry_     = lberry
  lelfield_   = lelfield
  efield_     = efield
  nberrycyc_  = nberrycyc
  efield_cart_ = efield_cart
  tqr_        = tqr
  real_space_ = real_space
  !
  title_      = title
  lkpoint_dir_=lkpoint_dir
  dt_         = dt
  tefield_    = tefield
  dipfield_   = dipfield
  prefix_     = trim( prefix )
  pseudo_dir_ = trimcheck( pseudo_dir )
  nstep_      = nstep
  iprint_     = iprint
  lecrpa_     = lecrpa
  !
  nat_     = nat
  ntyp_    = ntyp
  edir_    = edir
  emaxpos_ = emaxpos
  eopreg_  = eopreg
  eamp_    = eamp
  dfftp%nr1     = nr1
  dfftp%nr2     = nr2
  dfftp%nr3     = nr3
  ecutrho_ = ecutrho
  ecutwfc_ = ecutwfc
  ecfixed_ = ecfixed
  qcutz_   = qcutz
  q2sigma_ = q2sigma
  dffts%nr1    = nr1s
  dffts%nr2    = nr2s
  dffts%nr3    = nr3s
  degauss_ = degauss
  !
  tot_charge_        = tot_charge
  tot_magnetization_ = tot_magnetization
  !
  lspinorb_ = lspinorb
  starting_spin_angle_ = starting_spin_angle
  noncolin_ = noncolin
  angle1_   = angle1
  angle2_   = angle2
  report_   = report
  lambda_   = lambda
  one_atom_occupations_ = one_atom_occupations
  !
  no_t_rev_ = no_t_rev
  allfrac   = use_all_frac
  !
  spline_ps_ = spline_ps
  !
  Hubbard_U_(1:ntyp)      = hubbard_u(1:ntyp)
  Hubbard_alpha_(1:ntyp)  = hubbard_alpha(1:ntyp)
  lda_plus_u_             = lda_plus_u
  la2F_                   = la2F
  nspin_                  = nspin
  starting_magnetization_ = starting_magnetization
  starting_ns             = starting_ns_eigenvalue
  U_projection            = U_projection_type
  noinv_                  = noinv
  nosym_                  = nosym
  nosym_evc_              = nosym_evc
  nofrac                  = force_symmorphic
  nbnd_                   = nbnd
  !
#if defined (EXX)
  !
  x_gamma_extrapolation_ = x_gamma_extrapolation
  !
  nqx1_ = nqx1
  nqx2_ = nqx2
  nqx3_ = nqx3
  !
  exxdiv_treatment_ = trim(exxdiv_treatment)
  yukawa_   = yukawa
  ecutvcut_ = ecutvcut
  !
#endif
  !
  vdw_table_name_  = vdw_table_name
  !
  diago_full_acc_ = diago_full_acc
  starting_wfc    = startingwfc
  starting_pot    = startingpot
  mixing_beta_    = mixing_beta
  !
  remove_rigid_rot_ = remove_rigid_rot
  upscale_          = upscale
  refold_pos_       = refold_pos
  press_            = press
  cell_factor_      = cell_factor
  !
  ! ... for WANNIER_AC
  !
  use_wannier_ = use_wannier
  use_energy_int_ = use_energy_int
  nwan_ = nwan
  print_wannier_coeff_ = print_wannier_coeff
  !
  !
  ! ... BFGS specific
  !
  bfgs_ndim_        = bfgs_ndim
  trust_radius_max_ = trust_radius_max
  trust_radius_min_ = trust_radius_min
  trust_radius_ini_ = trust_radius_ini
  w_1_              = w_1
  w_2_              = w_2
  !
#ifdef __SOLVENT
  !
  ! ...  Solvent
  !
  do_solvent_ = do_solvent
  save_vltot = .false.
  IF ( do_solvent ) save_vltot = .true.
  verbose_  = verbose
  solvent_thr_  = solvent_thr
  !
  stype_    = stype
  rhomax_   = rhomax
  rhomin_   = rhomin
  tbeta_    = tbeta
  IF ( stype .EQ. 1 ) THEN
    tbeta_  = LOG( rhomax / rhomin )
  END IF
  !
  epszero_ = epszero
  eps_mode_ = eps_mode
  ALLOCATE( solvationrad_( ntyp ) )
  solvationrad_( 1:ntyp ) = solvationrad( 1:ntyp )
  ALLOCATE( atomicspread_( ntyp ) )
  atomicspread_( 1:ntyp ) = atomicspread( 1:ntyp )
  !
  ifdtype_ = ifdtype
  nfdpoint_ = nfdpoint
  !
  mixrhopol_ = mixrhopol
  tolrhopol_ = tolrhopol
  !
  gamma_      = gamma*1.D-3*bohr_radius_si**2/rydberg_si
  delta_      = delta
  !
  extpressure_ = extpressure*1.D9/rydberg_si*bohr_radius_si**3
  !
#endif
  !
  ! ... ESM
  !
  esm_bc_ = esm_bc
  esm_efield_ = esm_efield
  esm_w_ = esm_w
  esm_nfit_ = esm_nfit
  !
  IF (trim(occupations) /= 'from_input') one_atom_occupations_=.false.
  !
  llondon     = london
  lon_rcut    = london_rcut
  scal6       = london_s6
  !
#if defined __MS2
  !
  ! MS2 specific parameters
  !
  MS2_enabled_ = MS2_enabled
  MS2_handler_ = MS2_handler
#endif
  !
  SELECT CASE( trim( assume_isolated ) )
      !
    CASE( 'makov-payne', 'm-p', 'mp' )
      !
      do_makov_payne = .true.
      do_comp_mt     = .false.
      do_comp_esm    = .false.
      !
    CASE( 'dcc' )
      !
      CALL errore('iosys','density countercharge correction currently disabled',1)
      !
    CASE( 'martyna-tuckerman', 'm-t', 'mt' )
      !
      do_comp_mt     = .true.
      do_makov_payne = .false.
      do_comp_esm    = .false.
      !
    CASE( 'esm' )
      !
      do_comp_esm    = .true.
      do_comp_mt     = .false.
      do_makov_payne = .false.
      !
    CASE( 'none' )
      !
      do_makov_payne = .false.
      do_comp_mt     = .false.
      do_comp_esm    = .false.
      !
    CASE DEFAULT
      !
      call errore ('iosys','unrecognized value for assume_isolated',1)
  END SELECT
  !
  ! ... read following cards
  !
  ALLOCATE( ityp( nat_ ) )
  ALLOCATE( tau(    3, nat_ ) )
  ALLOCATE( force(  3, nat_ ) )
  ALLOCATE( if_pos( 3, nat_ ) )
  ALLOCATE( extfor( 3, nat_ ) )
  IF ( tfixed_occ ) THEN
     IF ( nspin_ == 4 ) THEN
        ALLOCATE( f_inp( nbnd_, 1 ) )
     ELSE
        ALLOCATE( f_inp( nbnd_, nspin_ ) )
     ENDIF
  ENDIF
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat_ ) )
  !
  ! ... note that read_cards_pw no longer reads cards!
  !
  CALL read_cards_pw ( psfile, tau_format )
  !
  ! ... set up atomic positions and crystal lattice
  !
  call cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                        trd_ht, rd_ht, cell_units )
  !
  ! ... set up k-points
  !
  CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, nkstot, xk, wk )
  gamma_only = ( k_points == 'gamma' )
  !
  CALL convert_tau ( tau_format, nat_, tau)
  !
  IF ( wmass == 0.D0 ) THEN
     !
     ! ... set default value of wmass
     !
#if defined __PGI
     DO ia = 1, nat_
        wmass = wmass + amass( ityp(ia) )
     ENDDO
#else
     wmass = sum( amass(ityp(:)) )
#endif
     !
     wmass = wmass * amconv
     IF ( calc == 'nd' .or. calc == 'nm' ) THEN
        wmass = 0.75D0 * wmass / pi / pi / omega**( 2.D0 / 3.D0 )
     ELSEIF ( calc == 'cd' .or. calc == 'cm' ) THEN
        wmass = 0.75D0 * wmass / pi / pi
     ENDIF
     !
     cmass  = wmass
     !
  ELSE
     !
     ! ... wmass is given in amu, Renata's dynamics uses masses in atomic units
     !
     cmass  = wmass * amconv
     !
  ENDIF
  !
  ! ... unit conversion for pressure
  !
  press_ = press_ / ry_kbar
  !
  ! ... set constraints for cell dynamics/optimization
  !
  CALL init_dofree ( cell_dofree )
  !
  ! ... read pseudopotentials (also sets DFT)
  !
  CALL readpp ( input_dft )
  !
#if defined(EXX)
    !
    ! Set variables for hybrid functional HSE
    !
    IF (exx_fraction >= 0.0_DP) CALL set_exx_fraction (exx_fraction)
    IF (screening_parameter >= 0.0_DP) &
        & CALL set_screening_parameter (screening_parameter)
    !
#endif
  !
  ! ... read the vdw kernel table if needed
  !
  inlc = get_inlc()
  if (inlc == 1 .or. inlc == 2) then
      call initialize_kernel_table()
  endif
  !
  ! ... if DFT finite size corrections are needed, define the appropriate volume
  !
  IF (dft_has_finite_size_correction()) &
      CALL set_finite_size_volume(REAL(omega*nk1*nk2*nk3))
  !
  ! ... In the case of variable cell dynamics save old cell variables
  ! ... and initialize a few other variables
  !
  IF ( lmovecell ) THEN
     !
     at_old    = at
     omega_old = omega
     IF ( cell_factor_ <= 0.D0 ) cell_factor_ = 1.2D0
     !
     IF ( cmass <= 0.D0 ) &
        CALL errore( 'iosys', &
                   & 'vcsmd: a positive value for cell mass is required', 1 )
     !
  ELSE
     !
     cell_factor_ = 1.D0
     !
  ENDIF
  !
  ! ... allocate arrays for dispersion correction
  !
  IF ( llondon) CALL init_london ( )
  !
  ! ... variables for constrained dynamics are set here
  !
  lconstrain = ( ncolvar_inp + nconstr_inp > 0 )
  !
  IF ( lconstrain ) THEN
     IF ( lbfgs .OR. lmovecell ) CALL errore( 'iosys', &
              'constraints only with fixed-cell dynamics', 1 )
     CALL init_constraint( nat, tau, ityp, alat )
  END IF
  !
  ! ... read atomic positions and unit cell from data file
  ! ... must be done before "verify_tmpdir" because the latter
  ! ... removes the data file in a run from scratch
  !
  IF ( startingconfig == 'file' ) CALL read_config_from_file()
  !
  CALL verify_tmpdir( tmp_dir )
  !
  IF ( .not. trim( wfcdir ) == 'undefined' ) THEN
     !
     wfc_dir = trimcheck ( wfcdir )
     !
     CALL verify_tmpdir( wfc_dir )
     !
  ENDIF
  !
  CALL restart_from_file()
  !
  ! ... Files
  !
  input_drho  = ' '
  output_drho = ' '
  !
  IF (real_space ) THEN
   WRITE( stdout, '(5x,"Real space treatment of Beta functions, &
         &V.1 (BE SURE TO CHECK MANUAL!)",/)' )
  ENDIF
  !
  ! Deallocation of temp input arrays
  !
  CALL deallocate_input_parameters ()  
  !
  RETURN
  !
END SUBROUTINE iosys
!
!----------------------------------------------------------------------------
SUBROUTINE read_cards_pw ( psfile, tau_format )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE input_parameters,   ONLY : atom_label, atom_pfile, atom_mass, taspc, &
                                 tapos, rd_pos, atomic_positions, if_pos,  &
                                 sp_pos, f_inp, rd_for, tavel, sp_vel, rd_vel
  USE dynamics_module,    ONLY : tavel_ => tavel, vel
  USE cell_base,          ONLY : at, ibrav
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp, tau, atm, extfor
  USE fixed_occ,          ONLY : tfixed_occ, f_inp_ => f_inp
  USE ions_base,          ONLY : if_pos_ =>  if_pos, amass, fixatom
  USE control_flags,      ONLY : lfixatom, textfor
  !
  IMPLICIT NONE
  !
  CHARACTER (len=256) :: psfile(ntyp)
  CHARACTER (len=80)  :: tau_format
  INTEGER, EXTERNAL :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  !
  INTEGER :: is, ia
  !
  !
  amass = 0
  !
  IF ( .not. taspc ) &
     CALL errore( 'read_cards_pw', 'atomic species info missing', 1 )
  IF ( .not. tapos ) &
     CALL errore( 'read_cards_pw', 'atomic position info missing', 1 )
  !
  DO is = 1, ntyp
     !
     amass(is)  = atom_mass(is)
     psfile(is) = atom_pfile(is)
     atm(is)    = atom_label(is)
     !
     IF ( amass(is) <= 0.0_DP ) amass(is)= &
              atom_weight(atomic_number(trim(atm(is))))

     IF ( amass(is) <= 0.D0 ) CALL errore( 'read_cards_pw', 'invalid  mass', is )
     !
  ENDDO
  !
  textfor = .false.
  IF( any( rd_for /= 0.0_DP ) ) textfor = .true.
  !
  DO ia = 1, nat
     !
     tau(:,ia) = rd_pos(:,ia)
     ityp(ia)  = sp_pos(ia)
     extfor(:,ia) = rd_for(:,ia)
     !
  ENDDO
  !
  ! ... check for initial velocities read from input file
  !
  IF ( tavel .AND. ANY ( sp_pos(:) /= sp_vel(:) ) ) &
      CALL errore("cards","list of species in block ATOMIC_VELOCITIES &
                 & must be identical to those in ATOMIC_POSITIONS",1)
  tavel_ = tavel
  IF ( tavel_ ) THEN
     ALLOCATE( vel(3, nat) )
     DO ia = 1, nat
        vel(:,ia) = rd_vel(:,ia)
     END DO
  END IF
  !
  ! ... The constrain on fixed coordinates is implemented using the array
  ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1
  ! ... otherwise. 
  !
  if_pos_(:,:) = if_pos(:,1:nat)
  fixatom = COUNT( if_pos_(1,:)==0 .AND. if_pos_(2,:)==0 .AND. if_pos_(3,:)==0 )
  lfixatom = ANY ( if_pos_ == 0 )
  !
  tau_format = trim( atomic_positions )
  !
  IF ( tfixed_occ ) THEN
     !
     f_inp_ = f_inp
     !
     DEALLOCATE ( f_inp )
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE read_cards_pw
!
!-----------------------------------------------------------------------
SUBROUTINE convert_tau (tau_format, nat_, tau)
!-----------------------------------------------------------------------
  !
  ! ... convert input atomic positions to internally used format:
  ! ... tau in a0 units
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : bohr_radius_angs
  USE cell_base,     ONLY : at, alat
  IMPLICIT NONE
  CHARACTER (len=*), INTENT(in)  :: tau_format
  INTEGER, INTENT(in)  :: nat_
  REAL (DP), INTENT(inout) :: tau(3,nat_)
  !
  SELECT CASE( tau_format )
  CASE( 'alat' )
     !
     ! ... input atomic positions are divided by a0: do nothing
     !
  CASE( 'bohr' )
     !
     ! ... input atomic positions are in a.u.: divide by alat
     !
     tau = tau / alat
     !
  CASE( 'crystal' )
     !
     ! ... input atomic positions are in crystal axis
     !
     CALL cryst_to_cart( nat_, tau, at, 1 )
     !
  CASE( 'angstrom' )
     !
     ! ... atomic positions in A: convert to a.u. and divide by alat
     !
     tau = tau / bohr_radius_angs / alat
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','tau_format=' // &
                & trim( tau_format ) // ' not implemented', 1 )
     !
  END SELECT
  !
END SUBROUTINE convert_tau
!-----------------------------------------------------------------------
SUBROUTINE verify_tmpdir( tmp_dir )
  !-----------------------------------------------------------------------
  !
  USE wrappers,         ONLY : f_mkdir
  USE input_parameters, ONLY : restart_mode
  USE control_flags,    ONLY : lbands
  USE io_files,         ONLY : prefix, xmlpun, &
                               delete_if_present, check_writable
  USE pw_restart,       ONLY : pw_readfile
  USE mp_global,        ONLY : mpime, nproc
  USE io_global,        ONLY : ionode
  USE mp,               ONLY : mp_barrier
  USE xml_io_base,      ONLY : copy_file
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(inout) :: tmp_dir
  !
  INTEGER             :: ios, image, proc, nofi
  LOGICAL             :: exst
  CHARACTER (len=256) :: file_path, filename
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  !
  file_path = trim( tmp_dir ) // trim( prefix )
  !
  !
  IF ( restart_mode == 'from_scratch' ) THEN
     !
     ! ... let us try to create the scratch directory
     !
     CALL parallel_mkdir ( tmp_dir )
     !
  ENDIF
  !
  !
  ! ... if starting from scratch all temporary files are removed
  ! ... from tmp_dir ( only by the master node )
  !
  IF ( restart_mode == 'from_scratch' ) THEN
     !
     ! ... xml data file in save directory is removed
     !     but, header is read anyway to store qexml version
     !
     CALL pw_readfile( 'header', ios )
     !
     IF ( ionode ) THEN
        !
        IF ( .not. lbands ) THEN
            !
            ! save a bck copy of datafile.xml (AF)
            !
            filename = trim( file_path ) // '.save/' // trim( xmlpun )
            INQUIRE( FILE = filename, EXIST = exst )
            !
            IF ( exst ) CALL copy_file( trim(filename), trim(filename) // '.bck' )
            !
            CALL delete_if_present( trim(filename) )
            !
        ENDIF
        !
        ! ... extrapolation file is removed
        !
        CALL delete_if_present( trim( file_path ) // '.update' )
        !
        ! ... MD restart file is removed
        !
        CALL delete_if_present( trim( file_path ) // '.md' )
        !
        ! ... BFGS restart file is removed
        !
        CALL delete_if_present( trim( file_path ) // '.bfgs' )
        !
     ENDIF
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE verify_tmpdir

!-----------------------------------------------------------------------
SUBROUTINE parallel_mkdir ( tmp_dir )
  !-----------------------------------------------------------------------
  !
  ! ... Safe creation of the scratch directory in the parallel case
  ! ... Works on both parallel and distributed file systems
  ! ... Not really a smart algorithm, though
  !
  USE wrappers,      ONLY : f_mkdir
  USE mp_global,     ONLY : mpime, nproc
  USE mp,            ONLY : mp_barrier, mp_sum
  USE io_files,      ONLY : check_writable
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: tmp_dir
  !
  INTEGER             :: ios, proc
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  ! ... the scratch directory is created sequentially by all the cpus
  !
  DO proc = 0, nproc - 1
     !
     IF ( proc == mpime ) ios = f_mkdir( trim( tmp_dir ) )
     CALL mp_barrier()
     !
  ENDDO
  !
  ! ... each job checks whether the scratch directory is writable
  ! ... note that tmp_dir should end by a "/"
  !
  IF ( ios /= 0 ) CALL errore( 'parallel_mkdir', trim( tmp_dir ) // &
                             & ' non existent or non writable', 1 )
  !
  RETURN
  !
END SUBROUTINE parallel_mkdir

