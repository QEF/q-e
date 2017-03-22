
! Copyright (C) 2002-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! included monopole related stuff, search for 'TB'
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE iosys()
  !-----------------------------------------------------------------------------
  !
  ! ...  Copy data read from input file (in subroutine "read_input_file") and
  ! ...  stored in modules input_parameters into internal modules
  ! ...  Note that many variables in internal modules, having the same name as
  ! ...  those in input_parameters, are locally renamed by adding a "_"
  !
  USE kinds,         ONLY : DP
  USE funct,         ONLY : dft_is_hybrid, dft_has_finite_size_correction, &
                            set_finite_size_volume, get_inlc 
  USE funct,         ONLY: set_exx_fraction, set_screening_parameter
  USE control_flags, ONLY: adapt_thr, tr2_init, tr2_multi
  USE constants,     ONLY : autoev, eV_to_kelvin, pi, rytoev, &
                            ry_kbar, amu_ry, bohr_radius_angs, eps8
  USE mp_pools,      ONLY : npool
  !
  USE io_global,     ONLY : stdout, ionode, ionode_id
  !
  USE kernel_table,  ONLY : initialize_kernel_table
  !
  USE bp,            ONLY : nppstr_    => nppstr, &
                            gdir_      => gdir, &
                            lberry_    => lberry, &
                            lelfield_  => lelfield, &
                            lorbm_     => lorbm, &
                            efield_    => efield, &
                            nberrycyc_ => nberrycyc, &
                            efield_cart_ => efield_cart, &
                            phase_control
  !
  USE cell_base,     ONLY : at, alat, omega, bg, &
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
  USE fcp_variables, ONLY : lfcpopt_ => lfcpopt, &
                            lfcpdyn_ => lfcpdyn, &
                            fcp_mu_ => fcp_mu, &
                            fcp_mass_ => fcp_mass, &
                            fcp_temperature, &
                            fcp_relax_step_ => fcp_relax_step, &
                            fcp_relax_crit_ => fcp_relax_crit
  !
  USE extfield,      ONLY : tefield_  => tefield, &
                            dipfield_ => dipfield, &
                            edir_     => edir, &
                            emaxpos_  => emaxpos, &
                            eopreg_   => eopreg, &
                            eamp_     => eamp, &
  ! TB added monopole related variables
                            zmon_     => zmon, &
                            monopole_ => monopole, &
                            relaxz_   => relaxz, &
                            block_   => block, &
                            block_1_   => block_1, &
                            block_2_   => block_2, &
                            block_height_   => block_height, &
                            forcefield, &
                            forcemono
  !
  USE io_files,      ONLY : input_drho, output_drho, &
                            psfile, tmp_dir, wfc_dir, &
                            prefix_     => prefix, &
                            pseudo_dir_ => pseudo_dir
  !
  USE force_mod,     ONLY : lforce, lstres, force
  !
  USE fft_base, ONLY : dfftp
  USE fft_base, ONLY : dffts
  !
  USE klist,         ONLY : ltetra, lgauss, ngauss, two_fermi_energies, &
                            smearing_          => smearing, &
                            degauss_           => degauss, &
                            tot_charge_        => tot_charge, &
                            tot_magnetization_ => tot_magnetization
  USE ktetra,        ONLY : tetra_type
  USE start_k,       ONLY : init_start_k
  !
  USE ldaU,          ONLY : Hubbard_U_     => hubbard_u, &
                            Hubbard_J0_ => hubbard_j0, &
                            Hubbard_J_ => hubbard_j, &
                            Hubbard_alpha_ => hubbard_alpha, &
                            Hubbard_beta_ => hubbard_beta, &
                            lda_plus_u_    => lda_plus_u, &
                            lda_plus_u_kind_    => lda_plus_u_kind, &
                            niter_with_fixed_ns, starting_ns, U_projection
  !
  USE martyna_tuckerman, ONLY: do_comp_mt
  !
  USE esm,           ONLY: do_comp_esm, &
                           esm_bc_ => esm_bc, &
                           esm_nfit_ => esm_nfit, &
                           esm_efield_ => esm_efield, &
                           esm_w_ => esm_w, &
                           esm_a_ => esm_a
  !
  USE a2F,           ONLY : la2F_ => la2F
  !
  USE exx,           ONLY : x_gamma_extrapolation_ => x_gamma_extrapolation, &
                            nqx1_ => nq1, &
                            nqx2_ => nq2, &
                            nqx3_ => nq3, &
                            exxdiv_treatment_ => exxdiv_treatment, &
                            yukawa_           => yukawa, &
                            ecutvcut_         => ecutvcut, &
                            ecutfock_         => ecutfock
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
  USE extrapolation, ONLY : pot_order, wfc_order
  USE control_flags, ONLY : isolve, max_cg_iter, david, tr2, imix, gamma_only,&
                            nmix, iverbosity, niter, &
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
                            tq_smoothing_     => tq_smoothing, &
                            tbeta_smoothing_  => tbeta_smoothing, &
                            io_level, ethr, lscf, lbfgs, lmd, &
                            lbands, lconstrain, restart, twfcollect, &
                            llondon, do_makov_payne, lxdm, &
                            ts_vdw_           => ts_vdw, &
                            lecrpa_           => lecrpa, &
                            smallmem
  USE control_flags, ONLY: scf_must_converge_ => scf_must_converge
  !
  USE wvfct,         ONLY : nbnd_ => nbnd
  USE gvecw,         ONLY : ecfixed_ => ecfixed, &
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
                       lforcet_ => lforcet,    &
                       starting_spin_angle_ => starting_spin_angle

  !
  USE symm_base, ONLY : no_t_rev_ => no_t_rev, nofrac, allfrac, &
                        nosym_ => nosym, nosym_evc_=> nosym_evc, spacegroup
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

  USE qmmm, ONLY : qmmm_config

  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : title, calculation, verbosity, restart_mode,    &
                               nstep, iprint, tstress, tprnfor, dt, outdir,    &
                               wfcdir, prefix, etot_conv_thr, forc_conv_thr,   &
                               pseudo_dir, disk_io, tefield, dipfield, lberry, &
                               gdir, nppstr, wf_collect,lelfield,lorbm,efield, &
                               nberrycyc, lkpoint_dir, efield_cart, lecrpa,    &
                               vdw_table_name, memory, tqmmm,                  &
                               efield_phase, monopole

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
                               lda_plus_U_kind, Hubbard_U, Hubbard_J,       &
                               Hubbard_J0, Hubbard_beta,                    &
                               Hubbard_alpha, input_dft, la2F,              &
                               starting_ns_eigenvalue, U_projection_type,   &
                               x_gamma_extrapolation, nqx1, nqx2, nqx3,     &
                               exxdiv_treatment, yukawa, ecutvcut,          &
                               exx_fraction, screening_parameter, ecutfock, &
                               gau_parameter,                               &
                               edir, emaxpos, eopreg, eamp, noncolin, lambda, &
                               angle1, angle2, constrained_magnetization,     &
                               B_field, fixed_magnetization, report, lspinorb,&
                               starting_spin_angle, assume_isolated,spline_ps,&
                               vdw_corr, london, london_s6, london_rcut, london_c6, &
                               london_rvdw, &
                               ts_vdw, ts_vdw_isolated, ts_vdw_econv_thr,     &
                               xdm, xdm_a1, xdm_a2, lforcet,                  &
                               one_atom_occupations,                          &
                               esm_bc, esm_efield, esm_w, esm_nfit, esm_a,    &
                               lfcpopt, lfcpdyn, fcp_mu, fcp_mass, fcp_tempw, & 
                               fcp_relax_step, fcp_relax_crit,                &
                               space_group, uniqueb, origin_choice,           &
                               rhombohedral, zmon, relaxz, block, block_1,    &
                               block_2, block_height
  !
  ! ... ELECTRONS namelist
  !
  USE input_parameters, ONLY : electron_maxstep, mixing_mode, mixing_beta, &
                               mixing_ndim, mixing_fixed_ns, conv_thr,     &
                               tqr, tq_smoothing, tbeta_smoothing,         &
                               diago_thr_init, diago_cg_maxiter,           &
                               diago_david_ndim, diagonalization,          &
                               diago_full_acc, startingwfc, startingpot,   &
                               real_space, scf_must_converge
  USE input_parameters, ONLY : adaptive_thr, conv_thr_init, conv_thr_multi
  !
  ! ... IONS namelist
  !
  USE input_parameters, ONLY : ion_dynamics, ion_positions, tolp, &
                               tempw, delta_t, nraise, ion_temperature,        &
                               refold_pos, remove_rigid_rot, upscale,          &
                               pot_extrapolation,  wfc_extrapolation,          &
                               w_1, w_2, trust_radius_max, trust_radius_min,   &
                               trust_radius_ini, bfgs_ndim, rd_pos, sp_pos, &
                               rd_for, rd_if_pos => if_pos, lsg
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
  USE input_parameters, ONLY : nconstr_inp, trd_ht, rd_ht, cell_units
  !
  USE constraints_module,    ONLY : init_constraint
  USE read_namelists_module, ONLY : read_namelists, sm_not_set
  USE london_module,         ONLY : init_london, lon_rcut, scal6, in_c6, in_rvdw
  USE xdm_module,            ONLY : init_xdm, a1i, a2i
  USE tsvdw_module,          ONLY : vdw_isolated, vdw_econv_thr
  USE us,                    ONLY : spline_ps_ => spline_ps
  !
  USE input_parameters,      ONLY : deallocate_input_parameters
  USE wyckoff,               ONLY : nattot, sup_spacegroup
  USE qexsd_module,          ONLY : qexsd_input_obj
  USE qes_types_module,      ONLY: input_type
  ! 
  IMPLICIT NONE
  !
  INTERFACE  
     SUBROUTINE   pw_init_qexsd_input(obj,obj_tagname)
     IMPORT                       :: input_type
     TYPE(input_type)             :: obj
     CHARACTER(LEN=*),INTENT(IN)  :: obj_tagname
     END SUBROUTINE
  END INTERFACE
!!!!  
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER, EXTERNAL :: read_config_from_file
  !
  INTEGER  :: ia, nt, inlc, ibrav_sg, ierr
  LOGICAL  :: exst, parallelfs
  REAL(DP) :: theta, phi, ecutwfc_pp, ecutrho_pp
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
        calc    = 'vm'
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
        calc        = 'vd'
        !
     CASE( 'langevin' )
        !
        calc        = 'ld'
        temperature = tempw
        !
     CASE( 'langevin-smc', 'langevin+smc' )
        !
        calc        = 'ls'
        temperature = tempw
        !
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
        !
     CASE DEFAULT
        !
        CALL errore( 'iosys', 'calculation=' // trim( calculation ) // &
                   & ': cell_dynamics=' // trim( cell_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( lbfgs .AND. TRIM(ion_dynamics) /= 'bfgs' ) &
        CALL infomsg( 'iosys', 'calculation='// trim( calculation ) // &
                   & ': ion dynamics ' // trim( ion_dynamics )// &
                   & " ignored, 'bfgs' assumed" )
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
                   & ': cell_dynamics=' // trim( cell_dynamics ) // &
                   & ' not supported', 1 )
        !
     END SELECT
     !
     IF ( trim( ion_dynamics ) /= 'beeman' ) &
        CALL infomsg( 'iosys', 'calculation=' // trim( calculation ) // &
                   & ': ion_dynamics=' // trim( ion_dynamics ) // &
                   & " ignored, assuming 'beeman'" )
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
  ! TB
  ! IF ( tefield .and. ( .not. nosym ) ) THEN
  IF ( tefield .and. ( .not. nosym ) .and. ( .not. monopole)) THEN
     nosym = .true.
     WRITE( stdout, &
            '(5x,"Presently no symmetry can be used with electric field",/)' )
  ENDIF
  !TB begin some checks on input
  IF ( (monopole) .AND. ( .NOT. nosym )) THEN
     WRITE( stdout,'(/,5x,"Presently symmetry can be used with monopole field",/)' )
     WRITE( stdout,'(5x,"setting verbosity to high",/)' )
     WRITE( stdout,'(5x,"CAREFULLY CHECK ALL SYMMETRIES",/)' )
     verbosity='high'
  ENDIF
  IF ((zmon>1.0).OR.(zmon<0.0)) &
     CALL errore( 'iosys', 'Position of the monopole has to be between within ]0,1[' , 1 )
  IF ( (monopole) .AND. ((tefield).OR.(dipfield)) ) THEN
     IF (edir .ne. 3) &
        CALL errore( 'iosys','Using monopole and tefield/dipfield, edir must be 3', 1)
     IF ((zmon>=emaxpos).AND.(zmon<=(emaxpos+eopreg))) &
        CALL errore( 'iosys', 'Monopole between the 2 dipole planes not allowed' , 1 )
     IF ((block).AND.(block_1.NE.emaxpos).AND.(block_2.NE.(emaxpos+eopreg))) THEN
        WRITE( stdout,'(/,5x,"Neither block_1=emaxpos, nor block_2=emaxpos+eopreg, CHECK IF THIS IS WHAT YOU WANT",/)' )
        WRITE( stdout,'(/,5x,"eopreg is nevertheless used for the smooth increase of the barrier",/)' )
  ENDIF
  IF ( (monopole) .AND. (block) ) THEN
     IF ((block_1<0.0) .OR. (block_1>1.0) .OR. (block_2<0.0) .OR. (block_2>1.0)) &
        CALL errore( 'iosys', 'Both block_1, block_2 have to be between wihtin ]0,1[' , 1 )
     IF (block_1>=block_2) &
        CALL errore( 'iosys', 'Wrong order of block_1, block_2, should be block_1<block_2' , 1 )
     ENDIF
  ENDIF
  !TB end
  IF ( (tefield.or.monopole) .and. tstress ) THEN !TB no stress with monopole
     lstres = .false.
     WRITE( stdout, &
            '(5x,"Presently stress not available with electric field and monopole",/)' )
  ENDIF
  !TB Why no E-field with SOC?
  !IF ( tefield .and. ( nspin > 2 ) ) THEN
  IF ( (tefield .and. ( nspin > 2 )) .and. (.not.monopole) ) THEN
     CALL errore( 'iosys', 'LSDA not available with electric field' , 1 )
  ENDIF
  !
  ! ... define memory- and disk-related internal switches
  !
  smallmem = ( TRIM( memory ) == 'small' )
  twfcollect = wf_collect
  !
  ! ... Set Values for electron and bands
  !
  tfixed_occ = .false.
  ltetra     = .false.
  lgauss     = .false.
  ngauss     = 0
  !
  SELECT CASE( trim( occupations ) )
  CASE( 'fixed' )
     !
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
     ltetra = .true.
     tetra_type = 0
     !
  CASE( 'tetrahedra_lin', 'tetrahedra-lin')
     !
     ltetra = .true.
     tetra_type = 1
     !
  CASE('tetrahedra_opt', 'tetrahedra-opt')
     !
     ltetra = .true.
     tetra_type = 2
     !
  CASE( 'from_input' )
     !
     ngauss     = 0
     tfixed_occ = .true.
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys','occupations ' // trim( occupations ) // &
                & ' not implemented', 1 )
     !
  END SELECT
  !
  IF( ltetra ) THEN
     IF( lforce ) CALL infomsg( 'iosys', &
       'BEWARE:  force calculation with tetrahedra (not recommanded)')
     IF( lstres ) CALL infomsg( 'iosys', &
       'BEWARE: stress calculation with tetrahedra (not recommanded)')
  END IF
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
     IF ( noncolin ) CALL errore( 'iosys', &
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
  IF ( lda_plus_u .AND. lda_plus_u_kind == 0 .AND. noncolin ) THEN
     CALL errore('iosys', 'simplified LDA+U not implemented with &
                          &noncol. magnetism, use lda_plus_u_kind = 1', 1)
  END IF
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
     IF ( TRIM(startingwfc) /= 'file' ) THEN
        CALL infomsg('input','WARNING: "startingwfc" set to '//TRIM(startingwfc)//' may spoil restart')
     END IF
     IF ( TRIM(startingpot) /= 'file' ) THEN
        CALL infomsg('input','WARNING: "startingpot" set to '//TRIM(startingpot)//' may spoil restart')
        startingpot = 'file'
     END IF
     IF ( trim( ion_positions ) == 'from_input' ) THEN
        CALL infomsg('input','WARNING: restarting from positions as given in input')
        startingconfig = 'input'
     ELSE
        startingconfig = 'file'
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
  CASE ( 'medium' )
     !
     io_level = 1
     !
  CASE ( 'low' )
     !
     io_level = 0
     !
  CASE ( 'none' )
     !
     io_level = -1
     IF ( twfcollect ) THEN
        CALL infomsg('iosys', 'minimal I/O required, wf_collect reset to FALSE')
        twfcollect= .false.
     ENDIF
     !
  CASE DEFAULT
     !
     ! In the scf case, it is usually convenient to write to RAM;
     ! otherwise it is preferrable to write to disk, since the number
     ! of k-points can be large, leading to large RAM requirements
     !
     IF ( lscf ) THEN
        io_level = 0
     ELSE
        io_level = 1
     END IF
     !
  END SELECT
  !
  Hubbard_U(:)    = Hubbard_U(:) / rytoev
  Hubbard_J0(:)   = Hubbard_J0(:) / rytoev
  Hubbard_J(:,:)  = Hubbard_J(:,:) / rytoev
  Hubbard_alpha(:)= Hubbard_alpha(:) / rytoev
  Hubbard_beta(:) = Hubbard_beta(:) / rytoev
  !
  ethr = diago_thr_init
  !
  IF ( startingpot /= 'atomic' .and. startingpot /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingpot: use default (1)' )
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
     startingpot = 'file'
     !
  ENDIF
  !
  IF (      startingwfc /= 'atomic' .and. &
            startingwfc /= 'random' .and. &
            startingwfc /= 'atomic+random' .and. &
            startingwfc /= 'file' ) THEN
     !
     CALL infomsg( 'iosys', 'wrong startingwfc: use default (atomic+random)' )
     startingwfc = 'atomic+random'
     !
  ENDIF
  ! 
  IF (one_atom_occupations .and. startingwfc /= 'atomic' ) THEN
     CALL infomsg( 'iosys', 'one_atom_occupations requires startingwfc atomic' )
     startingwfc = 'atomic'
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
  adapt_thr = adaptive_thr
  tr2_init  = conv_thr_init
  tr2_multi = conv_thr_multi
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
     imix = 0
  CASE( 'TF' )
     imix = 1
  CASE( 'local-TF' )
     imix = 2
  CASE( 'potential' )
     CALL errore( 'iosys', 'potential mixing no longer implemented', 1 )
  CASE DEFAULT
     CALL errore( 'iosys', 'unknown mixing ' // trim( mixing_mode ), 1 )
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
     iverbosity = 1
  CASE( 'low', 'default', 'minimal' )
     iverbosity = 0 
  CASE DEFAULT
     iverbosity = 0
  END SELECT
  !
  IF ( lberry .OR. lelfield .OR. lorbm ) THEN
     IF ( npool > 1 ) CALL errore( 'iosys', &
          'Berry Phase/electric fields not implemented with pools', 1 )
     IF ( lgauss .OR. ltetra ) CALL errore( 'iosys', &
          'Berry Phase/electric fields only for insulators!', 1 )
     IF ( lmovecell ) CALL errore( 'iosys', &
          'Berry Phase/electric fields not implemented with variable cell', 1 )
  END IF
  !
  ! ... Copy values from input module to PW internals
  !
  nppstr_     = nppstr
  gdir_       = gdir
  lberry_     = lberry
  lelfield_   = lelfield
  lorbm_      = lorbm
  efield_     = efield
  nberrycyc_  = nberrycyc
  efield_cart_ = efield_cart
  SELECT CASE(efield_phase)
     CASE( 'none' )
        phase_control=0
     CASE ('write')
        phase_control=1
     CASE ('read')
        phase_control=2
     CASE DEFAULT
        CALL errore( 'iosys', &
          'Unknown efield_phase', 1 )
  END SELECT
  tqr_        = tqr
  real_space_ = real_space
  !
  tq_smoothing_ = tq_smoothing
  tbeta_smoothing_ = tbeta_smoothing
  !
  title_      = title
  lkpoint_dir_=lkpoint_dir
  dt_         = dt
  tefield_    = tefield
  dipfield_   = dipfield
  !TB start
  monopole_   = monopole
  zmon_    = zmon
  relaxz_  = relaxz
  block_   = block
  block_1_ = block_1
  block_2_ = block_2
  block_height_ = block_height
  !TB end
  prefix_     = trim( prefix )
  pseudo_dir_ = trimcheck( pseudo_dir )
  nstep_      = nstep
  iprint_     = iprint
  lecrpa_     = lecrpa
  scf_must_converge_ = scf_must_converge
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
  lforcet_ = lforcet
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
  Hubbard_J_(1:3,1:ntyp)  = hubbard_j(1:3,1:ntyp)
  Hubbard_J0_(1:ntyp)     = hubbard_j0(1:ntyp)
  Hubbard_alpha_(1:ntyp)  = hubbard_alpha(1:ntyp)
  Hubbard_beta_(1:ntyp)   = hubbard_beta(1:ntyp)
  lda_plus_u_             = lda_plus_u
  lda_plus_u_kind_        = lda_plus_u_kind
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
  IF (trim(occupations) /= 'from_input') one_atom_occupations_=.false.
  !
  !  ... initialize variables for vdW (dispersions) corrections
  !
  SELECT CASE( TRIM( vdw_corr ) )
    !
    CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
      !
      llondon= .TRUE.
      ts_vdw_= .FALSE.
      lxdm   = .FALSE.
      !
    CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
      !
      llondon= .FALSE.
      ts_vdw_= .TRUE.
      lxdm   = .FALSE.
      !
    CASE( 'XDM', 'xdm' )
       !
      llondon= .FALSE.
      ts_vdw_= .FALSE.
      lxdm   = .TRUE.
      !
    CASE DEFAULT
      !
      llondon= .FALSE.
      ts_vdw_= .FALSE.
      lxdm   = .FALSE.
      !
  END SELECT
  IF ( london ) THEN
     CALL infomsg("iosys","london is obsolete, use ""vdw_corr='grimme-d2'"" instead")
     vdw_corr='grimme-d2'
     llondon = .TRUE.
  END IF
  IF ( xdm ) THEN
     CALL infomsg("iosys","xdm is obsolete, use ""vdw_corr='xdm'"" instead")
     vdw_corr='xdm'
     lxdm = .TRUE.
  END IF
  IF ( ts_vdw ) THEN
     CALL infomsg("iosys","ts_vdw is obsolete, use ""vdw_corr='TS'"" instead")
     vdw_corr='TS'
     ts_vdw_ = .TRUE.
  END IF
  IF ( llondon.AND.lxdm .OR. llondon.AND.ts_vdw_ .OR. lxdm.AND.ts_vdw_ ) &
     CALL errore("iosys","must choose a unique vdW correction!", 1)
  !
  IF ( llondon) THEN
     lon_rcut    = london_rcut
     scal6       = london_s6
     in_c6(:)    = london_c6(:)
     in_rvdw(:)  = london_rvdw(:)
  END IF
  IF ( lxdm ) THEN
     a1i = xdm_a1
     a2i = xdm_a2
  END IF
  IF ( ts_vdw_ ) THEN
     vdw_isolated = ts_vdw_isolated
     vdw_econv_thr= ts_vdw_econv_thr
  END IF
  !
  !  calculate all the atomic positions if only the inequivalent ones
  !  have been given.
  !  NB: ibrav is an output of this routine
  !
  IF (space_group /= 0 .AND. .NOT. lsg ) &
     CALL errore('input','space_group requires crystal_sg atomic &
                                                   &coordinates',1 )
  IF (lsg) THEN
     IF (space_group==0) &
        CALL errore('input','The option crystal_sg requires the space group &
                                                   &number',1 )
     CALL sup_spacegroup(rd_pos,sp_pos,rd_for,rd_if_pos,space_group,nat,&
              uniqueb,rhombohedral,origin_choice,ibrav_sg)
     spacegroup = space_group
     IF (ibrav==-1) THEN
        ibrav=ibrav_sg
     ELSEIF (ibrav /= ibrav_sg) THEN
        CALL errore ('input','Input ibrav not compatible with space group &
                                                   &number',1 )
     ENDIF
     nat_=nattot
  ENDIF
  !
  ! QM/MM specific parameters
  !
  IF (.NOT. tqmmm) CALL qmmm_config( mode=-1 )
  !
  do_makov_payne  = .false.
  do_comp_mt      = .false.
  do_comp_esm     = .false.
  !
  SELECT CASE( trim( assume_isolated ) )
      !
    CASE( 'makov-payne', 'm-p', 'mp' )
      !
      do_makov_payne = .true.
      IF ( ibrav < 1 .OR. ibrav > 3 ) CALL errore(' iosys', &
              'Makov-Payne correction defined only for cubic lattices', 1)
      !
    CASE( 'dcc' )
      !
      CALL errore('iosys','density countercharge correction currently disabled',1)
      !
    CASE( 'martyna-tuckerman', 'm-t', 'mt' )
      !
      do_comp_mt     = .true.
      !
    CASE( 'esm' )
      !
      do_comp_esm    = .true.
      !
  END SELECT
  !
  IF ( do_comp_mt .AND. lstres ) THEN
     lstres = .false.
     WRITE( stdout, &
          '(5x,"Stress calculation not meaningful in isolated systems",/)' )
  END IF
  !
  ! ... ESM
  !
  esm_bc_ = esm_bc
  esm_efield_ = esm_efield
  esm_w_ = esm_w
  esm_nfit_ = esm_nfit 
  esm_a_ = esm_a
  !
  IF ( esm_bc .EQ. 'bc4' ) THEN
    IF ( ABS(esm_w) .LT. 1.D-8 ) THEN
      CALL errore ('iosys','esm_w too small',1)
    ELSEIF ( esm_w .GT. 0.D0 ) THEN
      CALL errore ('iosys','positive esm_w not allowed for bc4',1)
    ENDIF
    IF ( esm_a .LT. 1.D-4 ) THEN
      CALL errore ('iosys','smoothness parameter for bc4 too small',1)
    ELSEIF ( esm_a .GT. 10.D0 ) THEN
      CALL errore ('iosys','smoothness parameter for bc4 too big',1)
    ENDIF
  ENDIF
  !
  ! ... FCP
  !
  lfcpopt_        = lfcpopt
  lfcpdyn_        = lfcpdyn
  fcp_mu_         = fcp_mu
  fcp_mass_       = fcp_mass
  fcp_temperature = fcp_tempw
  !
  IF ( lfcpopt .or. lfcpdyn ) THEN
     IF ( .not. do_comp_esm ) THEN
        CALL errore ('iosys','FCP optimise/dynamics currently not available without ESM',1)
     ENDIF
     IF ( trim( calculation ).NE.'relax'.AND.trim( calculation ).NE.'md')THEN
        CALL errore ('iosys',"FCP optimise/dynamics only available with calculation = 'relax' and 'md'",1)
     ENDIF
  ENDIF
  !
  IF ( fcp_temperature == 0.0_DP ) &
     fcp_temperature = temperature
  fcp_relax_step_ = fcp_relax_step
  fcp_relax_crit_ = fcp_relax_crit
  !
  CALL plugin_read_input()
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
  IF ( monopole ) ALLOCATE( forcemono( 3, nat_ ) ) !TB monopole forces
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
  ! ... Files (for compatibility) and directories
  !     This stuff must be done before calling read_config_from_file!
  !
  input_drho  = ' '
  output_drho = ' '
  tmp_dir = trimcheck ( outdir )
  IF ( .not. trim( wfcdir ) == 'undefined' ) THEN
     wfc_dir = trimcheck ( wfcdir )
  ELSE
     wfc_dir = tmp_dir
  ENDIF
  !
  at_old    = at
  omega_old = omega
  !
  ! ... Read atomic positions and unit cell from data file, if needed,
  ! ... overwriting what has just been read before from input
  !
  ierr = 1
  IF ( startingconfig == 'file' .AND. .NOT. lforcet ) &
     ierr = read_config_from_file(nat, at_old, omega_old, lmovecell, &
                                       at, bg, omega, tau)
  !
  ! ... read_config_from_file returns 0 if structure successfully read
  ! ... Atomic positions (tau) must be converted to internal units
  ! ... only if they were read from input, not from file
  !
  IF ( ierr /= 0 ) CALL convert_tau ( tau_format, nat_, tau)
  !
  ! ... set up k-points
  !
  CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, nkstot, xk, wk )
  gamma_only = ( k_points == 'gamma' )
  !
!  IF ( real_space .AND. .NOT. gamma_only ) &
!     CALL errore ('iosys', 'Real space only with Gamma point', 1)
  IF ( lelfield .AND. gamma_only ) &
      CALL errore( 'iosys', 'electric fields not available for k=0 only', 1 )
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
     wmass = wmass * amu_ry
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
     cmass  = wmass * amu_ry
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
  ! ... read pseudopotentials (also sets DFT and a few more variables)
  ! ... returns values read from PP files into ecutwfc_pp, ecutrho_pp
  !
  CALL readpp ( input_dft, .FALSE., ecutwfc_pp, ecutrho_pp )
  CALL set_cutoff ( ecutwfc, ecutrho, ecutwfc_pp, ecutrho_pp )
  !
  ! ... set parameters of hybrid functionals
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
  IF(ecutfock <= 0.0_DP) THEN
     ! default case
     ecutfock_ = 4.0_DP*ecutwfc
  ELSE
     IF(ecutfock < ecutwfc .OR. ecutfock > ecutrho) CALL errore('iosys', &
          'ecutfock can not be < ecutwfc or > ecutrho!', 1) 
     ecutfock_ = ecutfock
  END IF
  IF ( lstres .AND. dft_is_hybrid() .AND. npool > 1 )  CALL errore('iosys', &
         'stress for hybrid functionals not available with pools', 1)
  IF ( lmovecell.AND. dft_is_hybrid() ) CALL errore('iosys',&
         'Variable cell and hybrid XC not tested',1)
  !
  ! ... must be done AFTER dft is read from PP files and initialized
  ! ... or else the two following parameters will be overwritten
  !
  IF (exx_fraction >= 0.0_DP) CALL set_exx_fraction (exx_fraction)
  IF (screening_parameter >= 0.0_DP) &
        & CALL set_screening_parameter (screening_parameter)
  !
  ! ... read the vdw kernel table if needed
  !
  vdw_table_name_  = vdw_table_name
  inlc = get_inlc()
  IF (inlc > 0) CALL initialize_kernel_table(inlc)
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
     ! The next two lines have been moved before the call to read_config_from_file:
     !      at_old    = at
     !      omega_old = omega
     IF ( cell_factor_ <= 0.0_dp ) cell_factor_ = 2.0_dp
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
  IF ( lxdm) CALL init_xdm ( )
  !
  ! ... variables for constrained dynamics are set here
  !
  lconstrain = ( nconstr_inp > 0 )
  !
  IF ( lconstrain ) THEN
     IF ( lbfgs .OR. lmovecell ) CALL errore( 'iosys', &
              'constraints only with fixed-cell dynamics', 1 )
     CALL init_constraint( nat, tau, ityp, alat )
  END IF
  !
  ! ... End of reading input parameters
  !
  CALL pw_init_qexsd_input(qexsd_input_obj, obj_tagname="input")
  CALL deallocate_input_parameters ()  
  !
  ! ... Initialize temporary directory(-ies)
  !
  CALL check_tempdir ( tmp_dir, exst, parallelfs )
  IF ( .NOT. exst .AND. restart ) THEN
     CALL infomsg('iosys', 'restart disabled: needed files not found')
     restart = .false.
  ELSE IF ( .NOT. exst .AND. (lbands .OR. .NOT. lscf) ) THEN
     CALL errore('iosys', 'bands or non-scf calculation not possible: ' // &
                          'needed files are missing', 1)
  ELSE IF ( exst .AND. .NOT.restart ) THEN
     CALL clean_tempdir ( tmp_dir )
  END IF
  IF ( TRIM(wfc_dir) /= TRIM(tmp_dir) ) &
     CALL check_tempdir( wfc_dir, exst, parallelfs )

  RETURN
  !
END SUBROUTINE iosys
!
SUBROUTINE set_cutoff ( ecutwfc_in, ecutrho_in, ecutwfc_pp, ecutrho_pp )
  !
  ! Copy to modules the cutoffs, either read from input or from PP files
  ! Values of ecutwfc and ecutrho are returned in ecutwfc_in, ecutrho_in
  !
  USE kinds, ONLY : dp
  USE gvecs, ONLY : dual
  USE gvect, ONLY : ecutrho
  USE gvecw, ONLY : ecutwfc
  !
  IMPLICIT NONE
  REAL(dp), INTENT(INOUT) :: ecutwfc_in, ecutrho_in
  REAL(dp), INTENT(IN)    :: ecutwfc_pp, ecutrho_pp
  !
  IF( ecutwfc_in <= 0.0_dp ) THEN
     IF( ecutwfc_pp <= 0.0_DP ) THEN
        CALL errore( 'set_cutoff' ,' ecutwfc not set ',1)
     ELSE
        ecutwfc = ecutwfc_pp
     END IF
  ELSE
     ecutwfc = ecutwfc_in
  END IF
  IF( ecutrho_in <= 0.0_dp ) THEN
     IF( ecutwfc_in > 0.0_dp ) THEN
        ecutrho = 4.0_dp*ecutwfc_in
     ELSE IF( ecutrho_pp > 0.0_dp ) THEN
        ecutrho = ecutrho_pp
     ELSE IF( ecutwfc_pp > 0.0_dp ) THEN
        ecutrho = 4.0_dp*ecutwfc_pp
     END IF
  ELSE
     ecutrho = ecutrho_in
  ENDIF
  ecutwfc_in = ecutwfc
  ecutrho_in = ecutrho
  dual = ecutrho / ecutwfc
  IF ( dual <= 1.0_dp ) CALL errore( 'set_cutoff', 'ecutrho <= ecutwfc?!?', 1 )
  !
END SUBROUTINE set_cutoff
!
!----------------------------------------------------------------------------
SUBROUTINE read_cards_pw ( psfile, tau_format )
  !----------------------------------------------------------------------------
  !
  USE kinds,              ONLY : DP
  USE input_parameters,   ONLY : atom_label, atom_pfile, atom_mass, taspc, &
                                 tapos, rd_pos, atomic_positions, if_pos,  &
                                 sp_pos, f_inp, rd_for, tavel, sp_vel, rd_vel, &
                                 lsg
  USE dynamics_module,    ONLY : vel
  USE cell_base,          ONLY : at, ibrav
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp, tau, atm, extfor
  USE fixed_occ,          ONLY : tfixed_occ, f_inp_ => f_inp
  USE ions_base,          ONLY : if_pos_ =>  if_pos, amass, fixatom
  USE control_flags,      ONLY : textfor, tv0rd
  USE wyckoff,            ONLY : nattot, tautot, ityptot, extfortot, &
                                 if_postot, clean_spacegroup
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
  IF (lsg) THEN
     tau(:,:)=tautot(:,:)
     ityp(:) = ityptot(:)
     extfor(:,:) = extfortot(:,:)
     if_pos_(:,:) = if_postot(:,:)
     CALL clean_spacegroup()
  ELSE 
     DO ia = 1, nat
        !
        tau(:,ia) = rd_pos(:,ia)
        ityp(ia)  = sp_pos(ia)
        extfor(:,ia) = rd_for(:,ia)
        if_pos_(:,ia) = if_pos(:,ia)
        !
     ENDDO
  ENDIF
  !
  ! ... check for initial velocities read from input file
  !
  IF ( tavel .AND. ANY ( sp_pos(:) /= sp_vel(:) ) ) &
      CALL errore("cards","list of species in block ATOMIC_VELOCITIES &
                 & must be identical to those in ATOMIC_POSITIONS",1)
  tv0rd = tavel
  IF ( tv0rd ) THEN
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
  fixatom = COUNT( if_pos_(1,:)==0 .AND. if_pos_(2,:)==0 .AND. if_pos_(3,:)==0 )
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
SUBROUTINE check_tempdir ( tmp_dir, exst, pfs )
  !-----------------------------------------------------------------------
  !
  ! ... Verify if tmp_dir exists, creates it if not
  ! ... On output:
  ! ...    exst= .t. if tmp_dir exists
  ! ...    pfs = .t. if tmp_dir visible from all procs of an image
  !
  USE wrappers,      ONLY : f_mkdir_safe
  USE io_global,     ONLY : ionode, ionode_id
  USE mp_images,     ONLY : intra_image_comm, nproc_image, me_image
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: tmp_dir
  LOGICAL, INTENT(out)         :: exst, pfs
  !
  INTEGER             :: ios, image, proc, nofi
  CHARACTER (len=256) :: file_path, filename
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  ! ... create tmp_dir on ionode
  ! ... f_mkdir_safe returns -1 if tmp_dir already exists
  ! ...                       0 if         created
  ! ...                       1 if         cannot be created
  !
  IF ( ionode ) ios = f_mkdir_safe( TRIM(tmp_dir) )
  CALL mp_bcast ( ios, ionode_id, intra_image_comm )
  exst = ( ios == -1 )
  IF ( ios > 0 ) CALL errore ('check_tempdir','tmp_dir cannot be opened',1)
  !
  ! ... let us check now if tmp_dir is visible on all nodes
  ! ... if not, a local tmp_dir is created on each node
  !
  ios = f_mkdir_safe( TRIM(tmp_dir) )
  CALL mp_sum ( ios, intra_image_comm )
  pfs = ( ios == -nproc_image ) ! actually this is true only if .not.exst 
  !
  RETURN
  !
END SUBROUTINE check_tempdir
!
!-----------------------------------------------------------------------
SUBROUTINE clean_tempdir( tmp_dir )
  !-----------------------------------------------------------------------
  !
  USE io_files,         ONLY : prefix, delete_if_present
  USE io_global,        ONLY : ionode
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: tmp_dir
  !
  CHARACTER (len=256) :: file_path, filename
  !
  ! ... remove temporary files from tmp_dir ( only by the master node )
  !
  file_path = trim( tmp_dir ) // trim( prefix )
  IF ( ionode ) THEN
     CALL delete_if_present( trim( file_path ) // '.update' )
     CALL delete_if_present( trim( file_path ) // '.md' )
     CALL delete_if_present( trim( file_path ) // '.bfgs' )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE clean_tempdir
