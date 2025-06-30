
! Copyright (C) 2002-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE iosys()
  !-----------------------------------------------------------------------------
  !! Copy data read from input file (subroutine \(\tex/temttt{read_input_file}\)
  !! and stored in modules input_parameters into internal modules.  
  !! Wrapper routine: the original "iosys" was too long and was split into
  !! more *iosys* routines, containing input_parameters module
  !
  USE input_parameters,      ONLY : deallocate_input_parameters
  !
  USE qexsd_input,           ONLY : qexsd_input_obj
  USE qes_types_module,      ONLY : input_type
#if defined (__ENVIRON)
  USE plugin_flags,          ONLY : use_environ
  USE environ_base_module,   ONLY : read_environ_input, init_environ_setup
#endif
#if defined (__OSCDFT)
  USE plugin_flags,          ONLY : use_oscdft
  USE oscdft_base,           ONLY : oscdft_ctx
  USE oscdft_input,          ONLY : oscdft_read_input
#endif
 USE lsda_mod,               ONLY : starting_magnetization
 USE ions_base,              ONLY : nsp
  INTERFACE
     SUBROUTINE   pw_init_qexsd_input(obj,obj_tagname)
     IMPORT                       :: input_type
     TYPE(input_type)             :: obj
     CHARACTER(LEN=*),INTENT(IN)  :: obj_tagname
     END SUBROUTINE
  END INTERFACE
  !
  !
  ! ... copy control variables to internal variables
  !
  CALL control_iosys()
  !
  ! ... set temporary directory, check existence if restart required
  !
  CALL tmpdir_iosys()
  !
  ! ... read PP files, set DFT and cutoffs
  !
  CALL pseudo_iosys()
  !
  ! ... set up magnetization variables; done here because we need zv values
  !     initialized in pseudo_iosys 
  !
  CALL magnetization_iosys() 
  ! 
  ! ... set atomic structure (read from xml file if required)
  !
  CALL structure_iosys()
  !
  ! ... read optional plugin input files
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_read_input('PW')
#endif 
#if defined (__ENVIRON)
  IF (use_environ) THEN
     CALL read_environ_input()
     CALL init_environ_setup('PW')
  END IF
#endif
#if defined (__OSCDFT)
  IF (use_oscdft) THEN
     CALL oscdft_read_input(oscdft_ctx%inp)
  END IF
#endif
  !
  ! ... final reading and setup of variables depending upon previous routines
  ! ... (e.g. depending upon control variables, cell, positions, cutoffs, ...)
  !
  CALL iosys_end()
  !
  ! ... End of processing input parameters - write to xml file and free memory
  !
#if ! defined (__INTEL_COMPILER) || __INTEL_COMPILER >= 1300
  CALL pw_init_qexsd_input(qexsd_input_obj, obj_tagname="input")
#endif
  CALL deallocate_input_parameters ()  
  !
END SUBROUTINE iosys
!
!----------------------------------------------------------------------------
SUBROUTINE tmpdir_iosys()
  !-------------------------------------------------------------------------
  !!  Set and check temporary directory
  !!  Must be called before any I/O operations are performed 
  !!  Input : outdir, wfc_dir, prefix, plus control flags restart, lbands, lscf
  !!  Output: module io_files, restart
  !
  USE input_parameters, ONLY : outdir, wfcdir, prefix_ => prefix
  USE io_files,         ONLY : tmp_dir, wfc_dir, prefix, &
                               check_tempdir, clean_tempdir, nd_nmbr, &
                               pseudo_dir_cur, restart_dir
  !
  USE control_flags,    ONLY : restart, lbands, lscf
  IMPLICIT NONE
  !
  LOGICAL :: exst, parallelfs
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  !
  !! outdir = TRIM(outdir) // TRIM(nd_nmbr)
  !!   if previous line is uncommented, each process sees a different directory
  !!   (the process number is added at the end of tmp_dir)
  !!   For testing purposes only; works only if outdir does not end with '/' 
  ! 
  tmp_dir = trimcheck ( outdir )
  IF ( .not. trim( wfcdir ) == 'undefined' ) THEN
     wfc_dir = trimcheck ( wfcdir )
  ELSE
     wfc_dir = tmp_dir
  ENDIF
  prefix = trim( prefix_ )
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
  !
  ! If restarting or making a non-scf calculation:
  ! read the PP files from the .save directory of the previous run
  !
  IF ( restart .OR. lbands .OR. .NOT. lscf ) pseudo_dir_cur = restart_dir()
  !
END SUBROUTINE tmpdir_iosys
!
!----------------------------------------------------------------------------
SUBROUTINE pseudo_iosys()
  !-----------------------------------------------------------------------------
  !
  !! ... read PP files, set DFT and cutoffs
  !! ... needed input variables are in module input_parameters
  !! ... variables nsp and psfile set here are needed by readpp
  !! ... variables ecutwfc and ecutrho in modules are set here (see set_cutoff)
  !
  USE input_parameters, ONLY : input_dft, atom_pfile, &
                               pseudo_dir_ => pseudo_dir, &
                               ecutwfc, ecutrho, ntyp, &
                               nr1, nr2, nr3, nr1s, nr2s, nr3s
  USE kinds,            ONLY : dp
  USE io_files,         ONLY : psfile, pseudo_dir
  USE ions_base,        ONLY : nsp
  USE read_pseudo_mod,  ONLY : readpp
  !
  IMPLICIT NONE
  REAL(dp) :: ecutwfc_pp, ecutrho_pp
  INTEGER :: is
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  nsp = ntyp
  pseudo_dir = trimcheck( pseudo_dir_ )
  DO is = 1, nsp
     psfile(is) = atom_pfile(is)
  ENDDO
  !
  ! ... read pseudopotentials (also sets DFT and a few more variables)
  ! ... returns values read from PP files into ecutwfc_pp, ecutrho_pp
  !
  CALL readpp ( input_dft, .FALSE., ecutwfc_pp, ecutrho_pp )
  !
  CALL set_cutoff ( ecutwfc, ecutrho, ecutwfc_pp, ecutrho_pp, &
       nr1,nr2,nr3, nr1s,nr2s,nr3s )
  !
END SUBROUTINE pseudo_iosys
!
!----------------------------------------------------------------------------
SUBROUTINE control_iosys()
  !-----------------------------------------------------------------------------
  !! Copy data read from input file (in subroutine \(\texttt{read_input_file}\)
  !! and stored in modules input_parameters into internal modules.  
  !! Note that many variables in internal modules, having the same name as
  !! those in input_parameters, are locally renamed by adding an underscore "_".
  !
  USE kinds,         ONLY : DP
  !  
  USE control_flags, ONLY : adapt_thr, tr2_init, tr2_multi  
  !
  USE constants,     ONLY : pi, ry_kbar, eps8
  !
  USE ldaU,          ONLY : niter_with_fixed_ns

  USE xc_lib,        ONLY:  xclib_dft_is
  !
  USE io_global,     ONLY : stdout
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
  USE cell_base,     ONLY : init_dofree, &
                            press_       => press, &
                            wmass_       => wmass
  !
  USE ions_base,     ONLY : if_pos, ityp, tau, extfor, atm, nat, nsp, &
                            amass
  !
  USE starting_scf,  ONLY : startingconfig, starting_wfc, starting_pot
  USE starting_scf,  ONLY : starting_charge_ => starting_charge
  !
  USE run_info,      ONLY : title_ => title
  !
  USE cellmd,        ONLY : ntcheck, calc, lmovecell, &
                            cell_factor_ => cell_factor
  !
  USE dynamics_module, ONLY : control_temp, temperature, thermostat, &
                              dt_         => dt, &
                              delta_t_    => delta_t, &
                              nraise_     => nraise, &
                              refold_pos_ => refold_pos, &
                              fire_nmin_ => fire_nmin, &
                              fire_f_inc_ => fire_f_inc, &
                              fire_f_dec_ => fire_f_dec,  &
                              fire_alpha_init_ => fire_alpha_init, &  
                              fire_falpha_ => fire_falpha, &
                              fire_dtmax_ => fire_dtmax
  !
  USE extfield,      ONLY : tefield_  => tefield, &
                            dipfield_ => dipfield, &
                            edir_     => edir, &
                            emaxpos_  => emaxpos, &
                            eopreg_   => eopreg, &
                            eamp_     => eamp, &
  ! TB added gate related variables
                            zgate_    => zgate, &
                            gate_     => gate, &
                            relaxz_   => relaxz, &
                            block_    => block, &
                            block_1_   => block_1, &
                            block_2_   => block_2, &
                            block_height_   => block_height
  !
  USE klist,         ONLY : ltetra, lgauss, ngauss, two_fermi_energies, &
                            smearing_          => smearing, &
                            degauss_           => degauss, &
                            tot_charge_        => tot_charge, &
                            tot_magnetization_ => tot_magnetization, &
                            degauss_cond_      => degauss_cond, &
                            nelec_cond_        => nelec_cond
  USE ktetra,        ONLY : tetra_type
  !
  USE add_dmft_occ,  ONLY : dmft_ => dmft, &
                            dmft_prefix_ => dmft_prefix
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
  USE lsda_mod,      ONLY : nspin_                  => nspin, &
                            lsda
  !
  USE relax,         ONLY : epse, epsf, epsp, starting_scf_threshold
  !
  USE control_flags, ONLY : sic, scissor
  USE sic_mod,       ONLY : pol_type_ => pol_type, sic_gamma_ => sic_gamma, &
                            sic_energy_ => sic_energy
  USE sci_mod,       ONLY : sci_vb_ => sci_vb, sci_cb_ => sci_cb
 
  !
  USE extrapolation, ONLY : pot_order, wfc_order
  USE control_flags, ONLY : isolve, max_cg_iter, david, &
                            rmm_ndim, rmm_conv, gs_nblock, rmm_with_davidson, &
                            tr2, imix, gamma_only, tnosep, tnoseh, &
                            nmix, iverbosity, smallmem, nexxiter, niter, &
                            io_level, ethr, lscf, lbfgs, lmd, &
                            lbands, lconstrain, restart, &
                            llondon, ldftd3, do_makov_payne, lxdm, &
                            lensemb, lforce   => tprnfor, &
                            tstress_          => tstress, &
                            remove_rigid_rot_ => remove_rigid_rot, &
                            diago_full_acc_   => diago_full_acc, &
                            tolp_             => tolp, &
                            upscale_          => upscale, &
                            mixing_beta_      => mixing_beta, &
                            nstep_            => nstep, &
                            iprint_           => iprint, &
                            noinv_            => noinv, &
                            tqr_              => tqr, &
                            tq_smoothing_     => tq_smoothing, &
                            tbeta_smoothing_  => tbeta_smoothing, &
                            ts_vdw_           => ts_vdw, &
                            mbd_vdw_          => mbd_vdw, &
                            lecrpa_           => lecrpa, &
                            scf_must_converge_=> scf_must_converge, & 
                            treinit_gvecs_    => treinit_gvecs, &  
                            max_xml_steps_    => max_xml_steps, & 
                            use_spinflip_    => use_spinflip, & 
                            symm_by_label 
  USE check_stop,    ONLY : max_seconds_ => max_seconds
  !
  USE wvfct,         ONLY : nbnd_ => nbnd, &
                            nbnd_cond_ => nbnd_cond              
  !
  USE two_chem,      ONLY : twochem_ => twochem
  !  
  USE gvecw,         ONLY : ecfixed_ => ecfixed, &
                            qcutz_   => qcutz, &
                            q2sigma_ => q2sigma
  !
  USE fixed_occ,     ONLY : tfixed_occ, f_inp_ => f_inp, &
                            one_atom_occupations_ => one_atom_occupations
  !
  USE noncollin_module, ONLY : i_cons, mcons, bfield, &
                               noncolin_  => noncolin, &
                               lambda_    => lambda, &
                               report_    => report, &
                               lspinorb_ => lspinorb,  &
                               lforcet_ => lforcet,    &
                               starting_spin_angle_ => starting_spin_angle
  !
  USE symm_base, ONLY : no_t_rev_ => no_t_rev, nofrac, allfrac, &
                        nosym_ => nosym, nosym_evc_=> nosym_evc
  !
  USE bfgs_module,   ONLY : init_bfgs
  !
  USE wannier_new, ONLY :   use_wannier_      => use_wannier, &
                            use_energy_int_   => use_energy_int, &
                            nwan_             => nwan, &
                            print_wannier_coeff_    => print_wannier_coeff

  USE Coul_cut_2D,  ONLY :  do_cutoff_2D 

  USE realus,                ONLY : real_space_ => real_space

  USE qmmm,                  ONLY : qmmm_config
  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : title, calculation, verbosity, restart_mode,    &
                               nstep, iprint, tstress, tprnfor, dt,            &
                               etot_conv_thr, forc_conv_thr,                   &
                               disk_io, tefield, dipfield, lberry,             &
                               gdir, nppstr, wf_collect,lelfield,lorbm,efield, &
                               nberrycyc, efield_cart, lecrpa,                 &
                               lfcp, vdw_table_name, memory, max_seconds,      &
                               tqmmm, efield_phase, gate, max_xml_steps,       &
                               trism, twochem, symmetry_with_labels,           & 
                               use_spinflip 

  !
  ! ... SYSTEM namelist
  !
  USE input_parameters, ONLY : ntyp, nbnd,tot_charge,tot_magnetization,     &
                               ecutwfc, ecutrho, noinv, nosym, nosym_evc,   &
                               no_t_rev, use_all_frac, force_symmorphic,    &
                               starting_charge, occupations, degauss,       &
                               smearing, nspin, ecfixed, qcutz, q2sigma,    &
                               la2F, dmft, dmft_prefix,                     &
                               pol_type, sic_gamma, sic_energy, sci_vb, sci_cb, &
                               edir, emaxpos, eopreg, eamp, noncolin, lambda, &
                               report, lspinorb, assume_isolated,        &
                               vdw_corr, london, london_s6, london_rcut, london_c6, &
                               london_rvdw, &
                               ts_vdw, ts_vdw_isolated, ts_vdw_econv_thr,     &
                               mbd_vdw,     &
                               xdm, xdm_a1, xdm_a2, lforcet,                  &
                               one_atom_occupations,                          &
                               esm_bc, esm_efield, esm_w, esm_nfit, esm_a,    &
                               lgcscf,                                        &
                               zgate, relaxz, block, block_1, block_2,        &
                               block_height, lgcscf, nbnd_cond, nelec_cond,   &
                               degauss_cond
  !
  ! ... ELECTRONS namelist
  !
  USE input_parameters, ONLY : exx_maxstep, electron_maxstep, mixing_mode, mixing_beta, &
                               mixing_ndim, mixing_fixed_ns, conv_thr,     &
                               tqr, tq_smoothing, tbeta_smoothing,         &
                               diago_thr_init,                             &
                               diago_cg_maxiter,                           &
                               diago_david_ndim, diago_rmm_ndim,           &
                               diago_rmm_conv, diago_gs_nblock,            &
                               diagonalization, diago_full_acc,            &
                               startingwfc, startingpot,                   &
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
                               trust_radius_ini, bfgs_ndim, tgdiis_step,       &
                               fire_nmin, fire_f_inc, fire_f_dec, &
                               fire_alpha_init, fire_falpha, fire_dtmax
!
! ... IONS NOSE HOOVER THERMOSTAT

  USE input_parameters, ONLY: fnosep, nhpcl, nhptyp, ndega, nhgrp, fnhscl  
                           
  !
  ! ... CELL namelist
  !
  USE input_parameters, ONLY : cell_parameters, cell_dynamics, press, wmass,  &
                               cell_temperature, cell_factor, press_conv_thr, &
                               cell_dofree, treinit_gvecs 
  !
  ! ... CELL NOSE HOOVER THERMOSTAT
  ! 
  USE input_parameters, ONLY: fnoseh, temph                                      
  !
  ! ... WANNIER_NEW namelist
  !
  USE input_parameters, ONLY : use_wannier, nwan, use_energy_int, print_wannier_coeff
  !
  ! ... CARDS
  !
  USE input_parameters,      ONLY : k_points, xk, wk, nk1, nk2, nk3,  &
                                    k1, k2, k3, nkstot
  USE input_parameters,      ONLY : nconstr_inp, f_inp
  !
  USE read_namelists_module, ONLY : sm_not_set
  USE london_module,         ONLY : lon_rcut, scal6, in_c6, in_rvdw
  USE xdm_module,            ONLY : a1i, a2i
  USE tsvdw_module,          ONLY : vdw_isolated, vdw_econv_thr
  !
  !
  IMPLICIT NONE
  !
  INTEGER  :: ia, nt, i, j
  LOGICAL  :: domag, sm_wasnt_set
  REAL(DP) :: theta, phi, V
  !
  ! MAIN CONTROL VARIABLES, MD AND RELAX
  !
  title_      = title
  lecrpa_     = lecrpa  
  !
  lforce    = tprnfor
  !
  SELECT CASE( trim( calculation ) )
  CASE( 'scf' )
     !
     lscf  = .true.
     nstep = min(1, nstep)
     !
  CASE( 'ensemble' )
     !
     lscf  = .true.
     lensemb = .true.
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
     CASE ( 'fire' )
        !
        lmd     = .true.
        calc    = 'fi'
        ! set fire variables
        fire_nmin_ = fire_nmin
        fire_f_inc_ = fire_f_inc
        fire_f_dec_ = fire_f_dec
        fire_alpha_init_ = fire_alpha_init
        fire_falpha_ = fire_falpha
        fire_dtmax_ = fire_dtmax
        !
        ntcheck = nstep + 1
        !
     CASE ( 'ipi' )
        !
        CONTINUE
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
     CASE ('velocity-verlet')
        !
        calc = 'wd'
        !
     CASE( 'langevin' )
        !
        calc        = 'ld'
        !
     CASE( 'langevin-smc', 'langevin+smc' )
        !
        calc        = 'ls'
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
     treinit_gvecs_ = treinit_gvecs 
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
     CASE ( 'ipi' )
        !
        CONTINUE
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
     CASE ( 'ipi' )
        !
        CONTINUE
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
  dt_    = dt
  nstep_ = nstep
  tstress_ = lmovecell .OR. ( tstress .and. lscf )
  !
  sic_gamma_ = sic_gamma
  sic_energy_ = sic_energy
  IF(sic_gamma /= 0.d0 ) sic = .true.
  pol_type_ = trim(pol_type)
  sci_vb_ = sci_vb
  sci_cb_ = sci_cb
  IF(sci_vb .NE. 0.d0 .or. sci_cb .NE. 0.d0 ) scissor = .true.
  IF(scissor .and. nspin .ne. 2) CALL errore('allocate_scissor', 'spin polarized calculation required',1)
  !
  ! ELECTRIC FIELDS (SAWTOOTH), GATE FIELDS
  !
  IF ( tefield .and. ( .not. nosym ) .and. ( .not. gate ) ) THEN
     nosym = .true.
     WRITE( stdout, &
            '(5x,"Presently no symmetry can be used with electric field",/)' )
  ENDIF
  !
  IF ( (gate) .AND. ( .NOT. nosym )) THEN
     WRITE( stdout,'(/,5x,"Presently symmetry can be used with gate field",/)' )
     WRITE( stdout,'(5x,"setting verbosity to high",/)' )
     WRITE( stdout,'(5x,"CAREFULLY CHECK ALL SYMMETRIES",/)' )
     verbosity='high'
  ENDIF
  IF ((zgate>1.0).OR.(zgate<0.0)) &
     CALL errore( 'iosys', 'Position of the charged plate representing the gate has to be between within ]0,1[' , 1 )
  IF ( (gate) .AND. ((tefield).OR.(dipfield)) ) THEN
     IF (edir .ne. 3) &
        CALL errore( 'iosys','Using gate and tefield/dipfield, edir must be 3', 1)
     IF ((zgate>=emaxpos).AND.(zgate<=(emaxpos+eopreg))) &
        CALL errore( 'iosys', 'Charged plate between the 2 dipole planes not allowed' , 1 )
     IF ((block).AND.(block_1.NE.emaxpos).AND.(block_2.NE.(emaxpos+eopreg))) THEN
        WRITE( stdout,'(/,5x,"Neither block_1=emaxpos, nor block_2=emaxpos+eopreg, CHECK IF THIS IS WHAT YOU WANT",/)' )
        WRITE( stdout,'(/,5x,"eopreg is nevertheless used for the smooth increase of the barrier",/)' )
  ENDIF
  IF ( (gate) .AND. (block) ) THEN
     IF ((block_1<0.0) .OR. (block_1>1.0) .OR. (block_2<0.0) .OR. (block_2>1.0)) &
        CALL errore( 'iosys', 'Both block_1, block_2 have to be between within ]0,1[' , 1 )
     IF (block_1>=block_2) &
        CALL errore( 'iosys', 'Wrong order of block_1, block_2: should be block_1<block_2' , 1 )
     ENDIF
  ENDIF
  !
  IF ( (tefield.or.gate) .and. tstress ) THEN
     tstress_ = .false.
     WRITE( stdout, &
            '(5x,"Presently stress not available with electric field and gates",/)' )
  ENDIF
  ! FIXME: is the following check correct?
  IF ( (tefield .and. ( nspin > 2 )) .and. (.not.gate) ) &
     CALL errore( 'iosys', 'LSDA not available with electric field' , 1 )
  tefield_ = tefield
  dipfield_= dipfield
  edir_    = edir
  emaxpos_ = emaxpos
  eopreg_  = eopreg
  eamp_    = eamp
  gate_    = gate
  zgate_   = zgate
  relaxz_  = relaxz
  block_   = block
  block_1_ = block_1
  block_2_ = block_2
  block_height_ = block_height
  !
  ! SPIN POLARIZATION
  !
  SELECT CASE( nspin )
  CASE( 1 )
     !
     lsda = .false.
     IF ( noncolin ) THEN 
       nspin = 4
       symm_by_label = symmetry_with_labels
     END IF 
     !
  CASE( 2 )
     !
     lsda = .true.
     IF ( noncolin ) CALL errore( 'iosys', &
                     'noncolin .and. nspin==2 are conflicting flags', 1 )
     symm_by_label = symmetry_with_labels 
     use_spinflip_ = use_spinflip
     !
  CASE( 4 )
     !
     lsda = .false.
     noncolin = .true.
     symm_by_label = symmetry_with_labels 
     !
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'wrong input value for nspin', 1 )
     !
  END SELECT
  nspin_  = nspin
  !
  ! OCCUPATIONS
  !
  CALL set_occupations( occupations, smearing, degauss, &
       tfixed_occ, ltetra, tetra_type, lgauss, ngauss )
  !
  degauss_ = degauss
  smearing_ = smearing
  degauss_cond_ = degauss_cond
  nelec_cond_ = nelec_cond
  !
  IF( ltetra ) THEN
     IF( lforce ) CALL infomsg( 'iosys', &
       'BEWARE:  force calculation with tetrahedra (not recommanded)')
     IF( tstress_ ) CALL infomsg( 'iosys', &
       'BEWARE: stress calculation with tetrahedra (not recommanded)')
  END IF
  IF( nbnd < 1 ) CALL errore( 'iosys', 'nbnd less than 1', nbnd ) 
  nbnd_    = nbnd
  nbnd_cond_ = nbnd_cond
  !
  two_fermi_energies = ( tot_magnetization /= -10000._DP)
  IF ( two_fermi_energies .and. tot_magnetization < -9999._DP) &
     CALL errore( 'iosys', 'tot_magnetization only larger than -9999 is allowed', 1 )
  IF ( two_fermi_energies .and. .not. lsda ) &
     CALL errore( 'iosys', 'tot_magnetization requires nspin=2', 1 )
  !
  IF ( TRIM(occupations) == 'fixed' .and. lsda  .and. lscf ) THEN
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
  tot_charge_        = tot_charge
  tot_magnetization_ = tot_magnetization
  !
  IF ( one_atom_occupations .and. trim(occupations) /= 'from_input' ) THEN
     CALL infomsg( 'iosys', 'one_atom_occupations requires occupations from input' )
     one_atom_occupations =.false.
  END IF
  IF ( one_atom_occupations .and. startingwfc /= 'atomic' ) THEN
     CALL infomsg( 'iosys', 'one_atom_occupations requires startingwfc atomic' )
     startingwfc = 'atomic'
  ENDIF
  one_atom_occupations_ = one_atom_occupations
  !
  IF ( tfixed_occ ) THEN
     IF ( nspin == 4 ) THEN
        ALLOCATE( f_inp_( nbnd, 1 ) )
     ELSE
        ALLOCATE( f_inp_( nbnd, nspin ) )
     ENDIF
     f_inp_ = f_inp
  ENDIF
  !
  nsp = ntyp
  !
  !
  ! STARTING AND RESTARTING
  !
  SELECT CASE( trim( restart_mode ) )
     !
  CASE( 'from_scratch' )
     !
     restart        = .false.
     ! ... non-scf calculation: read atomic positions and cell from file
     ! ... so that they are consistent.  FIXME: lforcet?
     IF ( trim( ion_positions ) == 'from_file' .OR. &
          (.NOT. lscf .AND. .NOT. lforcet) ) THEN
        startingconfig = 'file'
     ELSE
        startingconfig = 'input'
     END IF
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
  starting_charge_ = starting_charge
  starting_wfc     = startingwfc
  starting_pot     = startingpot
  !
  ! MEMORY AND DISK USAGE, VERBOSITY
  !
  smallmem = ( TRIM( memory ) == 'small' )
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
  CASE ( 'nowf' )
     !
     io_level = -1
     !
  CASE ( 'minimal' )
     !
     io_level = -2
     !
  CASE ( 'none' )
     !
     io_level = -3
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
  SELECT CASE( trim( verbosity ) )
  CASE( 'debug', 'high', 'medium' )
     iverbosity = 1
  CASE( 'low', 'default', 'minimal' )
     iverbosity = 0 
  CASE DEFAULT
     iverbosity = 0
  END SELECT
  iprint_ = iprint
  max_xml_steps_ = max_xml_steps
  !
  ! DIAGONALIZATION
  !
  SELECT CASE( trim( diagonalization ) )
  CASE ( 'david', 'davidson' )
     !
     isolve = 0
     david = diago_david_ndim
     !
  CASE ( 'cg' )
     !
     isolve = 1
     max_cg_iter = diago_cg_maxiter
     !
  CASE ( 'ppcg', 'PPCG' )
     !
     isolve = 2
     CALL errore( 'iosys', 'PPCG diagonalization not supported anymore (Dec. 2024)', 1 )
     !
  CASE ( 'paro', 'ParO' )
     !
     isolve = 3
     !
  CASE ( 'rmm', 'rmm-diis', 'rmm-davidson' )
     !
     isolve = 4
     rmm_ndim  = diago_rmm_ndim
     rmm_conv  = diago_rmm_conv
     gs_nblock = diago_gs_nblock
     rmm_with_davidson = .TRUE. 
     !
  CASE  ( 'rmm-paro')
     !
     isolve = 4
     rmm_ndim = diago_rmm_ndim 
     rmm_conv = diago_rmm_conv 
     gs_nblock = diago_gs_nblock 
     rmm_with_davidson = .FALSE.  
  CASE DEFAULT
     !
     CALL errore( 'iosys', 'diagonalization ' // &
                & trim( diagonalization ) // ' not implemented', 1 )
     !
  END SELECT
  !
  ethr = diago_thr_init
  tr2   = conv_thr
  nexxiter = exx_maxstep
  niter = electron_maxstep
  adapt_thr = adaptive_thr
  tr2_init  = conv_thr_init
  tr2_multi = conv_thr_multi
  diago_full_acc_ = diago_full_acc
  !
  ! EXTRAPOLATION
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
  ! TEMPERATURE AND THERMOSTATS
  !
  SELECT CASE( trim( ion_temperature ) )
  CASE( 'not_controlled', 'not-controlled', 'not controlled' )
     !
     control_temp = .false.
     IF ( ion_dynamics(1:8) == 'langevin' ) THEN
        temperature  = tempw
     ELSE
        temperature  = 0.0_dp
     END IF
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
   CASE( 'svr', 'Svr', 'SVR' )
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
  CASE ('nose')
     thermostat = trim(ion_temperature)
     temperature = tempw 
     tnosep = .true. 
     control_temp = .true.


  CASE DEFAULT
     !
     CALL errore( 'iosys', &
                & 'unknown ion_temperature ' // trim( ion_temperature ), 1 )
     !
  END SELECT
  !
  ! CELL TEMPERATURE
  ! 
  SELECT CASE (trim(cell_temperature))
    CASE ('nose')
      IF (.not. (control_temp .or. tnosep)) temperature = temph 
      tnoseh = .true.
    CASE ('not_controlled', 'not-controlled', 'not controlled') 
      tnoseh = .false. 
    CASE DEFAULT 
      CALL errore('iosys', 'unsupported cell_temperature',1)
  END SELECT 
  !
  ! SELF-CONSISTENCY
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
  IF ( mixing_beta < 0.0_DP ) THEN
     !
     IF ( lgcscf .AND. trism ) THEN
        ! GC-SCF with ESM-RISM
        mixing_beta = 0.1_DP
     ELSE IF ( lgcscf ) THEN
        ! GC-SCF with ESM-BC2 or ESM-BC3
        mixing_beta = 0.2_DP
     ELSE IF ( trism ) THEN
        ! 3D-RISM or ESM-RISM
        mixing_beta = 0.2_DP
     ELSE
        ! default
        mixing_beta = 0.7_DP
     END IF
     !
  END IF
  !
  starting_scf_threshold = tr2
  nmix                   = mixing_ndim
  mixing_beta_           = mixing_beta
  niter_with_fixed_ns    = mixing_fixed_ns
  scf_must_converge_     = scf_must_converge
  !
  IF ( ion_dynamics == ' bfgs' .and. epse <= 20.D0 * ( tr2 / upscale ) ) &
       CALL errore( 'iosys', 'required etot_conv_thr is too small:' // &
                     & ' conv_thr must be reduced', 1 )
  !
  ! ELECTRIC FIELDS AND BERRY PHASE
  !
  IF ( lberry .OR. lelfield .OR. lorbm ) THEN
     IF ( lgauss .OR. ltetra ) CALL errore( 'iosys', &
          'Berry Phase/electric fields only for insulators!', 1 )
     IF ( lmovecell ) CALL errore( 'iosys', &
          'Berry Phase/electric fields not implemented with variable cell', 1 )
  END IF
  !
  nppstr_     = nppstr
  gdir_       = gdir
  lberry_     = lberry
  lelfield_   = lelfield
  lorbm_      = lorbm
  efield_     = efield
  nberrycyc_  = nberrycyc
  efield_cart_ = efield_cart
  twochem_    = twochem
  SELECT CASE(efield_phase)
     CASE( 'none' )
        phase_control=0
     CASE ('write')
        phase_control=1
     CASE ('read')
        phase_control=2
     CASE DEFAULT
        CALL errore( 'iosys', 'Unknown efield_phase', 1 )
  END SELECT
  !
  ! K-POINTS
  !
  gamma_only = ( k_points == 'gamma' )
  IF ( lelfield .AND. gamma_only ) &
      CALL errore( 'iosys', 'electric fields not available for k=0 only', 1 )
  !
  ! DMFT
  !
  dmft_             = dmft
  dmft_prefix_      = dmft_prefix
  !
#if defined __HDF5
  IF ( dmft) THEN
     IF ( nspin > 1 ) CALL errore( 'iosys', &
          'DMFT update not implemented with nspin > 1', 1 )
  ENDIF
#else
  IF ( dmft) THEN
      CALL errore( 'iosys', 'DMFT update not implemented without HDF5 library', 1 )
  ENDIF
#endif
  !
  ! DFT+U
  !
  CALL dftu_iosys ( nsp )
  !
  ! REAL-SPACE TREATMENT
  !
  tqr_        = tqr
  real_space_ = real_space
  tq_smoothing_ = tq_smoothing
  tbeta_smoothing_ = tbeta_smoothing
  !
  ! SYMMETRY
  !
  no_t_rev_ = no_t_rev
  allfrac   = use_all_frac
  noinv_    = noinv
  nosym_    = nosym
  nosym_evc_= nosym_evc
  nofrac    = force_symmorphic
  !
  ! MOLECULAR DYNAMICS AND VARIABLE-CELL MD
  ! 
  remove_rigid_rot_ = remove_rigid_rot
  upscale_          = upscale
  refold_pos_       = refold_pos
  ecfixed_ = ecfixed
  qcutz_   = qcutz
  q2sigma_ = q2sigma
  cell_factor_      = cell_factor
  IF ( lmovecell ) THEN
     IF ( cell_factor_ <= 0.0_dp ) cell_factor_ = 2.0_dp
  ELSE
     cell_factor_ = 1.D0
  ENDIF
  wmass_ = wmass
  !
  ! ... unit conversion for pressure
  !
  press_ = press / ry_kbar
  !
  ! CONSTRAINTS
  !
  CALL init_dofree ( cell_dofree )
  lconstrain = ( nconstr_inp > 0 )
  IF ( lconstrain .AND. ( lbfgs .OR. lmovecell ) ) &
     CALL errore( 'iosys', 'constraints only with fixed-cell dynamics', 1 )
  !
  ! MISCELLANEOUS VARIABLES
  !
  la2F_      = la2F
  max_seconds_ = max_seconds
  !
  ! ... for WANNIER_AC
  !
  use_wannier_ = use_wannier
  use_energy_int_ = use_energy_int
  nwan_ = nwan
  print_wannier_coeff_ = print_wannier_coeff
  !
  ! ... BFGS specific
  !
  CALL init_bfgs( stdout, bfgs_ndim, trust_radius_max, trust_radius_min, &
        trust_radius_ini, w_1, w_2, (bfgs_ndim .gt. 1) .and. tgdiis_step )
  !
  IF (trim(occupations) /= 'from_input') one_atom_occupations_=.false.
  !
  !  VdW CORRECTIONS (SEMI-EMPIRICAL)
  !
  CALL set_vdw_corr ( vdw_corr, llondon, ldftd3, ts_vdw_, mbd_vdw_, lxdm)
  !
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
  IF ( mbd_vdw ) THEN
     CALL infomsg("iosys","mbd_vdw is obsolete, use ""vdw_corr='MBD'"" instead")
     vdw_corr='MBD'
     mbd_vdw_ = .TRUE.
  END IF
  IF ( llondon.AND.lxdm .OR. llondon.AND.ts_vdw_ .OR. lxdm.AND.ts_vdw_ .OR. &
           ldftd3.AND.llondon .OR. ldftd3.AND.lxdm .OR. ldftd3.AND.ts_vdw ) &
     CALL errore("iosys","must choose a unique vdW correction!", 1)
  IF ( lxdm .AND. noncolin ) &
     CALL errore("iosys","XDM not implemented for noncolinear case", 1)
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
  ! BOUNDARY CONDITIONS
  !
  do_makov_payne  = .false.
  do_comp_mt      = .false.
  do_comp_esm     = .false.
  do_cutoff_2D    = .false.
  !
  SELECT CASE( trim( assume_isolated ) )
      !
    CASE( 'makov-payne', 'm-p', 'mp' )
      !
      do_makov_payne = .true.
      !
    CASE( 'martyna-tuckerman', 'm-t', 'mt' )
      !
      do_comp_mt     = .true.
      !
    CASE( 'esm' )
      !
      do_comp_esm    = .true.
      !
    CASE( '2D' )
      !
      do_cutoff_2D   = .true.
      !
    CASE ( 'none' )
      !
      CONTINUE
      !
    CASE DEFAULT
      !
      CALL errore('iosys','unknown value assume_isolated="' // &
              & TRIM(assume_isolated) // '"',1)
      !
  END SELECT
  !
  IF ( do_comp_mt .AND. tstress_ ) THEN
     tstress_ = .false.
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
  ! QM/MM specific parameters
  !
  IF (.NOT. tqmmm) CALL qmmm_config( mode=-1 )
  !
END SUBROUTINE control_iosys
!
!
SUBROUTINE magnetization_iosys()
  USE kinds, ONLY: DP
  USE constants, ONLY: pi
  USE control_flags, ONLY: lscf
  USE input_parameters, ONLY : tot_magnetization, starting_magnetization,          &
                              lambda, angle1, angle2, constrained_magnetization,  &
                              B_field, fixed_magnetization, report, lspinorb,&
                              starting_spin_angle, noncolin, lforcet
 
  USE lsda_mod,      ONLY : nspin, &
                            starting_magnetization_ => starting_magnetization, &
                            lsda
  !
  USE noncollin_module, ONLY: i_cons, mcons, bfield, &
                              noncolin_  => noncolin, &
                              lambda_    => lambda, &
                              angle1_    => angle1, &
                              angle2_    => angle2, &
                              report_    => report, &
                              lspinorb_ => lspinorb,  &
                              lforcet_ => lforcet,    &
                              starting_spin_angle_ => starting_spin_angle
  !
  USE fixed_occ, ONLY: tfixed_occ
  USE klist, ONLY: tot_magnetization_ => tot_magnetization, &
                   two_fermi_energies 
  USE ions_base, ONLY: nsp, zv
  USE read_namelists_module, ONLY: sm_not_set
  IMPLICIT NONE
  !
  LOGICAL  :: sm_wasnt_set, domag
  REAL(DP) :: theta, phi
  INTEGER ::  nt  
  !
  noncolin_ = noncolin
  lspinorb_ = lspinorb
  lforcet_ = lforcet
  !
  ! ... starting_magnetization(nt) = sm_not_set means "not set"
  ! ... take notice and set to the default (zero)
  !
  sm_wasnt_set = ALL (starting_magnetization(1:nsp) == sm_not_set)
  DO nt = 1, nsp
     IF ( starting_magnetization(nt) == sm_not_set ) &
          starting_magnetization(nt) = 0.0_dp
  END DO
  IF (ANY(ABS(starting_magnetization(1:nsp)) .ge. 1._DP)) &
    starting_magnetization(1:nsp) = starting_magnetization(1:nsp) / zv(1:nsp)
  !
  !
  ! NONCOLLINEAR MAGNETISM, MAGNETIC CONSTRAINTS
  !
  IF (noncolin) THEN
     DO nt = 1, nsp
        !
        angle1(nt) = pi * angle1(nt) / 180.D0
        angle2(nt) = pi * angle2(nt) / 180.D0
        !
     ENDDO
  ELSE
     angle1=0.d0
     angle2=0.d0
  ENDIF

  SELECT CASE( trim( constrained_magnetization ) )
  CASE( 'none' )
     !
     ! ... if no constraints are imposed on the magnetization, 
     ! ... starting_magnetization must be set for at least one atomic type
     !
     IF ( lscf .AND. lsda .AND. ( .NOT. tfixed_occ ) .AND. &
          ( .not. two_fermi_energies )  .AND. sm_wasnt_set ) &
        CALL errore('iosys','some starting_magnetization MUST be set', 1 )
     !
     ! ... bring starting_magnetization between -1 and 1
     !
     DO nt = 1, nsp
        starting_magnetization(nt) = MIN( 1.0_dp,starting_magnetization(nt))
        starting_magnetization(nt) = MAX(-1.0_dp,starting_magnetization(nt))
     ENDDO
     !
     i_cons = 0
     !
  CASE( 'atomic' )
     !
     ! ... if "atomic" constraints are imposed on the magnetization, 
     ! ... starting_magnetization must be set for at least one atomic type
     !
     IF ( nspin == 1 ) &
        CALL errore( 'iosys','constrained atomic magnetizations ' // &
                   & 'require nspin=2 or 4 ', 1 )
     IF ( sm_wasnt_set ) &
        CALL errore( 'iosys','constrained atomic magnetizations ' // &
                   & 'require that some starting_magnetization is set', 1 )
     !
     i_cons = 1
     !
     IF (nspin == 4) THEN
        ! non-collinear case
        DO nt = 1, nsp
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
        DO nt = 1, nsp
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
     DO nt = 1, nsp
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
  starting_magnetization_ = starting_magnetization
  starting_spin_angle_ = starting_spin_angle
  angle1_   = angle1
  angle2_   = angle2
  lambda_   = lambda
  domag     = ANY ( ABS( starting_magnetization(1:nsp) ) > 1.D-6 )
  !
  IF ( (i_cons == 1 .OR. nspin == 2) .AND. (report /= 0) ) THEN
     report_ = -1
  ELSE IF ( (i_cons /= 0 .OR. report /= 0) .AND. ( domag .AND. noncolin) ) THEN
     report_ = report
  ELSE
     report_ = 0
  END IF

END SUBROUTINE magnetization_iosys
!----------------------------------------------------------------------------
SUBROUTINE structure_iosys ( )
  !----------------------------------------------------------------------------
  !
  !! Set cell, lattice, atomic positions, space group from input data
  !
  USE input_parameters, ONLY : ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                               trd_ht, rd_ht, cell_units
  USE cell_base,        ONLY : cell_base_init, at, alat, omega, bg
  USE starting_scf,     ONLY : startingconfig
  USE control_flags,    ONLY : restart, lscf
  USE ions_base,        ONLY : nat, nsp, tau
  !
  IMPLICIT NONE
  LOGICAL :: stop_on_error, is_tau_read
  !
  !
  ! SPACE GROUP (sets lattice from space group, must precede cell_base_init)
  !
  CALL spacegroup_iosys ( )
  !
  ! CRYSTAL LATTICE (sets at, bg and other relevant cell-related variables)
  !
  call cell_base_init ( ibrav, celldm, a, b, c, cosab, cosac, cosbc, &
                        trd_ht, rd_ht, cell_units )
  !
  ! ATOMIC SPECIES AND POSITIONS (mat require lattice vectors at)
  !
  CALL pos_iosys ( )
  !
  ! IF REQUIRED: read atomic positions from file
  !
  IF ( .NOT. restart .AND. startingconfig=='file' ) THEN
     !
     ! If this is not an nscf run don't stop on error
     !
     stop_on_error = .NOT. lscf
     !
     CALL read_conf_from_file( stop_on_error, nat, nsp, tau, alat, at, &
                               is_tau_read )
     !
     ! Update reciprocal lattice and volume (may be updated if coming from a vc run)
     !
     CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
     CALL volume (alat, at(:,1), at(:,2), at(:,3), omega)
     !
  END IF
  !
END SUBROUTINE structure_iosys
!
!-------------------------------------------------------------------------------
SUBROUTINE iosys_end ( )
  !-----------------------------------------------------------------------------
  !
  !! Final copying of input variables and related settings
  !! Contains all initializations that must be done AFTER
  !! control variables, cell and atomic structure have been read and set
  !
  USE input_parameters,      ONLY : trism, lfcp, lgcscf
  USE input_parameters,      ONLY : k_points, xk, wk, nk1, nk2, nk3,  &
                                    k1, k2, k3, nkstot
  USE input_parameters,      ONLY : ibrav, ecutwfc, ecutrho
  !
  USE control_flags,         ONLY : lconstrain, lxdm, llondon, ldftd3, &
                                    do_makov_payne
  USE cell_base,             ONLY : alat, at
  USE ions_base,             ONLY : nat, ityp, tau
  USE constraints_module,    ONLY : init_constraint
  USE start_k,               ONLY : init_start_k
  USE extfield,              ONLY : tefield, forcefield, gate, forcegate
  USE fcp_module,            ONLY : fcp_iosys
  USE gcscf_module,          ONLY : gcscf_iosys
  USE rism_module,           ONLY : rism_iosys
  USE xdm_module,            ONLY : init_xdm
  USE london_module,         ONLY : init_london
  !
  IMPLICIT NONE
  INTEGER :: ibrav_mp
  INTEGER, EXTERNAL :: at2ibrav
  !
  ! ... set up k-points (may require reciprocal lattice vectors bg)
  !
  CALL init_start_k ( nk1, nk2, nk3, k1, k2, k3, k_points, nkstot, xk, wk )
  !
  ! ... set cell mass for variable-cell MD if not set in input (requires amass)
  !
  CALL set_wmass ( )
  !
  ! ... check Bravais lattice for Makov-Payne if ibrav = 0
  !     (done here to avoid annoying crash at the end of calculation)
  !
  IF ( do_makov_payne ) THEN
     ibrav_mp = ibrav
     IF ( ibrav_mp .EQ. 0 ) ibrav_mp = at2ibrav(at(:, 1), at(:, 2), at(:, 3))
     IF ( ibrav_mp < 1 .OR. ibrav_mp > 3 ) CALL errore(' iosys', &
          'Makov-Payne correction defined only for cubic lattices', 1)
  END IF
  !
  ! ... FIXME: is this the proper place for the two lines below?
  !
  IF ( tefield ) ALLOCATE( forcefield( 3, nat ) )
  IF ( gate )    ALLOCATE( forcegate ( 3, nat ) ) 
  !
  ! ... initialize arrays for XDM,DFT-D2, DFT-D3 dispersion corrections
  !
  IF ( lxdm )   CALL init_xdm ( )
  IF ( llondon) CALL init_london ( )
  IF ( ldftd3)  CALL dftd3_iosys ()
  !
  ! ... set parameters for hybrid functionals (requires cutoffs)
  !
  CALL exx_iosys ( ecutwfc, ecutrho )
  !
  ! ... set constraints (requires atomic positions)
  !
  IF ( lconstrain ) CALL init_constraint( nat, tau, ityp, alat )
  !
  ! ... set variables for RISM
  !
  CALL rism_iosys(trism)
  !
  ! ... set variables for FCP (this must be after RISM, to check condition)
  !
  CALL fcp_iosys(lfcp)
  !
  ! ... set variables for GC-SCF (this must be after RISM and FCP, to check condition)
  !
  CALL gcscf_iosys(lgcscf)
  !
END SUBROUTINE iosys_end
!
!-------------------------------------------------------------------------------
SUBROUTINE set_cutoff( ecutwfc_in, ecutrho_in, ecutwfc_pp, ecutrho_pp, &
                       nr1,nr2,nr3, nr1s,nr2s,nr3s )
  !-----------------------------------------------------------------------------
  !! Copy to modules the cutoffs, either read from input or from PP files,
  !! and the FFT dimensions for smooth and dense grids, if set from input.  
  !! Values of \(\text{ecutwfc}\) and \(\text{ecutrho}\) are returned in 
  !! \(\text{ecutwfc_in}\), \(\text{ecutrho_in}\).
  !
  USE kinds, ONLY : dp
  USE gvecs, ONLY : dual
  USE gvect, ONLY : ecutrho
  USE gvecw, ONLY : ecutwfc
  USE constants, ONLY : eps8
  USE fft_base,  ONLY : dffts
  USE fft_base,  ONLY : dfftp
  !
  !
  IMPLICIT NONE
  !
  REAL(dp), INTENT(INOUT) :: ecutwfc_in, ecutrho_in
  REAL(dp), INTENT(IN)    :: ecutwfc_pp, ecutrho_pp
  INTEGER, INTENT(IN)     :: nr1, nr2, nr3, nr1s, nr2s, nr3s
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
  IF ( dual < 4.0_dp - eps8 ) CALL infomsg( 'set_cutoff', &
          'ecutrho < 4*ecutwfc, are you sure?' )
  !
  ! ... ensure that smooth and dense grid coincide when ecutrho=4*ecutwfc
  ! ... even when the dense grid is set from input and the smooth grid is not
  !
  dfftp%nr1    = nr1
  dfftp%nr2    = nr2
  dfftp%nr3    = nr3
  IF ( ( nr1 /= 0 .AND. nr2 /= 0 .AND. nr3 /= 0 ) .AND. &
       ( nr1s== 0 .AND. nr2s== 0 .AND. nr3s== 0 ) .AND. &
       ( ABS(dual-4.0_dp) < eps8) ) THEN
     dffts%nr1 = nr1
     dffts%nr2 = nr2
     dffts%nr3 = nr3
  ELSE
     dffts%nr1 = nr1s
     dffts%nr2 = nr2s
     dffts%nr3 = nr3s
  END IF
  !
END SUBROUTINE set_cutoff
!
!----------------------------------------------------------------------------
SUBROUTINE spacegroup_iosys ( )
  !----------------------------------------------------------------------------
  !
  USE input_parameters,   ONLY : lsg, space_group, uniqueb, origin_choice,   &
                                 ibrav, nat_ => nat, ntyp, rhombohedral,     &
                                 sp_pos, rd_for, tavel, sp_vel, rd_vel,      &
                                 rd_pos, rd_if_pos
  USE wyckoff,            ONLY : sup_spacegroup
  USE symm_base,          ONLY : spacegroup
  !
  IMPLICIT NONE
  INTEGER :: ibrav_sg
  !
  IF (lsg) THEN
     IF (space_group==0) &
        CALL errore('input','The option crystal_sg requires the space group &
                                                   &number',1 )
     CALL sup_spacegroup( rd_pos, sp_pos, rd_for, rd_if_pos, space_group, &
          nat_, uniqueb, rhombohedral, origin_choice, ibrav_sg )
     spacegroup = space_group
     IF (ibrav==-1 .OR. ibrav == ibrav_sg) THEN
        ibrav = ibrav_sg
     ELSEIF (ibrav /= ibrav_sg) THEN
        CALL errore ('input','Input ibrav not compatible with space group &
                                                   &number',1 )
     ENDIF
  ELSE
     IF (space_group /= 0) &
          CALL errore('input','space_group requires crystal_sg atomic &
                                                   & coordinates',1 )
  END IF
  !
END SUBROUTINE spacegroup_iosys
!
!----------------------------------------------------------------------------
SUBROUTINE pos_iosys ( )
  !----------------------------------------------------------------------------
  !
  USE input_parameters,   ONLY : taspc, tapos, rd_pos, atomic_positions,  &
                                 rd_if_pos, ibrav, nat_ => nat, ntyp,     &
                                 sp_pos, rd_for, tavel, sp_vel, rd_vel,   &
                                 atom_mass, atom_label, lsg, tempw,&
                                 fnosep, nhpcl, nhptyp, ndega, nhgrp, fnhscl
  USE kinds,              ONLY : DP
  USE dynamics_module,    ONLY : vel, get_ndof
  USE force_mod,          ONLY : force
  USE ions_base,          ONLY : nat, nsp, ityp, tau, atm, &
                                 extfor, if_pos, amass, fixatom, tau_format, &
                                 tions_base_init
  USE ions_nose,          ONLY:  ions_nose_init 
  USE control_flags,      ONLY : textfor, tv0rd, tnosep
  USE wyckoff,            ONLY : nattot, tautot, ityptot, extfortot, &
                                 if_postot, clean_spacegroup
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  !
  INTEGER :: is, ia
  !
  amass = 0
  nsp = ntyp
  !
  IF ( .not. taspc ) &
     CALL errore( 'pos_iosys', 'atomic species info missing', 1 )
  IF ( .not. tapos ) &
     CALL errore( 'pos_iosys', 'atomic position info missing', 1 )
  !
  DO is = 1, nsp
     !
     amass(is)  = atom_mass(is)
     atm(is)    = atom_label(is)
     !
     IF ( amass(is) <= 0.0_DP ) amass(is)= &
              atom_weight(atomic_number(trim(atm(is))))

     IF ( amass(is) <= 0.D0 ) CALL errore( 'pos_iosys', 'invalid  mass', is )
     !
  ENDDO
  !
  textfor = .false.
  IF( any( rd_for /= 0.0_DP ) ) textfor = .true.
  !
  ! Beware: when Wyckoff positions are read, the number of atoms read from
  !          input is the number of independent atoms, not of all atoms
  !
  IF (lsg) THEN
     nat = nattot
  ELSE
     nat = nat_
  END IF
  !
  ALLOCATE( ityp( nat ) )
  ALLOCATE( tau(    3, nat ) )
  ALLOCATE( force(  3, nat ) )
  ALLOCATE( if_pos( 3, nat ) )
  ALLOCATE( extfor( 3, nat ) )

  IF (lsg) THEN
     tau(:,:)=tautot(:,:)
     ityp(:) = ityptot(:)
     extfor(:,:) = extfortot(:,:)
     if_pos(:,:) = if_postot(:,:)
     CALL clean_spacegroup()
     !
  ELSE
     !
     DO ia = 1, nat
        tau(:,ia) = rd_pos(:,ia)
        ityp(ia)  = sp_pos(ia)
        extfor(:,ia) = rd_for(:,ia)
        if_pos(:,ia) = rd_if_pos(:,ia)
     ENDDO
     !
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
     DO ia = 1, nat_
        vel(:,ia) = rd_vel(:,ia)
     END DO
  END IF
  !
  ! ... The constrain on fixed coordinates is implemented using the array
  ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1
  ! ... otherwise. 
  !
  fixatom = COUNT( if_pos(1,:)==0 .AND. if_pos(2,:)==0 .AND. if_pos(3,:)==0 )
  !
  ! ... Convert atomic positions (tau) to internal units
  !
  tau_format = trim( atomic_positions )
  CALL convert_tau ( tau_format, nat, tau )
  IF (tnosep) THEN 
     tions_base_init = .TRUE. 
     IF (ndega == 0) ndega = get_ndof() 
     call ions_nose_init(tempw, fnosep, nhpcl, nhptyp, ndega, nhgrp, fnhscl)
  END IF 
  !
END SUBROUTINE pos_iosys
!
!----------------------------------------------------------------------------
SUBROUTINE dftu_iosys ( ntyp )
  !----------------------------------------------------------------------------
  !
  ! Hubbard parameters: input
  !
  USE input_parameters, ONLY : Hubbard_beta, Hubbard_occ,&
                               Hubbard_alpha, Hubbard_alpha_back,           &
                               hub_pot_fix, reserv, reserv_back, backall,   &
                               starting_ns_eigenvalue,     &
                               Hubbard_U, Hubbard_J, Hubbard_J0, Hubbard_V, Hubbard_U2, &
                               Hubbard_n, Hubbard_l, Hubbard_projectors, &
                               Hubbard_n2, Hubbard_l2, Hubbard_n3, Hubbard_l3, &
                               lda_plus_u, lda_plus_u_kind, Hubbard_Um, Hubbard_alpha_m, &
                               Hubbard_Um_nc, Hubbard_alpha_m_nc
  !
  USE constants,     ONLY : rytoev
  !
  ! Hubbard parameters: output
  !
  USE ldaU,          ONLY : Hubbard_U_  => hubbard_u, &
                            Hubbard_Um_  => hubbard_um, &
                            Hubbard_Um_nc_  => hubbard_um_nc, &
                            Hubbard_J0_ => hubbard_j0, &
                            Hubbard_J_ => hubbard_j, &
                            Hubbard_n_ => hubbard_n, &
                            Hubbard_l_ => hubbard_l, &
                            Hubbard_alpha_ => hubbard_alpha, &
                            Hubbard_alpha_m_ => hubbard_alpha_m, &
                            Hubbard_alpha_m_nc_ => hubbard_alpha_m_nc, &
                            Hubbard_beta_ => hubbard_beta, &
                            lda_plus_u_    => lda_plus_u, &
                            lda_plus_u_kind_    => lda_plus_u_kind, &
                            Hubbard_projectors_ => hubbard_projectors, &
                            iso_sys_    => iso_sys, &
                            starting_ns, &
                            Hubbard_U2_ => hubbard_u2, &
                            Hubbard_n2_ => hubbard_n2, &
                            Hubbard_l2_ => hubbard_l2, &
                            Hubbard_n3_ => hubbard_n3, &
                            Hubbard_l3_ => hubbard_l3, &
                            Hubbard_alpha_back_ => hubbard_alpha_back, &
                            Hubbard_V_ => hubbard_v , &
                            Hubbard_occ_ => hubbard_occ, &
                            hub_pot_fix_ => hub_pot_fix, &
                            reserv_ => reserv, &
                            backall_ => backall, &
                            reserv_back_ => reserv_back
  !
  IMPLICIT NONE
  !
  ! Needed input variables
  INTEGER, INTENT(IN) :: ntyp
  !
  !
  lda_plus_u_      = lda_plus_u
  lda_plus_u_kind_ = lda_plus_u_kind
  !
  !
  Hubbard_U_(1:ntyp)          = hubbard_u(1:ntyp) / rytoev
  Hubbard_Um_(:,:,:)     = hubbard_um(:,:,:) / rytoev
  Hubbard_Um_nc_(:,:)    = hubbard_um_nc(:,:) / rytoev
  Hubbard_J_(1:3,1:ntyp)      = hubbard_j(1:3,1:ntyp) / rytoev
  Hubbard_J0_(1:ntyp)         = hubbard_j0(1:ntyp) / rytoev
  Hubbard_V_(:,:,:)           = hubbard_V(:,:,:) / rytoev
  Hubbard_U2_(:)              = hubbard_U2(:) / rytoev
  Hubbard_n_(1:ntyp)          = hubbard_n(1:ntyp)
  Hubbard_l_(1:ntyp)          = hubbard_l(1:ntyp)
  Hubbard_n2_(1:ntyp)         = hubbard_n2(1:ntyp)
  Hubbard_l2_(1:ntyp)         = hubbard_l2(1:ntyp)
  Hubbard_n3_(1:ntyp)         = hubbard_n3(1:ntyp)
  Hubbard_l3_(1:ntyp)         = hubbard_l3(1:ntyp) 
  Hubbard_projectors_         = hubbard_projectors
  Hubbard_alpha_(:)      = hubbard_alpha(:) / rytoev
  Hubbard_alpha_m_(:,:,:) = hubbard_alpha_m(:,:,:) / rytoev
  Hubbard_alpha_m_nc_(:,:) = hubbard_alpha_m_nc(:,:) / rytoev
  Hubbard_beta_(1:ntyp)       = hubbard_beta(1:ntyp) / rytoev
  Hubbard_occ_(1:ntyp,1:3)    = hubbard_occ(1:ntyp,1:3)
  Hubbard_alpha_back_(1:ntyp) = hubbard_alpha_back(1:ntyp) / rytoev
  starting_ns                 = starting_ns_eigenvalue
  backall_(1:ntyp)            = backall(1:ntyp)
  hub_pot_fix_                = hub_pot_fix
  reserv_                     = reserv
  reserv_back_                = reserv_back
  !
END SUBROUTINE dftu_iosys
!
!----------------------------------------------------------------------------
SUBROUTINE dftd3_iosys ( )
  !----------------------------------------------------------------------------
  !
  ! Setting DFT-D3 functional dependent parameters
  ! Must be call AFTER dft is read from PP files and initialized
  !
  USE input_parameters, ONLY : dftd3_threebody, dftd3_version
  USE dftd3_api,        ONLY : dftd3_init, dftd3_set_functional
  USE dftd3_qe,         ONLY : dftd3_xc, dftd3, dftd3_in
  USE funct,            ONLY : get_dft_short
  IMPLICIT NONE
  CHARACTER(LEN=256):: dft
  !
  if (dftd3_version==2) dftd3_threebody=.false.
  dftd3_in%threebody = dftd3_threebody
  CALL dftd3_init(dftd3, dftd3_in)
  dft = get_dft_short( )
  dft = dftd3_xc ( dft )
  CALL dftd3_set_functional(dftd3, func=dft, version=dftd3_version,tz=.false.)
  !
END SUBROUTINE dftd3_iosys 
!
!----------------------------------------------------------------------------
SUBROUTINE exx_iosys ( ecutwfc, ecutrho )
  !----------------------------------------------------------------------------
  !
  ! ... set parameters for hybrid functionals
  ! ... must be call AFTER dft is read from PP files and initialized
  ! ... (in readpp) and AFTER cutoffs have been set (in set_cutoff)
  !
  USE kinds, ONLY: dp
  !
  USE input_parameters, ONLY: x_gamma_extrapolation, nqx1, nqx2, nqx3,     &
                              exxdiv_treatment, yukawa, ecutvcut,          &
                              gau_parameter, localization_thr, scdm, ace,  &
                              scdmden, scdmgrd, nscdm, n_proj,             & 
                              exx_fraction, screening_parameter, ecutfock 
  USE io_global,     ONLY : stdout
  USE klist,         ONLY : tot_charge
  USE ions_base,     ONLY : nat, ityp, zv
  USE xc_lib,        ONLY:  xclib_dft_is
  USE xc_lib,        ONLY : xclib_set_exx_fraction, set_screening_parameter
  USE exx_base,      ONLY : x_gamma_extrapolation_ => x_gamma_extrapolation, &
                            nq1, nq2, nq3, &
                            exxdiv_treatment_ => exxdiv_treatment, &
                            yukawa_           => yukawa, &
                            ecutvcut_         => ecutvcut
  USE exx,          ONLY :  ecutfock_         => ecutfock, &
                            use_ace, nbndproj, local_thr 
  USE loc_scdm,      ONLY : use_scdm, scdm_den, scdm_grd, n_scdm
  !
  IMPLICIT NONE
  REAL(dp), INTENT(IN):: ecutwfc, ecutrho
  REAL(dp) :: nelec_
  !
  !
  x_gamma_extrapolation_ = x_gamma_extrapolation
  !
  nq1 = nqx1
  nq2 = nqx2
  nq3 = nqx3
  !
  exxdiv_treatment_ = trim(exxdiv_treatment)
  yukawa_   = yukawa
  ecutvcut_ = ecutvcut
  use_ace   = ace
  nbndproj  = n_proj
  local_thr = localization_thr
  use_scdm  = scdm
  scdm_den = scdmden
  scdm_grd = scdmgrd
  n_scdm   = nscdm   
  IF ( local_thr > 0.0_dp .AND. .NOT. use_ace ) &
     CALL errore('input','localization without ACE not implemented',1)
  IF ( use_scdm ) CALL errore('input','use_scdm not yet implemented',1)
  !
  nelec_ = SUM( zv(ityp(1:nat)) ) - tot_charge
  IF (nelec_ <= 1) THEN
    use_ace = .false.
    IF (xclib_dft_is('hybrid')) WRITE(stdout, &
    '(5x,"ACE is turned off because number of occupied orbitals <= 1",/)' )
  END IF
  !
  IF(ecutfock <= 0.0_DP) THEN
     ! default case
     ecutfock_ = MIN ( ecutrho, 4.0_DP*ecutwfc)
  ELSE
     IF(ecutfock < ecutwfc .OR. ecutfock > ecutrho) CALL errore('iosys', &
          'ecutfock can not be < ecutwfc or > ecutrho!', 1) 
     ecutfock_ = ecutfock
  END IF
  !
  IF (exx_fraction >= 0.0_DP) CALL xclib_set_exx_fraction (exx_fraction)
  !
  IF (screening_parameter >= 0.0_DP) &
        & CALL set_screening_parameter(screening_parameter)
  !
END SUBROUTINE exx_iosys

SUBROUTINE set_wmass ( )
  !
  ! ... set default value of wmass if no value is set on input
  !
  USE constants,     ONLY : pi
  USE cell_base,     ONLY : wmass, omega
  USE ions_base,     ONLY : ityp, nat, amass
  USE cell_nose,     ONLY : cell_nose_init
  USE input_parameters, ONLY: fnoseh, temph
  USE cellmd,        ONLY : calc, lmovecell
  USE control_flags, ONLY : tnoseh
  !
  IMPLICIT NONE
  INTEGER :: ia
  !
  IF ( .NOT. lmovecell ) RETURN
  IF ( wmass == 0.D0 ) THEN
     !
#if defined __PGI
     DO ia = 1, nat
        wmass = wmass + amass( ityp(ia) )
     ENDDO
#else
     wmass = sum( amass(ityp(:)) )
#endif
     IF ( calc == 'nd' .or. calc == 'nm' ) THEN
        wmass = 0.75D0 * wmass / pi / pi / omega**( 2.D0 / 3.D0 )
     ELSEIF ( calc == 'cd' .or. calc == 'cm' ) THEN
        wmass = 0.75D0 * wmass / pi / pi
     ENDIF
     !
  ENDIF
  IF ( wmass <= 0.D0 ) CALL errore( 'set_wmass', &
            & 'vcsmd: a positive value for cell mass is required', 1 )
  IF (tnoseh) CALL cell_nose_init(temph, fnoseh)
  !
END SUBROUTINE set_wmass
