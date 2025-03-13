;; pwtk-mode.el
;;
;; Copyright (C) 2016 Quantum ESPRESSO group
;; This file is distributed under the terms of the
;; GNU General Public License. See the file `License'
;; in the root directory of the present distribution,
;; or http://www.gnu.org/copyleft/gpl.txt .
;;
;; Author: Anton Kokalj (tone.kokalj at ijs.si)
;;
;; Acknowledgments:
;;
;; The implementation of qe-modes package was made possible by several
;; useful and helpful resources that are gratefully acknowledged, in
;; particular: "Mode Tutorial" of Scott Andrew Borton
;; (https://www.emacswiki.org/emacs/ModeTutorial, for indentation
;; code), "Derived Mode" and "Sample Mode" pages
;; (https://www.emacswiki.org/emacs/DerivedMode,
;; https://www.emacswiki.org/emacs/SampleMode) as well as the very
;; useful resources of Xah Lee
;; (http://ergoemacs.org/emacs/elisp_syntax_coloring.html). Sebastijan
;; Peljhan is acknowledged for his work on `xsf-mode' that inspired
;; the idea of writing the qe-modes. Last but not the least,
;; Hongyi Zhao contributed the ido-completion-read snippet of
;; code for selecting the values for the card's flags.


;; This file is not part of GNU Emacs.


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.
;;
;; This lisp script is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
;;
;; Permission is granted to distribute copies of this lisp script
;; provided the copyright notice and this permission are preserved in
;; all copies.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program; if not, you can either send email to this
;; program's maintainer or write to: The Free Software Foundation,
;; Inc.; 675 Massachusetts Avenue; Cambridge, MA 02139, USA.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; send bug reports to the author (tone.kokalj at ijs.si)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; Commentary:

;; This is the `pwtk-mode', a major mode for composing the PWTK scripts.
;; For the installation and usage, see the user_guide.pdf in the Doc/
;; subdirectory of the original package (quick installation
;; instructions are also available in the README file of the original
;; package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar pwtk-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; QE's supercards (if any)
(defvar qe-open-supercards   (list "FIRST_IMAGE" "INTERMEDIATE_IMAGE" "LAST_IMAGE" ))
(defvar qe-closed-supercards (list "AUTOPILOT" "BEGIN" "BEGIN_ENGINE_INPUT" "BEGIN_PATH_INPUT" "BEGIN_POSITIONS" "END" "ENDRULES" "END_ENGINE_INPUT" "END_PATH_INPUT" "END_POSITIONS" ))
  
;; QE's namelists
(defvar pwtk-qe-namelists (list "BANDS" "CELL" "CONTROL" "DOS" "ELECTRONS" "ENERGY_CURRENT" "FCP" "HAM" "INPUT" "INPUTCOND" "INPUTHP" "INPUTMOPDOS" "INPUTP" "INPUTPH" "INPUTPP" "INPUT_BGW2PW" "INPUT_PW2BGW" "INTERPOLATION" "IONS" "LR_CONTROL" "LR_DAV" "LR_INPUT" "LR_POST" "OSCDFT" "OSCDFT_ET_NAMELIST" "OSCDFT_PP_NAMELIST" "PATH" "PLOT" "PPACF" "PRESS_AI" "PROJWFC" "RISM" "SCREEN" "SYSTEM" "TEST" "WANNIER" ))

;; QE's variables
(defvar qe-vars (list "A" "B" "C" "CI_scheme" "DeltaE" "Emax" "Emin" "Hubbard_U" "Hubbard_alpha" "Hubbard_beta" "Hubbard_occ" "P_ext" "P_fin" "P_in" "Surf_t" "abisur" "abivol" "ace" "adapt" "adaptive_thr" "add_i_current_b" "ahc_dir" "ahc_nbnd" "ahc_nbndskip" "alpha_mix" "amass" "amass_amu" "ampre" "amprp" "angle1" "angle2" "approximation" "array_convergence_func" "asr" "assume_isolated" "at" "atom" "atom_proj" "atom_proj_dir" "atom_proj_exclude" "atom_proj_ext" "atom_proj_ortho" "atomic_number" "author" "axis" "band_file" "bdl" "bdr" "bds" "beta" "beta_gamma_z_prefix" "bfgs_ndim" "block" "block_1" "block_2" "block_height" "broadening" "bz_sum" "calculation" "calculator" "calwf" "cau_fact" "cell_damping" "cell_dofree" "cell_dynamics" "cell_factor" "cell_parameters" "cell_temperature" "cell_velocities" "celldm" "charge_response" "check_ks" "check_periodicity" "check_spread" "closure" "code_num" "compute_hp" "config" "configts" "constrained_magnetization" "continue_not_converged" "conv_thr" "conv_thr_chi" "conv_thr_init" "conv_thr_multi" "conv_thr_multiplier" "convergence_type" "cosAB" "cosAC" "cosBC" "d0psi_rs" "debug_print" "decut" "degauss" "degauss_cond" "degauss_ldos" "dek" "deld" "deltaE" "delta_e" "delta_t" "denergy" "determine_num_pert_only" "determine_q_mesh_only" "dft" "dftd3_hess" "dftd3_threebody" "dftd3_version" "diag_basis" "diago_cg_maxiter" "diago_david_ndim" "diago_full_acc" "diago_gs_nblock" "diago_rmm_conv" "diago_rmm_ndim" "diago_thr_init" "diagonalization" "dipfield" "disk_io" "dist_thr" "dmft" "dmft_prefix" "do_bands" "do_charge_neutral" "do_long_range" "docc_thr" "dos" "drho_star" "ds" "dt" "dthr" "dvscf_star" "dx" "e1" "e2" "e3" "eamp" "ecfixed" "ecut2d" "ecutfock" "ecutmax" "ecutmin" "ecutrho" "ecutsolv" "ecutvcut" "ecutwfc" "edir" "eels" "efermi" "efield" "efield_cart" "efield_phase" "efx0" "efx1" "efy0" "efy1" "efz0" "efz1" "eigen_similarity" "eign_file" "ekin_conv_thr" "ekincw" "el_ph_nsig" "el_ph_nsigma" "el_ph_sigma" "electron_damping" "electron_dynamics" "electron_maxstep" "electron_phonon" "electron_temperature" "electron_velocities" "elop" "emass" "emass_cutoff" "emax" "emaxld" "emaxpos" "emin" "eminld" "end" "energy0" "ensemble_energies" "eopreg" "epol" "eps_inf" "epsil" "epsproj" "equiv_type" "esm_bc" "esm_efield" "esm_nfit" "esm_w" "eta" "eth_ns" "eth_rps" "ethr_big_step" "ethr_nscf" "ethr_small_step" "etot_conv_thr" "ewind" "extrapolation" "exx_dis_cutoff" "exx_fraction" "exx_maxstep" "exx_me_rcut_pair" "exx_me_rcut_self" "exx_neigh" "exx_poisson_eps" "exx_ps_rcut_pair" "exx_ps_rcut_self" "exx_use_cube_domain" "exxdiv_treatment" "fcp_conv_thr" "fcp_delta_t" "fcp_dynamics" "fcp_mass" "fcp_mu" "fcp_ndiis" "fcp_nraise" "fcp_scheme" "fcp_temperature" "fcp_tempw" "fcp_thr" "fcp_tolp" "fcp_velocity" "fd" "fil_loc" "filband" "fildos" "fildrho" "fildvscf" "fildyn" "file_beta" "file_charge" "file_chi" "file_core" "file_output" "file_pseudo" "file_pseudopw" "file_qvan" "file_recon" "file_screen" "file_wfcaegen" "file_wfcncgen" "file_wfcusgen" "fileig" "fileout" "filepp" "filhess" "filmol" "filout" "filp" "filpdos" "filplot" "filproj" "filxsf" "final_conv_thr" "final_dir" "final_prefix" "find_atpert" "finish" "fire_alpha_init" "fire_dtmax" "fire_f_dec" "fire_f_inc" "fire_falpha" "fire_nmin" "first_last_opt" "first_step" "firstk" "fixed_magnetization" "fldos" "fldyn" "fleig" "flfrc" "flfrq" "fltau" "flvec" "fnhscl" "fnosee" "fnoseh" "fnosep" "forc_conv_thr" "force_symmorphic" "fpol" "freeze_all_atoms" "frozen_core" "gate" "gcscf_beta" "gcscf_conv_thr" "gcscf_mu" "gdir" "get_ground_state_first" "grease" "greash" "greasp" "has_disentangle" "has_max_multiplier" "has_min_multiplier" "have_empty" "homo_only" "huang" "i_atmwfc_beg_full" "i_atmwfc_beg_part" "i_atmwfc_end_full" "i_atmwfc_end_part" "i_bnd_beg_full" "i_bnd_beg_part" "i_bnd_end_full" "i_bnd_end_part" "i_orb" "ibrav" "iesr" "if_dft_spectrum" "if_random_init" "iflag" "ikind" "increment" "initial_dir" "initial_prefix" "input_dft" "interpolation" "iofspin" "ion_damping" "ion_dynamics" "ion_nstepe" "ion_positions" "ion_radius" "ion_temperature" "ion_velocities" "ipol" "iprint" "irmax" "irmin" "irr_bz" "isave" "isic" "iswitch" "iteration_type" "itermax" "itermax0" "iverbosity" "k1" "k2" "k3" "k_max" "k_min" "kband" "kcw_at_ks" "kcw_iverbosity" "kpoint" "kresolveddos" "l1" "l2" "l3" "l_vcut" "la2F" "lambda" "lambda_cold" "last_e" "last_irr" "last_k" "last_q" "last_step" "lastk" "latt" "laue_both_hands" "laue_buffer_left" "laue_buffer_right" "laue_expand_left" "laue_expand_right" "laue_nfit" "laue_starting_left" "laue_starting_right" "laue_wall" "laue_wall_epsilon" "laue_wall_lj6" "laue_wall_rho" "laue_wall_sigma" "laue_wall_z" "lberry" "lbinary_data" "lcharge" "lda_plus_u" "ldiag" "ldisp" "ldvscf_interpolate" "ldynamics" "lebedev" "lelfield" "lfcp" "lfock" "lforces" "lforcet" "lgcscf" "lgipaw_reconstruction" "lkpoint_dir" "lloc" "llocal" "lmin" "lnoloc" "localization_thr" "london" "london_c6" "london_rcut" "london_rvdw" "london_s6" "loop_ek" "lorbm" "loto_2d" "loto_disable" "low_directory_check" "lp" "lpaw" "lpdb" "lperm" "lplasma" "lplot" "lplot_drho" "lpunch" "lqdir" "lr_verbosity" "lraman" "lread_cond" "lread_loc" "lrotation" "lrpa" "lsave_wfc" "lsd" "lsdts" "lshift_d0psi" "lshift_q" "lsigma" "lsign" "lsmall" "lspinorb" "lsym" "ltammd" "ltks" "lwrite_cond" "lwrite_loc" "lwrite_overlaps" "magnons" "max_conv_thr" "max_iter" "max_multiplier" "max_out_wfc" "max_seconds" "maxiter" "maxwfdt" "mdiis1d_size" "mdiis1d_step" "mdiis3d_size" "mdiis3d_step" "memory" "method" "miller_max" "min_conv_thr" "min_gamma_n" "min_multiplier" "minimum_image" "miniter" "mixing_beta" "mixing_fixed_ns" "mixing_mode" "mixing_ndim" "modenum" "mp1" "mp2" "mp3" "n_inner" "n_ipol" "n_lambda" "n_max" "n_oscdft" "n_proj_boxes" "n_repeat_every_step" "n_workers" "na_ifc" "nat" "nat_todo" "nberrycyc" "nbnd" "nbnd_cond" "nconf" "ndega" "ndos" "ndr" "ndw" "nelec_cond" "nenergy" "new_core_ps" "nextffield" "nfile" "nframes" "ngauss" "nhgrp" "nhpcl" "nhptyp" "nit" "niter" "niter_cg_restart" "niter_cold_restart" "niter_max" "niter_ph" "nk" "nk1" "nk2" "nk3" "nlcc" "nld" "nmix" "nmix_ph" "no_hxc" "no_overlap" "no_t_rev" "nogg" "noinv" "noncolin" "normalize_swfc" "noscf" "nosym" "nosym_evc" "np1" "np2" "np3" "nppstr" "nq" "nq1" "nq2" "nq3" "nqx1" "nqx2" "nqx3" "nr1" "nr1b" "nr1s" "nr2" "nr2b" "nr2s" "nr3" "nr3b" "nr3s" "nraise" "nsd" "nsolv" "nspin" "nstep" "nstep_path" "nsteps" "ntyp" "num_basis_max" "num_eign" "num_init" "num_neigh" "num_of_images" "num_wann_emp" "num_wann_occ" "nwf" "nx" "ny" "nz" "nz1" "occupations" "omeg" "on_site_only" "one_atom_occupations" "only_init" "opt_scheme" "optimization_method" "orbj_fin" "orbj_in" "origin_choice" "ortho_eps" "ortho_max" "ortho_para" "orthogonalization" "orthogonalize_swfc" "outdir" "output" "output_format" "p_metric" "p_nbnd_occ" "p_nbnd_virt" "passop" "path_thr" "pawproj" "perturb_only_atom" "plot_2d" "plot_num" "plot_type" "plotboxes" "pol_type" "poor_of_ram" "poor_of_ram2" "pot_extrapolation" "pre_state" "precondition" "prefix" "prefixl" "prefixr" "prefixs" "prefixt" "press" "press_conv_thr" "print_debug" "print_eigvect" "print_matrix" "print_occupation_eigenvectors" "print_occupation_matrix" "pseudo_dir" "pseudo_hermitian" "pseudotype" "pvar" "q" "q1" "q2" "q2d" "q2sigma" "q3" "q_in_band_form" "q_in_cryst_coord" "qcutz" "qplda" "qplot" "radius" "rcloc" "rcore" "rcutv" "re_init_wfc_1" "re_init_wfc_2" "re_init_wfc_3" "read_dns_bare" "read_lr" "read_sym" "read_unitary_matrix" "readtau" "real_or_complex" "real_space" "recover" "reduce_io" "reduce_unk" "reduce_unk_factor" "reference" "refold_pos" "rel" "rel_dist" "relaxz" "relpert" "remove_interaction_blocks" "remove_rigid_rot" "report" "residue_conv_thr" "restart" "restart_mode" "restart_step" "rho0" "rho_thr" "rhog_file" "rhog_flag" "rhog_nvmax" "rhog_nvmin" "rhombohedral" "rism1d_bond_width" "rism1d_conv_thr" "rism1d_dielectric" "rism1d_maxstep" "rism1d_molesize" "rism1d_nproc" "rism3d_conv_level" "rism3d_conv_thr" "rism3d_maxstep" "rism3d_planar_average" "rlderiv" "rm" "rmatch_augfun" "rmatch_augfun_nc" "rmax" "rpwe" "rytoev_fact" "sHu_formatted" "sIu_formatted" "sample_bias" "save_dvpsi" "save_file" "saverho" "scale_sphere" "scdm_entanglement" "scdm_mu" "scdm_proj" "scdm_sigma" "scf_must_converge" "sci_cb" "sci_vb" "scissor" "screening_parameter" "search_sym" "seedname" "sic_energy" "sic_gamma" "single_pole" "skip_dw" "skip_equivalence_q" "skip_type" "skip_upperfan" "smear1d" "smear3d" "smearing" "solute_epsilon" "solute_lj" "solute_sigma" "space_group" "spin_component" "spn_formatted" "spread_thr" "start" "start_e" "start_irr" "start_k" "start_q" "starting1d" "starting3d" "starting_charge" "starting_magnetization" "starting_ns_eigenvalue" "starting_spin_angle" "startingpot" "startingwfc" "step" "step_mul" "step_rem" "string_method" "subtract_cm_vel" "sum_pertq" "sw_len" "swapping_technique" "symm_type" "tabps" "tcg" "td" "tdosinboxes" "tefield" "temp_kelvin" "temp_req" "temph" "tempv" "tempw" "three_point_derivative" "thresh_init" "title" "tk_plot" "tm" "tolp" "tolw" "tot_charge" "tot_magnetization" "tprnfor" "tqr" "tr2" "tr2_ph" "trajdir" "tran_file" "tran_prefix" "tranp" "trans" "trism" "trust_radius_ini" "trust_radius_max" "trust_radius_min" "ts_vdw" "ts_vdw_econv_thr" "ts_vdw_isolated" "tstress" "twochem" "uHu_formatted" "uIu_formatted" "uniqueb" "units" "upscale" "use_ace" "use_all_frac" "use_freezing" "use_gauss_ldos" "use_masses" "use_paw_as_gipaw" "use_ws_distance" "vdW_analysis" "vdw" "vdw_corr" "vel_input_units" "verbosity" "vkb" "vkbg_file" "vkbg_flag" "vscg_file" "vscg_flag" "vxc0_file" "vxc0_flag" "vxc_diag_nmax" "vxc_diag_nmin" "vxc_file" "vxc_flag" "vxc_integral" "vxc_offdiag_nmax" "vxc_offdiag_nmin" "vxc_zero_rho_core" "vxcdiag" "vxcg_file" "vxcg_flag" "w_1" "w_2" "w_T_npol" "wan_mode" "warm_up_niter" "weight" "wf_collect" "wf_efield" "wf_friction" "wf_q" "wf_switch" "wfc_extrapolation" "wfcdir" "wfdt" "wffort" "wfng_dk1" "wfng_dk2" "wfng_dk3" "wfng_file" "wfng_flag" "wfng_kgrid" "wfng_nband" "wfng_nk1" "wfng_nk2" "wfng_nk3" "wfng_nvmax" "wfng_nvmin" "wfng_occupation" "wfsd" "what" "which_augfun" "wmass" "worker_id" "wpot_dir" "write_amn" "write_coulomb" "write_dmn" "write_frc" "write_hr" "write_lr" "write_mmn" "write_sHu" "write_sIu" "write_spn" "write_uHu" "write_uIu" "write_unk" "write_unkg" "writev" "wvfn_formatted" "x0" "x_gamma_extrapolation" "xdm" "xdm_a1" "xdm_a2" "xmin" "xmlfile_full" "xmlfile_part" "zasr" "zed" "zeu" "zgate" "zue" "zval" ))

;; QE's cards & keywords
(defvar qe-cards (list "ADDITIONAL_K_POINTS" "ATOMIC_FORCES" "ATOMIC_POSITIONS" "ATOMIC_SPECIES" "ATOMIC_VELOCITIES" "CELL_PARAMETERS" "CLIMBING_IMAGES" "CONSTRAINTS" "GAMMA_VAL" "HUBBARD" "J0" "K_POINTS" "OCCUPATIONS" "PLOT_WANNIER" "REF_CELL_PARAMETERS" "ROUGHNESS" "SOLVENTS" "TARGET_OCCUPATION_NUMBERS" "TOTAL_CHARGE" "U" "USER_STARS" "V" "ON_STEP" ))

;; QE's flags
(defvar qe-flags (list "1/cell" "a.u" "alat" "angstrom" "atomic" "automatic" "bohr" "crystal" "crystal_b" "crystal_c" "crystal_sg" "g/cm^3" "gamma" "mol/L" "norm-atomic" "ortho-atomic" "pseudo" "tpiba" "tpiba_b" "tpiba_c" "wf" ))

;; PWTK's cmds
(defvar pwtk-cmds (list "ARTN_PARAMETERS" "ATOMIC_POSITIONS_fromCRD" "ATOMIC_POSITIONS_fromPWO" "ATOMIC_POSITIONS_fromXSF" "AUTOPILOT" "BOUNDARY" "CELL_PARAMETERS_and_ATOMIC_POSITIONS_fromPWO" "CELL_PARAMETERS_and_ATOMIC_POSITIONS_fromXSF" "CELL_PARAMETERS_fromPWO" "CELL_PARAMETERS_fromXSF" "CPPP" "D3HESS" "DIELECTRIC_REGIONS" "DIFDEN" "DYNMAT" "ELECTROSTATIC" "ENVIRON" "EOS" "EXTERNAL_CHARGES" "FASTNEB" "INPUT_PW2GW" "INPUT_PW2WAN" "KCW_CONTROL" "LD1_INPUT" "LD1_INPUTP" "LD1_TEST" "LL" "LOBSTER" "LSF" "MATDYN" "MOPDOS" "OSCDFT_ET" "OSCDFT_PP" "PBS" "POSITIONS" "POSITIONS_fromAXSF" "POSITIONS_fromCRD" "POSITIONS_fromPATH" "POSTAHC" "Q2R" "RISM_INPUTPP" "RISM_PLOT" "SBCO" "SH" "SLURM" "SPINNEB" "abort" "addImage" "addImageCoor_fromCRD" "addImageCoor_fromPWO" "addImageCoor_fromXSF" "afterDone" "all_curri_clear" "all_curri_fprint" "all_curri_get" "artn" "artn_clear" "artn_fprint" "artn_get" "auto_neb" "backup_io" "band_inti_clear" "band_inti_fprint" "band_inti_get" "bands" "bands_plot" "bands_run" "bg" "bg_wait" "bgw2pwi_clear" "bgw2pwi_fprint" "bgw2pwi_get" "bi_clear" "bi_fprint" "bi_get" "bin_dir" "bin_query" "calcmol" "card" "cardAppend" "cardClear" "cardContent" "cardFlags" "cardGet" "cardGetContent" "cardGetFlags" "cardGetPWTK" "cardPrepend" "cardPrint" "catch_error" "container" "cp" "cpi_clear" "cpi_fprint" "cpi_get" "cpi_visualize" "cpppi_clear" "cpppi_fprint" "cpppi_get" "d3hessi_clear" "d3hessi_fprint" "d3hessi_get" "davi_clear" "davi_fprint" "davi_get" "delayed_eval" "delayed_eval_finish" "deleteAtoms" "deleteAtoms1st" "deleteAtomsLast" "deleteAtomsRange" "deleteImageCharges" "deleteImages" "di_clear" "di_fprint" "di_get" "difden_run" "difden_segmentSpec" "displaceAtoms" "display" "dmi_clear" "dmi_fprint" "dmi_get" "dmo_vibSpectrum" "dos" "dos_plot" "dos_run" "dryrunCP" "dryrunNEB" "dryrunPW" "dummy" "dynmat" "eelsi_clear" "eelsi_fprint" "eelsi_get" "environ" "environ_clear" "environ_fprint" "environ_get" "eos_run" "eval_in_dir" "execute" "fast_neb" "finish" "fixAtoms" "fixAtoms1st" "fixAtomsLast" "fixAtomsRange" "flo_run" "getImage" "getImageAtmPos" "getImageCharge" "getImageUnit" "getPositions" "heatmap" "hp" "hpi_clear" "hpi_fprint" "hpi_get" "ifset" "ihandle" "imageTile" "import" "infoMsg" "input_clear" "input_handle" "input_peek" "input_pop" "input_push" "input_pushpop" "input_stack_level" "insertAtoms" "kcwi_clear" "kcwi_fprint" "kcwi_get" "lani_clear" "lani_fprint" "lani_get" "ld1i_clear" "ld1i_fprint" "ld1i_get" "ldos_fullplot" "ldos_multiplot" "ldos_plot" "ll_clear" "ll_fprint" "ll_get" "ll_head" "ll_profile" "ll_profileDefault" "ll_pwtk_profile" "ll_pwtk_profileDefault" "ll_tail" "load_artn.in" "load_environ.in" "load_from" "load_fromALL_CURRI" "load_fromBAND_INTI" "load_fromBGW2PWI" "load_fromBI" "load_fromCPI" "load_fromCPPPI" "load_fromD3HESSI" "load_fromDAVI" "load_fromDI" "load_fromDMI" "load_fromEELSI" "load_fromHPI" "load_fromKCWI" "load_fromLANI" "load_fromLD1I" "load_fromMAGI" "load_fromMDI" "load_fromMPDI" "load_fromNEBI" "load_fromOSCDFT_ETI" "load_fromOSCDFT_PPI" "load_fromPHI" "load_fromPOSTAHCI" "load_fromPPACFI" "load_fromPPI" "load_fromPPRISMI" "load_fromPRI" "load_fromPW2BGWI" "load_fromPW2GWI" "load_fromPW2WAN90I" "load_fromPWCONDI" "load_fromPWI" "load_fromQ2RI" "load_fromSPECI" "load_oscdft.in" "lobster_datafy" "lobster_getinfo" "lobster_multiplot" "lobster_plot" "lobster_run" "lobster_sumdos" "lsf_clear" "lsf_fprint" "lsf_get" "lsf_head" "lsf_profile" "lsf_profileDefault" "lsf_pwtk_profile" "lsf_pwtk_profileDefault" "lsf_tail" "magi_clear" "magi_fprint" "magi_get" "manual_neb" "matdyn" "mdi_clear" "mdi_fprint" "mdi_get" "mdi_getAffix" "measure_time" "molecularpdos" "mopdos_plot" "mopdos_run" "mpdi_clear" "mpdi_fprint" "mpdi_get" "multiplot" "namelist" "namelist.affix" "namelist.affixGet" "namelist.affixGetAffix" "namelist.affixPrint" "namelist.affixPrintAffix" "namelistClear" "namelistGet" "namelistGetPWTK" "namelistGetVarValue" "namelistPrint" "neb" "neb_plot" "neb_refine" "neb_refine_auto" "neb_refine_merge" "nebi_clear" "nebi_fprint" "nebi_get" "nebi_visualize" "nice" "orbitalGroup" "oscdft" "oscdft_clear" "oscdft_eti_clear" "oscdft_eti_fprint" "oscdft_eti_get" "oscdft_fprint" "oscdft_get" "oscdft_ppi_clear" "oscdft_ppi_fprint" "oscdft_ppi_get" "outdir" "outdir_backup" "outdir_clean" "outdir_create" "outdir_postfix" "outdir_prefix" "outdir_prefix_append" "outdir_query" "pbs_clear" "pbs_fprint" "pbs_get" "pbs_head" "pbs_profile" "pbs_profileDefault" "pbs_pwtk_profile" "pbs_pwtk_profileDefault" "pbs_tail" "pdos_atm_files" "pdos_fullplot" "pdos_multiplot" "pdos_plot" "pdos_run" "ph" "phi_clear" "phi_fprint" "phi_get" "phi_getAffix" "phi_getTitle" "plot" "postahci_clear" "postahci_fprint" "postahci_get" "postfix" "postrun" "pp" "pp_bader" "ppacfi_clear" "ppacfi_fprint" "ppacfi_get" "ppi_clear" "ppi_fprint" "ppi_get" "pprismi_clear" "pprismi_fprint" "pprismi_get" "prefix" "prerun" "pri_clear" "pri_fprint" "pri_get" "print" "printAll" "printTitle" "prog" "projwfc" "propagate" "propagate_clear" "pseudo_auto_download" "pseudo_dir" "pseudo_download" "pseudo_etot_test" "pseudo_etot_test_many" "pw" "pw2bgwi_clear" "pw2bgwi_fprint" "pw2bgwi_get" "pw2gwi_clear" "pw2gwi_fprint" "pw2gwi_get" "pw2wan90i_clear" "pw2wan90i_fprint" "pw2wan90i_get" "pw_pp_bader" "pwcondi_clear" "pwcondi_fprint" "pwcondi_get" "pwi_clear" "pwi_fprint" "pwi_get" "pwi_visualize" "pwo_alat" "pwo_conv_thr" "pwo_dipole" "pwo_ecutrho" "pwo_ecutwfc" "pwo_efermi" "pwo_etot" "pwo_etot_conv_thr" "pwo_forc_conv_thr" "pwo_ibrav" "pwo_mixing_beta" "pwo_nat" "pwo_nbnd" "pwo_nelec" "pwo_ntyp" "pwo_press" "pwo_totene" "pwo_totfor" "pwo_totmag" "pwo_volume" "pwo_xc" "q2r" "q2ri_clear" "q2ri_fprint" "q2ri_get" "q2ri_getAffix" "rcp" "re_spin_neb" "readFile" "relax" "relax_fromPWO" "relax_fromXSF" "remedy" "remote" "remote_get" "remote_id" "remote_ids" "remote_pidwait" "remote_wait" "replaceCoor" "replaceImage" "report" "rerunCP" "rerunDAVID" "rerunDAVIDSON" "rerunEELS" "rerunLANCZOS" "rerunNEB" "rerunPH" "rerunPW" "rerunPWCOND" "rerunXX" "restart" "restart_fromIO" "rexec" "rmkdir" "rsync" "rtclsh" "run" "runALL_CURRENTS" "runBANDS" "runBAND_INTERPOLATION" "runBGW2PW" "runCP" "runCPPP" "runD3HESS" "runDAVID" "runDAVIDSON" "runDOS" "runDYNMAT" "runEELS" "runHP" "runKCW" "runLANCZOS" "runLD1" "runMATDYN" "runMOLECULARPDOS" "runMOPDOS" "runNEB" "runOSCDFT_ET" "runOSCDFT_PP" "runPH" "runPOSTAHC" "runPP" "runPPACF" "runPPRISM" "runPROJWFC" "runPW" "runPW2BGW" "runPW2GW" "runPW2WANNIER90" "runPWCOND" "runPW_remedy" "runQ2R" "runSPECTRUM" "runXX" "run_fromPWO" "run_fromXSF" "run_id" "run_ids" "run_wait" "save_artn.in" "save_environ.in" "save_oscdft.in" "sbco_postNeb" "sbco_run" "scanpar" "script" "scriptAppend" "scriptClear" "scriptGet" "scriptPrepend" "scriptPrint" "seq" "serial_postfix" "serial_prefix" "setImageCharge" "setImageCharges" "sh_clear" "sh_fprint" "sh_get" "sh_head" "sh_profile" "sh_profileDefault" "sh_pwtk_profile" "sh_pwtk_profileDefault" "sh_tail" "slurm_clear" "slurm_fprint" "slurm_get" "slurm_head" "slurm_profile" "slurm_profileDefault" "slurm_pwtk_profile" "slurm_pwtk_profileDefault" "slurm_tail" "smi2coor" "smi2png" "smi2svg" "smi2xsf" "smilesToCoor" "speci_clear" "speci_fprint" "speci_get" "splot" "stopOnError" "substituteAtmType" "substituteAtom" "substituteAtoms" "sumldos" "sumldosFiles" "sumpdos" "sumpdosFiles" "thread" "thread_id" "thread_ids" "thread_wait" "time2ms" "time2s" "time2unit" "tpool" "treatrun" "try_exec" "turbo_davidson" "turbo_eels" "turbo_lanczos" "turbo_spectrum" "varvalue" "vc-relax" "vc-relax_fromPWO" "vc-relax_fromXSF" "warning" "wfcdir" "wfcdir_backup" "wfcdir_clean" "wfcdir_create" "wfcdir_postfix" "wfcdir_prefix" "wfcdir_prefix_append" "wfcdir_query" "write" "writeFile" ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar qe-open-supercards-regexp   (regexp-opt qe-open-supercards   'symbols)) ; may not exists
(defvar qe-closed-supercards-regexp (regexp-opt qe-closed-supercards 'symbols)) ; may not exists

(defvar qe-cards-regexp (regexp-opt
			    (append qe-cards qe-open-supercards) 'symbols))
(defvar qe-flags-regexp (regexp-opt qe-flags 'symbols))

(defvar pwtk-qe-namelist-face (cons (regexp-opt (append pwtk-qe-namelists) 'symbols) font-lock-function-name-face))
(defvar pwtk-cmds-face (cons (regexp-opt (append pwtk-cmds) 'symbols) font-lock-function-name-face))
(defvar qe-variable-face (cons (regexp-opt qe-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun qe-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward qe-closed-supercards-regexp limit 'no-error)))

(defun qe-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward qe-cards-regexp limit 'no-error)))

(defun qe-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward qe-flags-regexp limit 'no-error)))


(font-lock-add-keywords 'pwtk-mode (list
				     pwtk-qe-namelist-face 
				     pwtk-cmds-face 
				     qe-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'pwtk-mode '(
				      (qe-closed-supercards-matcher 1 font-lock-preprocessor-face t)
                                      (qe-cards-matcher 1 font-lock-function-name-face t)
				      (qe-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar qe-keywords '(pwtk-qe-namelist-face qe-variable-face))
(defvar qe-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar pwtk-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `pwtk-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the pwtk-mode as derived-mode

(define-derived-mode pwtk-mode tcl-mode 
  "PWTK-mode"  
  "Major mode for editing PWTK scripts"    
  )

;; free memory

(setq pwtk-qe-namelists nil)
(setq pwtk-cmds nil)
(setq qe-vars nil)
(setq qe-cards nil)
(setq qe-flags nil)
(setq qe-open-supercards   nil)
(setq qe-closed-supercards nil)


(require 'qe-funcs)
(provide 'pwtk-mode)

