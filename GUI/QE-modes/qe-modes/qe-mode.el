;; qe-mode.el
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
;; (http://ergoemacs.org/emacs/elisp_syntax_coloring.html).  Last but
;; not the least Sebastijan Peljhan is acknowledged for his work on
;; `xsf-mode' that inspired the idea of writing the qe-modes.


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

;; This is the `qe-mode', a major mode for composing the Quantum ESPRESSO
;; QE-generic input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar qe-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; qe's supercards (if any)
(defvar qe-open-supercards   (list "FIRST_IMAGE" "INTERMEDIATE_IMAGE" "LAST_IMAGE" ))
(defvar qe-closed-supercards (list "AUTOPILOT" "BEGIN" "BEGIN_ENGINE_INPUT" "BEGIN_PATH_INPUT" "BEGIN_POSITIONS" "END" "ENDRULES" "END_ENGINE_INPUT" "END_PATH_INPUT" "END_POSITIONS" ))
  
;; qe's namelists
(defvar qe-namelists (list "&BANDS" "&CELL" "&CONTROL" "&DOS" "&ELECTRONS" "&FCP" "&INPUT" "&INPUTCOND" "&INPUTHP" "&INPUTMOPDOS" "&INPUTP" "&INPUTPH" "&INPUTPP" "&INPUT_BGW2PW" "&INPUT_PW2BGW" "&IONS" "&LR_CONTROL" "&LR_DAV" "&LR_INPUT" "&LR_POST" "&PATH" "&PLOT" "&PPACF" "&PRESS_AI" "&PROJWFC" "&SYSTEM" "&TEST" "&WANNIER" ))
(defvar qe-end-namelist (list "&END" "/"))

;; qe's variables
(defvar qe-vars (list "A" "B" "C" "CI_scheme" "DeltaE" "Emax" "Emin" "Hubbard_J" "Hubbard_J0" "Hubbard_U" "Hubbard_V" "Hubbard_alpha" "Hubbard_beta" "Hubbard_parameters" "P_ext" "P_fin" "P_in" "Surf_t" "U_projection_type" "abisur" "abivol" "ace" "adapt" "adaptive_thr" "ahc_dir" "ahc_nbnd" "ahc_nbndskip" "alpha_mix" "amass" "amass_amu" "ampre" "amprp" "angle1" "angle2" "approximation" "asr" "assume_isolated" "at" "atom" "atomic_number" "author" "axis" "band_file" "bdl" "bdr" "bds" "beta" "beta_gamma_z_prefix" "bfgs_ndim" "block" "block_1" "block_2" "block_height" "broadening" "bz_sum" "calculation" "calculator" "calwf" "cau_fact" "cell_damping" "cell_dofree" "cell_dynamics" "cell_factor" "cell_parameters" "cell_temperature" "cell_velocities" "celldm" "charge_response" "code_num" "compute_hp" "config" "configts" "constrained_magnetization" "conv_thr" "conv_thr_chi" "conv_thr_init" "conv_thr_multi" "cosAB" "cosAC" "cosBC" "d0psi_rs" "decut" "degauss" "degauss_ldos" "dek" "deld" "deltaE" "delta_e" "delta_t" "denergy" "determine_num_pert_only" "dft" "dftd3_threebody" "dftd3_version" "diago_cg_maxiter" "diago_david_ndim" "diago_full_acc" "diago_thr_init" "diagonalization" "dipfield" "direction" "disk_io" "do_charge_neutral" "do_long_range" "docc_thr" "dos" "drho_star" "ds" "dt" "dthr" "dvscf_star" "dx" "e1" "e2" "e3" "eamp" "ecfixed" "ecut2d" "ecutfock" "ecutmax" "ecutmin" "ecutrho" "ecutvcut" "ecutwfc" "edir" "eels" "efermi" "efield" "efield_cart" "efield_phase" "efx0" "efx1" "efy0" "efy1" "efz0" "efz1" "eigen_similarity" "eign_file" "ekin_conv_thr" "ekincw" "el_ph_nsig" "el_ph_nsigma" "el_ph_sigma" "electron_damping" "electron_dynamics" "electron_maxstep" "electron_phonon" "electron_temperature" "electron_velocities" "elop" "emass" "emass_cutoff" "emax" "emaxld" "emaxpos" "emin" "eminld" "end" "energy0" "ensemble_energies" "eopreg" "epol" "epsil" "epsproj" "equiv_type" "esm_bc" "esm_efield" "esm_nfit" "esm_w" "eta" "eth_ns" "eth_rps" "ethr_nscf" "etot_conv_thr" "ewind" "extrapolation" "exx_dis_cutoff" "exx_fraction" "exx_me_rcut_pair" "exx_me_rcut_self" "exx_neigh" "exx_poisson_eps" "exx_ps_rcut_pair" "exx_ps_rcut_self" "exx_use_cube_domain" "exxdiv_treatment" "fcp_conv_thr" "fcp_delta_t" "fcp_dynamics" "fcp_mass" "fcp_mu" "fcp_ndiis" "fcp_nraise" "fcp_scheme" "fcp_temperature" "fcp_tempw" "fcp_thr" "fcp_tolp" "fcp_velocity" "fd" "fil_loc" "filband" "fildos" "fildrho" "fildvscf" "fildyn" "file_beta" "file_charge" "file_chi" "file_core" "file_pseudo" "file_pseudopw" "file_qvan" "file_recon" "file_screen" "file_wfcaegen" "file_wfcncgen" "file_wfcusgen" "fileig" "fileout" "filepp" "filmol" "filout" "filp" "filpdos" "filplot" "filproj" "filxsf" "find_atpert" "finish" "first_last_opt" "firstk" "fixed_magnetization" "fldos" "fldyn" "fleig" "flfrc" "flfrq" "fltau" "flvec" "fnhscl" "fnosee" "fnoseh" "fnosep" "forc_conv_thr" "force_symmorphic" "fpol" "freeze_all_atoms" "frozen_core" "gate" "gcscf_beta" "gcscf_conv_thr" "gcscf_mu" "gdir" "grease" "greash" "greasp" "i_atmwfc_beg_full" "i_atmwfc_beg_part" "i_atmwfc_end_full" "i_atmwfc_end_part" "i_bnd_beg_full" "i_bnd_beg_part" "i_bnd_end_full" "i_bnd_end_part" "ibrav" "iesr" "if_dft_spectrum" "if_random_init" "iflag" "ikind" "increment" "input_dft" "interpolation" "iofspin" "ion_damping" "ion_dynamics" "ion_nstepe" "ion_positions" "ion_radius" "ion_temperature" "ion_velocities" "ipol" "iprint" "irmax" "irmin" "isave" "isic" "iswitch" "itermax" "itermax0" "iverbosity" "k1" "k2" "k3" "k_max" "k_min" "kband" "kpoint" "kresolveddos" "l1" "l2" "l3" "la2F" "lambda" "lambda_cold" "last_e" "last_irr" "last_k" "last_q" "lastk" "latt" "lberry" "lbinary_data" "lcharge" "lda_plus_u" "lda_plus_u_kind" "ldiag" "ldisp" "ldvscf_interpolate" "ldynamics" "lelfield" "lfcp" "lfock" "lforces" "lforcet" "lgcscf" "lgipaw_reconstruction" "lkpoint_dir" "lloc" "llocal" "lmin" "lnoloc" "localization_thr" "london" "london_c6" "london_rcut" "london_rvdw" "london_s6" "loop_ek" "lorbm" "loto_2d" "loto_disable" "low_directory_check" "lp" "lpaw" "lpdb" "lperm" "lplasma" "lplot" "lplot_drho" "lqdir" "lr_verbosity" "lraman" "lread_cond" "lread_loc" "lrotation" "lrpa" "lsave_wfc" "lsd" "lsdts" "lshift_d0psi" "lshift_q" "lsigma" "lsign" "lsmall" "lspinorb" "lsym" "ltammd" "ltks" "lwrite_cond" "lwrite_loc" "lwrite_overlaps" "max_iter" "max_out_wfc" "max_seconds" "maxiter" "maxwfdt" "memory" "minimum_image" "mixing_beta" "mixing_fixed_ns" "mixing_mode" "mixing_ndim" "modenum" "n_inner" "n_ipol" "n_lambda" "n_proj_boxes" "na_ifc" "nat" "nat_todo" "nberrycyc" "nbnd" "nconf" "ndega" "ndos" "ndr" "ndw" "nenergy" "new_core_ps" "newoutdir" "nfile" "nframes" "ngauss" "nhgrp" "nhpcl" "nhptyp" "ninter_cold_restart" "nit" "niter_cg_restart" "niter_max" "niter_ph" "nk" "nk1" "nk2" "nk3" "nlcc" "nld" "nmix" "nmix_ph" "no_hxc" "no_overlap" "no_t_rev" "nogg" "noinv" "noncolin" "noscf" "nosym" "nosym_evc" "np1" "np2" "np3" "nppstr" "nq" "nq1" "nq2" "nq3" "nqx1" "nqx2" "nqx3" "nr1" "nr1b" "nr1s" "nr2" "nr2b" "nr2s" "nr3" "nr3b" "nr3s" "nraise" "nsd" "nspin" "nstep" "nstep_path" "nsteps" "ntyp" "num_basis_max" "num_eign" "num_init" "num_neigh" "num_of_images" "nwf" "nx" "ny" "nz" "nz1" "occupations" "omeg" "one_atom_occupations" "only_init" "opt_scheme" "orbj_fin" "orbj_in" "origin_choice" "ortho_eps" "ortho_max" "ortho_para" "orthogonalization" "outdir" "output" "output_format" "p_nbnd_occ" "p_nbnd_virt" "passop" "path_thr" "pawproj" "perturb_only_atom" "plot_2d" "plot_num" "plot_type" "plotboxes" "poor_of_ram" "poor_of_ram2" "pot_extrapolation" "precondition" "prefix" "prefixl" "prefixr" "prefixs" "prefixt" "press" "press_conv_thr" "pseudo_dir" "pseudo_hermitian" "pseudotype" "pvar" "q" "q1" "q2" "q2d" "q2sigma" "q3" "q_in_band_form" "q_in_cryst_coord" "qcutz" "qplda" "qplot" "radius" "rcloc" "rcore" "rcutv" "read_dns_bare" "readtau" "real_or_complex" "real_space" "recover" "reduce_io" "reference" "refold_pos" "rel" "rel_dist" "relaxz" "relpert" "remove_rigid_rot" "report" "residue_conv_thr" "restart" "restart_mode" "restart_step" "rho0" "rho_thr" "rhog_file" "rhog_flag" "rhog_nvmax" "rhog_nvmin" "rhombohedral" "rlderiv" "rm" "rmatch_augfun" "rmatch_augfun_nc" "rmax" "rpwe" "rytoev_fact" "sample_bias" "save_file" "saverho" "scf_must_converge" "scissor" "screening_parameter" "search_sym" "single_pole" "skip_dw" "skip_equivalence_q" "skip_type" "skip_upperfan" "smearing" "space_group" "spin_component" "start" "start_e" "start_irr" "start_k" "start_q" "starting_charge" "starting_magnetization" "starting_ns_eigenvalue" "starting_spin_angle" "startingpot" "startingwfc" "step" "string_method" "sum_pertq" "sw_len" "symm_type" "tabps" "tcg" "td" "tdosinboxes" "tefield" "temp_kelvin" "temp_req" "temph" "tempw" "thresh_init" "title" "tk_plot" "tm" "tolp" "tolw" "tot_charge" "tot_magnetization" "tprnfor" "tqr" "tr2" "tr2_ph" "tran_file" "tran_prefix" "tranp" "trans" "trust_radius_ini" "trust_radius_max" "trust_radius_min" "ts_vdw" "ts_vdw_econv_thr" "ts_vdw_isolated" "tstress" "uniqueb" "units" "upscale" "use_ace" "use_all_frac" "use_freezing" "use_masses" "use_paw_as_gipaw" "vdW_analysis" "vdw" "vdw_corr" "verbosity" "vkb" "vkbg_file" "vkbg_flag" "vscg_file" "vscg_flag" "vxc0_file" "vxc0_flag" "vxc_diag_nmax" "vxc_diag_nmin" "vxc_file" "vxc_flag" "vxc_integral" "vxc_offdiag_nmax" "vxc_offdiag_nmin" "vxc_zero_rho_core" "vxcdiag" "vxcg_file" "vxcg_flag" "w_1" "w_2" "w_T_npol" "weight" "wf_collect" "wf_efield" "wf_friction" "wf_q" "wf_switch" "wfc_extrapolation" "wfcdir" "wfdt" "wffort" "wfng_dk1" "wfng_dk2" "wfng_dk3" "wfng_file" "wfng_flag" "wfng_kgrid" "wfng_nband" "wfng_nk1" "wfng_nk2" "wfng_nk3" "wfng_nvmax" "wfng_nvmin" "wfng_occupation" "wfsd" "what" "which_augfun" "wmass" "wpot_dir" "write_coulomb" "writev" "x0" "x_gamma_extrapolation" "xdm" "xdm_a1" "xdm_a2" "xmin" "xmlfile_full" "xmlfile_part" "zasr" "zed" "zeu" "zgate" "zue" "zval" ))

;; qe's cards & keywords
(defvar qe-cards (list "ADDITIONAL_K_POINTS" "ATOMIC_FORCES" "ATOMIC_POSITIONS" "ATOMIC_SPECIES" "ATOMIC_VELOCITIES" "CELL_PARAMETERS" "CLIMBING_IMAGES" "CONSTRAINTS" "K_POINTS" "OCCUPATIONS" "PLOT_WANNIER" "REF_CELL_PARAMETERS" "TOTAL_CHARGE" "on_step" ))

;; qe's flags
(defvar qe-flags (list "a.u" "alat" "angstrom" "automatic" "bohr" "crystal" "crystal_b" "crystal_c" "crystal_sg" "gamma" "tpiba" "tpiba_b" "tpiba_c" ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar qe-open-supercards-regexp   (regexp-opt qe-open-supercards   'symbols)) ; may not exists
(defvar qe-closed-supercards-regexp (regexp-opt qe-closed-supercards 'symbols)) ; may not exists

(defvar qe-cards-regexp (regexp-opt
			    (append qe-cards qe-open-supercards) 'symbols))
(defvar qe-flags-regexp (regexp-opt qe-flags 'symbols))

(defvar qe-namelist-face (cons (regexp-opt (append qe-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar qe-variable-face (cons (regexp-opt qe-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar qe-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar qe-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append qe-cards qe-open-supercards qe-closed-supercards) t)))
;;
(defvar qe-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar qe-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt qe-namelists t)))
(defvar qe-indent-re             (concat "^[ \t]*" (regexp-opt qe-cards     t)))


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


(font-lock-add-keywords 'qe-mode (list
				     qe-namelist-face 
				     qe-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'qe-mode '(
				      (qe-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      (qe-cards-matcher 1 font-lock-keyword-face t)
				      (qe-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar qe-keywords '(qe-namelist-face qe-variable-face))
(defvar qe-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar qe-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `qe-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun qe-indent-line ()
  "Indent current line according to qe input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at qe-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at qe-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at qe-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at qe-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at qe-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the qe-mode as derived-mode

(define-derived-mode qe-mode prog-mode 
  "QE-generic"  
  "Major mode for editing Qunatum-ESPRESSO input files (QE-generic mode)"  
  (setq font-lock-defaults '((qe-keywords) nil t))
  (set (make-local-variable 'indent-line-function) 'qe-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" qe-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" qe-mode-syntax-table)
  (modify-syntax-entry ?=  " " qe-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq qe-namelists nil)
(setq qe-vars nil)
(setq qe-cards nil)
(setq qe-flags nil)
(setq qe-open-supercards   nil)
(setq qe-closed-supercards nil)


(require 'qe-funcs)

(provide 'qe-mode)

