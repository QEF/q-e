source commands.tcl

set ::guib::settings(filename_only_tail) 1

module PH\#auto -title "PWSCF GUI: module PH.x" -script {

    readfilter  ::pwscf::phReadFilter
    #writefilter ::pwscf::phWriteFilter
    
    line job_title -name "Job Title" {
	var title_line {
	    -label    "Job title:"
	    -fmt      %S
	}
    }

    page inputph -name "INPUTPH" {
	
	namelist inputpp -name "INPUTPH" {

	    page files -name "Files/Diretories" {
		var outdir {
		    -label    "Temporary directory where punch file resides (outdir):"
		    -widget   entrydirselectquote
		    -fmt      %S -validate string
		}

		var prefix -label "Prefix of data file saved by PW.X (prefix):" \
		    -widget   [list entrybutton "Prefix ..." [list ::pwscf::phSelectPunchFile $this prefix]] \
		    -fmt      %S -validate string
		

		separator -label "--- Output Data Files ---"

		var fildyn {
		    -label    "File containing the dynamical matrix (fildyn):" 
		    -widget   entryfileselectquote
		    -fmt      %S -validate string
		}
		var fildrho {
		    -label    "File containing the charge density variations (fildrho):"
		    -widget   entryfileselectquote
		    -fmt      %S -validate string
		}
		var fildvscf {
		    -label    "File containing the potential variation (fildvscf):"
		    -widget   entryfileselectquote
		    -fmt      %S -validate string
		}
	    }
	    
	    page calcs -name "What to Compute" {

		var ldisp {
		    -label     "Compute phonon dispersions (ldisp):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		var nogg {
		    -label     "Disable gamma_gamma tricks (nogg):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		var ldiag {
		    -label     "Force diagonalization of the incomplete dynamical matrix (ldiag):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		#group nq -name q-grid {
		#    group nq_grid {
		#	packwidgets left
		#	var nq1 {
		#	    -label     "nq\#1:"
		#	    -validate posint
		#	}
		#	var nq2 {
		#	    -label     "nq\#2:"
		#	    -validate posint
		#	}
		#	var nq3 {
		#	    -label     "nq\#3:"
		#	    -validate posint
		#	}
		#    }
		#    group iq_grid {
		#	packwidgets left
		#	var iq1 {
		#	    -label     "iq\#1:"
		#	    -validate posint
		#	}
		#	var iq2 {
		#	    -label     "iq\#2:"
		#	    -validate posint
		#	}
		#	var iq3 {
		#	    -label     "iq\#3:"
		#	    -validate posint
		#	}
		#    }
		#}

		var trans {
		    -label     "Compute phonons for a single q vector (trans):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}
		
		var epsil {
		    -label     "Compute the macroscopic dielectric constant for q=0 (epsil):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		var electron_phonon {
		    -label     "Compute electron-phonon (electron_phonon):"
		    -textvalue {
			"electron-phonon lambda coefficients <simple>"
			"electron-phonon by interpolation <interpolated>"
			"el.-phon. lambda coeff. via opt.tetrahedron <lambda_tetra>"
			"phonon linewidth via opt.tetrahedron <gamma_tetra>"
		    }
		    -value     {'simple' 'interpolated' 'lambda_tetra' 'gamma_tetra'}
		    -widget    optionmenu
		}
		
		var lrpa {
		    -label     "Compute dielectric constant with RPA and dV_xc=0 (lrpa):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}

		var lnoloc {
		    -label     "Compute dielectric constant with dV_H=0 and  dV_xc=0 (lnoloc):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}
		
		var fpol {
		    -label     "Compute dynamic polarizabilities (fpol):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}

		var lshift_q {
		    -label "Use a wave-vector half-step grid displaced (lshift_q):"
		    -widget    radiobox
		    -textvalue { Yes No }	      
		    -value     { .true. .false. }
		}

		var zeu {
		    -label     "Compute effective charges from the dielectric responses (zeu):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		var zue {
		    -label     "Compute effective charges from the phonon density responses (zue):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		var elop {
		    -label     "Compute electro-optic coefficients (elop):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		var lraman {
		    -label     "Compute Raman coefficients (lraman):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}

		var search_sym {
		    -label     "Perform mode symmetry analysis (search_sym):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}

		group ramanthreshold -name "Thresholds for Raman" -decor normal {
		    var eth_rps {
			-label    "Threshold for calculation of x|Psi> (eth_rps):"
			-validate fortranreal
		    }
		    var eth_ns {
			-label    "Threshold for non-scf wavefunction calculation (eth_ns):"
			-validate fortranreal
		    }
		    var dek {
			-label    "Delta k used for wavefunction derivation wtr k (dek)::"
			-validate fortranreal
		    }
		}

		##separator -label "--- Potential variations ---"

		group dvscf -name "Potential variations" -decor normal {
		    # to trick the chek_gui only
		    var dvscf_star
		    #/
		    var dvscf_star_open {
			-variable  dvscf_star%open
			-label     "Compute potential variations at all q-points in the star (dvscf_star%open):"
			-textvalue {Yes No}
			-value     {.true. .false.}
			-widget    radiobox
		    }
		    group dvscf_star_specs -decor none {
			var dvscf_star_dir {
			    -variable dvscf_star%dir
			    -label    "Directory where to store rotated dvscf files (dvscf_star%dir):" 
			    -widget   entrydirselectquote
			    -fmt      %S -validate string
			}
			var dvscf_star_ext {
			    -variable dvscf_star%ext
			    -label    "Extension to use for the name of dvscf files (dvscf_star%ext):" 
			    -fmt      %S -validate string
			}
			var dvscf_star_basis {
			    -variable dvscf_star%basis
			    -label    "Basis used for the dvscf files (dvscf_star%basis):" 
			    -widget   radiobox
			    -fmt      %s
			    -value {
				'cartesian'
				'modes'
			    }
			    -textvalue {
				"basis of cartesian 1-atom displacements <cartesian>"
				"basis of the modes at the rotated q-point <modes>"
			    }
			}
			var dvscf_star_pat {
			    -variable  dvscf_star%pat
			    -label     "Save displacement patterns and q vector (dvscf_star%pat):"
			    -textvalue {Yes No}
			    -value     {.true. .false.}
			    -widget    radiobox
			}
		    }
		}

		##separator -label "--- Perturbation of the charge density ---"

		group drho -name "Perturbation of the charge density" -decor normal {
		    # to trich the chek_gui only
		    var drho_star
		    #/
		    var drho_star_open {
			-variable  drho_star%open
			-label     "Compute charge density perturbation at all q-points in the star (drho_star%open):"
			-textvalue {Yes No}
			-value     {.true. .false.}
			-widget    radiobox
		    }
		    group drho_star_specs -decor none {
			var drho_star_dir {
			    -variable drho_star%dir
			    -label    "Directory where to store rotated drho files (drho_star%dir):" 
			    -widget   entrydirselectquote
			    -fmt      %S -validate string
			}
			var drho_star_ext {
			    -variable drho_star%ext
			    -label    "Extension to use for the name of drho files (drho_star%ext):" 
			    -fmt      %S -validate string
			}
			var drho_star_basis {
			    -variable drho_star%basis
			    -label    "Basis used for the drho files (drho_star%basis):" 
			    -widget   radiobox
			    -fmt      %s
			    -value {
				'cartesian'
				'modes'
			    }
			    -textvalue {
				"basis of cartesian 1-atom displacements <cartesian>"
				"basis of the modes at the rotated q-point <modes>"
			    }
			}
			var drho_star_pat {
			    -variable  drho_star%pat
			    -label     "Save displacement patterns and q vector (drho_star%pat):"
			    -textvalue {Yes No}
			    -value     {.true. .false.}
			    -widget    radiobox
			}
		    }
		}
	    }

	    page misc -name "Control options" {

		var recover {
		    -label "Restart from an interrupted run (recover):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}

		var low_directory_check {
		    -label "(low_directory_check):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}
		
		var only_init {
		    -label "Calculate only bands and other initialization quantities (only_init):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}
		
		var qplot {
		    -label "Read a list of q-points from input (qplot):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}
		
		var q2d {
		    -label "(q2d):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}
		
		var q_in_band_form {
		    -label "(q_in_band_form):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}
		
		var lqdir {
		    -label     "Creates inside outdir a separate subdirectory for each q vector (lqdir):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		}

		auxilvar reps_type {
		    -label "How to specify irreducible representations:"
		    -textvalue {"with start_irr/last_irr" "with nat_todo"}
		    -value { start_last_irr nat_todo }
		    -widget radiobox
		}

		group irrep_spec -name "Specification of irreducible representation(s)" -decor normal {

		    group start_last_irr -name "Range of irreducible representations" -decor normal {
			var start_irr  {
			    -validate posint 
			    -label "First irreducible representations in range (first_irr):"
			}
			var last_irr  {
			    -validate posint
			    -label "Last irreducible representations in range (last_irr):"
			}
		    }
		    
		    var nat_todo {
			-label "Number of atom to be displaced (nat_todo):"
			-validate nonnegint
			-widget spinint
		    }
		    
		    var modenum {
			-validate nonnegint
			-label "Index of the irreducible representation for single-mode calculation (modenum):"
		    }
		}

		group q_spec -name "q-point specification" -decor normal {
		    var start_q  {
			-validate posint
			-label "First q-point in range (start_q):"
		    }
		    var last_q  {
			-validate posint
			-label "Last q-point in range (last_q):"
		    }
		}

		separator -label "--- Atomic Masses ---"
		
		auxilvar ntyp {
		    -label   "Number of types of atoms in the unit cell (ntyp):"
		    -validate posint
		    -fmt      %d
		    -default  1
		    -widget   spinint  
		}			    

		dimension amass {
		    -label     "Atomic mass [amu] of each atomic type"
		    -validate  fortranreal
		    -start     1
		    -end       1
		}

		separator -label "--- Misc control options ---"

		var verbosity {
		    -label     "Verbosity of output (verbosity):"
		    -textvalue {"low" "high"}
		    -value     {'low' 'high'}
		    -widget    optionmenu
		}
		var reduce_io {
		    -label    "Reduce I/O to the strict minimum (reduce_io):"
		    -textvalue {Yes No}
		    -value     {.true. .false.}
		    -widget    radiobox
		    -fmt       %s
		}
		var max_seconds {
		    -label    "Maximum allowed CPU run-time [in seconds] (max_seconds):"
		    -validate posint
		    -widget   spinint
		    -fmt      %d
		}

		separator -label "--- SCF settings ---"

		var niter_ph {
		    -label    "Maximum number of iterations in an SCF step (niter_ph):"
		    -widget   spinint
		    -fmt      %d
		}
		var tr2_ph {
		    -label    "Threshold for selfconsistency (tr2_ph):"
		    -validate fortranreal
		}
		var alpha_mix {
		    -variable alpha_mix(1) 
		    -label    "Mixing factor for updating the SCF potential (alpha_mix(1)):"
		    -validate fortranreal
		}
		var nmix_ph {
		    -label    "Number of iterations used in mixing of potential (nmix_ph):"
		    -widget   spinint
		    -fmt      %d
		}	
	    }	

	    page q_k_points -name "q- and k-points" {
	    
		separator -label "--- K-point mesh + shift ---" 
		
		group k_points {
		    packwidgets left
		    var nk1 -label "nk1:"  -widget spinint  -validate posint
		    var nk2 -label "nk2:"  -widget spinint  -validate posint
		    var nk3 -label "nk3:"  -widget spinint  -validate posint
		}
		group k_shift {
		    packwidgets left
		    var k1 -label "k1:"  -widget spinint  -validate binary
		    var k2 -label "k2:"  -widget spinint  -validate binary
		    var k3 -label "k3:"  -widget spinint  -validate binary
		}
		
		separator -label "--- Q-point mesh ---" 
		
		group q_points {
		    packwidgets left
		    var nq1 -label "nq1:"  -widget spinint  -validate posint
		    var nq2 -label "nq2:"  -widget spinint  -validate posint
		    var nq3 -label "nq3:"  -widget spinint  -validate posint
		}    
	    }
	}
    }

    page suffixCards -name "Suffix cards" {
	line xq_list -name "The phonon wavevector" {
	    packwidgets left
	    var xq1 {
		-variable xq(1)
		-label    "xq(1):"
		-validate fortranreal
	    }
	    var xq2 {
		-variable xq(2)
		-label    "xq(2):"
		-validate fortranreal
	    }
	    var xq3 {
		-variable xq(3)
		-label    "xq(3):"
		-validate fortranreal
	    }
	}

	group qPointsSpec -name "q-points specification" -decor normal {
	    line nqs_line -decor none {
		var nqs -label "Number of q-points:" -validate posint -widget spinint -default 1 -outfmt %3d
	    }
	    table qPoints {
		-caption "Enter q-points coordinates and weights:"
		-head     {xq1 xq2 xq3 nq}
		-validate {fortranreal fortranreal fortranreal int}
		-cols     4
		-rows     1
		-outfmt   {%13.8f  %13.8f  %13.8f  %3d}
	    }
	}

	line atom_disp_line -name "List of atoms to displace:" {
	    var nat_todo_list {
		-label "Indices of atoms (comma or whitespace separated):"
	    }
	}
    } 
    # ----------------------------------------------------------------------
    # take care of specialities
    # ----------------------------------------------------------------------
    source ph-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source ph-help.tcl
}
