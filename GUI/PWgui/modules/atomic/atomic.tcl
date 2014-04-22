source commands.tcl

set ::guib::settings(filename_only_tail) 1

# maximum number of nconf
set ::pwscf::atomic_max_nconf 20

module LD1\#auto -title "PWSCF GUI: module LD1.x" -script {

    readfilter  ::pwscf::atomicReadFilter
    writefilter ::pwscf::atomicDFTFilter

    #
    # PAGE: General
    #
    page general -name "General" {
	#
	namelist input -name "INPUT" {
	    #
	    required {
		var iswitch {
		    -label     "Type of calculation (iswitch):"
		    -widget    radiobox
		    -textvalue {
			"All-Electron"
			"PseudoPotential Generation"
			"PseudoPotential Test"
			"LDA-1/2 correction"
		    }
		    -value     {1 3 2 4}
		    -default "All-Electron" 
		}
		var zed {
		    -label    "Nuclear charge (zed):"
		    -validate fortranreal
		}
		var atom  {
		    -label    "Atomic symbol (atom):"
		    -fmt      %S
		    -validate string
		}
		var config  {
		    -label    "Electronic configuration (config):"
		    -fmt      %S
		    -validate string
		}
		var relpert {
		    -label    "Compute relativistic corrections to non-relativistic\nKohn-Sham energy levels (relpert):"
		    -widget   radiobox
		    -textvalue {Yes No} -value {.true. .false.}
		}
		var rel_dist {
		    -label "How to fill the electronic states (rel_dist):"
		    -value {'energy' 'average'} 
		    -textvalue {
			"by increasing energy of states  <energy>" 
			"distrubute electrons uniformly  <average>"
		    }
		    -widget optionmenu
		}
		var write_coulomb {
		    -label "Write a fake pseuopotential (write_coulomb):"
		    -textvalue {Yes No} -value {.true. .false.} -widget radiobox
		}
		var rel  {
		    -label    "Relativistic Effects (rel):"
		    -widget    radiobox
		    -textvalue {"Non Relativistic (Schroedinger)"
			"Scalar Relativistic"
			"Full Relativistic (Dirac)"}
		    -value     {0 1 2}
		}
		var lsmall {
		    -label    "Write on files the small component (lsmall):"
		    -widget   radiobox
		    -textvalue {Yes No} -value {.true. .false.}
		}
		var max_out_wfc {
		    -label    "Maximum number of atomic wavefunctions written (max_out_wfc):"
		    -validate posint 
		    -widget spinint
		}

		
		var noscf {
		    -label    "Skip charge calculation and ignore the occupations (noscf):"
		    -widget    radiobox
		    -textvalue {Yes No} -value {.true. .false.}
		}
		var lsd  {
		    -label    "LSDA Spin Polarization (lsd):"
		    -widget    radiobox
		    -textvalue {"No" "Yes"}
		    -value     {0 1}
		    -default  0
		}
		var dft  {
		    -label    "Exchange-Correlation (dft):"
		    -widget    radiobox
		    -textvalue {"Ceperley-Alder LDA, Perdew-Zunger data (PZ)"
			"Perdew-Wang  GGA (PW91)"
			"Becke-Perdew GGA (BP)"
			"Perdew-Becke-Ernzerhof (PBE)"
			"Becke-Lee-Yang-Parr (BLYP)"
			"Other"
		    }
		    -value     {'PZ' 'PW91' 'BP' 'PBE' 'BLYP' 'REPLACE_ME'}
		}
		var dft_  {
		    -label    "Enter Exchange-Correlation Functional (dft_):"
		    -fmt      %S -validate string
		}
	    }
	    #
	    optional {
		var title {
		    -label    "Job name or Comment (optional):"
		    -fmt      %S -validate string
		}
		var prefix  {
		    -label    "Prefix for output file names (prefix):"
		    -fmt      %S -validate string
		}
		var verbosity {
		    -label "Verbosity of output (verbosity):"
		    -textvalue {low high} -value {'low' 'high'}
		    -widget radiobox
		}
		var file_charge {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the all-electron total charge (file_charge):"
		    -fmt      %S -validate string
		}
		var beta  {
		    -label    "Mixing parameter for self-consistency (beta):"
		    -validate  fortranposreal
		}
		var tr2  {
		    -label    "Convergence Threshold for self-consistency (tr2):"
		    -validate  fortranposreal
		}
		
		var latt  {
		    -label    "Latter Correction (latt):"
		    -widget    radiobox
		    -textvalue {"No" "Yes"}
		    -value     {0 1}
		    -default  0
		}
		var isic  {
		    -label    "Self-Interaction Correction (isic):"
		    -widget    radiobox
		    -textvalue {"No" "Yes"}
		    -value     {0 1}
		    -default  0
		}
		var rytoev_fact {
		    -label "Conversion factor from Ry to eV (rytoev_fact):"
		    -validate fortranposreal
		}
		var cau_fact {
		    -label "Speed of light in a.u. (cau_fact):"
		    -validate fortranposreal
		}
		var vdw {
		    -label    "Calculation of van der Waals coefficients:"
		    -widget    radiobox
		    -textvalue {"No" "Yes"}
		    -value     {.true. .false.}
		    -default  .false.
		}

		separator -label "--- Grid: r(i)= exp (xmin + (i-1)*dx) / Z ---"

		var xmin  {
		    -label    "Grid parameter (xmin):"
		    -validate  fortrannegreal
		}
		var dx  {
		    -label    "Grid parameter (dx):"
		    -validate  fortranposreal
		}
		var rmax  {
		    -label    "Grid parameter (rmax):"
		    -validate  fortranposreal
		}

		separator -label "--- Parameters for Logarithmic derivative calculation ---"
		var nld  {
		    -label    "Number of logarithmic derivatives to calculate (nld):"
		    -widget   optionmenu
		    -textvalue {
			"None" "1" "2" "3" "4"
		    }
		    -value { 0 1 2 3 4 }
		    -default  "None"
		}
		var rlderiv  {
		    -label    "Radius at which logarithmic derivatives are calculated (rlderiv):"
		    -validate  fortranposreal
		}
		var eminld  {
		    -label    "Minimum energy [Ry] for Plotting (eminld):"
		    -validate  fortranreal
		}
		var emaxld  {
		    -label    "Maximum energy [Ry] for Plotting (emaxld):"
		    -validate  fortranreal
		}
		var deld  {
		    -label    "Plotting in steps of Delta E [Ry] (deld):"
		    -validate  fortranposreal
		}
		var rpwe {
		    -label "Radius [a.u.] for partial wave expansions (rpwe):"
		    -validate  fortranposreal
		}
		#
	    }
	    #
	}
    }
    
    # all-electron cards

    group AE_cards -name "All-electron cards" -decor normal {
	line nwf_line -decor none {
	    var nwf -label "Number of states:" -validate posint -widget spinint -default 1 -outfmt %3d
	}
	table AE_wfs {
	    -caption "Wavefunction specifications:"
	    -head    {Label N L Occupancy "Spin index"}
	    -cols 5
	    -rows 1
	    -validate {string int int fortranreal fortranreal}
	    -outfmt   {"  %s  "  %3d %3d %f " %s"}
	}
    }

    #
    # PAGE: pseudo
    #
    page pseudoPotential -name "PseudoPotential Generation" {
	#
	namelist inputp -name "InputP" {
	    #
	    required {
		var pseudotype {
		    -label     "Type of PseudoPotential (pseudotype):"
		    -widget    radiobox
		    -textvalue {
			"Norm Conserving, one channel per angular momentum"
			"Norm Conserving, more than one channel per angular momentum"
			"Ultrasoft PP or PAW"}
		    -value     {1 2 3}
		}
		var file_pseudopw {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output PP (file_pseudopw):"
		    -fmt      %S -validate string
		}
		var lloc {
		    -label    "Local Potential channel (lloc):"
		    -widget   optionmenu
		    -textvalue {
			"all-electron potential" "L=0" "L=1" "L=2" "L=3" "L=4"
		    }
		    -value { -1 0 1 2 3 4 }
		    -default  "all-electron potential"
		}
		var rcloc {
		    -label    "Matching Radius for Local Potential\n  [optional if an L channel was specified] (rcloc):"
		    -validate  fortranposreal
		}
		var new_core_ps {
		    -label "Pseudize the core charge with bessel functions (new_core_ps):"
		    -textvalue {Yes No} -value {.true. .false.} -widget radiobox
		}
		var tm {
		    -label    "Type of pseudization procedure (tm):"
		    -widget    radiobox
		    -textvalue {
			"Troullier-Martins"
			"Rabe-Rappe-Kaxiras-Joannopoulos"
		    }
		    -value     {.true. .false.}
		}
		var rho0 {
		    -label    "Charge at r=0 ( rho0):"
		    -validate fortrannonnegreal
		    -default  0.0
		}

	    }
	    #
	    optional {
		var nlcc {
		    -label    "Nonlinear Core Correction (nlcc):"
		    -widget   radiobox
		    -textvalue {"No" "Yes"}
		    -value     {.false. .true.}
		    -default  "No"
		}
		var rcore {
		    -label    "Matching Radius for Nonlinear Core Correction (rcore):"
		    -validate  fortranposreal
		}
		var lpaw {
		    -label    "Generate PAW dataset [experimental feature] (lpaw):"
		    -widget   radiobox
		    -textvalue {"No" "Yes"}
		    -value     {.false. .true.}
		    -default  "No"
		}
		var which_augfun {
		    -label "Which pseudization augmentation functions (which_augfun):"
		    -value {'AE' 'BESSEL' 'GAUSS' 'BG' 'PSQ'}
		    -widget optionmenu
		}
		var rmatch_augfun {
		    -label "Pseudization radius for the augmentation functions (rmatch_augfun):"
		    -validate fortranposreal
		}
		var rmatch_augfun_nc {
		    -label "Pseudize augmentation functions from origin to min(rcut(ns),rcut(ns1) (rmatch_augfun_nc):"
		    -widget   radiobox
		    -textvalue {"No" "Yes"}
		    -value     {.false. .true.}
		}
		
		var lsave_wfc {
		    -label "Save all-electron and pseudo wavefunctions (lsave_wfc):"
		    -textvalue {Yes No} -value {.true. .false.} -widget radiobox	
 		}

		var lgipaw_reconstruction {
		    -label "Generate all-electron data for GIPAW (lgipaw_reconstruction):"
		    -textvalue {Yes No} -value {.true. .false.} -widget radiobox	
 		}

		var use_paw_as_gipaw {
		    -label "Use PAW as GIPAW (use_paw_as_gipaw):"
		    -textvalue {Yes No} -value {.true. .false.} -widget radiobox	
 		}
		
		var author {
		    -label "Name of the author (author):"
		    -validate string
		}
		var file_recon {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing data needed for PAW reconstruction (file_recon):"
		    -fmt      %S -validate string
		}
		var file_chi {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output pseudo-orbitals (file_chi):"
		    -fmt      %S -validate string
		}
		var file_beta {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output beta functions (file_beta):"
		    -fmt      %S -validate string
		}
		var file_qvan {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output Q functions (file_qvan):"
		    -fmt      %S -validate string
		}
		var file_screen {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output screeing potential (file_screen):"
		    -fmt      %S -validate string
		}
		var file_core {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output total and core charge (file_core):"
		    -fmt      %S -validate string
		}
		var file_wfcaegen {
		    -label    "Name of the output file with all-electron wfc (file_wfcaegen):"
		    -validate string
		    -fmt      %S -validate string
		}
		var file_wfcncgen {
		    -label    "Name of the output file with norm-conserving wfc (file_wfcncgen):"
		    -validate string
		    -fmt      %S -validate string
		}
		var file_wfcusgen {
		    -label    "Name of the output file with ultra-soft wfc (file_wfcusgen):"
		    -validate string
		    -fmt      %S -validate string
		}
		var zval {
		    -label    "Valence charge [for special cases only!] (zval):"
		    -validate  fortranreal
		}
	    }
	    #
	}
		
	group PP_cards -name "All-electron cards" -decor normal {
	    line PP_nwf_line -decor none {
		var nwfs -label "Number of states:" -validate posint -widget spinint -default 1 -outfmt %3d
	    }
	    table PP_wfs {
		-caption  "Enter Wavefunctions to be Pseudized:"
		-head     {Label N L Occupancy Energy Rcut "US Rcut" "Tot.ang.moment"}
		-validate {whatever int int fortranreal fortranreal fortranreal fortranreal fortranreal}
		-cols     8
		-rows     1		    
		-outfmt   {"  %s  " %3d %3d %8.3f %8.3f %6.2f %6.2f}
	    }	
	}
    }
    #
    # PAGE: test
    #
    page testing -name "PseudoPotential Test" {
	#
	namelist test -name "TEST" {
	    scriptvar old_nconf
	    var nconf {
		-label "Number of testing configurations (nconf):" 
		-widget spinint 
		-validate posint
	    }
	    var file_pseudo {
		-widget   entryfileselectquote
		-label    "Name of the file containing the input PP (file_pseudo):"
		-fmt      %S -validate string
	    }
	    dimension configts {
		-label     "Test electronic configurations:"
		-start     1
		-end       1
		-validate  string
		-fmt       %S -validate string
	    }
	    dimension lsdts {
		-label     "Use LSDA in test of electronic configurations:"
		-start     1
		-end       1
		-widget    radiobox
		-textvalue {"No" "Yes"} -value {0 1}
	    }
	    var frozen_core {
		-label "Calculate only the core wavefunctions of 1st configuration (frozen_core):"
		-value {.false. .true.} -textvalue {"Yes" "No"}
		-widget radiobox
	    } 
	    var rm  {
		-label    "Box radius for spherical Bessel basis set (rm):"
		    -validate  fortranposreal
	    }
	    var ecutmin  {
		-label    "Minimum energy cutoff [Ry] for ghost/convergence test (ecutmin):"
		-validate  fortranposreal
	    }
	    var ecutmax {
		-label    "Maximum energy cutoff [Ry] for ghost/convergence test (ecutmax):"
		-validate  fortranposreal
	    }
	    var decut  {
		-label    "Step [Ry] for cutoff in ghost/convergence test (decut):"
		-validate  fortranposreal
	    }
	    var rcutv {
		-label "Cutoff distance for the inclusion of LDA-1/2 potential (rcutv):"
		-validate  fortranreal
	    }
	}

	# testing configuration cards

	group test_cards -name "Testing configurations not specified by configts(*)" -decor normal {
	    for {set ic 1} {$ic <= $::pwscf::atomic_max_nconf} {incr ic} {
		group test_conf_$ic -name "Testing configuration \#.$ic" -decor normal [subst -nocommands {
		    line nwfts_line_$ic -decor none [list var nwfts_$ic -label "Number of states:" -validate posint -widget spinint -default 1 -outfmt %3d]
		    
		    table test_wfs_$ic {
			-caption  "Wavefunctions specifications:"
			-head     {Label N L Occupancy enerts rcutts rcutusts "Spin index"}
			-cols     8
			-rows     1
			-validate {string int int fortranreal}
			-outfmt   {"  %s  "  %3d %3d %8.3f %8.3f %6.2f %6.2f " %s"}
		    }
		}]								]
	    }
	}
    }

    # ----------------------------------------------------------------------
    # take care of specialities
    # ----------------------------------------------------------------------
    source atomic-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source atomic-help.tcl
}
