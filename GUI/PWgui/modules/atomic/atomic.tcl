source commands.tcl

set ::guib::settings(filename_only_tail) 1

module LD1\#auto -title "PWSCF GUI: module LD1.x" -script {

    readfilter  ::pwscf::atomicReadFilter
    writefilter ::pwscf::atomicDFTFilter

    #
    # PAGE: General
    #
    page general -name "General" {
	namelist input -name "INPUT" {
	    #
	    required {
		var iswitch {
		    -label     "Type of calculation:"
		    -widget    radiobox
		    -textvalue {"All-Electron"
			"PseudoPotential Generation"
			"PseudoPotential Test"}
		    -value     {1 3 2}
		    -default "All-Electron" 
		}
		var atom  {
		    -label    "Atomic Symbol:"
		    -fmt      %S
		}
		var config  {
		    -label    "Electronic Configuration:"
		    -fmt      %S
		}
		var rel  {
		    -label    "Relativistic Effects:"
		    -widget    radiobox
		    -textvalue {"Non Relativistic (Schroedinger)"
			"Scalar Relativistic"
			"Full Relativistic (Dirac)"}
		    -value     {0 1 2}
		}
		var lsd  {
		    -label    "LSDA Spin Polarization:"
		    -widget    radiobox
		    -textvalue {"No" "Yes"}
		    -value     {0 1}
		    -default  0
		}
		var dft  {
		    -label    "Exchange-Correlation:"
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
		    -label    "Enter Exchange-Correlation Functional:"
		    -fmt      %S
		}
	    }
	    #
	    optional {
		var title {
		    -label    "Job name or Comment (optional):"
		    -fmt      %S
		}
		var prefix  {
		    -label    "Prefix for output file names:"
		    -fmt      %S
		}
		var beta  {
		    -label    "Mixing parameter for self-consistency:"
		    -validate  fortranposreal
		}
		var tr2  {
		    -label    "Convergence Threshold for self-consistency:"
		    -validate  fortranposreal
		}
		
		var latt  {
		    -label    "Latter Correction:"
		    -widget    radiobox
		    -textvalue {"No" "Yes"}
		    -value     {0 1}
		    -default  0
		}
		var isic  {
		    -label    "Self-Interaction Correction:"
		    -widget    radiobox
		    -textvalue {"No" "Yes"}
		    -value     {0 1}
		    -default  0
		}
		separator -label "--- Grid: r(i)= exp (xmin + (i-1)*dx) / Z ---"
		var xmin  {
		    -label    "Grid parameter xmin:"
		    -validate  fortrannegreal
		}
		var dx  {
		    -label    "Grid parameter dx:"
		    -validate  fortranposreal
		}
		var rmax  {
		    -label    "Grid parameter rmax:"
		    -validate  fortranposreal
		}
		separator -label "--- Parameters for Logarithmic derivative calculation ---"
		var nld  {
		    -label    "Number of logarithmic derivatives to calculate:"
		    -widget   optionmenu
		    -textvalue {
			"None" "1" "2" "3" "4"
		    }
		    -value { 0 1 2 3 4 }
		    -default  "None"
		}
		var rlderiv  {
		    -label    "Radius at which logarithmic derivatives are calculated:"
		    -validate  fortranposreal
		}
		var eminld  {
		    -label    "Minimum energy (Ry) for Plotting:"
		    -validate  fortranreal
		}
		var emaxld  {
		    -label    "Maximum energy (Ry) for Plotting:"
		    -validate  fortranreal
		}
		var deld  {
		    -label    "Plotting in steps of Delta E (Ry):"
		    -validate  fortranposreal
		}
		#
	    }
	    #
	}
    }

    #
    # Page: pseudo
    #
    page pseudoPotential -name "PseudoPotential Generation" {
	namelist inputp -name "InputP" {
	    #
	    required {
		var pseudotype {
		    -label     "Type of PseudoPotential:"
		    -widget    radiobox
		    -textvalue {
			"Norm Conserving, one channel per angular momentum"
			"Norm Conserving, more than one channel per angular momentum"
			"UltraSoft"}
		    -value     {1 2 3}
		}
		var file_pseudopw {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output PP:"
		    -fmt      %S
		}
		var lloc {
		    -label    "Local Potential channel:"
		    -widget   optionmenu
		    -textvalue {
			"all-electron potential" "L=0" "L=1" "L=2" "L=3" "L=4"
		    }
		    -value { -1 0 1 2 3 4 }
		    -default  "all-electron potential"
		}
		var rcloc {
		    -label    "Matching Radius for Local Potential (optional if a L channel was specified):"
		    -validate  fortranposreal
		}
		var tm {
		    -label    "Type of pseudization procedure:"
		    -widget    radiobox
		    -textvalue {
			"Troullier-Martins"
			"Rabe-Rappe-Kaxiras-Joannopoulos"
		    }
		    -value     {.true. .false.}
		}
		var rho0 {
		    -label    "Charge at r=0:"
		    -validate fortrannonnegreal
		    -default  0.0
		}

	    }
	    #
	    optional {
		var nlcc {
		    -label    "Nonlinear Core Correction:"
		    -widget   radiobox
		    -textvalue {"No" "Yes"}
		    -value     {.false. .true.}
		    -default  "No"
		}
		var rcore {
		    -label    "Matching Radius for Nonlinear Core Correction:"
		    -validate  fortranposreal
		}
		var lpaw {
		    -label    "Generate PAW dataset (experimental feature):"
		    -widget   radiobox
		    -textvalue {"No" "Yes"}
		    -value     {.false. .true.}
		    -default  "No"
		}
		var file_recon {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing data needed for PAW reconstruction:"
		    -fmt      %S
		}
		var file_chi {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output pseudo-orbitals:"
		    -fmt      %S
		}
		var file_beta {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output beta functions:"
		    -fmt      %S
		}
		var file_qvan {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output Q functions:"
		    -fmt      %S
		}
		var file_screen {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output screeing potential:"
		    -fmt      %S
		}
		var file_core {
		    -widget   entryfileselectquote
		    -label    "Name of the file containing the output total and core charge:"
		    -fmt      %S
		}
	    }
	    #
	}
	#
	# This section is not yet working
	#
	group pseudization -name "Specify states to be pseudized" -decor normal {
	    line Nwfs -name "Number of states" -decor none {
		var nwfs -label "Number of states:" -widget spinint -validate posint -default 1
	    }

	    line Wfs -name "Wavefunctions" -decor none {
		table wfs {
		    -caption  "Enter Wavefunctions to be Pseudized:"
		    -head     {Label N L Occupancy Energy Rcut "US Rcut"}
		    -validate {whatever int int fortranreal fortranreal fortranreal fortranreal}
		    -cols     7
		    -rows     1
		    -outfmt   {" %2s " %1d %1d %8.3f %8.3f %6.2f %6.2f}
		}		
	    }
	}
    }
    page testing -name "PseudoPotential Test" {
	#
	#
	#
	namelist test -name "TEST" {
	    var nconf {
		-label "Number of testing configurations:" 
		-widget spinint 
		-validate posint
		-default 1
	    }
	    var file_pseudo {
		-widget   entryfileselectquote
		-label    "Name of the file containing the input PP:"
		-fmt      %S
	    }
	    dimension configts {
		-label     "Test electronic configurations"
		-start     1
		-end       1
	    }
	}
    }

    # ----------------------------------------------------------------------
    # take care of specialties
    # ----------------------------------------------------------------------
    source atomic-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source atomic-help.tcl
}
