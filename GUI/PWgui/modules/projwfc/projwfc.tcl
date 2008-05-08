source commands.tcl

module ProjWfc\#auto -title "PWSCF GUI: module ProjWfc.x" -script {

    readfilter ::pwscf::projwfcReadFilter

    namelist inputpp -name "INPUTPP" {
	optional {
	    var prefix {
		-label    "Prefix of punch file saved by program PW.X (prefix):" 
		-fmt      %S -validate string
	    }
	    
	    var outdir {
		-label    "Temporary directory where PW.X files resides (outdir):"
		-widget   entrydirselectquote
		-fmt      %S -validate string
	    }

	    var filpdos {
		-label "Prefix for output files containing PDOS(E) (filpdos):"
		-validate string
	    }

	    var filproj {
		-label "File containing the projections (filproj):"
		-validate string
	    }
	    
	    separator -label "--- PDOS ploting options ---"

	    var ngauss {
		-label   "Type of gaussian broadening (ngauss):"
		-widget  optionmenu
		-value   {0 1 -1 99}
		-textvalue {
		    "Simple Gaussian (default)"
		    "Methfessel-Paxton of order 1"
		    " Marzari-Vanderbilt \"cold smearing\""
		    "Fermi-Dirac function"
		}
	    }

	    var degauss {
		-label     "Gaussian broadening \[in Ry\] (degauss):"
		-validate  fortranreal
	    }

	    var DeltaE {
		-label    "Resolution of PDOS plots \[in eV\] (DeltaE):"
		-validate fortranreal
		-default  0.01
	    }
	    
	    var lsym {
		-label "Symmetrize projections (lsym):"
		-value { 1 0 }
		-textvalue { Yes No }
		-widget radiobox
	    }	

	    separator -label "--- Energy window for PDOS ---"

	    var Emin {
		-label    "Minimum energy \[in eV\] (Emin):"
		-validate fortranreal
	    }

	    var Emax {
		-label    "Maximum energy \[in eV\] (Emin):"
		-validate fortranreal
	    }
	}
    }

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source projwfc-help.tcl
}
