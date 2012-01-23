source commands.tcl

module Dos\#auto -title "PWSCF GUI: module Dos.x" -script {

    readfilter ::pwscf::dosReadFilter

    namelist dos -name "DOS" {
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

	    var fildos {
		-label "Prefix for output files containing DOS(E) (fildos):"
		-validate string
	    }
	    
	    separator -label "--- DOS ploting options ---"

	    var ngauss {
		-label   "Type of gaussian broadening (ngauss):"
		-widget  radiobox
		-value   {0 1 -1 99}
		-textvalue {
		    "Simple Gaussian (default)"
		    "Methfessel-Paxton of order 1"
		    "Marzari-Vanderbilt \"cold smearing\""
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
    source dos-help.tcl
}
