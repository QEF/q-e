#source commands.tcl

module ProjWfc\#auto -title "PWSCF GUI: module ProjWfc.x" -script {
    namelist inputpp -name "INPUTPP" {
	optional {
	    var prefix {
		-label    "Prefix of punch file saved by program PW.X (prefix):" 
		-fmt      %S
	    }
	    
	    var outdir {
		-label    "Temporary directory where PW.X files resides (outdir):"
		-widget   entrydirselectquote
		-fmt      %S
	    }
	    
	    var io_choice {
		-label     "Type of the output (io_choice):"
		-widget    optionmenu
		-value     {'standard' 'files' 'both'}
		-textvalue {
		    "write output to stdout only"
		    "write output to files only"
		    "write output to stdout and files"
		}
	    }

	    separator -label "--- PDOS ploting options ---"

	    var DeltaE {
		-label    "Resolution of PDOS plots \[in eV\] (DeltaE):"
		-validate fortranreal
		-default  0.01
	    }
	    
	    var smoothing {
		-label    "Gaussian smoothing for PDOS \[in eV\] (smoothing):"
		-validate fortranreal
		-default  0.15
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
