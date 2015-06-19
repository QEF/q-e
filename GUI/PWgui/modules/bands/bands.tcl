source commands.tcl

module Bands\#auto -title "PWSCF GUI: module Bands.x" -script {

    readfilter ::pwscf::bandsReadFilter

    namelist bands -name "BANDS" {

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
	    
	    var filband {
                -label "Output file containing the bands (filband):"
                -validate string
            }
	    
	    var spin_component {
		-label    "Spin component (spin_component):"
		-widget   radiobox
		-textvalue {
		    "spin up"
		    "spin down"
		}
		-value { 1 2 }
	    }
	    
	    group sigma {
		packwidgets left
		dimension lsigma { 
		    -start 1
		    -end   3
		    -label "Expectation values of the spin operator on the spinor wave-functions (lsigma)"
		    -textvalue { Yes No }
		    -value     { .true. .false. }
		    -widget    radiobox
		}
	    }
	    
	    var lp {
		-label "Write matrix elements of the momentum operator p to a file (lp):"
		-textvalue { Yes No }
		-value     { .true. .false. }
		-widget    radiobox		
	    }

	    var filp {
		-label "Output file containing matrix elements (filp):"
                -validate string	
		file name for matrix elements of p
	    }
	   
	    var lsym {
		-label "Classify bands according to irreducible representations of small group of k (lsym):"
		-textvalue { Yes No }
		-value     { .true. .false. }
		-widget    radiobox
	    }
	 
	    var no_overlap { 
		-label "Don't change the order of eigenvalues in output file (no_overlap:)"
		-textvalue { Yes No }
		-value     { .true. .false. }
		-widget    radiobox
	    }

	    var plot_2d { 
		-label "Print the eigenvalues in 2D gnuplot's format (plot_2d):"
		-textvalue { Yes No }
		-value     { .true. .false. }
		-widget    radiobox
	    }
	    
	    separator -label "--- Range of k-points for symmetry analysis ---"

	    group kpoints {
		packwidgets left

		var firstk -label "First k-point (firstk):" -widget spinint -validate posint
		var lastk  -label "Last k-point (firstk):" -widget spinint -validate posint
	    }
	}
    }

    # ----------------------------------------------------------------------
    # take care of specialities
    # ----------------------------------------------------------------------
    source bands-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source bands-help.tcl
}
