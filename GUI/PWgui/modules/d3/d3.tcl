source commands.tcl

module d3\#auto -title "PWSCF GUI: module D3.x" -script {

    readfilter ::pwscf::d3ReadFilter

    namelist inputph -name "INPUTPH" {
	required {
	    var fildrho {
		-label "File containing the variation of the charge at q (fildrho):"
		-widget   entryfileselectquote
		-validate string
	    }
	    var fild0rho {
		-label "File containing the variation of the charge at q=0 (fild0rho):"
		-widget   entryfileselectquote
		-validate string
	    }
	    auxilvar ntyp {
		-label   "Number of types of atoms:"
		-validate posint
		-fmt      %d
		-default  1
		-widget   spinint  
	    }		
	    dimension amass {
		-label     "Atomic mass of each atomic type:"
		-validate  fortranreal
		-start     1
		-end       1
	    }
	}
	optional {
	    var prefix \
		-label  "Prefix for file names (prefix):" \
		-widget [list entrybutton "Prefix ..." [list ::pwscf::selectFileRoot $this prefix]] \
		-fmt    %S -validate string
	    
	    var outdir {
		-label  "Temporary directory (outdir):"
		-widget entrydirselectquote
		-validate string
	    }
	    var iverbosity {
		-label "Verbosity of output (iverbosity):"
		-textvalue {high low}
		-value     {1 0}
		-widget    radiobox
	    }
	    var fildyn -label "Output file with the derivative of the dynamical matrix (fildyn):" -validate string
	    var ethr_ph {
		-label "Threshold for iterative diagonalization (ethr_ph):" 
		-validate fortranposreal
	    }
	    auxilvar nmode {
		-label    "Number of q=0 modes to compute:"
		-validate posint
		-default  1
		-widget   spinint
	    }
	    dimension q0mode_todo {
		-label    "The q=0 modes to compute"
		-validate nonnegint
		-widget   spinint
		-start     1
		-end       1
	    }
	    var wraux {
		-label "Write different terms of the matrix on different files (wraux):"
		-textvalue {Yes No}
		-value     {.true. .false.}
		-widget     radiobox
	    }
	    var recv {
		-label "Is this recover run (rcev):"
		-textvalue {Yes No}
		-value     {.true. .false.}
		-widget     radiobox
	    }
	    var istop  {
		-label    "Where to stop calculation, istop=0 means do not stop (istop):"
		-widget   spinint
		-validate nonnegint
	    }
	}
    }

    # ----------------------------------------------------------------------
    # take care of specialities
    # ----------------------------------------------------------------------
    source d3-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source d3-help.tcl
}
