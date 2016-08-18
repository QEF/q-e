source commands.tcl

module ProjWfc\#auto -title "PWSCF GUI: module ProjWfc.x" -script {

    readfilter ::pwscf::projwfcReadFilter

    namelist projwfc -name "PROJWFC" {
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
	    
	    var lsym {
		-label "Symmetrize projections (lsym):"
		-value { .true. .false. }
		-textvalue { Yes No }
		-widget radiobox
	    }	

	    var kresolveddos {
		-label "Compute k-resolved DOS (kresolveddos):"
		-value { .true. .false. }
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

	    separator -label "--- PROJECTIONS options ---"
	    
	    var lwrite_overlaps {
		-label "Print the overlap matrix of atomic orbitals (lwrite_overlaps):"
		-value { .true. .false. }
		-textvalue { Yes No }
		-widget radiobox
	    }
	    
	    var lbinary_data {
		-label "Write atomic_proj datafile in binary format (lbinary_data):"
		-value { .true. .false. }
		-textvalue { Yes No }
		-widget radiobox
	    }	    

	    separator -label "--- PAW option ---"

	    var pawproj {
		-label "use PAW projectors and all-electron PAW basis (pawproj):"
		-value { .true. .false. }
		-textvalue { Yes No }
		-widget radiobox
	    }
	    separator -label "--- Local DOS options ---"
	    
	    var tdosinboxes {
		-label "Compute the local DOS computed in volumes (tdosinboxes):"
		-value { .true. .false. }
		-textvalue { Yes No }
		-widget radiobox
	    }

	    group local_dos -decor normal {
		var n_proj_boxes {
		    -label   "Number of boxes where the local DOS is computed (n_proj_boxes):"
		    -widget   spinint
		    -validate posint
		    -default  1
		}
		
		var plotboxes {
		    -label "Write the boxes into XSF 3D datagrid file (plotboxes):"
		    -value { .true. .false. }
		    -textvalue { Yes No }
		    -widget radiobox
		}
		
		table irmin {
		    -caption  "First point to be included in the given box (irmin):"
		    -head     {(1,*) (2,*) (3,*)}
		    -validate {int int int}
		    -cols     3
		    -rows     1
		}	    
		
		table irmax {
		    -caption  "Last point to be included in the given box (irmax):"
		    -head     {(1,*) (2,*) (3,*)}
		    -validate {int int int}
		    -cols     3
		    -rows     1
		}	    
	    }
	}
    }

    # ----------------------------------------------------------------------
    # take care of specialties
    # ----------------------------------------------------------------------
    source projwfc-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source projwfc-help.tcl
}
