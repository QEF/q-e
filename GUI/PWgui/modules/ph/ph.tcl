source commands.tcl

set ::guib::settings(filename_only_tail) 1

module PH\#auto -title "PWSCF GUI: module PH.x" -script {

    readfilter  ::pwscf::phReadFilter
    #writefilter ::pwscf::phWriteFilter

    line job_title -name "Job Title" {
	var title_ph {
	    -label    "Job title:"
	    -fmt      %S
	}
    }
    namelist inputpp -name "INPUTPH" {
	required {
	    auxilvar ntyp {
		-label   "Number of types of atoms in the unit cell (ntyp):"
		-validate posint
		-fmt      %d
		-default  1
		-widget   spinint  
	    }		
	    
	    dimension amass {
		-label     "Atomic mass of each atomic type"
		-validate  fortranreal
		-start     1
		-end       1
	    }
	}
	
	optional {
	    var outdir {
	    	-label    "Temporary directory where punch file resides (outdir):"
	    	-widget   entrydirselectquote
	    	-fmt      %S
	    }
	    var prefix -label "Prefix of punch file saved by PW.X (prefix):" \
		-widget   [list entrybutton "Prefix ..." [list ::pwscf::phSelectPunchFile $this prefix]] \
		-fmt      %S

	    var iverbosity {
	    	-label     "Verbosity of output (iverbosity):"
	    	-textvalue {"short output" "verbose output"}
	    	-value     {0 1}
	    	-widget    optionmenu
	    }
	    var time_max {
	    	-label    "Maximum allowed CPU run-time [in seconds] (time_max):"
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
	    var alpha_mix1 {
	    	-variable alpha_mix(1) 
	    	-label    "Mixing factor for updating the SCF potential (alpha_mix(1)):"
	    	-validate fortranreal
	    }
	    var nmix_ph {
	    	-label    "Number of iterations used in mixing of potential (nmix_ph):"
	    	-widget   spinint
	    	-fmt      %d
	    }

	    separator -label "--- Files ---"

	    var fildyn {
	    	-label    "File containing the dynamical matrix (fildyn):" 
	    	-widget   entryfileselectquote
	    	-fmt      %S
	    }
	    var fildrho {
	    	-label    "File containing the charge density variations (fildrho):"
	    	-widget   entryfileselectquote
	    	-fmt      %S
	    }
	    var filelph {
	    	-label    "File containing the electron-phonon matrix elements (filelph):"
	    	-widget   entryfileselectquote
	    	-fmt      %S
	    }
	    var fildvscf {
	    	-label    "File containing the potential variation (fildvscf):"
	    	-widget   entryfileselectquote
	    	-fmt      %S
	    }

	    separator -label "--- What to compute ---"

	    var maxirr {
	    	-label    "Maximum number of irreducible representation (maxirr):"
	    	-widget   spinint
	    	-fmt      %d
	    }

	    var epsil {
	    	-label     "Compute the macroscopic dielectric constant (epsil):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
	    }
	    var trans {
	    	-label     "Compute phonons (trans):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
	    }
	    var elph {
	    	-label     "Compute electron-phonon lambda coefficients (elph):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
	    }
	    var zue {
	    	-label     "Computed effective charges from the phonon density responses (zue):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
	    }
	}
    }

    line xq_line -name "The phonon wavevector" {
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

    # ----------------------------------------------------------------------
    # take care of specialties
    # ----------------------------------------------------------------------
    source ph-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source ph-help.tcl
}
