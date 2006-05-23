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
	    var outdir {
	    	-label    "Temporary directory where punch file resides (outdir):"
	    	-widget   entrydirselectquote
	    	-fmt      %S
	    }
	    var recover {
		-label "Restart from an interrupted run (recover):"
		-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    }
	    var prefix -label "Prefix of data file saved by PW.X (prefix):" \
		-widget   [list entrybutton "Prefix ..." [list ::pwscf::phSelectPunchFile $this prefix]] \
		-fmt      %S

	    separator -label "--- What to Calculate ---"

	    var ldisp {
	    	-label     "Compute phonon dispersions (ldisp):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
	    }
	    group nq -name q-grid {
		var nq1 {
		    -label     "NQ\#1:"
		    -validate posint
		}
		var nq2 {
		    -label     "NQ\#2:"
		    -validate posint
		}
		var nq3 {
		    -label     "NQ\#3:"
		    -validate posint
		}
	    }

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

	}
	
	optional {


	    separator -label "--- More options on what to compute ---"

	    var elph {
	    	-label     "Compute electron-phonon lambda coefficients (elph):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
		-default   .false.
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
	    var elop {
	    	-label     "Compute electro-optic coefficients (elop):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
	    }

	    var fpol {
		-label     "Compute dynamic polarizabilities (fpol):"
		-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    }

	    var zue {
	    	-label     "Compute effective charges from the phonon density responses (zue):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
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

	    var maxirr {
	    	-label    "Maximum number of irreducible representation (maxirr):"
	    	-widget   spinint
	    	-fmt      %d
	    }

	    var modenum {
		-label    "Mode number for single-mode calculation (modenum):"
		-validate integer
		-widget   spinint
	    }
	    
	    var nrapp {
		-label    "Number of representations to do (nrapp):"
		-validate nonnegint
		-widget   spinint
	    }

	    var nat_todo {
		-label "Number of atom to be displaced (nat_todo):"
		-validate nonnegint
		-widget spinint
	    }
	    
	    separator -label "--- Misc control options ---"

	    var lnscf {
	    	-label    "Perform the non-scf calculation first (lnscf):"
	    	-textvalue {Yes No}
	    	-value     {.true. .false.}
	    	-widget    radiobox
	    	-fmt       %s
	    }
	    var iverbosity {
	    	-label     "Verbosity of output (iverbosity):"
	    	-textvalue {"short output" "verbose output"}
	    	-value     {0 1}
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

	    separator -label "--- Output Data Files ---"

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

    line representation_line -name "List of the representation to do:" {
	var representation_list {
	    -label "Indices of representations (comma or whitespace separated):"
	}
    }

    line atom_disp_line -name "List of atoms to displace:" {
	var atom_disp_list {
	    -label "Indices of atoms (comma or whitespace separated):"
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
