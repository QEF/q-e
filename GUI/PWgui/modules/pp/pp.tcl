source commands.tcl

module PP\#auto -title "PWSCF GUI: module PP.x" -script {

    readfilter  ::pwscf::ppReadFilter

    namelist inputpp -name "INPUTPP" {

	var prefix \
	    -label    "Prefix of punch file saved by program PW.X (prefix):" \
	    -widget [list entrybutton "Prefix ..." [list ::pwscf::selectFileRoot $this prefix]] \
	    -fmt      %S

	var outdir {
	    -label    "Temporary directory where PW.X files resides (outdir):"
	    -widget   entrydirselectquote
	    -fmt      %S
	}
	var filplot {
	    -label    "Output file that will contain the calculated quantity (filplot):"
	    -fmt      %S
	}
	var plot_num {
	    -label    "What to calculate (plot_num):"
	    -widget   radiobox
	    -textvalue {
		"charge density"
		"total potential (= V_bare + V_H + V_xc)"
		"local ionic potential"
		"local density of states at E_fermi" 
		"local density of electronic entropy"
		"STM images"
		"spin polarization (= rho(up) - rho(down))"
		"|psi|^2"
		"electron localization function (ELF)"
		"planar average of all |psi|^2"
		"integrated local density of states (ILDOS)"
		"the V_bare + V_H potential"
		"the electric field potential"
		"the noncolinear magnetization"
	    }
	    -value { 0 1 2 3 4 5 6 7 8 9 10 11 12 13 }
	    -fmt %d
	}
	var spin_component {
	    -label    "Charge/potential/magnetization spin component (spin_component):"
	    -widget   optionmenu
	    -textvalue {
		"total charge/potential"
		"spin up charge/potential"
		"spin down charge/potential"
		"absolute value"
		"x component of the magnetization"
		"y component of the magnetization"
		"z component of the magnetization"
	    }
	    -value { 0 1 2  0 1 2 3 }
	}	

	group stm -name "STM" {
	    var sample_bias {
		-label    "For STM: the bias of the sample [in Ryd] in STM images (sample_bias):"
		-validate fortranreal
	    }
	    var stm_wfc_matching {
		-label     "For STM: wave-function matching (stm_wfc_matching):"
		-widget    radiobox
		-textvalue {Yes No}
		-value     {.true. .false.}
	    }
	    var z {
		-label    "For STM: height of matching [in celldm(3) units] (z):"
		-validate fortranreal
	    }
	    var dz {
		-label    "For STM: distance of next STM image calculation (dz):"
		-validate fortranreal
	    }
	}

	group psi2 -name "Psi2" {
	    var kpoint {
		-label    "For |psi^2|: which k-point (kpoint):"
		-widget    spinint
		-validate  posint
		-fmt       %d
	    }	
	    var kband {
		-label    "For |psi^2|: which band (kband):"
		-widget    spinint
		-validate  posint
		-fmt       %d
	    }
	    var lsign {
		-label    "For |psi^2| & Gamma: save the sign(psi) (lsign):"
		-widget    radiobox
		-textvalue {Yes No}
		-value     {.true. .false.}
	    }
	}

	group ildos -name "ILDOS" {
	    var emin {
		-label    "For ILDOS: miminum energy [in eV] (emin):"
		-validate  fortranreal
	    }	
	    var emax {
		-label    "FOR ILDOS: maximum energy [in eV] (emax):"
		-validate  fortranreal
	    }	
	}
    }

    # ----------------------------------------------------------------------
    # take care of specialties
    # ----------------------------------------------------------------------
    source pp-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source pp-help.tcl
}
