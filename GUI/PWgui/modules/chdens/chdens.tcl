# --------------------------------------------------
# enable DEBUGGING::
#set ::tclu::DEBUG 1
# --------------------------------------------------


source commands.tcl
#set ::guib::settings(filename_only_tail) 1

module ChDens\#auto -title "PWSCF GUI: module ChDens.x" -script {

    readfilter ::pwscf::chdensReadFilter

    namelist input -name "INPUT" {

	var nfile {
	    -label    "Number of data files (nfile):"
	    -widget   spinint
	    -validate posint
	    -default  1
	}

	dimension filepp {
	    -label    "Name of the data file"
	    -start    1
	    -end      1
	    -widget   entryfileselectquote
	}

	dimension weight {
	    -label    "Weight of the charge"
	    -start    1
	    -end      1
	    -widget   entry
	    -validate fortranreal
	    -default  1.0
	}

	separator -label "--- Plot info ---"

	var fileout -label "Name of output file (fileout):"

	var iflag {
	    -label     "Dimensionality of plot (iflag):"
	    -textvalue {
		"1D plot"
		"2D plot"
		"3D plot"
		"2D polar plot"
	    }
	    -value  { 1 2 3 4 }
	    -widget optionmenu
	}

	var plot_out {
	    -label    "What to plot (plot_out):"
	    -textvalue {
		"normal plot"
		"spherical averaged plot"
		"induced polarization along x"
		"induced polarization along y"
		"induced polarization along z"
	    }
	    -value     { 1 0 2 3 4 }
	    -widget    optionmenu
	}

	var output_format {
	    -label     "Format of the output (output_format):"
	    -textvalue {
		"XCRYSDEN's XSF format"
		"XCRYSDEN's XSF format (whole unit cell)"
		"format suitable for gnuplot"
		"format suitable for contour.x"
		"format suitable for plotrho"
		"format suitable for gOpenMol"
		"Gaussian cube-file format"
	    }
	    -value     { 3 5 0 1 2 4 6 }
	    -widget    optionmenu
	}

	separator -label "--- Polarization input ---"

	var epsilon {
	    -label    "Dielectric constant for polarization calculation (epsilon):"
	    -validate fortranreal
	}

	var filepol -label "Name of output file for induced polarization (filepol):"
	
	separator -label "--- Spanning vectors & origin ---"

	dimension e1 {
	    -label    "1st spanning vector:"
	    -validate fortranreal
	    -start    1
	    -end      3
	    -pack     left
	}

	dimension e2 {
	    -label    "2nd spanning vector"
	    -validate fortranreal
	    -start    1
	    -end      3
	    -pack     left
	}

	dimension e3 {
	    -label    "3rd spanning vector"
	    -validate fortranreal
	    -start    1
	    -end      3
	    -pack     left
	}

	dimension x0 {
	    -label    "Origin of the plot"
	    -validate fortranreal
	    -start    1
	    -end      3
	    -pack     left
	}


	separator -label "--- Number of points in each direction ---" 

	group nxnynz -name nxnynz {
	    packwidgets left
	    var nx -label "nx:" -validate posint -widget spinint
	    var ny -label "ny:" -validate posint -widget spinint
	    var nz -label "nz:" -validate posint -widget spinint
	}

	separator -label "--- Polar plot ---"
	var radius -label "Radius of the sphere (radius):" -validate real

	separator -label "--- Ionic and electronic dipole moment ---"
	var idpol {
	    -label "Which dipole moment to compute (idpol):" 
	    -textvalue {
		"ionic + electronic dipole moment"
		"electronic-only dipole moment"
	    }
	    -value {1 2} 
	    -widget optionmenu
	}

	separator -label "--- Makov-Payne (MP) correction for charged supercells ---"
	var makov {
	    -label "Compute the 1st and 2d order MP corrections (makov):"
	    -textvalue { Yes No }
	    -value     { .true. .false. }
	    -widget    radiobox
	}
    }
    
    # ----------------------------------------------------------------------
    # take care of specialties
    # ----------------------------------------------------------------------
    source chdens-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source chdens-help.tcl
}
