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
		"1D plot, spherical average"
		"1D plot"
		"2D plot"
		"3D plot"
		"2D polar plot"
	    }
	    -value  { 0 1 2 3 4 }
	    -widget optionmenu
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
