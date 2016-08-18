source commands.tcl

module PP\#auto -title "PWSCF GUI: module PP.x" -script {

    readfilter  ::pwscf::ppReadFilter

    #
    # Namelist: INPUTPP
    #

    page extract -name "Specify property to calculate" {
	namelist inputpp -name "INPUTPP" {
	    
	    var prefix \
		-label    "Prefix of punch file saved by program PW.X (prefix):" \
		-widget   [list entrybutton "Prefix ..." [list ::pwscf::selectFileRoot $this prefix]] \
		-fmt      %S -validate string \
		-validate string

	    var outdir {
		-label    "Temporary directory where PW.X files resides (outdir):"
		-widget   entrydirselectquote
		-fmt      %S -validate string
		-validate string
	    }
	    var filplot {
		-label    "Output file that will contain the calculated quantity (filplot):"
		-fmt      %S -validate string
		-validate string
	    }
	    var plot_num {
		-label    "What to calculate (plot_num):"
		-widget   radiobox
		-textvalue {
		    "charge density"
		    "total potential (= V_bare + V_H + V_xc)"
		    "local ionic potential (= V_bare)"
		    "local density of states at E_fermi" 
		    "local density of electronic entropy"
		    "STM images"
		    "spin polarization (= rho(up) - rho(down))"
		    "|psi|^2"
		    "|psi|^2 (noncollinear case)"
		    "electron localization function (ELF)"
		    "charge density minus superposition of atomic densities"
		    "integrated local density of states (ILDOS)"
		    "electrostatic potential (= V_bare + V_H)"
		    "sawtooth electric field potential"
		    "noncolinear magnetization"
		    "all-electron valence charge density (for PAW)"
		    "exchange-correlation magnetic field (for noncollinear)"
		    "reduced density gradient"
		    "product of density and 2nd-eigenvalue of density Hessian matrix"
		    "all-electron charge density (valence+core, for PAW)"
		}
		-value { 0 1 2 3 4 5 6 7 7 8 9 10 11 12 13 17 18 19 20 21 }
		-fmt %d
	    }
	    var spin_component {
		-label    "Spin component (spin_component):"
		-widget   optionmenu
		-textvalue {
		    "total charge/potential"
		    "spin up charge/potential"
		    "spin down charge/potential"
		    "charge"
		    "absolute value"
		    "x component of the magnetization"
		    "y component of the magnetization"
		    "z component of the magnetization"
		}
		-value { 0 1 2  0  0 1 2 3 }
	    }	
	    
	    separator -label "--- Options for STM images ---"

	    group stm -name "STM" {
		var sample_bias {
		    -label    "For STM: the bias of the sample [in Ryd] in STM images (sample_bias):"
		    -validate fortranreal
		}
	    }
	    
	    separator -label "--- Options for |psi|^2 ---"

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
	    
	    separator -label "--- Options for ILDOS ---"

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
    }

    #
    # Namelist: PLOT
    #

    page chdens -name "Specify Plot " {

	namelist plot -name "PLOT" {
	    
	    var nfile {
		-label    "Number of data files (nfile):"
		-widget   spinint
		-validate posint
		-default  1
	    }
	    
	    dimension filepp {
		-label    "Filenames of data files:"
		-start    1
		-end      1
		-widget   entryfileselectquote
		-validate string
		-fmt      %S
	    }
	    
	    dimension weight {
		-label    "Weighting factors:"
		-start    1
		-end      1
		-widget   entry
		-validate fortranreal
		-default  1.0
	    }
	    
	    separator -label "--- Plot info ---"
	    
	    var fileout -label "Name of output file (fileout):" -validate string
	    
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
		    "XCRYSDEN's XSF format (2D or 3D)" 
		    "XCRYSDEN's XSF format (whole unit cell) (3D)"
		    "format suitable for gnuplot (1D)"
		    "format suitable for contour.x (2D)"
		    "format suitable for plotrho (2D)"
		    "format suitable for gOpenMol (3D)"
		    "Gaussian cube-file format (3D)"
		    "format suitable for gnuplot (2D)"
		}
		-value     { 3 5 0 1 2 4 6 7 }
		-widget    optionmenu
	    }
	    var interpolation { 
		-label     "Interpolation (interpolation):"
		-value     { 'fourier' 'bspline' }
		-textvalue { Fourier bspline }
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
