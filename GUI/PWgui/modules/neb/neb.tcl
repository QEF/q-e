source commands.tcl

module NEB -title "PWSCF GUI: module NEB.x" -script {

    readfilter ::pwscf::nebReadFilter

    #
    # PATH
    # 

    namelist inputpp -name "PATH" {

	var string_method {
	    -label     "Type of path calculation (string_method):"
	    -widget    radiobox
	    -textvalue {
		"Nudged Elastic Band  <neb>"
		"String Method Dynamics  <smd>"
	    }
	    -value {
		'neb'
		'smd'
	    }
	    -default "Nudged Elastic Band  <neb>"
	}
	
	var restart_mode {
	    -label    "Restart mode (restart_mode):"
	    -widget   radiobox
	    -textvalue {
		"from scratch  <from_scratch>"
		"from previous interrupted run  <restart>"
	    }
	    -value {
		'from_scratch'
		'restart'
	    }
	}
	
	var nstep_path {
	    -label    "Number of ionic steps (nstep_path):"
	    -widget   spinint
	    -validate posint
	}


	var num_of_images {
	    -label   "Number of images used to discretize the path (num_of_images):"
	    -widget   spinint
	    -validate posint
	}
	
	var opt_scheme {
	    -label "Type of optimization scheme (opt_scheme):"
	    -value {
		'quick-min' 
		'broyden'
		'broyden2'
		'sd'
		'langevin'
	    }
	    -textvalue {
		"optimization algorithm based on molecular dynamics  <quick-min>"
		"Broyden method  <broyden>"
		"Alternate Broyden method  <broyden2>"
		"steepest descent  <sd>"
		"finite temperature langevin dynamics  <langevin>"
	    }
	    -widget optionmenu
	}
	
	var CI_scheme {
	    -label "Type of climbing image scheme (CI_scheme):"
	    -textvalue {
		"do not use climbing image  <no-CI>"
		"image highest in energy is allowed to climb  <auto>"
		"climbing images are manually selected  <manual>"
	    }
	    -value {
		'no-CI'
		'auto'
		'manual'
	    }
	    -widget optionmenu
	}

	var first_last_opt {
	    -label "Optimize also the first and the last configurations (first_last_opt):"
	    -textvalue { Yes No }
	    -value     { .TRUE. .FALSE. }
	    -widget    radiobox
	}

	var minimum_image {
	    -label     "Use a minimum image criterion (minimum_image):"
	    -textvalue { Yes No }
	    -value     { .TRUE. .FALSE. }
	    -widget    radiobox
	}

	var temp_req {
	    -label    "Temperature used for langevin dynamics of the string (temp_req):"
	    -validate fortranposreal
	}

	var ds {
	    -label    "Optimization step length (ds):"
	    -validate fortranposreal
	}
	
	var path_thr {
	    -label "Convergence threshold for path optimization (path_thr):"
	    -validate fortranposreal
	}
	
	var use_freezing {
	    -label "Only the images with larger errors are optimised (use_freezing):"
	    -textvalue { Yes No }
	    -value     { .TRUE. .FALSE. }
	    -widget    radiobox
	}

	var use_masses {
	    -label "The optimisation is done with mass-weighted coordinates (use_masses):"
	    -textvalue { Yes No }
	    -value     { .TRUE. .FALSE. }
	    -widget    radiobox
	}	    
	
	group elastic_constants -name "Elastic Constants for NEB spring:" -decor normal {
	    packwidgets left
	    var k_max -label "k_max:" -validate fortranposreal			
	    var k_min -label "k_min:" -validate fortranposreal			
	}

	group fcp -decor normal {
	    separator -label "--- Fictitious Charge Particle (FCP) options ---"

	    var lfcpopt {
		-label "Perform a constant bias potential calculation using ESM (lfcpopt):"
		-textvalue { Yes No }
		-value     { .TRUE. .FALSE. }
		-widget    radiobox
	    }

	    group fcp_specs -decor none {	
		var fcp_mu {
		    -label "target Fermi energy in Ry (fcp_mu):"	    
		    -validate fortranreal
		}
		
		var fcp_tot_charge_first {
		    -label "Total charge for the first image (fcp_tot_charge_first):"
		    -validate fortranreal
		}    
		
		var fcp_tot_charge_last {
		    -label "Total charge for the last image (fcp_tot_charge_last):"
		    -validate fortranreal
		}
	    }
	}
    }

    # 
    # CLIMBING_IMAGES
    # 
    group climbing_images -name "Card: CLIMBING_IMAGES" -decor normal {
	keyword climbing_images_key CLIMBING_IMAGES\n
	line climbing_images_line -decor none {
	    var climbing_images_list {
		-label "List of climbing images, separated by a comma:"
		-infmt %S
	    }
	}
    }

    # ----------------------------------------------------------------------
    # take care of specialities
    # ----------------------------------------------------------------------
    source neb-event.tcl

    # ------------------------------------------------------------------------
    # source the HELP file
    # ------------------------------------------------------------------------
    source neb-help.tcl
}
