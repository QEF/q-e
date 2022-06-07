tracevar plot_num w {

    switch -exact -- [vartextvalue plot_num] {
	"charge density" -
	"total potential (= V_bare + V_H + V_xc)" -
        "charge density minus superposition of atomic densities" -
	"all-electron valence charge density (for PAW)" -
	"all-electron charge density (valence+core, for PAW)" {
	    widget spin_component enable
	    widget spin_component(2) disable
	    widgetconfigure spin_component(1) -textvalues {
		"spin up + spin down"
		"spin up only"
		"spin down only"		
	    }
            switch -exact [varvalue spin_component(1)] {
                0 { varset spin_component(1) -textvalue "spin up + spin down" }
                1 { varset spin_component(1) -textvalue "spin up only" }
                2 { varset spin_component(1) -textvalue "spin down only" }
            }
	    groupwidget stm   disable 
	    groupwidget psi2  disable 
	    groupwidget ildos disable
	    groupwidget ldos  disable
	}

	"local density of states at specific energies (LDOS)" {
	    widget spin_component enable
	    widget spin_component(2) disable
	    widgetconfigure spin_component(1) -textvalues {
		"spin up + spin down"
		"spin up only"
		"spin down only"		
	    }            
            switch -exact [varvalue spin_component(1)] {
                0 { varset spin_component(1) -textvalue "spin up + spin down" }
                1 { varset spin_component(1) -textvalue "spin up only" }
                2 { varset spin_component(1) -textvalue "spin down only" }
            }
            groupwidget stm   disable 
	    groupwidget psi2  disable  
	    groupwidget ildos enable
	    groupwidget ldos  enable
	}
	
	"STM images" {
	    widget spin_component disable
	    groupwidget stm   enable  
	    groupwidget psi2  disable 
	    groupwidget ildos disable
	    groupwidget ldos  disable
	}

	"|psi|^2" {
	    widget spin_component disable
	    groupwidget stm   disable 
	    groupwidget psi2  enable  
	    groupwidget ildos disable
	    groupwidget ldos  disable
	}

	"|psi|^2 (noncollinear case)" {
	    widget spin_component enable
	    widget spin_component(2) enable
	    foreach i {1 2} {
                widgetconfigure spin_component($i) -textvalues {
		    "charge"
		    "x component of the magnetization"
		    "y component of the magnetization"
		    "z component of the magnetization"
		}
                switch -exact [varvalue spin_component($i)] {
                    0 { varset spin_component($i) -textvalue "charge" }
                    1 { varset spin_component($i) -textvalue "x component of the magnetization" }
                    2 { varset spin_component($i) -textvalue "y component of the magnetization" }
                    3 { varset spin_component($i) -textvalue "z component of the magnetization" }
                }
            }
	    groupwidget stm   disable 
	    groupwidget psi2  enable  
	    groupwidget ildos disable
	    groupwidget ldos  disable
	}
	
	"integrated local density of states (ILDOS)" {
	    widget spin_component enable
	    widget spin_component(2) disable
	    widgetconfigure spin_component(1) -textvalues {
		"spin up + spin down"
		"spin up only"
		"spin down only"		
	    }
            switch -exact [varvalue spin_component(1)] {
                0 { varset spin_component(1) -textvalue "spin up + spin down" }
                1 { varset spin_component(1) -textvalue "spin up only" }
                2 { varset spin_component(1) -textvalue "spin down only" }
            }
	    groupwidget stm   disable 
	    groupwidget psi2  disable  
	    groupwidget ildos enable
	    groupwidget ldos  disable
	}

	"noncolinear magnetization" {
	    widget spin_component enable
	    widget spin_component(2) disable
	    widgetconfigure spin_component(1) -textvalues {
		"absolute value"
		"x component of the magnetization"
		"y component of the magnetization"
		"z component of the magnetization"
	    }	
            switch -exact [varvalue spin_component(1)] {
                0 { varset spin_component(1) -textvalue "absolute value" }
                1 { varset spin_component(1) -textvalue "x component of the magnetization" }
                2 { varset spin_component(1) -textvalue "y component of the magnetization" }
                2 { varset spin_component(1) -textvalue "z component of the magnetization" }
            }
	    groupwidget stm   disable 
	    groupwidget psi2  disable  
	    groupwidget ildos disable
	    groupwidget ldos  disable
	}

	default {
	    widget spin_component disable 
	    groupwidget stm   disable 
	    groupwidget ildos disable
	    groupwidget psi2  disable
	    groupwidget ldos  disable
	}	    
    }
}

tracevar nfile w {
    widgetconfigure filepp -end [varvalue nfile]
    widgetconfigure weight -end [varvalue nfile]
}

tracevar iflag w {
    switch -exact -- [vartextvalue iflag] {
	"1D plot, spherical average" { 
	    widget x0 enable
	    widget e1 enable
	    widget nx enable

	    widget e2 disable
	    widget e3 disable
	    widget ny disable
	    widget nz disable
	    widget radius disable
	    if { [string match "*2D*" [vartextvalue output_format]] || [string match "*3D*" [vartextvalue output_format]] } {
		varset output_format -textvalue "format suitable for gnuplot (1D)"
	    }
	}

	"1D plot" { 
	    widget x0 enable
	    widget e1 enable
	    widget nx enable

	    widget e2 disable
	    widget e3 disable
	    widget ny disable
	    widget nz disable
	    widget radius disable
	    if { [string match "*2D*" [vartextvalue output_format]] || [string match "*3D*" [vartextvalue output_format]] } {
		varset output_format -textvalue "format suitable for gnuplot (1D)"
	    }
	}
	
	"2D plot" {
	    widget x0 enable
	    widget e1 enable
	    widget e2 enable
	    widget nx enable
	    widget ny enable

	    widget e3 disable
	    widget nz disable
	    widget radius disable
	    if { [string match "*3D*" [vartextvalue output_format]] || [string match "*1D*" [vartextvalue output_format]] } {
		varset output_format -textvalue "XCRYSDEN's XSF format - slow (2D or 3D)"
	    }
	}
	
	"3D plot" {
            set fmt [vartextvalue output_format]

	    if { ([string match "*1D*" $fmt] || [string match "*2D*" $fmt]) && ![string match "*XSF*" $fmt] } {
		varset output_format -textvalue "XCRYSDEN's XSF format - fast (whole unit cell, 3D)"
	    }
	    
	    switch -exact -- [vartextvalue output_format] {
		"XCRYSDEN's XSF format - fast (whole unit cell, 3D)" -
		"Gaussian cube-file format (3D)" {
		    widget e1 disable
		    widget e2 disable
		    widget e3 disable
		    widget x0 disable
		    widget nx disable
		    widget ny disable
		    widget nz disable
		    widget radius disable
		}
		default {	        
		    widget e1 enable
		    widget e2 enable
		    widget e3 enable
		    widget x0 enable
		    widget nx enable
		    widget ny enable
		    widget nz enable
		    widget radius disable
		}
	    }
	}
	
	"2D polar plot" {
	    widget e1 enable
	    widget e2 enable
	    widget e3 disable
	    widget x0 enable
	    widget nx enable
	    widget ny enable
	    widget nz disable
	    widget radius enable

	    if { [string match "*1D*" [vartextvalue output_format]] || [string match "*3D*" [vartextvalue output_format]] } {
		varset output_format -textvalue ""
	    }	    
	}
    }
}


tracevar output_format w {
    # also take output_format into account
    set iflag [varvalue iflag]
    switch -exact -- [vartextvalue output_format] {
	"XCRYSDEN's XSF format - fast (whole unit cell, 3D)" -
	"Gaussian cube-file format (3D)" {
	    widget e1 disable
	    widget e2 disable
	    widget e3 disable
	    widget x0 disable
	    widget nx disable
	    widget ny disable
	    widget nz disable
	    widget radius disable
	    varset iflag -textvalue "3D plot"
	}
        "XCRYSDEN's XSF format - slow (2D or 3D)"  {
            widget e1 enable
            widget e2 enable
            widget x0 enable
            widget nx enable
            widget ny enable
            widget radius disable
            if { [string is integer $iflag] && $iflag > 2 } {
                widget e3 enable
                widget nz enable
            }
            if { $iflag != 2 && $iflag != 3 } {
                varset iflag -textvalue "3D plot"
            }
        }
        "format suitable for gnuplot (1D)" {
            if { $iflag != 0 && $iflag != 1 } {
                varset iflag -textvalue "1D plot"
            }
        }
        "format suitable for plotrho (2D)" -
        "format suitable for gnuplot (2D)" {
            if { $iflag != 2 && $iflag != 4 } {
                varset iflag -textvalue "2D plot"
            }            
        }
	default {
	    # induce a trace-action on iflag variable
	    varset iflag -textvalue [vartextvalue iflag]
	}
    }
}

tracevar spin_component(2) w {
    varset plot_num -textvalue "|psi|^2 (noncollinear case)"
}

postprocess {
    varset iflag         -textvalue "3D plot"
    varset output_format -textvalue "XCRYSDEN's XSF format - fast (whole unit cell, 3D)"
    varset plot_num      -textvalue ""
}
