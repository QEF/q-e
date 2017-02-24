tracevar plot_num w {

    set spin_component_text [vartextvalue spin_component(1)]

    ::tclu::DEBUG "Plot_Num ..."

    switch -exact -- [vartextvalue plot_num] {
	"charge density" -
	"total potential (= V_bare + V_H + V_xc)" -
	"electrostatic potential (= V_bare + V_H)" -
	"all-electron valence charge density (for PAW)" -
	"all-electron charge density (valence+core, for PAW)" {
	    widget spin_component enable
	    widget spin_component(2) disable
	    widgetconfigure spin_component(1) -textvalues {
		"total charge/potential"
		"spin up charge/potential"
		"spin down charge/potential"		
	    }
	    groupwidget stm   disable 
	    groupwidget psi2  disable 
	    groupwidget ildos disable
	    groupwidget ldos  disable

	    if { ! [regexp {charge/potential} $spin_component_text] } {
		varset spin_component(1) -value {}
	    }		
	}

	"local density of states at specific energies (LDOS)" {
	    widget spin_component enable
	    widget spin_component(2) disable
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
	    }
	    groupwidget stm   disable 
	    groupwidget psi2  enable  
	    groupwidget ildos disable
	    groupwidget ldos  disable

	    if { [regexp {charge/potential} $spin_component_text] } {
		varset spin_component(1) -textvalue "charge"
	    }
	    if { [regexp {charge/potential} [vartextvalue spin_component(2)]] } {
		varset spin_component(2) -textvalue "charge"
	    }
	}
	
	"integrated local density of states (ILDOS)" {
	    widget spin_component enable
	    widget spin_component(2) disable
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
	    groupwidget stm   disable 
	    groupwidget psi2  disable  
	    groupwidget ildos disable
	    groupwidget ldos  disable

	    if { [regexp {charge} $spin_component_text] } {
		varset spin_component(1) -textvalue "absolute value"
	    }
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
		varset output_format -textvalue "XCRYSDEN's XSF format (2D or 3D)"
	    }
	}
	
	"3D plot" {
	    
	    if { [string match "*1D*" [vartextvalue output_format]] || [string match "*2D*" [vartextvalue output_format]] } {
		varset output_format -textvalue "XCRYSDEN's XSF format (whole unit cell) (3D)"
	    }
	    
	    switch -exact -- [vartextvalue output_format] {
		"XCRYSDEN's XSF format (whole unit cell) (3D)" -
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
    switch -exact -- [vartextvalue output_format] {
	"XCRYSDEN's XSF format (whole unit cell)" -
	"Gaussian cube-file format" {
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
	default {
	    # induce a trace-action on iflag variable
	    set iflag [vartextvalue iflag]
	    varset iflag -textvalue $iflag
	}
    }
}

postprocess {
    varset iflag         -textvalue "3D plot"
    varset output_format -textvalue "XCRYSDEN's XSF format (whole unit cell)"
}

postprocess {
    varset plot_num -textvalue ""
}
