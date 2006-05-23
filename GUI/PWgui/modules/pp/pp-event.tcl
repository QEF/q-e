tracevar plot_num w {

    ::tclu::DEBUG "Plot_Num ..."

    switch -exact -- [vartextvalue plot_num] {
	"charge density" -
	"total potential (= V_bare + V_H + V_xc)" -
	"the V_bare + V_H potential" {
	    widget spin_component enable
	    widgetconfigure spin_component -textvalues {
		"total charge/potential"
		"spin up charge/potential"
		"spin down charge/potential"		
	    }
	    groupwidget stm   disable 
	    groupwidget psi2  disable 
	    groupwidget ildos disable
	}

	"STM images" {
	    widget spin_component disable
	    groupwidget stm   enable  
	    groupwidget psi2  disable 
	    groupwidget ildos disable
	}

	"|psi|^2" {
	    widget spin_component disable 
	    groupwidget stm   disable 
	    groupwidget psi2  enable  
	    groupwidget ildos disable
	}

	"|psi|^2 (noncollinear case)" {
	    widget spin_component enable
	    widgetconfigure spin_component -textvalues {
		"absolute value"
		"x component of the magnetization"
		"y component of the magnetization"
		"z component of the magnetization"
	    }	
	    groupwidget stm   disable 
	    groupwidget psi2  enable  
	    groupwidget ildos disable
	}
	
	"integrated local density of states (ILDOS)" {
	    widget spin_component disable 
	    groupwidget stm   disable 
	    groupwidget psi2  disable  
	    groupwidget ildos enable
	}
	"the noncolinear magnetization" {
	    widget spin_component enable
	    widgetconfigure spin_component -textvalues {
		"absolute value"
		"x component of the magnetization"
		"y component of the magnetization"
		"z component of the magnetization"
	    }	
	    groupwidget stm   disable 
	    groupwidget psi2  disable  
	    groupwidget ildos disable
	    
	}
	default {
	    widget spin_component disable 
	    groupwidget stm   disable 
	    groupwidget ildos disable
	    groupwidget psi2  disable 
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
	    widget e1 enable
	    widget e2 disable
	    widget e3 disable
	    widget x0 enable
	    widget nx enable
	    widget ny disable
	    widget nz disable
	    widget radius disable
	    if { [string match "XCRYSDEN*" [vartextvalue output_format]] } {
		varset output_format -textvalue ""
	    }
	}

	"1D plot" { 
	    widget e1 enable
	    widget e2 disable
	    widget e3 disable
	    widget x0 enable
	    widget nx enable
	    widget ny disable
	    widget nz disable
	    widget radius disable
	    if { [string match "XCRYSDEN*" [vartextvalue output_format]] } {
		varset output_format -textvalue ""
	    }
	}
	
	"2D plot" {
	    widget e1 enable
	    widget e2 enable
	    widget e3 disable
	    widget x0 enable
	    widget nx enable
	    widget ny enable
	    widget nz disable
	    widget radius disable
	    if { [string match "XCRYSDEN's XSF format (whole unit cell)" [vartextvalue output_format]] } {
		varset output_format -textvalue "XCRYSDEN's XSF format"
	    }
	}
	
	"3D plot" {
	    switch -exact -- [vartextvalue output_format] {
		"XCRYSDEN's XSF format (whole unit cell)" {
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
	    if { [string match "XCRYSDEN*" [vartextvalue output_format]] } {
		varset output_format -textvalue ""
	    }	    
	}
    }
}


tracevar output_format w {
    # also take output_format into account
    switch -exact -- [vartextvalue output_format] {
	"XCRYSDEN's XSF format (whole unit cell)" {
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