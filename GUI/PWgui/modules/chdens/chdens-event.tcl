tracevar nfile w {
    widgetconfigure filepp -end [varvalue nfile]
    widgetconfigure weight -end [varvalue nfile]
}

tracevar iflag w {
    switch -exact -- [vartextvalue iflag] {
	"1D plot" { 
	    widget e1 enable
	    widget e2 disable
	    widget e3 disable
	    widget x0 enable
	    widget nx enable
	    widget ny disable
	    widget nz disable
	    widget radius disable
	    widget idpol disable
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
	    widget idpol disable
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
	    if { [varvalue plot_out] == 1 } {
		widget idpol enable
	    } else {
		widget idpol disable
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
	    widget idpol disable
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

tracevar plot_out w {
    switch -exact -- [vartextvalue plot_out] {
	"induced polarization along x" -
	"induced polarization along y" -
	"induced polarization along z" {
	    widget epsilon enable
	    widget filepol enable
	}
	"normal plot" {
	    if { [varvalue iflag] == 3 } {
		widget idpol enable
	    } else {
		widget idpol disable
	    }
	}
	default {
	    widget epsilon disable
	    widget filepol disable
	    widget idpol disable
	}
    }
}

postprocess {
    varset iflag         -textvalue "3D plot"
    varset plot_out      -textvalue "normal plot"
    varset output_format -textvalue "XCRYSDEN's XSF format (whole unit cell)"
}