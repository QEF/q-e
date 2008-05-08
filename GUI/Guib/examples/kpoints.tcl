module kpoints -title "An example GUI" -script {
    
    line ktype -name "K-point input" {        
        keyword kpoints K_POINTS        
        var kpoint_type {
            -label     "K-Point input" 
            -textvalue { "Automatic generation" "Gamma point only" }
            -value     { automatic gamma }            
            -widget    radiobox
        }
    }

    packwidgets left
    line kmesh -name "K-point mesh" {
        var nk1  -label "nk1:"  -widget entry  -validate posint  -default 1
        var nk2  -label "nk2:"  -widget entry  -validate posint  -default 1
        var nk3  -label "nk3:"  -widget entry  -validate posint  -default 1
    }
    #line kshift -name "K-point mesh shift" {
    #    var s1 -label "s1:"  -widget entry  -validate posint  -default 1
    #    var s2 -label "s2:"  -widget entry  -validate posint  -default 1
    #    var s3 -label "s3:"  -widget entry  -validate posint  -default 1
    #}
    line kshift -name "K-point mesh shift" {
        var s1 -label "s1:"  -widget scale -default 1
        var s2 -label "s2:"  -widget scale -default 1
        var s3 -label "s3:"  -widget scale -default 1
    }

    tracevar kpoint_type w {
    	set value [vartextvalue kpoint_type]    
    	if { $value == "Automatic generation" } {
    	    set status enable
    	} else {
    	    set status disable
    	}
    	
    	groupwidget kmesh $status
    	groupwidget kshift $status
    }

    postprocess {
	foreach scale {s1 s2 s3} {
	    widgetconfigure $scale -from 0 -to 1 -resolution 1
	}
    }

    help s1 -helptext {This is shift in first direction}
}
