module two_pages -title "Testing GUIB implemantation: testing pages" -script {
    
    # Page #.2
    
    page p2 -name "Page.2" {
	namelist n1 -name first {
	    puts stderr "*** nml"
	    foreach {var} {v1 v2 v11 v12} label {"First var:" "Second var:" "3rd var:" "4th var:"} {
		puts stderr "*** var $var"
		var $var -label $label
	    }
	}
    }

    # Page #.1
    
    page p1 -name "Page.1" {
    	line l1 -name "choose" {
    	    var choo {
    		-label     "What to do:"
    		-variable  choose
    		-widget    optionmenu
    		-textvalue {
    		    "Enable Page.2" 
    		    "Enable Page.3" 
    		    "Enable Pages 2&3"
    		    "Disable Pages 2&3"
    		}
    		-value     {0 1 2 3}
    	    }
    	}
    }
}