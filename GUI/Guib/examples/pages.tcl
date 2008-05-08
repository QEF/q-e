module pages -title "Testing GUIB implemantation: test No.1" -script {

    # support for undefined variables in the input files
    set ::guib::settings(NAMELIST.variable_support_undefined) 1

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
	    page p1.1 -name "Inside Page: feel " {
		var feel {
		    -label     "How do you feel:"
		    -widget    radiobox
		    -value     {0 1 2}
		    -textvalue {bad good excellent}
		    -default   bad
		}
	    }
	    page p1.2 -name "Inside Page: day" {
		var day {
		    -label     "Which day is today:"
		    -widget    optionmenu
		    -value     {0 1 2 3 4 5 6 7}
		    -textvalue {Mon Tue Wed Thu Fri Sat Sun}
		    -default   Wed
		}
	    }
	}
    }

    # Page #.2

    page p2 -name "Page.2" {
	namelist n1 -name first {
	    foreach {var} {v1 v2 v11 v12} label {"First var:" "Second var:" "3rd var:" "4th var:"} {
		var $var -label $label
	    }
	    table t1 -caption "Test table" -head "X Y Z N" -rows 2 -cols 4
	    packwidgets left
	    dimension d1 -variable dim1 -label "Dimension \#.1" -start 1 -end 3	    
	}
	line l2 -name "Table" {
	    table t_line -caption "Simple Table" -head "c1 c2 c3" -rows 3 -cols 3	    
	}
	table t12 -caption "Simple Table \#.2" -head "c1 c2 c3" -rows 2 -cols 3 	
    }
    
    
    # Page #.3
    
    page p3 -name "Page.3" {
	namelist n2 -name second {
	    var v3 -label "3rd var:"
	    var v4 -label "4th var:"
	    table t2 -caption "Test table" -head "X Y Z N" -rows 2 -cols 4	    
	    dimension d2 -variable dim1 -label "Dimension \#.3" -start 1 -end 3	    	
	}
    }
    
    #
    # specialities
    #
    tracevar choo w {
	switch -- [varvalue choo] {
	    "Enable Page.2" {
		groupwidget p2 enable
		groupwidget p3 disable
	    }
	    "Enable Page.3" {
		groupwidget p2 disable
		groupwidget p3 enable
	    }
	    "Enable Pages 2&3" {
		groupwidget p2 enable
		groupwidget p3 enable
	    }
	    "Disable Pages 2&3" {
		groupwidget p2 disable
		groupwidget p3 disable
	    }
	}
    }

    postprocess {
	varset choo -textvalue "Select a value here"
    }

    help feel -helpfmt txt2html -helptext {
	Describe how you weel.
    }
    help day -helpfmt txt2html -helptext {
	Tell me which day it is.
    }
}