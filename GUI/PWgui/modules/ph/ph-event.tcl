tracevar ntyp w {    
    widgetconfigure amass -end [varvalue ntyp]
}

tracevar ldisp w {
    groupwidget ramanthreshold disable
    if {  [varvalue ldisp] == ".true." } {
	widget lraman  disable
	widget elop    disable
	widget trans   disable
	groupwidget xq_line disable
	groupwidget nq enable
	groupwidget q_spec enable
    } else {
	widget lraman  enable
	widget elop    enable
	widget trans   enable
	groupwidget xq_line enable
	groupwidget nq disable
	groupwidget q_spec disable
    }
}

tracevar lraman w {
    if {  [varvalue lraman] == ".true." || [varvalue elop] == ".true." } {
	groupwidget ramanthreshold enable
    } else {
	groupwidget ramanthreshold disable
    }
}

tracevar elop w {
    if {  [varvalue lraman] == ".true." || [varvalue elop] == ".true." } {
	groupwidget ramanthreshold enable
    } else {
	groupwidget ramanthreshold disable
    }
}

tracevar reps_type w {
    
    widget nat_todo disable	        
    groupwidget start_last_irr      disable
    groupwidget representation_line disable 
    groupwidget atom_disp_line      disable 
    
    switch -- [varvalue reps_type] {
	start_last_irr {
	    groupwidget start_last_irr enable
	}
	nat_todo {
	    widget nat_todo enable
	    groupwidget atom_disp_line enable
	}
    }
}

tracevar nat_todo w {
    if { [varvalue nat_todo] < 1 } {
	groupwidget atom_disp_line disable 
    } else {
	groupwidget atom_disp_line enable
    }
}

postprocess {
    varset ldisp  -value  .false.
    varset lraman -value  .false.
    varset elop   -value  .false.
    varset trans  -value  .true.
    varset epsil  -value  .false.
    varset elph   -value  .false.
    varset recover -value .false.
    varset fpol   -value  .false.
    #varset nat_todo -value 0
    varset reps_type -value {}
}