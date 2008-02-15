tracevar ntyp w {    
    widgetconfigure amass -end [varvalue ntyp]
}

tracevar ldisp w {
    groupwidget ramanthreshold disable
    if {  [varvalue ldisp] == ".true." } {
	widget lraman  disable
	widget elop    disable
	widget trans   disable
	widget lnscf   disable
	groupwidget xq_line disable
	groupwidget nq enable
    } else {
	widget lraman  enable
	widget elop    enable
	widget trans   enable
	widget lnscf   enable
	groupwidget xq_line enable
	groupwidget nq disable
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
    
    widget maxirr   disable
    widget nrapp    disable
    widget nat_todo disable	        
    groupwidget representation_line disable 
    groupwidget atom_disp_line      disable 
    
    switch -- [varvalue reps_type] {
	maxirr {
	    widget maxirr enable	    
	}
	nrapp {
	    widget nrapp enable
	    groupwidget representation_line enable
	}
	nat_todo {
	    widget nat_todo enable
	    groupwidget atom_disp_line enable
	}
    }
}

tracevar nrapp w {
    if { [varvalue nrapp] < 1 } {
	groupwidget representation_line disable 
    } else {
	groupwidget representation_line enable
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
    #varset nrapp  -value   0
    #varset nat_todo -value 0
    varset reps_type -value {}
}