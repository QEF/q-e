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
    varset nrapp  -value   0
    varset nat_todo -value 0
}