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
    if {  [varvalue lraman] == ".true." } {
	groupwidget ramanthreshold enable
    }
}
tracevar elop w {
    if {  [varvalue elop] == ".true." } {
	groupwidget ramanthreshold enable
    }
}
postprocess {
    varset ldisp  -value .false.
    varset lraman -value .false.
    varset elop   -value .false.
    varset trans  -value .true.
    varset epsil  -value .false.
    varset elph   -value .false.
}