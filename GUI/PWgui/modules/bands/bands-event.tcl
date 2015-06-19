tracevar lsym w {
    if { [vartextvalue lsym] == "Yes" } {
	groupwidget kpoints enable
    } else {
	groupwidget kpoints disable
    }
}

tracevar lp w {
    if { [vartextvalue lp] == "Yes" } {
	widget filp enable
    } else {
	widget filp disable
    }
}

postprocess {
    varset lsym  -value {}
    varset lp    -value {}
}
