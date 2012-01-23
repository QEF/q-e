tracevar lsym w {
    if { [vartextvalue lsym] == "Yes" } {
	groupwidget kpoints enable
    } else {
	groupwidget kpoints disable
    }
}

postprocess {
    varset lsym  -value {}
}