tracevar tdosinboxes w {
    if { [vartextvalue tdosinboxes] == "Yes" } {
	groupwidget local_dos enable
    } else {
	groupwidget local_dos disable
    }
}

postprocess {
    varset tdosinboxes  -value {}
}