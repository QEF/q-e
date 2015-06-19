tracevar tdosinboxes w {
    if { [vartextvalue tdosinboxes] == "Yes" } {
	groupwidget local_dos enable
    } else {
	groupwidget local_dos disable
    }
}

tracevar n_proj_boxes w {
    set nbox [varvalue  n_proj_boxes]

    widgetconfigure irmin -rows $nbox
    widgetconfigure irmax -rows $nbox
}

postprocess {
    varset tdosinboxes  -value {}
}
