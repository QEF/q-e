tracevar ntyp w {
    set ntyp [varvalue ntyp]
    foreach w {skip_type equiv_type perturb_only_atom} {
        widgetconfigure $w -end $ntyp
    }
}

tracevar only_q w {
    if { [vartextvalue only_q] == "Yes" } {
        groupwidget start_last_q enable
        widget sum_pertq enable
    } else {
        groupwidget start_last_q disable
        widget sum_pertq disable
    }     
}

tracevar start_q w {
    if { [varvalue start_q] != "" } {
        varset only_q -value .true.
        groupwidget start_last_q enable
        widget sum_pertq enable
    }
}
tracevar last_q w {
    if { [varvalue last_q] != "" } {
        varset only_q -value .true.
        groupwidget start_last_q enable
        widget sum_pertq enable
    }
}
