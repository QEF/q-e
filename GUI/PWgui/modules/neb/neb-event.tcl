tracevar CI_scheme w {
    if { [varvalue string_method] == "'neb'" && [varvalue CI_scheme] == "'manual'" } {
	groupwidget climbing_images enable
    } else {
	groupwidget climbing_images disable
    }
}

tracevar lfcpopt w {
    if { [vartextvalue lfcpopt] == "Yes" } {
	groupwidget fcp_specs enable
    } else {
	groupwidget fcp_specs disable
    }
}

postprocess {
    varset string_method   -value 'neb'
    varset CI_scheme       -value {}
    varset lfcpopt         -value {}
}
