tracevar CI_scheme w {
    if { [varvalue string_method] == "'neb'" && [varvalue CI_scheme] == "'manual'" } {
	groupwidget climbing_images enable
    } else {
	groupwidget climbing_images disable
    }
}

postprocess {
    varset string_method   -value 'neb'
    varset CI_scheme       -value {}
}
