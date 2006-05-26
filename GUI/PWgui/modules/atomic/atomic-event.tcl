tracevar iswitch w {
    switch -exact -- [varvalue iswitch] {
        1 {
            groupwidget inputp disable
            groupwidget test disable
            groupwidget pseudization disable
            widget sic enable
        }
        2 {
            groupwidget inputp disable
            groupwidget test enable
            groupwidget pseudization disable
            widget sic disable
        }
        3 {
            groupwidget inputp enable
            groupwidget test enable
            groupwidget pseudization enable
            widget sic disable
        }
    }
}

tracevar dft w {    
    if {  [varvalue dft] == "'REPLACE_ME'" } {
        widget dft_ enable
    } else {
        widget dft_ disable
    }
}

tracevar nld w {    
    if {  [varvalue nld] == 0 } {
        widget rlderiv disable
        widget eminld  disable
        widget emaxld  disable
        widget deld    disable
    } else {
        widget rlderiv enable
        widget eminld  enable
        widget emaxld  enable
        widget deld    enable
    }
}

tracevar lloc w {    
    if {  [varvalue lloc] == -1 } {
        widget rcloc enable
    } else {
        widget rcloc disable
    }
}

tracevar nlcc w {    
    if {  [varvalue nlcc] == ".true." } {
        widget rcore enable
    } else {
        widget rcore disable
    }
}
tracevar lpaw w {    
    if {  [varvalue lpaw] == ".true." } {
        widget file_recon enable
    } else {
        widget file_recon disable
    }
}
tracevar tm w {    
    if {  [varvalue tm] == ".false." } {
        widget rho0 enable
    } else {
        widget rho0 disable
    }
}

tracevar rel w {    
    if {  [varvalue rel] == 2 } {
        widget lsd disable
    } else {
        widget lsd enable
    }
}

tracevar nwfs w {
    # wfc is table
    widgetconfigure wfs -rows [varvalue nwfs]
}


tracevar nconf w {
    # configts is dimension
    widgetconfigure configts -end [varvalue nconf]
}

postprocess {
    varset iswitch -value 1
    varset rel -value 1
    varset dft -value 'PZ'
    varset nld -value 0
    varset lloc -value -1
    varset nlcc -value .false.
    varset lpaw -value .false.
    varset tm   -value .false.
    varset zval  -value 0.0
    varset ecutmin -value 0.0
    varset ecutmax -value 0.0
}