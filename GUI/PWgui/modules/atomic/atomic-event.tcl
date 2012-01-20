tracevar iswitch w {
    switch -exact -- [varvalue iswitch] {
        1 {
	    # all-electron calculation
            groupwidget inputp     disable
            groupwidget test       disable
            groupwidget PP_cards   disable	    
	    groupwidget test_cards disable

	    if { [string trim [varvalue config]] == {} } {
		groupwidget AE_cards enable
	    }
	    widget sic enable	    
        }
        2 {
	    # PP test calculation
            groupwidget inputp   disable
            groupwidget AE_cards disable	    
            groupwidget PP_cards disable

            groupwidget test enable
	    groupwidget test_cards enable; # this enables all test_conf_#; disable those with indices > nconf
	    for {set i [expr [varvalue nconf] + 1]} {$i <= $::pwscf::atomic_max_nconf} {incr i} {
		groupwidget test_conf_$i forget
	    }
    
            widget sic disable
        }
        3 {
	    # PP generation
            groupwidget AE_cards disable	    

            groupwidget inputp     enable
            groupwidget PP_cards   enable
            groupwidget test       disable
	    groupwidget test_cards disable; # this enables all test_conf_#; disable those with indices > nconf
	    for {set i [expr [varvalue nconf] + 1]} {$i <= $::pwscf::atomic_max_nconf} {incr i} {
		groupwidget test_conf_$i forget
	    }

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
    set rel [varvalue rel]

    # AE cards and lsd

    if { $rel == 2 } {
        widget lsd disable
	
	widgetconfigure AE_wfs \
	    -head {Label N L Occupancy "Tot.ang.moment"} \
	    -cols 5
	
	for {set ic 1} {$ic <= $::pwscf::atomic_max_nconf} {incr ic} {
	    widgetconfigure test_wfs_$ic \
		-cols 8 \
		-head {Label N L Occupancy enerts rcutts rcutusts "Tot.ang.moment"}
	}
    } else {
        widget lsd enable

	widgetconfigure AE_wfs \
	    -head   {Label N L Occupancy "Spin index"}
	
	# force the correct number of columns
	varset lsd -value [varvalue lsd]
    }
	
    # PP cards

    if { $rel == 0 || $rel == 2 } {
	widgetconfigure PP_wfs \
	    -head {Label N L Occupancy Energy Rcut "US Rcut" "Tot.ang.moment"} \
	    -cols 8	
    } else {
	widgetconfigure PP_wfs \
	    -head {Label N L Occupancy Energy Rcut "US Rcut"} \
	    -cols 7
    }
}

tracevar nwfs w {
    # wfc is table
    widgetconfigure PP_wfs -rows [varvalue nwfs]
}


tracevar nconf w {
    set nconf     [varvalue nconf]
    set old_nconf [varvalue old_nconf]

    if { $nconf == {} || $old_nconf == {} } { return }

    widgetconfigure configts -end $nconf
    widgetconfigure lsdts    -end $nconf
    
    if { $old_nconf > $nconf } {
	# delete ...

	for {set i $old_nconf} { $i > $nconf } {incr i -1} {
	    groupwidget test_conf_$i forget
	}
    } else {
	# create ...
	for {set i $old_nconf} { $i <= $nconf } {incr i} {
	    groupwidget test_conf_$i create
	}
    }

    #
    varset old_nconf -value $nconf
}

for {set ic 1} {$ic <= $::pwscf::atomic_max_nconf} {incr ic} {
    tracevar nwfts_$ic w [subst -nocommands {
	widgetconfigure test_wfs_$ic -rows [varvalue nwfts_$ic]
    }]
}

# manage correctly cards 1.1
tracevar config w {
    set config [string trim [varvalue config]]    
    if { $config == {} && [varvalue iswitch] == 1 } {
	groupwidget AE_cards enable 
    } else {
	groupwidget AE_cards disable 
    }
}
tracevar nwf w {
    widgetconfigure AE_wfs -rows [varvalue nwf]
}
tracevar lsd w {
    if { [varvalue lsd] == 1 } {

	widgetconfigure AE_wfs -cols 5	

	for {set ic 1} {$ic <= $::pwscf::atomic_max_nconf} {incr ic} {
	    widgetconfigure test_wfs_$ic \
		-cols 8 \
		-head {Label N L Occupancy enerts rcutts rcutusts "Spin index"}
	}
    } else {
	widgetconfigure AE_wfs -cols 4
	for {set ic 1} {$ic <= $::pwscf::atomic_max_nconf} {incr ic} {
	    widgetconfigure test_wfs_$ic -cols 7
	}
    }
}

tracevar rmatch_augfun_nc w {
    if { [vartextvalue rmatch_augfun_nc] == "Yes" } {
	widget rmatch_augfun disable
    } else {
	if { [varvalue iswitch] == 3 } {
	    widget rmatch_augfun enable
	}
    }
}

postprocess {
    varset old_nconf -value $::pwscf::atomic_max_nconf; # this is dirty
    varset nconf -value 1

    varset iswitch -value 1
    varset rel -value 1
    varset dft -value 'PZ'
    varset nld -value 0
    varset lloc -value -1
    varset nlcc -value .false.
    varset lpaw -value .false.
    varset tm   -value .false.

    varset rmatch_augfun_nc -value {}
}