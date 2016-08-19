tracevar ntyp w {    
    widgetconfigure amass -end [varvalue ntyp]
}

tracevar ldisp w {
    groupwidget ramanthreshold disable
    if {  [varvalue ldisp] == ".true." } {
	widget lraman  disable
	widget elop    disable
	widget trans   disable
	groupwidget xq_line disable
	groupwidget nq enable
	groupwidget q_spec enable
    } else {
	widget lraman  enable
	widget elop    enable
	widget trans   enable
	groupwidget xq_line enable
	groupwidget nq disable
	groupwidget q_spec disable
    }
}

tracevar lraman w {
    if {  [varvalue lraman] == ".true." || [varvalue elop] == ".true." } {
	groupwidget ramanthreshold enable
    } else {
	groupwidget ramanthreshold disable
    }
}

tracevar elop w {
    if {  [varvalue lraman] == ".true." || [varvalue elop] == ".true." } {
	groupwidget ramanthreshold enable
    } else {
	groupwidget ramanthreshold disable
    }
}

tracevar reps_type w {
    
    widget nat_todo disable	        
    groupwidget start_last_irr      disable
    groupwidget representation_line disable 
    groupwidget atom_disp_line      disable 
    
    switch -- [varvalue reps_type] {
	start_last_irr {
	    groupwidget start_last_irr enable
	}
	nat_todo {
	    widget nat_todo enable
	    groupwidget atom_disp_line enable
	}
    }
}

tracevar nat_todo w {
    if { [varvalue nat_todo] == "" } {
	widget nat_todo_list disable
    } else {
	if { [varvalue nat_todo] < 1 } {
	    widget nat_todo_list disable 
	} else {
	    widget nat_todo_list enable
	}
    }
}


tracevar qplot w {
    if { [varvalue qplot] == ".true." } {
	groupwidget qPointsSpec enable
	groupwidget xq_list     disable
    } else {
	groupwidget qPointsSpec disable
	groupwidget xq_list     enable
    }
}

tracevar nqs w {
    set nqs [varvalue nqs]
    widgetconfigure qPoints  -rows $nqs
}

tracevar dvscf_star_open w {
    if { [varvalue dvscf_star_open] == ".true." } {
	groupwidget dvscf_star_specs enable
    } else {
	groupwidget dvscf_star_specs disable
    }
}
tracevar drho_star_open w {
    if { [varvalue drho_star_open] == ".true." } {
	groupwidget drho_star_specs enable
    } else {
	groupwidget drho_star_specs disable
    }
}

## help postproccessing (hack for help of dvscf_star & drho_star structures)
#    
#foreach ident {dvscf_star drho_star} {
#    
#    set obj      [_getObjFromVarident $ident]
#    set id       [$obj getIdFromVarident $ident]
#    set helptext [$obj getOptionValue $id helptext]
#    
#    foreach elem {open dir ext basis pat} {
#	help ${ident}_$elem -helpfmt helpdoc -helptext $helptext
#    }
#}


postprocess {
    widget dvscf_star forget
    widget drho_star  forget

    varset ldisp     -value  {}; #.false.
    varset lraman    -value  {}; #.false.
    varset elop      -value  {}; #.false.
    varset trans     -value  .true.
    varset epsil     -value  {}; #.false.
    varset recover   -value  {}; #.false.
    varset fpol      -value  {}; #.false.
    varset reps_type -value  {}
    varset nat_todo  -value  {}
    varset qplot     -value  {}
    varset dvscf_star_open -value {}
    varset drho_star_open -value {}
}
