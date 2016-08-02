namespace eval ::helpdoc {}

# Purpose: definition of QE modules that PWgui knows about
proc ::helpdoc::core_modules {} {
    variable mod_def_mappings

    # this is the list of "core" mod_def_mappings, i.e., modules (QE
    # programs) which the PWgui knows
    
    set mod_def_mappings {
	pw        PW/Doc     INPUT_PW 
	neb       NEB/Doc    INPUT_NEB
	ph        PHonon/Doc INPUT_PH
	pp        PP/Doc     INPUT_PP
	projwfc   PP/Doc     INPUT_PROJWFC
	bands     PP/Doc     INPUT_BANDS
	dos       PP/Doc     INPUT_DOS
	atomic    atomic/Doc INPUT_LD1
	d3        PHonon/Doc INPUT_D3
    }
}


###
# Purpose: print the usage message for utility programs that requires "$0 module" execution
proc ::helpdoc::prog_module_usage {} {    
    global argv0
    
    puts stderr "
  Usage: $argv0 module
  
  Where module is one of:
  
      [::helpdoc::get_supported_modules]
    "
    exit 1
}


# Purpose: get a list of all modules that PWgui supports
proc ::helpdoc::get_supported_modules {} {
    variable mod_def_mappings

    if { ! [info exists mod_def_mappings] } {
	core_modules
    }

    set modules ""
    foreach {mod subdir def_prefix} $mod_def_mappings {
	lappend modules $mod
    }
    
    return $modules
}


##
# Purpose: get the pathname of INPUT_XX.def file from the module name
proc ::helpdoc::get_def_filename {module {must_exists 1}} {
    variable mod_def_mappings
    global topdir

    if { ! [info exists topdir] } {
	puts stderr "variable \"topdir\" does not exist; define it before calling get_def_filename"
	exit 1
    }
    if { ! [info exists mod_def_mappings] } {
	core_modules
    }

    set deffile ""
    
    foreach {mod subdir def_prefix} $mod_def_mappings {
	if { $mod == $module } {
	    set deffile [file join $topdir $subdir $def_prefix.def]
	}
    }

    if { $must_exists == 1 && "$deffile" eq {} } {
	puts stderr "\[failed\]

module \"$module\" is not defined by \"mod_def_mappings\" variable,
list of defined modules:   
                          [get_supported_modules]
"
	exit 1	
    }
    
    return $deffile
}

##
# Purpose: get the pathname of PWgui's module file from the module name
proc ::helpdoc::get_gui_module_filename {module {must_exists 1}} {
    variable mod_def_mappings
    global moduledir
    
    if { ! [info exists moduledir] } {
	puts stderr "variable \"moduledir\" does not exist; define it before calling get_def_filename"
	exit 1
    }
    if { ! [info exists mod_def_mappings] } {
	core_modules
    }
    
    set modfile ""
    
    foreach {mod subdir def_prefix} $mod_def_mappings {
	if { $mod == $module } {
	    set modfile [file join  $moduledir $mod $mod.tcl]
	}
    }
    
    if { $must_exists == 1 && "$modfile" eq {} } {
	puts stderr "\[failed\]

module \"$module\" is not defined by \"mod_def_mappings\" variable,
list of defined modules:   
                         [get_supported_modules]
"
	exit 1	
    }
    
    return $modfile
}
